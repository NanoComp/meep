/* Copyright (C) 2005-2015 Massachusetts Institute of Technology
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2, or (at your option)
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex>

#include "meep.hpp"
#include "meep_internals.hpp"

using namespace std;

namespace meep {

fields::fields(structure *s, double m, double beta,
	       bool zero_fields_near_cylorigin) :
  S(s->S), gv(s->gv), user_volume(s->user_volume), v(s->v), m(m), beta(beta)
{
  shared_chunks = s->shared_chunks;
  verbosity = 0;
  components_allocated = false;
  synchronized_magnetic_fields = 0;
  outdir = new char[strlen(s->outdir) + 1]; strcpy(outdir, s->outdir);
  if (gv.dim == Dcyl)
    S = S + r_to_minus_r_symmetry(m);
  phasein_time = 0;
  for (int d=0;d<5;d++) { k[d] = 0.0; eikna[d] = 1.0; }
  is_real = 0;
  a = gv.a;
  dt = s->dt;
  t = 0;
  sources = NULL;
  fluxes = NULL;
  // Time stuff:
  for (int i = 0; i < MEEP_TIMING_STACK_SZ; ++i) was_working_on[i] = Other;
  working_on = Other;
  for (int i=0;i<=Other;i++) times_spent[i] = 0.0;
  last_wall_time = last_step_output_wall_time = -1;
  am_now_working_on(Other);

  num_chunks = s->num_chunks;
  typedef fields_chunk *fields_chunk_ptr;
  chunks = new fields_chunk_ptr[num_chunks];
  for (int i=0;i<num_chunks;i++)
    chunks[i] = new fields_chunk(s->chunks[i], outdir, m,
				 beta, zero_fields_near_cylorigin);
  FOR_FIELD_TYPES(ft) {
    for (int ip=0;ip<3;ip++) {
      comm_sizes[ft][ip] = new size_t[num_chunks*num_chunks];
      for (int i=0;i<num_chunks*num_chunks;i++) comm_sizes[ft][ip][i] = 0;
    }
    typedef realnum *realnum_ptr;
    comm_blocks[ft] = new realnum_ptr[num_chunks*num_chunks];
    for (int i=0;i<num_chunks*num_chunks;i++)
      comm_blocks[ft][i] = 0;
  }
  for (int b=0;b<2;b++) FOR_DIRECTIONS(d)
    if (gv.has_boundary((boundary_side)b, d)) boundaries[b][d] = Metallic;
    else boundaries[b][d] = None;
  chunk_connections_valid = false;

  // unit directions are periodic by default:
  FOR_DIRECTIONS(d)
    if (gv.has_boundary(High, d) && gv.has_boundary(Low, d) && d != R
	&& s->user_volume.num_direction(d) == 1)
      use_bloch(d, 0.0);
}

fields::fields(const fields &thef) :
  S(thef.S), gv(thef.gv), user_volume(thef.user_volume), v(thef.v)
{
  shared_chunks = thef.shared_chunks;
  verbosity = 0;
  components_allocated = thef.components_allocated;
  synchronized_magnetic_fields = thef.synchronized_magnetic_fields;
  outdir = new char[strlen(thef.outdir) + 1]; strcpy(outdir, thef.outdir);
  m = thef.m;
  beta = thef.beta;
  phasein_time = thef.phasein_time;
  for (int d=0;d<5;d++) { k[d] = thef.k[d]; eikna[d] = thef.eikna[d]; }
  is_real = thef.is_real;
  a = thef.a;
  dt = thef.dt;
  t = thef.t;
  sources = NULL;
  fluxes = NULL;
  // Time stuff:
  for (int i = 0; i < MEEP_TIMING_STACK_SZ; ++i) was_working_on[i] = Other;
  working_on = Other;
  for (int i=0;i<=Other;i++) times_spent[i] = 0.0;
  last_wall_time = -1;
  am_now_working_on(Other);

  num_chunks = thef.num_chunks;
  typedef fields_chunk *fields_chunk_ptr;
  chunks = new fields_chunk_ptr[num_chunks];
  for (int i=0;i<num_chunks;i++)
    chunks[i] = new fields_chunk(*thef.chunks[i]);
  FOR_FIELD_TYPES(ft) {
    for (int ip=0;ip<3;ip++) {
      comm_sizes[ft][ip] = new size_t[num_chunks*num_chunks];
      for (int i=0;i<num_chunks*num_chunks;i++) comm_sizes[ft][ip][i] = 0;
    }
    typedef realnum *realnum_ptr;
    comm_blocks[ft] = new realnum_ptr[num_chunks*num_chunks];
    for (int i=0;i<num_chunks*num_chunks;i++)
      comm_blocks[ft][i] = 0;
  }
  for (int b=0;b<2;b++) FOR_DIRECTIONS(d)
    boundaries[b][d] = thef.boundaries[b][d];
  chunk_connections_valid = false;
}

fields::~fields() {
  for (int i=0;i<num_chunks;i++) delete chunks[i];
  delete[] chunks;
  FOR_FIELD_TYPES(ft) {
    for (int i=0;i<num_chunks*num_chunks;i++)
      delete[] comm_blocks[ft][i];
    delete[] comm_blocks[ft];
    for (int ip=0;ip<3;ip++)
      delete[] comm_sizes[ft][ip];
  }
  delete sources;
  delete fluxes;
  delete[] outdir;
  if (!quiet) print_times();
}

void fields::verbose(int gv) {
  verbosity = gv;
  for (int i=0;i<num_chunks;i++) chunks[i]->verbose(gv);
}

void fields::use_real_fields() {
  LOOP_OVER_DIRECTIONS(gv.dim, d)
    if (boundaries[High][d] == Periodic && k[d] != 0.0)
      abort("Can't use real fields with bloch boundary conditions!\n");
  is_real = 1;
  for (int i=0;i<num_chunks;i++) chunks[i]->use_real_fields();
  chunk_connections_valid = false;
}

bool fields::have_component(component c) {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->f[c][0])
      return true;
  return false;
}

fields_chunk::~fields_chunk() {
  is_real = 0; // So that we can make sure to delete everything...
  // for mu=1 non-PML regions, H==B to save space/time - don't delete twice!
  DOCMP2 FOR_H_AND_B(hc,bc) if (f[hc][cmp] == f[bc][cmp]) f[bc][cmp] = NULL;
  DOCMP2 FOR_COMPONENTS(c) {
    delete[] f[c][cmp];
    delete[] f_u[c][cmp];
    delete[] f_w[c][cmp];
    delete[] f_cond[c][cmp];
    delete[] f_minus_p[c][cmp];
    delete[] f_w_prev[c][cmp];
    delete[] f_backup[c][cmp];
    delete[] f_u_backup[c][cmp];
    delete[] f_w_backup[c][cmp];
    delete[] f_cond_backup[c][cmp];
  }
  delete[] f_rderiv_int;
  FOR_FIELD_TYPES(ft)
    for (int ip=0;ip<3;ip++)
      for (int io=0;io<2;io++)
	delete[] connections[ft][ip][io];
  FOR_FIELD_TYPES(ft) delete[] connection_phases[ft];
  while (dft_chunks) {
    dft_chunk *nxt = dft_chunks->next_in_chunk;
    delete dft_chunks;
    dft_chunks = nxt;
  }
  FOR_FIELD_TYPES(ft) {
    delete sources[ft];
    delete[] zeroes[ft];
  }
  FOR_FIELD_TYPES(ft) for (polarization_state *cur = pol[ft]; cur; ) {
    polarization_state *p = cur;
    cur = cur->next;
    p->s->delete_internal_data(p->data);
    delete p;
  }
  if (s->refcount-- <= 1) delete s; // delete if not shared
  if (new_s && new_s->refcount-- <= 1) delete new_s; // delete if not shared
}

fields_chunk::fields_chunk(structure_chunk *the_s, const char *od,
			   double m, double beta,
			   bool zero_fields_near_cylorigin) : gv(the_s->gv), v(the_s->v), m(m), zero_fields_near_cylorigin(zero_fields_near_cylorigin), beta(beta) {
  s = the_s; s->refcount++;
  verbosity = 0;
  outdir = od;
  new_s = NULL;
  is_real = 0;
  a = s->a;
  Courant = s->Courant;
  dt = s->dt;
  dft_chunks = NULL;
  FOR_FIELD_TYPES(ft) {
    polarization_state *cur = NULL;
    pol[ft] = NULL;
    for (susceptibility *chiP = the_s->chiP[ft]; chiP; chiP = chiP->next) {
      polarization_state *p = new polarization_state;
      // P and data lazily allocated in update_pols
      p->data = NULL;
      p->s = chiP;
      p->next = NULL;
      if (cur) { cur->next = p; cur = p; }
      else { pol[ft] = cur = p; }
    }
  }
  doing_solve_cw = false;
  solve_cw_omega = 0.0;
  FOR_FIELD_TYPES(ft) sources[ft] = NULL;
  FOR_COMPONENTS(c) DOCMP2 {
    f[c][cmp] = NULL;
    f_u[c][cmp] = NULL;
    f_w[c][cmp] = NULL;
    f_cond[c][cmp] = NULL;
    f_minus_p[c][cmp] = NULL;
    f_w_prev[c][cmp] = NULL;
    f_backup[c][cmp] = NULL;
    f_u_backup[c][cmp] = NULL;
    f_w_backup[c][cmp] = NULL;
    f_cond_backup[c][cmp] = NULL;
  }
  f_rderiv_int = NULL;
  FOR_FIELD_TYPES(ft) {
    for (int ip=0;ip<3;ip++)
      num_connections[ft][ip][Incoming]
	= num_connections[ft][ip][Outgoing] = 0;
    connection_phases[ft] = 0;
    for (int ip=0;ip<3;ip++) for (int io=0;io<2;io++)
      connections[ft][ip][io] = NULL;
    zeroes[ft] = NULL;
    num_zeroes[ft] = 0;
  }
  figure_out_step_plan();
}

fields_chunk::fields_chunk(const fields_chunk &thef)
  : gv(thef.gv), v(thef.v) {
  s = thef.s; s->refcount++;
  verbosity = thef.verbosity;
  outdir = thef.outdir;
  m = thef.m;
  zero_fields_near_cylorigin = thef.zero_fields_near_cylorigin;
  beta = thef.beta;
  new_s = thef.new_s; new_s->refcount++;
  is_real = thef.is_real;
  a = thef.a;
  Courant = thef.Courant;
  dt = thef.dt;
  dft_chunks = NULL;
  FOR_FIELD_TYPES(ft) {
    polarization_state *cur = NULL;
    for (polarization_state *ocur = thef.pol[ft]; ocur; ocur = ocur->next) {
      polarization_state *p = new polarization_state;
      p->data = NULL;
      p->s = ocur->s;
      p->next = NULL;
      pol[ft] = NULL;
      if (ocur->data) p->data = p->s->copy_internal_data(p->data);
      if (cur) { cur->next = p; cur = p; }
      else { pol[ft] = cur = p; }
    }
  }
  doing_solve_cw = thef.doing_solve_cw;
  solve_cw_omega = thef.solve_cw_omega;
  FOR_FIELD_TYPES(ft) sources[ft] = NULL;
  FOR_COMPONENTS(c) DOCMP2 {
    f[c][cmp] = NULL;
    f_u[c][cmp] = NULL;
    f_w[c][cmp] = NULL;
    f_cond[c][cmp] = NULL;
    f_backup[c][cmp] = NULL;
    f_u_backup[c][cmp] = NULL;
    f_w_backup[c][cmp] = NULL;
    f_cond_backup[c][cmp] = NULL;
  }
  FOR_COMPONENTS(c) DOCMP {
    if (!is_magnetic(c) && thef.f[c][cmp]) {
      f[c][cmp] = new realnum[gv.ntot()];
      memcpy(f[c][cmp], thef.f[c][cmp], sizeof(realnum) * gv.ntot());
    }
    if (thef.f_u[c][cmp]) {
      f_u[c][cmp] = new realnum[gv.ntot()];
      memcpy(f_u[c][cmp], thef.f_u[c][cmp], sizeof(realnum) * gv.ntot());
    }
    if (thef.f_w[c][cmp]) {
      f_w[c][cmp] = new realnum[gv.ntot()];
      memcpy(f_w[c][cmp], thef.f_w[c][cmp], sizeof(realnum) * gv.ntot());
    }
    if (thef.f_cond[c][cmp]) {
      f_cond[c][cmp] = new realnum[gv.ntot()];
      memcpy(f_cond[c][cmp], thef.f_cond[c][cmp], sizeof(realnum) * gv.ntot());
    }
  }
  FOR_MAGNETIC_COMPONENTS(c) DOCMP {
    if (thef.f[c][cmp] == thef.f[c-Hx+Bx][cmp])
      f[c][cmp] = f[c-Hx+Bx][cmp];
    else if (thef.f[c][cmp]) {
      f[c][cmp] = new realnum[gv.ntot()];
      memcpy(f[c][cmp], thef.f[c][cmp], sizeof(realnum) * gv.ntot());
    }
  }
  FOR_FIELD_TYPES(ft) {
    for (int ip=0;ip<3;ip++)
      num_connections[ft][ip][Incoming]
	= num_connections[ft][ip][Outgoing] = 0;
    connection_phases[ft] = 0;
    for (int ip=0;ip<3;ip++) for (int io=0;io<2;io++)
      connections[ft][ip][io] = NULL;
    zeroes[ft] = NULL;
    num_zeroes[ft] = 0;
  }
  FOR_COMPONENTS(c) DOCMP2 {
    if (thef.f_minus_p[c][cmp]) {
      f_minus_p[c][cmp] = new realnum[gv.ntot()];
      memcpy(f_minus_p[c][cmp], thef.f_minus_p[c][cmp],
	     sizeof(realnum) * gv.ntot());
    }
    if (thef.f_w_prev[c][cmp]) {
      f_w_prev[c][cmp] = new realnum[gv.ntot()];
      memcpy(f_w_prev[c][cmp], thef.f_w_prev[c][cmp],
	     sizeof(realnum) * gv.ntot());
    }
  }
  f_rderiv_int = NULL;
  figure_out_step_plan();
}

static inline bool cross_negative(direction a, direction b) {
  if (a >= R) a = direction(a - 3);
  if (b >= R) b = direction(b - 3);
  return ((3+b-a)%3) == 2;
}

static inline direction cross(direction a, direction b) {
  if (a == b) abort("bug - cross expects different directions");
  bool dcyl = a >= R || b >= R;
  if (a >= R) a = direction(a - 3);
  if (b >= R) b = direction(b - 3);
  direction c = direction((3+2*a-b)%3);
  if (dcyl && c < Z) return direction(c + 3);
  return c;
}

/* Call this whenever we modify the structure_chunk (fields_chunk::s) to
   implement copy-on-write semantics.  See also structure::changing_chunks. */
void fields_chunk::changing_structure() {
  if (s->refcount > 1) { // this chunk is shared, so make a copy
    s->refcount--;
    s = new structure_chunk(s);
  }
}

void fields::figure_out_step_plan() {
  for (int i = 0; i < num_chunks; ++i)
    if (chunks[i]->is_mine()) chunks[i]->figure_out_step_plan();
}

void fields_chunk::figure_out_step_plan() {
  FOR_COMPONENTS(cc)
    have_minus_deriv[cc] = have_plus_deriv[cc] = false;
  FOR_COMPONENTS(c1)
    if (f[c1][0]) {
      const direction dc1 = component_direction(c1);
      // Figure out which field components contribute.
      FOR_COMPONENTS(c2)
        if ((is_electric(c1) && is_magnetic(c2)) ||
            (is_D(c1) && is_magnetic(c2)) ||
            (is_magnetic(c1) && is_electric(c2)) ||
	    (is_B(c1) && is_electric(c2))) {
          const direction dc2 = component_direction(c2);
          if (dc1 != dc2 && gv.has_field(c2) && gv.has_field(c1) &&
              (has_direction(gv.dim,cross(dc1,dc2)) ||
	       (gv.dim == Dcyl && has_field_direction(gv.dim,cross(dc1,dc2))))) {
            direction d_deriv = cross(dc1,dc2);
            if (cross_negative(dc2, dc1)) {
              minus_component[c1] = c2;
              have_minus_deriv[c1] = true;
              minus_deriv_direction[c1] = d_deriv;
            } else {
              plus_component[c1] = c2;
              have_plus_deriv[c1] = true;
              plus_deriv_direction[c1] = d_deriv;
            }
          }
        }
    }
}

bool is_tm(component c) {
  switch (c) {
  case Hx: case Hy: case Bx: case By: case Ez: case Dz: return true;
  default: return false;
  }
  return false;
}

static bool is_like(ndim d, component c1, component c2) {
  if (d != D2) return true;
  return !(is_tm(c1) ^ is_tm(c2));
}

// this function should ordinarily not be called directly;
// instead it should be called via require_component,
// since only require_component knows what other field components
// need to be allocated in addition to c
bool fields_chunk::alloc_f(component c) {
  bool changed = false;
  if (is_mine())
    DOCMP {
      if (!f[c][cmp]) {
	changed = true;
	if (is_magnetic(c)) {
	  /* initially, we just set H == B ... later on, we lazily allocate
	     H fields if needed (if mu != 1 or in PML) in update_eh */
	  component bc = direction_component(Bx, component_direction(c));
	  if (!f[bc][cmp]) {
	    f[bc][cmp] = new realnum[gv.ntot()];
	    for (size_t i=0;i<gv.ntot();i++) f[bc][cmp][i] = 0.0;
	  }
	  f[c][cmp] = f[bc][cmp];
	}
	else {
	  f[c][cmp] = new realnum[gv.ntot()];
	  for (size_t i=0;i<gv.ntot();i++) f[c][cmp][i] = 0.0;
	}
      }
    }
  return changed;
}

void fields::require_component(component c) {
  if (!gv.has_field(c))
    abort("cannot require a %s component in a %s grid",
	  component_name(c), dimension_name(gv.dim));

  if (beta != 0 && gv.dim != D2)
    abort("Nonzero beta unsupported in dimensions other than 2.");

  components_allocated = true;

  // check if we are in 2d but anisotropy couples xy with z
  bool aniso2d = false;
  if (gv.dim == D2) {
    int i;
    for (i = 0; i < num_chunks; ++i)
      if (chunks[i]->s->has_chi(Ex, Z) ||
	  chunks[i]->s->has_chi(Ey, Z) ||
	  chunks[i]->s->has_chi(Ez, X) ||
	  chunks[i]->s->has_chi(Ez, Y) ||
	  chunks[i]->s->has_chi(Hx, Z) ||
	  chunks[i]->s->has_chi(Hy, Z) ||
	  chunks[i]->s->has_chi(Hz, X) ||
	  chunks[i]->s->has_chi(Hz, Y))
	break;
    aniso2d = or_to_all(i < num_chunks);
  }
  if (aniso2d && beta != 0 && is_real)
    abort("Nonzero beta need complex fields when mu/epsilon couple TE and TM");
  aniso2d = aniso2d || (beta != 0); // beta couples TE/TM

  // allocate fields if they haven't been allocated yet for this component
  int need_to_reconnect = 0;
  FOR_COMPONENTS(c_alloc)
    if (gv.has_field(c_alloc) && (is_like(gv.dim, c, c_alloc) || aniso2d))
      for (int i = 0; i < num_chunks; ++i)
      	if (chunks[i]->alloc_f(c_alloc))
      	  need_to_reconnect++;

  if (need_to_reconnect) figure_out_step_plan();
  if (sum_to_all(need_to_reconnect)) chunk_connections_valid = false;
}

void fields_chunk::remove_sources() {
  FOR_FIELD_TYPES(ft) { delete sources[ft]; sources[ft] = NULL; }
}

void fields::remove_sources() {
  delete sources;
  sources = NULL;
  for (int i=0;i<num_chunks;i++)
    chunks[i]->remove_sources();
}

void fields_chunk::remove_susceptibilities(bool shared_chunks) {
  FOR_FIELD_TYPES(ft) {
    for (polarization_state *cur = pol[ft]; cur; ) {
      polarization_state *p = cur;
      cur = cur->next;
      p->s->delete_internal_data(p->data);
      delete p;
    }
    pol[ft] = NULL;
  }

  if (!shared_chunks) {
    changing_structure();
  }
  s->remove_susceptibilities();
}

void fields::remove_susceptibilities() {
  for (int i=0;i<num_chunks;i++)
    chunks[i]->remove_susceptibilities(shared_chunks);
}

void fields::remove_fluxes() {
  delete fluxes;
  fluxes = NULL;
}

void fields_chunk::zero_fields() {
  FOR_COMPONENTS(c) DOCMP {
#define ZERO(array) if (array) memset(array, 0, sizeof(realnum) * gv.ntot())
    ZERO(f[c][cmp]);
    ZERO(f_u[c][cmp]);
    ZERO(f_w[c][cmp]);
    ZERO(f_cond[c][cmp]);
    ZERO(f_backup[c][cmp]);
    ZERO(f_u_backup[c][cmp]);
    ZERO(f_w_backup[c][cmp]);
    ZERO(f_cond_backup[c][cmp]);
#undef ZERO
  }
  if (is_mine()) FOR_FIELD_TYPES(ft)
      for (polarization_state *p = pol[ft]; p; p = p->next) {
	if (p->data) p->s->init_internal_data(f, dt, gv, p->data);
      }
}

void fields::zero_fields() {
  for (int i=0;i<num_chunks;i++)
    chunks[i]->zero_fields();
}

void fields::reset() {
  remove_sources();
  remove_fluxes();
  zero_fields();
  t = 0;
}

void fields_chunk::use_real_fields() {
  is_real = 1;
  // for mu=1 non-PML regions, H==B to save space/time - don't delete twice!
  FOR_H_AND_B(hc,bc) if (f[hc][1] == f[bc][1]) f[bc][1] = NULL;
  FOR_COMPONENTS(c) if (f[c][1]) {
    delete[] f[c][1];
    f[c][1] = 0;
  }
  if (is_mine()) FOR_FIELD_TYPES(ft)
    for (polarization_state *p = pol[ft]; p; p = p->next) {
      if (p->data) { // TODO: print an error message in this case?
	p->s->delete_internal_data(p->data);
	p->data = p->s->new_internal_data(f, gv);
	p->s->init_internal_data(f, dt, gv, p->data);
      }
    }
}

int fields::phase_in_material(const structure *snew, double time) {
  if (snew->num_chunks != num_chunks)
    abort("Can only phase in similar sets of chunks: %d vs %d\n",
	  snew->num_chunks, num_chunks);
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->phase_in_material(snew->chunks[i]);
  phasein_time = (int) (time/dt);
  // FIXME: how to handle changes in susceptibilities?
  return phasein_time;
}

void fields_chunk::phase_in_material(structure_chunk *snew) {
  new_s = snew; new_s->refcount++;
}

int fields::is_phasing() {
  return phasein_time > 0;
}

bool fields::equal_layout(const fields &f) const {
  if (a != f.a ||
      num_chunks != f.num_chunks ||
      v != f.v ||
      S != f.S)
    return false;
  for (int d=0;d<5;d++)
    if (k[d] != f.k[d])
      return false;
  for (int i = 0; i < num_chunks; ++i)
    if (chunks[i]->a != f.chunks[i]->a ||
	chunks[i]->v != f.chunks[i]->v)
      return false;
  return true;
}

// total computational grid_volume, including regions redundant by symmetry
volume fields::total_volume(void) const {
  volume gv0 = gv.interior();
  volume v = gv0;
  for (int n = 1; n < S.multiplicity(); ++n)
    v = v | S.transform(gv0, n);
  if (v.dim == Dcyl && v.in_direction_min(R) < 0)
    v.set_direction_min(R, 0);
  return v;
}

/* One-pixel periodic dimensions are used almost exclusively to
   emulate lower-dimensional computations, so if the user passes an
   empty size in that direction, they probably really intended to
   specify that whole dimension.  This function detects that case. */
bool fields::nosize_direction(direction d) const {
  return (gv.has_boundary(Low, d) && gv.has_boundary(High, d) &&
	  boundaries[Low][d] == Periodic && boundaries[High][d] == Periodic
	  && gv.num_direction(d) == 1);
}

void fields::set_solve_cw_omega(complex<double> omega) {
  for (int i = 0; i < num_chunks; ++i)
    chunks[i]->set_solve_cw_omega(omega);
}

void fields::unset_solve_cw_omega() {
  for (int i = 0; i < num_chunks; ++i)
    chunks[i]->unset_solve_cw_omega();
}

} // namespace meep
