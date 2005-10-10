/* Copyright (C) 2003 Massachusetts Institute of Technology
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
#include <stdlib.h>
#include <math.h>
#include <complex>

#include "meep.hpp"
#include "meep_internals.hpp"

namespace meep {

fields::fields(const structure *s, double m) :
  S(s->S), v(s->v), user_volume(s->user_volume), gv(s->gv), m(m)
{
  verbosity = 0;
  outdir = s->outdir;
  if (v.dim == Dcyl && m == int(m))
    S = S + r_to_minus_r_symmetry(int(m));
  phasein_time = 0;
  bands = NULL;
  for (int d=0;d<5;d++) k[d] = 0.0;
  is_real = 0;
  a = v.a;
  dt = s->dt;
  t = 0;
  sources = NULL;
  disable_sources = false;
  fluxes = NULL;
  // Time stuff:
  was_working_on = working_on = Other;
  for (int i=0;i<=Other;i++) times_spent[i] = 0.0;
  last_wall_time = last_step_output_wall_time = -1;
  am_now_working_on(Other);

  num_chunks = s->num_chunks;
  typedef fields_chunk *fields_chunk_ptr;
  chunks = new fields_chunk_ptr[num_chunks];
  for (int i=0;i<num_chunks;i++)
    chunks[i] = new fields_chunk(s->chunks[i], outdir, m);
  FOR_FIELD_TYPES(ft) {
    for (int ip=0;ip<3;ip++) {
      comm_sizes[ft][ip] = new int[num_chunks*num_chunks];
      for (int i=0;i<num_chunks*num_chunks;i++) comm_sizes[ft][ip][i] = 0;
    }
    typedef double *double_ptr;
    comm_blocks[ft] = new double_ptr[num_chunks*num_chunks];
    for (int i=0;i<num_chunks*num_chunks;i++)
      comm_blocks[ft][i] = 0;
  }
  for (int b=0;b<2;b++) FOR_DIRECTIONS(d)
    if (v.has_boundary((boundary_side)b, d)) boundaries[b][d] = Metallic;
    else boundaries[b][d] = None;
  chunk_connections_valid = false;
}

fields::fields(const fields &thef) :
  S(thef.S), v(thef.v), user_volume(thef.user_volume), gv(thef.gv)
{
  verbosity = 0;
  outdir = thef.outdir;
  m = thef.m;
  phasein_time = thef.phasein_time;
  bands = NULL;
  for (int d=0;d<5;d++) k[d] = thef.k[d];
  is_real = thef.is_real;
  a = thef.a;
  dt = thef.dt;
  t = thef.t;
  sources = NULL;
  disable_sources = thef.disable_sources;
  fluxes = NULL;
  // Time stuff:
  was_working_on = working_on = Other;
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
      comm_sizes[ft][ip] = new int[num_chunks*num_chunks];
      for (int i=0;i<num_chunks*num_chunks;i++) comm_sizes[ft][ip][i] = 0;
    }
    typedef double *double_ptr;
    comm_blocks[ft] = new double_ptr[num_chunks*num_chunks];
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
  delete bands;
  if (!quiet) print_times();
}

void fields::verbose(int v) {
  verbosity = v;
  for (int i=0;i<num_chunks;i++) chunks[i]->verbose(v);
}

void fields::use_real_fields() {
  LOOP_OVER_DIRECTIONS(v.dim, d)
    if (boundaries[High][d] == Periodic && k[d] != 0.0)
      abort("Can't use real fields with bloch boundary conditions!\n");
  is_real = 1;
  for (int i=0;i<num_chunks;i++) chunks[i]->use_real_fields();
  chunk_connections_valid = false;
}

bool fields::have_component(component c) {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      return chunks[i]->f[c][0] != NULL;
  return false;
}

fields_chunk::~fields_chunk() {
  delete s;
  is_real = 0; // So that we can make sure to delete everything...
  DOCMP2 FOR_COMPONENTS(c) {
    delete[] f[c][cmp];
    delete[] f_backup[c][cmp];
    delete[] f_p_pml[c][cmp];
    delete[] f_m_pml[c][cmp];
    delete[] f_backup_p_pml[c][cmp];
    delete[] f_backup_m_pml[c][cmp];
  }
  FOR_FIELD_TYPES(ft)
    for (int ip=0;ip<3;ip++)
      for (int io=0;io<2;io++)
	delete[] connections[ft][ip][io];
  FOR_FIELD_TYPES(ft) delete[] connection_phases[ft];
  FOR_ELECTRIC_COMPONENTS(ec) DOCMP2 delete[] d_minus_p[ec][cmp];
  delete h_sources;
  delete e_sources;
  delete pol;
  delete olpol;
  FOR_FIELD_TYPES(ft) delete[] zeroes[ft];
}

fields_chunk::fields_chunk(const structure_chunk *the_s, const char *od,
			   double m) : v(the_s->v), gv(the_s->gv), m(m) {
  s = new structure_chunk(the_s);
  rshift = 0;
  verbosity = 0;
  outdir = od;
  new_s = NULL;
  bands = NULL;
  is_real = 0;
  a = s->a;
  Courant = s->Courant;
  dt = s->dt;
  dft_chunks = NULL;
  pol = polarization::set_up_polarizations(s, is_real);
  olpol = polarization::set_up_polarizations(s, is_real);
  h_sources = e_sources = NULL;
  FOR_COMPONENTS(c) DOCMP2 {
    f[c][cmp] = NULL;
    f_backup[c][cmp] = NULL;
    f_p_pml[c][cmp] = NULL;
    f_m_pml[c][cmp] = NULL;
    f_backup_p_pml[c][cmp] = NULL;
    f_backup_m_pml[c][cmp] = NULL;
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
  have_d_minus_p = false;
  FOR_ELECTRIC_COMPONENTS(ec) DOCMP2 d_minus_p[ec][cmp] = NULL;
  figure_out_step_plan();
}

fields_chunk::fields_chunk(const fields_chunk &thef)
  : v(thef.v), gv(thef.gv) {
  s = new structure_chunk(thef.s);
  rshift = thef.rshift;
  verbosity = thef.verbosity;
  outdir = thef.outdir;
  m = thef.m;
  new_s = NULL;
  bands = NULL;
  is_real = thef.is_real;
  a = thef.a;
  Courant = thef.Courant;
  dt = thef.dt;
  dft_chunks = NULL;
  pol = polarization::set_up_polarizations(s, is_real);
  olpol = polarization::set_up_polarizations(s, is_real);
  h_sources = e_sources = NULL;
  FOR_COMPONENTS(c) DOCMP2 {
    f[c][cmp] = NULL;
    f_backup[c][cmp] = NULL;
    f_p_pml[c][cmp] = NULL;
    f_m_pml[c][cmp] = NULL;
    f_backup_p_pml[c][cmp] = NULL;
    f_backup_m_pml[c][cmp] = NULL;
  }
  FOR_COMPONENTS(c) DOCMP {
    if (thef.f[c][cmp])
      f[c][cmp] = new double[v.ntot()];
    if (thef.f_p_pml[c][cmp]) {
      f_p_pml[c][cmp] = new double[v.ntot()];
      f_m_pml[c][cmp] = new double[v.ntot()];
    }
    if (f[c][cmp])
      for (int i=0;i<v.ntot();i++)
	f[c][cmp][i] = thef.f[c][cmp][i];
    // Now for pml extra fields_chunk...
    if (f_p_pml[c][cmp])
      for (int i=0;i<v.ntot();i++)
	f_p_pml[c][cmp][i] = thef.f_p_pml[c][cmp][i];
    if (f_m_pml[c][cmp])
      for (int i=0;i<v.ntot();i++)
	f_m_pml[c][cmp][i] = thef.f_m_pml[c][cmp][i];
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
  have_d_minus_p = thef.have_d_minus_p;
  FOR_ELECTRIC_COMPONENTS(ec) DOCMP2 
    if (thef.d_minus_p[ec][cmp]) {
      d_minus_p[ec][cmp] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++)
	d_minus_p[ec][cmp][i] = thef.d_minus_p[ec][cmp][i];
    }
  figure_out_step_plan();
}

static inline bool cross_negative(direction a, direction b) {
  return ((3+b-a)%3) == 2;
}

static inline direction cross(direction a, direction b) {
  if (a < R && b < R) return (direction)((3+2*a-b)%3);
  return (direction) (2 + (3+2*(a-2)-(b-2))%3);
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
            (is_magnetic(c1) && is_electric(c2))) {
          const direction dc2 = component_direction(c2);
          if (dc1 != dc2 && v.has_field(c2) && v.has_field(c1) &&
              has_direction(v.dim,cross(dc1,dc2))) {
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
  for (int i=0;i<3;i++) {
    num_each_direction[i] = v.yucky_num(i);
    stride_each_direction[i] = v.stride(v.yucky_direction(i));
  }
  FOR_DIRECTIONS(d) {
    num_any_direction[d] = 1;
    stride_any_direction[d] = 0;
    for (int i=0;i<3;i++)
      if (d == v.yucky_direction(i)) {
        num_any_direction[d] = v.yucky_num(i);
        stride_any_direction[d] = v.stride(v.yucky_direction(i));
      }
  }
}

static bool is_tm(component c) {
  switch (c) {
  case Hx: case Hy: case Ez: case Dz: return true;
  default: return false;
  }
  return false;
}

static bool is_like(ndim d, component c1, component c2) {
  if (d != D2) return true;
  return !(is_tm(c1) ^ is_tm(c2));
}

void fields_chunk::alloc_f(component the_c) {
  FOR_COMPONENTS(c)
    if (is_mine() && v.has_field(c) && is_like(v.dim, the_c, c))
      DOCMP {
        if (!f[c][cmp]) {
          f[c][cmp] = new double[v.ntot()];
          for (int i=0;i<v.ntot();i++) f[c][cmp][i] = 0.0;
          if (!f_p_pml[c][cmp] && !is_electric(c)) {
	    bool need_pml = false;
	    LOOP_OVER_DIRECTIONS(v.dim, d)
	      if (s->C[d][c]) { need_pml = true; break; }
	    if (need_pml) {
	      // FIXME: we don't necessarily need both f_p_pml and f_m_pml
	      f_p_pml[c][cmp] = new double[v.ntot()];
	      f_m_pml[c][cmp] = new double[v.ntot()];
	      for (int i=0;i<v.ntot();i++) {
		f_p_pml[c][cmp][i] = 0.0;
		f_m_pml[c][cmp][i] = 0.0;
	      }
	    }
          }
	}
    }
  figure_out_step_plan();
}

void fields_chunk::remove_sources() {
  delete h_sources;
  delete e_sources;
  e_sources = h_sources = NULL;
}

void fields::remove_sources() {
  delete sources;
  sources = NULL;
  for (int i=0;i<num_chunks;i++) 
    chunks[i]->remove_sources();
}

void fields::remove_fluxes() {
  delete fluxes;
}

void fields_chunk::zero_fields() {
  FOR_COMPONENTS(c) DOCMP {
    if (f[c][cmp]) for (int i=0;i<v.ntot();i++) f[c][cmp][i] = 0.0;
    if (f_p_pml[c][cmp])
      for (int i=0;i<v.ntot();i++) f_p_pml[c][cmp][i] = 0.0;
    if (f_m_pml[c][cmp])
      for (int i=0;i<v.ntot();i++) f_m_pml[c][cmp][i] = 0.0;
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
  FOR_COMPONENTS(c) if (f[c][1]) {
    delete[] f[c][1];
    f[c][1] = 0;
  }
  if (is_mine() && pol) pol->use_real_fields();
  if (is_mine() && olpol) olpol->use_real_fields();
}

int fields::phase_in_material(const structure *snew, double time) {
  if (snew->num_chunks != num_chunks)
    abort("Can only phase in similar sets of chunks: %d vs %d\n", 
	  snew->num_chunks, num_chunks);
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->phase_in_material(snew->chunks[i]);
  phasein_time = (int) (time/dt);
  return phasein_time;
}

void fields_chunk::phase_in_material(const structure_chunk *snew) {
  new_s = snew;
}

int fields::is_phasing() {
  return phasein_time > 0;
}

// This is used for phasing the *radial origin* of a cylindrical structure
void fields::set_rshift(double rshift) {
  if (v.dim != Dcyl) abort("set_rshift is only for cylindrical coords");
  if (gv.in_direction_min(R) <= 0 && gv.in_direction_max(R) >= 0)
    abort("set_rshift is invalid if volume contains r=0");
  for (int i = 0; i < num_chunks; ++i)
    chunks[i]->rshift = rshift;
}

bool fields::equal_layout(const fields &f) const {
  if (a != f.a || 
      num_chunks != f.num_chunks ||
      gv != f.gv ||
      S != f.S)
    return false;
  for (int d=0;d<5;d++)
    if (k[d] != f.k[d])
      return false;
  for (int i = 0; i < num_chunks; ++i)
    if (chunks[i]->a != f.chunks[i]->a ||
	chunks[i]->gv != f.chunks[i]->gv)
      return false;
  return true;
}

} // namespace meep
