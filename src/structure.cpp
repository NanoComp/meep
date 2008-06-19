/* Copyright (C) 2005-2007 Massachusetts Institute of Technology
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
#include <string.h>

#include "meep.hpp"
#include "meep_internals.hpp"

namespace meep {

structure::structure()
  : Courant(0.5), gv(D1) // Aaack, this is very hokey.
{
  num_chunks = 0;
  num_effort_volumes = 0;
  effort_volumes = NULL;
  effort = NULL;
  outdir = ".";
  S = identity();
  a = 1; dt = Courant/a;
}

typedef structure_chunk *structure_chunk_ptr;

structure::structure(const volume &thev, material_function &eps,
		     const boundary_region &br,
		     const symmetry &s,
		     int num, double Courant, bool use_anisotropic_averaging,
		     double tol, int maxeval) :
  Courant(Courant), gv(D1) // Aaack, this is very hokey.
{
  outdir = ".";
  choose_chunkdivision(thev, num == 0 ? count_processors() : num, br, s);
  set_materials(eps, use_anisotropic_averaging, tol, maxeval);
}

structure::structure(const volume &thev, double eps(const vec &),
		     const boundary_region &br,
		     const symmetry &s,
		     int num, double Courant, bool use_anisotropic_averaging,
		     double tol, int maxeval) :
  Courant(Courant), gv(D1) // Aaack, this is very hokey.
{
  outdir = ".";
  choose_chunkdivision(thev, num == 0 ? count_processors() : num, br, s);
  simple_material_function epsilon(eps);
  set_materials(epsilon, use_anisotropic_averaging, tol, maxeval);
}

void structure::choose_chunkdivision(const volume &thev, 
				     int desired_num_chunks, 
				     const boundary_region &br,
				     const symmetry &s) {
  user_volume = thev;
  if (thev.dim == Dcyl && thev.get_origin().r() < 0)
    abort("r < 0 origins are not supported");
  v = thev;
  gv = v.surroundings();
  S = s;
  a = v.a;
  dt = Courant/a;

  // First, reduce overall volume v by symmetries:
  if (S.multiplicity() > 1) {
    bool break_this[3];
    for (int dd=0;dd<3;dd++) {
      const direction d = (direction) dd;
      break_this[dd] = false;
      for (int n=0;n<S.multiplicity();n++)
        if (has_direction(thev.dim,d) &&
            (S.transform(d,n).d != d || S.transform(d,n).flipped)) {
          if (thev.num_direction(d) & 1 && !break_this[d])
            master_printf("Padding %s to even number of grid points.\n",
			  direction_name(d));
          break_this[dd] = true;
        }
    }
    int break_mult = 1;
    for (int d=0;d<3;d++) {
      if (break_mult == S.multiplicity()) break_this[d] = false;
      if (break_this[d]) {
	break_mult *= 2;
	master_printf("Halving computational cell along direction %s\n",
		      direction_name(direction(d)));
	v = v.halve((direction)d);
      }
    }
    // Before padding, find the corresponding geometric volume.
    gv = v.surroundings();
    // Pad the little cell in any direction that we've shrunk:
    for (int d=0;d<3;d++)
      if (break_this[d]) v = v.pad((direction)d);
  }

  // initialize effort volumes
  num_effort_volumes = 1;
  effort_volumes = new volume[num_effort_volumes];
  effort_volumes[0] = v;
  effort = new double[num_effort_volumes];
  effort[0] = 1.0;

  // Next, add effort volumes for PML boundary regions:
  br.apply(this);

  // Finally, create the chunks:
  num_chunks = 0;
  chunks = new structure_chunk_ptr[desired_num_chunks * num_effort_volumes];
  for (int i = 0; i < desired_num_chunks; i++) {
    const int proc = i * count_processors() / desired_num_chunks;
    volume vi = v.split_by_effort(desired_num_chunks, i,
                                  num_effort_volumes, effort_volumes, effort);
    for (int j = 0; j < num_effort_volumes; j++) {
      volume vc;
      if (vi.intersect_with(effort_volumes[j], &vc)) {
	chunks[num_chunks] = new structure_chunk(vc, gv, Courant, proc);
	br.apply(this, chunks[num_chunks++]);
      }
    }
  }

  check_chunks();
}

void boundary_region::apply(structure *s) const {
  if (has_direction(s->v.dim, d) && s->user_volume.has_boundary(side, d)
      && s->user_volume.num_direction(d) > 1) {
    switch (kind) {
    case NOTHING_SPECIAL: break;
    case PML: s->use_pml(d, side, thickness, strength); break;
    default: abort("unknown boundary region kind");
    }
  }
  if (next)
    next->apply(s);
}

void boundary_region::apply(const structure *s, structure_chunk *sc) const {
  if (has_direction(s->v.dim, d) && s->user_volume.has_boundary(side, d)
      && s->user_volume.num_direction(d) > 1) {    
    switch (kind) {
    case NOTHING_SPECIAL: break;
    case PML: 
      sc->use_pml(d, thickness, s->user_volume.boundary_location(side, d),
		  strength); 
      break;
    default: abort("unknown boundary region kind");
    }
  }
  if (next)
    next->apply(s, sc);
}

boundary_region pml(double thickness, direction d, boundary_side side) {
  return boundary_region(boundary_region::PML, thickness, 1.0, d, side, NULL);
}
boundary_region pml(double thickness, direction d) {
  return (pml(thickness, d, Low) + pml(thickness, d, High));
}
boundary_region pml(double thickness) {
  boundary_region r;
  for (int id = 0; id < 5; ++id)
    r = r + pml(thickness, (direction) id);
  return r;
}

// First check that the chunk volumes do not intersect and that they add
// up to the total volume
void structure::check_chunks() {
  volume vol_intersection;
  for (int i=0; i<num_chunks; i++)
    for (int j=i+1; j<num_chunks; j++)
      if (chunks[i]->v.intersect_with(chunks[j]->v, &vol_intersection))
        abort("chunks[%d] intersects with chunks[%d]\n", i, j);
  // FIXME: should use 'long long' else will fail if grid > 2e9 points
  int sum = 0;
  for (int i=0; i<num_chunks; i++) {
    int grid_points = 1;
    LOOP_OVER_DIRECTIONS(chunks[i]->v.dim, d)
      grid_points *= chunks[i]->v.num_direction(d);
    sum += grid_points;
  }
  int v_grid_points = 1;
  LOOP_OVER_DIRECTIONS(v.dim, d) v_grid_points *= v.num_direction(d);
  if (sum != v_grid_points)
    abort("v_grid_points = %d, sum(chunks) = %d\n", v_grid_points, sum);
}

void structure::add_to_effort_volumes(const volume &new_effort_volume,
                                double extra_effort) {
  volume *temp_volumes =
    new volume[(2*number_of_directions(v.dim)+1)*num_effort_volumes]; 
  double *temp_effort =
    new double[(2*number_of_directions(v.dim)+1)*num_effort_volumes];
  // Intersect previous mat_volumes with this new_effort_volume
  int counter = 0;
  for (int j=0; j<num_effort_volumes; j++) {
    volume intersection, others[6];
    int num_others;
    if (effort_volumes[j].intersect_with(new_effort_volume, &intersection,
                                         others, &num_others)) {
      if (num_others > 1) {
        printf("effort_volumes[%d]  ", j);
        effort_volumes[j].print();
        printf("new_effort_volume  ");
        new_effort_volume.print();
        // NOTE: this may not be a bug if this function is used for
        // something other than PML.
        abort("Did not expect num_others > 1 in add_to_effort_volumes\n");
      }
      temp_effort[counter] = extra_effort + effort[j];
      temp_volumes[counter] = intersection;
      counter++;
      for (int k = 0; k<num_others; k++) {
        temp_effort[counter] = effort[j];	  
        temp_volumes[counter] = others[k];
        counter++;
      }
    } else {
      temp_effort[counter] = effort[j];	  
      temp_volumes[counter] = effort_volumes[j];
      counter++;
    }
  }

  delete[] effort_volumes;
  delete[] effort;
  effort_volumes = temp_volumes;
  effort = temp_effort;
  num_effort_volumes = counter;
}

structure::structure(const structure *s) : gv(s->gv) {
  num_chunks = s->num_chunks;
  outdir = s->outdir;
  v = s->v;
  S = s->S;
  user_volume = s->user_volume;
  chunks = new structure_chunk_ptr[num_chunks];
  for (int i=0;i<num_chunks;i++)
    chunks[i] = new structure_chunk(s->chunks[i]);
  num_effort_volumes = s->num_effort_volumes;
  effort_volumes = new volume[num_effort_volumes];
  effort = new double[num_effort_volumes];
  for (int i=0;i<num_effort_volumes;i++) {
    effort_volumes[i] = s->effort_volumes[i];
    effort[i] = s->effort[i];
  }
  a = s->a;
  Courant = s->Courant;
  dt = s->dt;
}

structure::structure(const structure &s) : gv(s.gv) {
  num_chunks = s.num_chunks;
  outdir = s.outdir;
  v = s.v;
  S = s.S;
  user_volume = s.user_volume;
  chunks = new structure_chunk_ptr[num_chunks];
  for (int i=0;i<num_chunks;i++) {
    chunks[i] = new structure_chunk(s.chunks[i]);
  }
  num_effort_volumes = s.num_effort_volumes;
  effort_volumes = new volume[num_effort_volumes];
  effort = new double[num_effort_volumes];
  for (int i=0;i<num_effort_volumes;i++) {
    effort_volumes[i] = s.effort_volumes[i];
    effort[i] = s.effort[i];
  }  
  a = s.a;
  Courant = s.Courant;
  dt = s.dt;
}

structure::~structure() {
  for (int i=0;i<num_chunks;i++) {
    if (chunks[i]->refcount-- <= 1) delete chunks[i];
    chunks[i] = NULL; // Just to be sure...
  }
  delete[] chunks;
  delete[] effort_volumes;
  delete[] effort;
}

/* To save memory, the structure chunks are shared with the
   fields_chunk objects instead of making a copy.  However, to
   preserve the illusion that the structure and fields are
   independent objects, we implement copy-on-write semantics. */
void structure::changing_chunks() { // call this whenever chunks are modified
  for (int i=0; i<num_chunks; i++)
    if (chunks[i]->refcount > 1) { // this chunk is shared, so make a copy
      chunks[i]->refcount--;
      chunks[i] = new structure_chunk(chunks[i]);
    }
}

void structure::set_materials(material_function &mat, 
			      bool use_anisotropic_averaging,
			      double tol, int maxeval) {
  set_epsilon(mat, use_anisotropic_averaging, tol, maxeval);
  if (mat.has_mu()) set_mu(mat);
  FOR_D_COMPONENTS(c) if (mat.has_conductivity(c)) set_conductivity(c, mat);
  FOR_B_COMPONENTS(c) if (mat.has_conductivity(c)) set_conductivity(c, mat);
  if (mat.has_chi3()) set_chi3(mat);
  if (mat.has_chi2()) set_chi2(mat);
}

void structure::set_epsilon(material_function &eps, 
			    bool use_anisotropic_averaging,
			    double tol, int maxeval) {
  double tstart = wall_time();
  changing_chunks();
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->set_epsilon(eps, use_anisotropic_averaging, tol, maxeval);
  if (!quiet)
    master_printf("time for set_epsilon = %g s\n", wall_time() - tstart);
}

void structure::set_epsilon(double eps(const vec &),
                            bool use_anisotropic_averaging,
			    double tol, int maxeval) {
  simple_material_function epsilon(eps);
  set_epsilon(epsilon, use_anisotropic_averaging, tol, maxeval);
}

void structure::set_mu(material_function &mu) {
  double tstart = wall_time();
  changing_chunks();
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->set_mu(mu);
  if (!quiet)
    master_printf("time for set_mu = %g s\n", wall_time() - tstart);
}

void structure::set_mu(double mufunc(const vec &)) {
  simple_material_function mu(mufunc);
  set_mu(mu);
}

void structure::set_conductivity(component c, material_function &C) {
  if (!v.has_field(c)) return;
  double tstart = wall_time();
  changing_chunks();
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->set_conductivity(c, C);
  if (!quiet)
    master_printf("time for set_conductivity = %g s\n", wall_time() - tstart);
}

void structure::set_conductivity(component c, double Cfunc(const vec &)) {
  simple_material_function conductivity(Cfunc);
  set_conductivity(c, conductivity);
}

void structure::set_chi3(material_function &eps) {
  changing_chunks();
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->set_chi3(eps);
}

void structure::set_chi3(double eps(const vec &)) {
  simple_material_function epsilon(eps);
  set_chi3(epsilon);
}

void structure::set_chi2(material_function &eps) {
  changing_chunks();
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->set_chi2(eps);
}

void structure::set_chi2(double eps(const vec &)) {
  simple_material_function epsilon(eps);
  set_chi2(epsilon);
}

void structure::use_pml(direction d, boundary_side b, double dx,
			double strength) {
  if (strength == 0.0 || dx <= 0.0) return;
  volume pml_volume = v;
  pml_volume.set_num_direction(d, int(dx*user_volume.a + 1 + 0.5)); //FIXME: exact value?
  if (b == High)
    pml_volume.set_origin(d, user_volume.big_corner().in_direction(d)
  			  - pml_volume.num_direction(d) * 2);
  const int v_to_user_shift = (user_volume.little_corner().in_direction(d) 
			       - v.little_corner().in_direction(d)) / 2;
  if (b == Low && v_to_user_shift != 0)
    pml_volume.set_num_direction(d, pml_volume.num_direction(d) + v_to_user_shift);
  add_to_effort_volumes(pml_volume, 0.60); // FIXME: manual value for pml effort
}

void structure::mix_with(const structure *oth, double f) {
  if (num_chunks != oth->num_chunks)
    abort("You can't phase materials with different chunk topologies...\n");
  changing_chunks();
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->mix_with(oth->chunks[i], f);
}

structure_chunk::~structure_chunk() {
  FOR_COMPONENTS(c) {
    FOR_DIRECTIONS(d) {
      delete[] inveps[c][d];
      delete[] invmu[c][d];
      delete[] conductivity[c][d];
      delete[] condinv[c][d];
    }
    delete[] chi2[c];
    delete[] chi3[c];
  }
  delete[] eps;
  FOR_DIRECTIONS(d) { 
    delete[] sig[d];
    delete[] siginv[d];
  }
  if (pb) delete pb;
}

void structure_chunk::mix_with(const structure_chunk *n, double f) {
  for (int i=0;i<v.ntot();i++)
    eps[i] = 1.0/(1.0/eps[i] + f*(1.0/n->eps[i]-1.0/eps[i]));
  FOR_COMPONENTS(c) FOR_DIRECTIONS(d) {
    if (!inveps[c][d] && n->inveps[c][d]) {
      inveps[c][d] = new double[v.ntot()];
      if (component_direction(c) == d) // diagonal components = 1 by default
	for (int i=0;i<v.ntot();i++) inveps[c][d][i] = 1.0;
      else
	for (int i=0;i<v.ntot();i++) inveps[c][d][i] = 0.0;
    }
    if (!invmu[c][d] && n->invmu[c][d]) {
      invmu[c][d] = new double[v.ntot()];
      if (component_direction(c) == d) // diagonal components = 1 by default
	for (int i=0;i<v.ntot();i++) invmu[c][d][i] = 1.0;
      else
	for (int i=0;i<v.ntot();i++) invmu[c][d][i] = 0.0;
    }
    if (!conductivity[c][d] && n->conductivity[c][d]) {
      conductivity[c][d] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) conductivity[c][d][i] = 0.0;
    }
    if (inveps[c][d]) {
      if (n->inveps[c][d])
	for (int i=0;i<v.ntot();i++)
	  inveps[c][d][i] += f*(n->inveps[c][d][i] - inveps[c][d][i]);
      else {
	double nval = component_direction(c) == d ? 1.0 : 0.0; // default
	for (int i=0;i<v.ntot();i++)
	  inveps[c][d][i] += f*(nval - inveps[c][d][i]);
      }
    }
    if (invmu[c][d]) {
      if (n->invmu[c][d])
	for (int i=0;i<v.ntot();i++)
	  invmu[c][d][i] += f*(n->invmu[c][d][i] - invmu[c][d][i]);
      else {
	double nval = component_direction(c) == d ? 1.0 : 0.0; // default
	for (int i=0;i<v.ntot();i++)
	  invmu[c][d][i] += f*(nval - invmu[c][d][i]);
      }
    }
    if (conductivity[c][d]) {
      if (n->conductivity[c][d])
	for (int i=0;i<v.ntot();i++)
	  conductivity[c][d][i] += f*(n->conductivity[c][d][i]
				      - conductivity[c][d][i]);
      else
	for (int i=0;i<v.ntot();i++)
	  conductivity[c][d][i] += f*(0.0 - conductivity[c][d][i]);
    }
    condinv_stale = true;
  }
  // Mix in the polarizability...
  polarizability *po = pb, *pn = n->pb;
  while (po && pn) {
    FOR_COMPONENTS(c)
      if (v.has_field(c) && is_electric(c))
        for (int i=0;i<v.ntot();i++)
          po->s[c][i] += f*(pn->s[c][i] - po->s[c][i]);
    po = po->next;
    pn = pn->next;
  }
}

// reflections from non-absorbed wave in continuum
const double Rabs = 1e-15;

void structure_chunk::use_pml(direction d, double dx, double bloc,
			      double strength) {
  if (strength == 0.0 || dx <= 0.0) return;
  const double prefac = strength * (-log(Rabs))/(4*dx*0.3333);
  // Don't bother with PML if we don't even overlap with the PML region...
  if (bloc > v.boundary_location(High,d) + dx + 1.0/a - 1e-10 ||
      bloc < v.boundary_location(Low,d) - dx - 1.0/a + 1e-10) return;
  if (is_mine()) {
    if (sig[d]) {
      delete[] sig[d]; 
      delete[] siginv[d];
      sig[d] = NULL;
      siginv[d] = NULL;
    }
    LOOP_OVER_DIRECTIONS(D3,dd) {
      if (!sig[dd]) {
	int spml = (dd==d)?(2*v.num_direction(d)+2):1;
	sigsize[dd] = spml;
	sig[dd] = new double[spml];
	siginv[dd] = new double[spml];
	for (int i=0;i<spml;++i) {
	  sig[dd][i] = 0.0;
	  siginv[dd][i] = 1.0;
	}
      }
    }
    for (int i=v.little_corner().in_direction(d);
	 i<=v.big_corner().in_direction(d)+1;++i) {
      int idx = i - v.little_corner().in_direction(d);
      double here = i * 0.5/a;
      const double x =
	0.5/a*((int)(dx*(2*a)+0.5) - (int)(fabs(bloc-here)*(2*a)+0.5));
      if (x > 0) {
	sig[d][idx]=0.5*dt*prefac*(x/dx)*(x/dx);
	siginv[d][idx] = 1/(1+sig[d][idx]);	
      }
    }
  }
  condinv_stale = true;
}

void structure_chunk::update_condinv() {
  if (!condinv_stale || !is_mine()) return;
  FOR_COMPONENTS(c) {
    direction d = component_direction(c);
    if (conductivity[c][d]) {
      if (!condinv[c][d]) condinv[c][d] = new double[v.ntot()];
      const direction dsig = direction((d+1)%3); // FIXME for Dcyl!
      const bool have_pml = sigsize[dsig] > 1;
      if (!have_pml) {
	LOOP_OVER_VOL(v, c, i)
	  condinv[c][d][i] = 1 / (1 + conductivity[c][d][i] * dt * 0.5);
      }
      else { // include PML conductivity in condinv
	int k0 = v.little_corner().in_direction(dsig);
	LOOP_OVER_VOL(v, c, i) {
	  IVEC_LOOP_ILOC(v, iloc);
          int k = iloc.in_direction(dsig) - k0;
	  condinv[c][d][i] = 1 / (1 + conductivity[c][d][i] * dt * 0.5
				  + sig[dsig][k]);
	}
      }
    }
    else if (condinv[c][d]) { // condinv not needed
      delete[] condinv[c][d];
      condinv[c][d] = NULL;
    }
  }
  condinv_stale = false;
}

structure_chunk::structure_chunk(const structure_chunk *o) : gv(o->gv) {
  refcount = 1;
  if (o->pb) pb = new polarizability(o->pb);
  else pb = NULL;
  a = o->a;
  Courant = o->Courant;
  dt = o->dt;
  v = o->v;
  the_proc = o->the_proc;
  the_is_mine = my_rank() == n_proc();
  if (is_mine() && o->eps) {
    eps = new double[v.ntot()];
    if (eps == NULL) abort("Out of memory!\n");
    for (int i=0;i<v.ntot();i++) eps[i] = o->eps[i];
  } else {
    eps = NULL;
  }
  FOR_COMPONENTS(c) {
    if (is_mine() && o->chi3[c]) {
      chi3[c] = new double[v.ntot()];
      if (chi3[c] == NULL) abort("Out of memory!\n");
      for (int i=0;i<v.ntot();i++) chi3[c][i] = o->chi3[c][i];
    } else {
      chi3[c] = NULL;
    }
    if (is_mine() && o->chi2[c]) {
      chi2[c] = new double[v.ntot()];
      if (chi2[c] == NULL) abort("Out of memory!\n");
      for (int i=0;i<v.ntot();i++) chi2[c][i] = o->chi2[c][i];
    } else {
      chi2[c] = NULL;
    }
  }
  FOR_COMPONENTS(c) FOR_DIRECTIONS(d) if (is_mine()) {
    if (o->inveps[c][d]) {
      inveps[c][d] = new double[v.ntot()];
      memcpy(inveps[c][d], o->inveps[c][d], v.ntot()*sizeof(double));
    } else inveps[c][d] = NULL;
    if (o->invmu[c][d]) {
      invmu[c][d] = new double[v.ntot()];
      memcpy(invmu[c][d], o->invmu[c][d], v.ntot()*sizeof(double));
    } else invmu[c][d] = NULL;
    if (o->conductivity[c][d]) {
      conductivity[c][d] = new double[v.ntot()];
      memcpy(conductivity[c][d], o->conductivity[c][d],
	     v.ntot()*sizeof(double));
      condinv[c][d] = new double[v.ntot()];
      memcpy(condinv[c][d], o->condinv[c][d], v.ntot()*sizeof(double));
    } else conductivity[c][d] = condinv[c][d] = NULL;
  }
  condinv_stale = o->condinv_stale;
  // Allocate the PML conductivity arrays:
  FOR_DIRECTIONS(d) { 
    sig[d] = NULL; 
    siginv[d] = NULL;
    sigsize[d] = 0;
  }
  for (int i=0;i<5;++i) sigsize[i] = 0;
  // Copy over the PML conductivity arrays:
  if (is_mine())
    FOR_DIRECTIONS(d) 
      if (o->sig[d]) {
	sig[d] = new double[2*v.num_direction(d)+1];
	siginv[d] = new double[2*v.num_direction(d)+1];
	sigsize[d] = o->sigsize[d];
	for (int i=0;i<2*v.num_direction(d)+1;i++) {
	  sig[d][i] = o->sig[d][i];
	  siginv[d][i] = o->sig[d][i];
	}
      }
}

void structure_chunk::set_chi3(material_function &epsilon) {
  if (!is_mine()) return;
  
  epsilon.set_volume(v.pad().surroundings());

  FOR_ELECTRIC_COMPONENTS(c)
    if (inveps[c][component_direction(c)]) {
      if (!chi3[c]) chi3[c] = new double[v.ntot()];
      bool trivial = true;
      LOOP_OVER_VOL(v, c, i) {
	IVEC_LOOP_LOC(v, here);
        chi3[c][i] = epsilon.chi3(here);
	trivial = trivial && (chi3[c][i] == 0.0);
      }
      
      /* currently, our update_e_from_d routine requires that
	 chi2 be present if chi3 is, and vice versa */
      if (!chi2[c]) {
	if (!trivial) {
	  chi2[c] = new double[v.ntot()]; 
	  memset(chi2[c], 0, v.ntot() * sizeof(double)); // chi2 = 0
	}
	else { // no chi3, and chi2 is trivial (== 0), so delete
	  delete[] chi3[c];
	  chi3[c] = NULL;
	}
      }
    }
}

void structure_chunk::set_chi2(material_function &epsilon) {
  if (!is_mine()) return;
  
  epsilon.set_volume(v.pad().surroundings());

  FOR_ELECTRIC_COMPONENTS(c)
    if (inveps[c][component_direction(c)]) {
      if (!chi2[c]) chi2[c] = new double[v.ntot()];
      bool trivial = true;
      LOOP_OVER_VOL(v, c, i) {
	IVEC_LOOP_LOC(v, here);
        chi2[c][i] = epsilon.chi2(here);
	trivial = trivial && (chi2[c][i] == 0.0);
      }

      /* currently, our update_e_from_d routine requires that
	 chi3 be present if chi2 is, and vice versa */
      if (!chi3[c]) {
	if (!trivial) {
	  chi3[c] = new double[v.ntot()]; 
	  memset(chi3[c], 0, v.ntot() * sizeof(double)); // chi3 = 0
	}
	else { // no chi2, and chi3 is trivial (== 0), so delete 
	  delete[] chi2[c];
          chi2[c] = NULL;
	}
      }
    }
}

void structure_chunk::set_conductivity(component c, material_function &C) {
  if (!is_mine() || !v.has_field(c)) return;

  C.set_volume(v.pad().surroundings());

  if (!is_electric(c) && !is_magnetic(c) && !is_D(c) && !is_B(c))
    abort("invalid component for conductivity");

  direction c_d = component_direction(c);
  component c_C = is_electric(c) ? direction_component(Dx, c_d) :
    (is_magnetic(c) ? direction_component(Bx, c_d) : c);
  double *multby = is_electric(c) ? inveps[c][c_d] :
    (is_magnetic(c) ? invmu[c][c_d] : NULL);
  if (!conductivity[c_C][c_d]) 
    conductivity[c_C][c_d] = new double[v.ntot()]; 
  if (!conductivity[c_C][c_d]) abort("Memory allocation error.\n");
  bool trivial = true;
  double *cnd = conductivity[c_C][c_d];
  if (multby) {
    LOOP_OVER_VOL(v, c_C, i) {
      IVEC_LOOP_LOC(v, here);
      cnd[i] = C.conductivity(c, here) * multby[i];
      trivial = trivial && (cnd[i] == 0.0);
    }
  }
  else {
    LOOP_OVER_VOL(v, c_C, i) {
      IVEC_LOOP_LOC(v, here);
      cnd[i] = C.conductivity(c, here);
      trivial = trivial && (cnd[i] == 0.0);
    }
  }
  if (trivial) { // skip conductivity computations if conductivity == 0
    delete[] conductivity[c_C][c_d];
    conductivity[c_C][c_d] = NULL;
  }
  condinv_stale = true;
}

structure_chunk::structure_chunk(const volume &thev, 
				 const geometric_volume &vol_limit, 
				 double Courant, int pr)
  : Courant(Courant), gv(thev.surroundings() & vol_limit) {
  refcount = 1;
  pml_fmin = 0.2;
  pb = NULL;
  v = thev;
  a = thev.a;
  dt = Courant/a;
  the_proc = pr;
  the_is_mine = n_proc() == my_rank();

  // initialize materials arrays to NULL
  eps = NULL;
  FOR_COMPONENTS(c) chi3[c] = NULL;
  FOR_COMPONENTS(c) chi2[c] = NULL;
  FOR_COMPONENTS(c) FOR_DIRECTIONS(d) {
    inveps[c][d] = NULL;
    invmu[c][d] = NULL;
    conductivity[c][d] = NULL;
    condinv[c][d] = NULL;
  }
  condinv_stale = false;
  FOR_DIRECTIONS(d) { 
    sig[d] = NULL; 
    siginv[d] = NULL;
    sigsize[d] = 0;
  }
}

double structure::max_eps() const {
  double themax = 0.0;
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      themax = max(themax,chunks[i]->max_eps());
  return max_to_all(themax);
}

double fields::max_eps() const {
  double themax = 0.0;
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      themax = max(themax,chunks[i]->s->max_eps());
  return max_to_all(themax);
}

double structure_chunk::max_eps() const {
  double themax = 0.0;
  for (int i=0;i<v.ntot();i++) themax = max(themax,eps[i]);
  return themax;
}

bool structure::equal_layout(const structure &s) const {
  if (a != s.a || 
      num_chunks != s.num_chunks ||
      gv != s.gv ||
      S != s.S)
    return false;
  for (int i = 0; i < num_chunks; ++i)
    if (chunks[i]->a != s.chunks[i]->a ||
	chunks[i]->gv != s.chunks[i]->gv)
      return false;
  return true;
}

void structure_chunk::remove_polarizabilities() {
  delete pb;
  pb = NULL;
}

void structure::remove_polarizabilities() {
  changing_chunks();
  for (int i=0;i<num_chunks;i++) 
    chunks[i]->remove_polarizabilities();
}

// for debugging, display the chunk layout
void structure::print_layout(void) const {
  direction d0 = v.yucky_direction(0);
  direction d1 = v.yucky_direction(1);
  direction d2 = v.yucky_direction(2);
  for (int i = 0; i < num_chunks; ++i) {
    master_printf("chunk[%d] on process %d, resolution %g (%s,%s,%s):"
		  " (%d,%d,%d) - (%d,%d,%d)\n",
		  i, chunks[i]->n_proc(), chunks[i]->a,
		  direction_name(d0),direction_name(d1),direction_name(d2),
		  chunks[i]->v.little_corner().yucky_val(0),
		  chunks[i]->v.little_corner().yucky_val(1),
		  chunks[i]->v.little_corner().yucky_val(2),
		  chunks[i]->v.big_corner().yucky_val(0),
		  chunks[i]->v.big_corner().yucky_val(1),
		  chunks[i]->v.big_corner().yucky_val(2));
  }
}

} // namespace meep
