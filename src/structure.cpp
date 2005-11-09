/* Copyright (C) 2005 Massachusetts Institute of Technology
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
		     int num, double Courant, bool use_anisotropic_averaging) :
  Courant(Courant), gv(D1) // Aaack, this is very hokey.
{
  outdir = ".";
  choose_chunkdivision(thev, num == 0 ? count_processors() : num, br, s);
  set_materials(eps, use_anisotropic_averaging);
}

structure::structure(const volume &thev, double eps(const vec &), 
		     const boundary_region &br,
		     const symmetry &s,
		     int num, double Courant, bool use_anisotropic_averaging) :
  Courant(Courant), gv(D1) // Aaack, this is very hokey.
{
  outdir = ".";
  choose_chunkdivision(thev, num == 0 ? count_processors() : num, br, s);
  simple_material_function epsilon(eps);
  set_materials(epsilon, use_anisotropic_averaging);
}

void structure::choose_chunkdivision(const volume &thev, 
				     int desired_num_chunks, 
				     const boundary_region &br,
				     const symmetry &s) {
  user_volume = thev;
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
      break_this[d] = false;
      for (int n=0;n<S.multiplicity();n++)
        if (has_direction(thev.dim,(direction)d) &&
            (S.transform(d,n).d != d || S.transform(d,n).flipped)) {
          break_this[d] = true;
          if (thev.num_direction(d) & 1)
            abort("Aaack, odd number of grid points in %s direction!\n",
		  direction_name(d));
        }
    }
    for (int d=0;d<3;d++)
      if (break_this[d]) v = v.split_specifically(2,0,(direction)d);
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
  if (has_direction(s->v.dim, d) && s->user_volume.has_boundary(side, d)) {
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
  if (has_direction(s->v.dim, d) && s->user_volume.has_boundary(side, d)) {
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
  for (int i=0;i<num_chunks;i++) chunks[i] = new structure_chunk(s->chunks[i]);
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
  for (int i=0;i<num_chunks;i++) chunks[i] = new structure_chunk(s.chunks[i]);
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
    delete chunks[i];
    chunks[i] = NULL; // Just to be sure...
  }
  delete[] chunks;
  delete[] effort_volumes;
  delete[] effort;
}

void structure::set_materials(material_function &mat, 
			      bool use_anisotropic_averaging, double minvol) {
  set_epsilon(mat, use_anisotropic_averaging, minvol);
  if (mat.has_kerr()) set_kerr(mat);
}

void structure::set_epsilon(material_function &eps, 
			    bool use_anisotropic_averaging, double minvol) {
  double tstart = wall_time();
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->set_epsilon(eps, use_anisotropic_averaging, minvol);
  if (!quiet)
    master_printf("time for set_epsilon = %g s\n", wall_time() - tstart);
}

void structure::set_epsilon(double eps(const vec &),
                            bool use_anisotropic_averaging, double minvol) {
  simple_material_function epsilon(eps);
  set_epsilon(epsilon, use_anisotropic_averaging, minvol);
}

void structure::set_kerr(material_function &eps) {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->set_kerr(eps);
}

void structure::set_kerr(double eps(const vec &)) {
  simple_material_function epsilon(eps);
  set_kerr(epsilon);
}

void structure::use_pml(direction d, boundary_side b, double dx,
			double strength) {
  if (strength == 0.0 || dx <= 0.0) return;
  volume pml_volume = v;
  pml_volume.set_num_direction(d, (int) (dx*user_volume.a + 1 + 0.5)); //FIXME: exact value?
  if ((boundary_side) b == High)
    pml_volume.origin_set_direction(d, (user_volume.num_direction(d)
                                        - pml_volume.num_direction(d))/user_volume.a);
  const int v_to_user_shift = (int)
    floor((user_volume.origin_in_direction(d) - v.origin_in_direction(d))*v.a + 0.5);
  if ((boundary_side) b == Low && v_to_user_shift != 0)
    pml_volume.set_num_direction(d, pml_volume.num_direction(d) + v_to_user_shift);
  add_to_effort_volumes(pml_volume, 0.60); // FIXME: manual value for pml effort
}

void structure::mix_with(const structure *oth, double f) {
  if (num_chunks != oth->num_chunks)
    abort("You can't phase materials with different chunk topologies...\n");
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->mix_with(oth->chunks[i], f);
}

structure_chunk::~structure_chunk() {
  FOR_ELECTRIC_COMPONENTS(c) {
    FOR_DIRECTIONS(d)
      delete[] inveps[c][d];
    delete[] kerr[c];
  }
  delete[] eps;

  FOR_COMPONENTS(c) FOR_DIRECTIONS(d) delete[] C[d][c];
  FOR_COMPONENTS(c)
    FOR_DIRECTIONS(d) FOR_DIRECTIONS(d2)
        delete[] Cdecay[d][c][d2];
  if (pb) delete pb;
}

void structure_chunk::mix_with(const structure_chunk *n, double f) {
  for (int i=0;i<v.ntot();i++)
    eps[i] = 1.0/(1.0/eps[i] + f*(1.0/n->eps[i]-1.0/eps[i]));
  FOR_ELECTRIC_COMPONENTS(c) FOR_DIRECTIONS(d)
    if (inveps[c][d])
      for (int i=0;i<v.ntot();i++)
        inveps[c][d][i] += f*(n->inveps[c][d][i] - inveps[c][d][i]);
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

const double Cmax = 0.5;

void structure_chunk::use_pml(direction d, double dx, double bloc,
			      double strength) {
  if (strength == 0.0 || dx <= 0.0) return;
  const double prefac = strength * Cmax/(dx*dx);
  // Don't bother with PML if we don't even overlap with the PML region...
  if (bloc > v.boundary_location(High,d) + dx + 1.0/a ||
      bloc < v.boundary_location(Low,d) - dx - 1.0/a) return;
  if (is_mine()) {
    FOR_COMPONENTS(c)
      if (v.has_field(c) && component_direction(c) != d) {
        if (!C[d][c]) {
          C[d][c] = new double[v.ntot()];
          for (int i=0;i<v.ntot();i++) C[d][c][i] = 0.0;
        }
	LOOP_OVER_VOL(v, (component) c, i) {
	  IVEC_LOOP_LOC(v, here);
          const double x =
            0.5/a*((int)(dx*(2*a)+0.5) -
                   (int)(2*a*fabs(bloc-here.in_direction(d))+0.5));
          if (x > 0) C[d][c][i] = prefac*x*x;
        }
      }
  }
}

void structure_chunk::update_pml_arrays() {
  if (!is_mine()) return;
  FOR_DIRECTIONS(d)
    FOR_COMPONENTS(c) 
      if (C[d][c] != NULL) {
	bool all_zeros = true;
	LOOP_OVER_VOL_OWNED(v, c, i) {
	  if (C[d][c][i] != 0.0) {
	    all_zeros = false;
	    goto done; // can't use 'break' since LOOP has nested loops
	  }
	}
      done:
	if (all_zeros) {
	  delete[] C[d][c];
	  C[d][c] = NULL;
	}
      }
  update_Cdecay();
}

void structure_chunk::update_Cdecay() {
  FOR_DIRECTIONS(d) FOR_COMPONENTS(c) FOR_DIRECTIONS(d2) {
    delete[] Cdecay[d][c][d2];
    Cdecay[d][c][d2] = NULL;
  }
  FOR_DIRECTIONS(d) FOR_COMPONENTS(c) 
    if (C[d][c] != NULL)
      FOR_DIRECTIONS(d2) 
	if ((inveps[c][d2] || d2 == component_direction(c)) && d2 != d) {
	  Cdecay[d][c][d2] = new double[v.ntot()];
	  for (int i=0;i<v.ntot();i++) {
	    if (is_magnetic(c)) Cdecay[d][c][d2][i] = 1.0/(1.0+0.5*C[d][c][i]);
	    else if (is_electric(c))
          Cdecay[d][c][d2][i] = inveps[c][d2][i]/(1.0+0.5*C[d][c][i]*inveps[c][d2][i]);
        else Cdecay[d][c][d2][i] = 1.0/(1.0+0.5*C[d][c][i]); // FIXME: Is this right?
	    /*
	      TODO: THIS IS NOT CORRECT FOR NON-DIAGONAL INVEPS...maybe
              above code is not correct either????
	    if (Cdecay[d][c][d2][i] == 0.0)
	      abort("In update_Cdecay: Cdecay == 0\n");
	    */
	  }
	}
}

structure_chunk::structure_chunk(const structure_chunk *o) : gv(o->gv) {
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
  FOR_COMPONENTS(c)
    if (is_mine() && o->kerr[c]) {
      kerr[c] = new double[v.ntot()];
      if (kerr[c] == NULL) abort("Out of memory!\n");
      for (int i=0;i<v.ntot();i++) kerr[c][i] = o->kerr[c][i];
    } else {
      kerr[c] = NULL;
    }
  FOR_COMPONENTS(c) FOR_DIRECTIONS(d)
    if (is_mine() && o->inveps[c][d]) {
      inveps[c][d] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) inveps[c][d][i] = o->inveps[c][d][i];
    } else {
      inveps[c][d] = NULL;
    }
  // Allocate the conductivity arrays:
  FOR_DIRECTIONS(d) FOR_COMPONENTS(c) C[d][c] = NULL;
  FOR_DIRECTIONS(d) FOR_COMPONENTS(c) FOR_DIRECTIONS(d2)
    Cdecay[d][c][d2] = NULL;
  // Copy over the conductivity arrays:
  if (is_mine())
    FOR_DIRECTIONS(d) FOR_COMPONENTS(c)
      if (o->C[d][c]) {
        C[d][c] = new double[v.ntot()];
        for (int i=0;i<v.ntot();i++) C[d][c][i] = o->C[d][c][i];
        FOR_DIRECTIONS(d2)
          if (o->Cdecay[d][c][d2]) {
            Cdecay[d][c][d2] = new double[v.ntot()];
            for (int i=0;i<v.ntot();i++)
              Cdecay[d][c][d2][i] = o->Cdecay[d][c][d2][i];
          }
      }
}

// The following is defined in anisotropic_averaging.cpp:
double anisoaverage(component ec, direction d, material_function &eps,
                    const geometric_volume &vol, double minvol);

void structure_chunk::set_kerr(material_function &epsilon) {
  if (!is_mine()) return;
  
  epsilon.set_volume(v.pad().surroundings());

  FOR_ELECTRIC_COMPONENTS(c)
    if (inveps[c][component_direction(c)]) {
      if (!kerr[c]) kerr[c] = new double[v.ntot()];
      LOOP_OVER_VOL(v, c, i) {
	IVEC_LOOP_LOC(v, here);
        kerr[c][i] = epsilon.kerr(here);
        LOOP_OVER_DIRECTIONS(v.dim,d)
          if (d != component_direction(c)) {
            vec dx = zero_vec(v.dim);
            dx.set_direction(d,0.5/a);
            // Following TD3D, we set kerr coefficient to zero if any of
            // the adjoining epsilon points is linear.
            kerr[c][i] = min(kerr[c][i], epsilon.kerr(here + dx));
            kerr[c][i] = min(kerr[c][i], epsilon.kerr(here - dx));
          }
      }
    }
}

void structure_chunk::set_epsilon(material_function &epsilon,
				  bool use_anisotropic_averaging,
				  double minvol) {
  if (!is_mine()) return;

  epsilon.set_volume(v.pad().surroundings());

  if (!eps) eps = new double[v.ntot()];
  LOOP_OVER_VOL(v, v.eps_component(), i) {
    IVEC_LOOP_LOC(v, here);
    eps[i] = epsilon.eps(here);
  }
  
  if (!use_anisotropic_averaging) {
    FOR_ELECTRIC_COMPONENTS(c)
      if (v.has_field(c)) {
	bool have_other_direction = false;
	vec dxa = zero_vec(v.dim);
	vec dxb = zero_vec(v.dim);
	direction c_d = component_direction(c);
	LOOP_OVER_DIRECTIONS(v.dim,da)
	  if (da != c_d) {
	    dxa.set_direction(da,0.5/a);
	    LOOP_OVER_DIRECTIONS(v.dim,db)
	      if (db != c_d && db != da) {
		dxb.set_direction(db,0.5/a);
		have_other_direction = true;
	      }
	    break;
	  }
	if (!inveps[c][c_d]) inveps[c][c_d] = new double[v.ntot()];
	LOOP_OVER_VOL(v, c, i) {
	  IVEC_LOOP_LOC(v, here);
	  if (!have_other_direction)
	    inveps[c][c_d][i] =
	      2.0/(epsilon.eps(here + dxa) + epsilon.eps(here - dxa));
	  else
	    inveps[c][c_d][i] = 4.0/(epsilon.eps(here + dxa + dxb) +
				     epsilon.eps(here + dxa - dxb) +
				     epsilon.eps(here - dxa + dxb) +
				     epsilon.eps(here - dxa - dxb));
	}
      }
  } else {
    if (minvol == 0.0) minvol = v.dV(zero_ivec(v.dim)).full_volume()/100.0;
    FOR_ELECTRIC_COMPONENTS(c)
      if (v.has_field(c))
	FOR_ELECTRIC_COMPONENTS(c2) if (v.has_field(c2)) {
	  direction d = component_direction(c2);
          if (!inveps[c][d]) inveps[c][d] = new double[v.ntot()];
          if (!inveps[c][d]) abort("Memory allocation error.\n");
	  LOOP_OVER_VOL(v, c, i) {
	    IVEC_LOOP_ILOC(v, here);
	    inveps[c][d][i] = anisoaverage(c, d, epsilon, v.dV(here), minvol);
	  }
      }
  }

  update_pml_arrays(); // PML stuff depends on epsilon
}

structure_chunk::structure_chunk(const volume &thev, 
				 const geometric_volume &vol_limit, 
				 double Courant, int pr)
  : Courant(Courant), gv(thev.surroundings() & vol_limit) {
  pml_fmin = 0.2;
  pb = NULL;
  v = thev;
  a = thev.a;
  dt = Courant/a;
  the_proc = pr;
  the_is_mine = n_proc() == my_rank();

  // initialize materials arrays to NULL
  eps = NULL;
  FOR_COMPONENTS(c) kerr[c] = NULL;
  FOR_COMPONENTS(c) FOR_DIRECTIONS(d) inveps[c][d] = NULL;
  FOR_DIRECTIONS(d) FOR_COMPONENTS(c) C[d][c] = NULL;
  FOR_DIRECTIONS(d) FOR_DIRECTIONS(d2) FOR_COMPONENTS(c) Cdecay[d][c][d2] = NULL;
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
  for (int i=0;i<num_chunks;i++) 
    chunks[i]->remove_polarizabilities();
}

} // namespace meep
