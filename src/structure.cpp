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

#include "meep.h"
#include "meep_internals.h"

namespace meep {

static double one_function(const vec &r) { return 1.0; }

structure::structure()
  : gv(D1) // Aaack, this is very hokey.
{
  num_chunks = 0;
  desired_num_chunks = 0;
  num_effort_volumes = 0;
  effort_volumes = NULL;
  effort = NULL;
  outdir = ".";
  S = identity();
}

typedef structure_chunk *structure_chunk_ptr;

structure::structure(const volume &thev, material_function &eps, int num, const symmetry &s)
  : gv(D1) // Aaack, this is very hokey.
{
  outdir = ".";
  if (num == 0) num = count_processors();
  desired_num_chunks = num;
  choose_chunkdivision(thev, eps, num, s);
  num_effort_volumes = 1;
  effort_volumes = new volume[num_effort_volumes];
  effort_volumes[0] = v;
  effort = new double[num_effort_volumes];
  effort[0] = 1.0;
}

structure::structure(const volume &thev, double eps(const vec &), 
	 int num, const symmetry &s)
  : gv(D1) // Aaack, this is very hokey.
{
  outdir = ".";
  if (num == 0) num = count_processors();
  desired_num_chunks = num;
  simple_material_function epsilon(eps);
  choose_chunkdivision(thev, epsilon, num, s);
  num_effort_volumes = 1;
  effort_volumes = new volume[num_effort_volumes];
  effort_volumes[0] = v;
  effort = new double[num_effort_volumes];
  effort[0] = 1.0;
}

void structure::optimize_volumes(int *Nv, volume *new_volumes, int *procs) {
  const int num_processors = count_processors();
  int counter = 0;
  for (int i=0;i<desired_num_chunks;i++) {
    const int proc = num_processors * i / desired_num_chunks;
    volume vi = v.split_by_effort(desired_num_chunks, i,
                                  num_effort_volumes, effort_volumes, effort);
    volume v_intersection;
    for (int j=0;j<num_effort_volumes;j++) {
      if (vi.intersect_with(effort_volumes[j], &v_intersection)) {
    	new_volumes[counter] = v_intersection;
        procs[counter] = proc;
        counter++;
      }
    }
  }
  *Nv = counter;
}

void structure::optimize_chunks() {
  int Nv;
  volume *new_volumes = new volume[333]; // fix me
  int *procs = new int[333]; // fix me
  optimize_volumes(&Nv, new_volumes, procs);
  redefine_chunks(Nv, new_volumes, procs);
  
  delete[] new_volumes;
  delete[] procs;
}

inline double zero_function(const vec &v) { return 0.0; }

static inline void copy_from(int from, int to,
                             double *f, double *t, int size=1) {
  double temp;
  if (my_rank() == from) temp = *f;
  send(from, to, &temp);
  if (my_rank() == to) *t = temp;
}

void structure::redefine_chunks(const int Nv, const volume *new_volumes,
                          const int *procs) {
  // First check that the new_volumes do not intersect and that they add
  // up to the total volume

  volume vol_intersection;
  for (int i=0; i<Nv; i++)
    for (int j=i+1; j<Nv; j++)
      if (new_volumes[i].intersect_with(new_volumes[j], &vol_intersection))
        abort("new_volumes[%d] intersects with new_volumes[%d]\n", i, j);
  int sum = 0;
  for (int i=0; i<Nv; i++) {
    int grid_points = 1;
    LOOP_OVER_DIRECTIONS(new_volumes[i].dim, d)
      grid_points *= new_volumes[i].num_direction(d);
    sum += grid_points;
  }
  int v_grid_points = 1;
  LOOP_OVER_DIRECTIONS(v.dim, d) v_grid_points *= v.num_direction(d);
  if (sum != v_grid_points)
    abort("v_grid_points = %d, sum(new_volumes) = %d\n", v_grid_points, sum);

  structure_chunk **new_chunks = new structure_chunk_ptr[Nv];
  for (int j=0; j<Nv; j++)
    new_chunks[j] = NULL;
  for (int j=0; j<Nv; j++) {
    for (int i=0; i<num_chunks; i++)
      if (chunks[i]->v.intersect_with(new_volumes[j], &vol_intersection)) {
        if (new_chunks[j] == NULL) {
	  simple_material_function eps_zero(zero_function);
          new_chunks[j] = new structure_chunk(new_volumes[j], eps_zero,
                                        gv, procs[j]);
          // the above happens even if chunk not owned by proc
        }
        // chunk "printing" is parallelized below

        // FIXME: The following code is *terribly* parallelized! The
        // boolean broadcasts should be removed in favor of sends from one
        // processor to another (and should only be done if their chunks
        // intersect).  The copy_froms should all be elimiated (they're a
        // crude kludge) in favor of sending the entire chunk's data all at
        // once, followed by a more leisurely copying and deletion of the
        // buffer.

        // eps
        for (int l=0; l<vol_intersection.ntot(); l++) {
          component c = vol_intersection.eps_component();
          ivec iv = vol_intersection.iloc(c, l);
          int index_old = chunks[i]->v.index(c, iv);
          int index_new = new_chunks[j]->v.index(c, iv);
          copy_from(chunks[i]->n_proc(), new_chunks[j]->n_proc(),
                    &chunks[i]->eps[index_old],
                    &new_chunks[j]->eps[index_new]);
        }
        // inveps
        FOR_COMPONENTS(c)
          FOR_DIRECTIONS(d)
            if (broadcast(chunks[i]->n_proc(),
                          chunks[i]->inveps[c][d] != NULL)) {
              if (new_chunks[j]->is_mine() && ! new_chunks[j]->inveps[c][d]) {
                new_chunks[j]->inveps[c][d] = new double[new_chunks[j]->v.ntot()];
                for (int i=0;i<new_chunks[j]->v.ntot();i++)
                  new_chunks[j]->inveps[c][d][i] = 0.0;
              }
              for (int l=0; l<vol_intersection.ntot(); l++) {
                ivec iv = vol_intersection.iloc(c, l);
                int index_old = chunks[i]->v.index(c, iv);
                int index_new = new_chunks[j]->v.index(c, iv);
                copy_from(chunks[i]->n_proc(), new_chunks[j]->n_proc(),
                          &chunks[i]->inveps[c][d][index_old],
                          &new_chunks[j]->inveps[c][d][index_new]);
              }
            }
        FOR_DIRECTIONS(d)
          FOR_COMPONENTS(c) {
            // C
            if (broadcast(chunks[i]->n_proc(),
                          chunks[i]->C[d][c] != NULL)) {
              if (new_chunks[j]->is_mine() &&
                  new_chunks[j]->C[d][c] == NULL) {
                new_chunks[j]->C[d][c] =
                  new double[new_chunks[j]->v.ntot()];
                for (int l=0; l<new_chunks[j]->v.ntot(); l++)
                  new_chunks[j]->C[d][c][l] = 0.0;
              }
              for (int l=0; l<vol_intersection.ntot(); l++) {
                ivec iv = vol_intersection.iloc(c, l);
                int index_old = chunks[i]->v.index(c, iv);
                int index_new = new_chunks[j]->v.index(c, iv);
                copy_from(chunks[i]->n_proc(), new_chunks[j]->n_proc(),
                          &chunks[i]->C[d][c][index_old],
                          &new_chunks[j]->C[d][c][index_new]);
              }
            }
            // polarization! FIXME
        }
      }
    if (new_chunks[j]->is_mine())
      new_chunks[j]->update_pml_arrays();
  }
  // Finally replace old chunks with new chunks
  delete[] chunks;
  chunks = new_chunks;
  num_chunks = Nv;
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

void structure::choose_chunkdivision(const volume &thev, material_function &eps,
                               int num, const symmetry &s) {
  num_chunks = num;
  user_volume = thev;
  v = thev;
  gv = v.surroundings();
  S = s;
  if (S.multiplicity() > 1) {
    // Have to work out the symmetry point and volume to use.
    bool break_this[3];
    for (int dd=0;dd<3;dd++) {
      const direction d = (direction) dd;
      break_this[d] = false;
      for (int n=0;n<S.multiplicity();n++)
        if (has_direction(thev.dim,(direction)d) &&
            (S.transform(d,n).d != d || S.transform(d,n).flipped)) {
          break_this[d] = true;
          if (thev.num_direction(d) & 1)
            abort("Aaack, odd number of grid points!\n");
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
  chunks = new structure_chunk_ptr[num_chunks];
  for (int i=0;i<num_chunks;i++) {
    const int proc = i*count_processors()/num_chunks;
    chunks[i] = new structure_chunk( v.split(num_chunks,i), eps, gv, proc);
  }
}

structure::structure(const structure *s) : gv(s->gv) {
  num_chunks = s->num_chunks;
  desired_num_chunks = s->desired_num_chunks;
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
}

structure::structure(const structure &s) : gv(s.gv) {
  num_chunks = s.num_chunks;
  desired_num_chunks = s.desired_num_chunks;
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

void structure::make_average_eps() {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->make_average_eps(); // FIXME
}

void structure::set_epsilon(material_function &eps, double minvol,
                      bool use_anisotropic_averaging) {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->set_epsilon(eps, minvol, use_anisotropic_averaging);
}

void structure::set_epsilon(double eps(const vec &), double minvol,
                            bool use_anisotropic_averaging) {
  simple_material_function epsilon(eps);
  set_epsilon(epsilon, minvol, use_anisotropic_averaging);
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

void structure::use_pml(direction d, boundary_side b, double dx, bool recalculate_chunks) {
  volume pml_volume = v;
  pml_volume.set_num_direction(d, (int) (dx*user_volume.a + 1 + 0.5)); //FIXME: exact value?
  if ((boundary_side) b == High)
    pml_volume.origin.set_direction(d, (user_volume.num_direction(d)
                                        - pml_volume.num_direction(d))/user_volume.a);
  const int v_to_user_shift = (int)
    floor((user_volume.origin.in_direction(d) - v.origin.in_direction(d))*v.a + 0.5);
  if ((boundary_side) b == Low && v_to_user_shift != 0)
    pml_volume.set_num_direction(d, pml_volume.num_direction(d) + v_to_user_shift);
  add_to_effort_volumes(pml_volume, 0.60); // FIXME: manual value for pml effort

  if (recalculate_chunks) optimize_chunks();  

  for (int i=0;i<num_chunks;i++)
    chunks[i]->use_pml(d, dx, user_volume.boundary_location(b,d));
}

void structure::use_pml_everywhere(double dx, bool recalculate_chunks) {
  for (int b=0;b<2;b++) FOR_DIRECTIONS(d)
    if (user_volume.has_boundary((boundary_side)b, d))
      use_pml(d, (boundary_side)b, dx, false);
  if (recalculate_chunks) optimize_chunks();
}

void structure::mix_with(const structure *oth, double f) {
  if (num_chunks != oth->num_chunks)
    abort("You can't phase materials with different chunk topologies...\n");
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->mix_with(oth->chunks[i], f);
}

structure_chunk::~structure_chunk() {
  FOR_ELECTRIC_COMPONENTS(c)
    FOR_DIRECTIONS(d)
      delete[] inveps[c][d];
  delete[] eps;

  FOR_COMPONENTS(c) FOR_DIRECTIONS(d) delete[] C[d][c];
  FOR_COMPONENTS(c)
    FOR_DIRECTIONS(d) FOR_DIRECTIONS(d2)
        delete[] Cdecay[d][c][d2];
  if (pb) delete pb;
}

// static double sig(double r, double power);

// static double minimize_badness(double sig[], int thickness, double eps, double fmin, int i);
// inline void reverse(double sig[], int l) {
//   for (int i=0;i<l/2;i++) {
//     double temp = sig[i];
//     sig[i] = sig[l-1-i];
//     sig[l-1-i] = temp;
//   }
// }

// static double badness(double sig[], int thickness, double epsilon, double fmin) {
//   if (thickness < 1) return 1;
//   const double A = .0001/fmin*.1/fmin, K = 6.0/epsilon*2.25/epsilon;
//   double sofar = 1.0;
//   for (int i=0;i<thickness-1;i++) {
//     double first_trans = exp(-K*sig[i+1]);
//     double refl = A*fabs(sig[i]-sig[i+1])*fabs(sig[i]-sig[i+1]);
//     double total_trans = exp(-K*sig[i])*first_trans;
//     sofar = refl + (1-refl)*total_trans*sofar;
//     if (sofar > 1.0) sofar = 1.0;
//   }
//   double last_refl = A*fabs(sig[thickness-1]);
//   sofar = last_refl + (1-last_refl)*sofar;
//   return sofar;
// }

// static double minimize_badness(double sig[], int thickness,
//                                double epsilon, double fmin, int i) {
//   double behind_reflection = badness(sig, i-1, epsilon, fmin);
  

//   double now = badness(sig, thickness, epsilon, fmin);
//   double tried = now;
//   do {
//     now = tried;
//     sig[i] *= 1.001;
//     tried = badness(sig, thickness, epsilon, fmin);
//   } while (tried < now);
//   sig[i] /= 1.001;
//   tried = now = badness(sig, thickness, epsilon, fmin);
//   do {
//     now = tried;
//     sig[i] /= 1.001;
//     tried = badness(sig, thickness, epsilon, fmin);
//   } while (tried < now);
//   sig[i] *= 1.001;
//   return badness(sig, thickness, epsilon, fmin);
// }

// static double sig(double r, double power) {
//   return pow(r, power);
// }

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

void structure_chunk::make_average_eps() {
  double meaneps = 0;
  for (int i=0;i<v.ntot();i++) {
    meaneps += eps[i]; // This is totally wrong, as it needs parallelization.
  }
  meaneps /= v.ntot();
  for (int i=0;i<v.ntot();i++)
    eps[i] = meaneps;
  FOR_ELECTRIC_COMPONENTS(c)
    if (v.has_field(c))
      for (int i=0;i<v.ntot();i++)
        inveps[c][component_direction(c)][i] = 1/meaneps;
}

const double Cmax = 0.5;

void structure_chunk::use_pml(direction d, double dx, double bloc) {
  const double prefac = Cmax/(dx*dx);
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
        for (int i=0;i<v.ntot();i++) {
          const double x =
            0.5/a*((int)(dx*(2*a)+0.5) -
                   (int)(2*a*fabs(bloc-v.loc((component)c,i).in_direction(d))+0.5));
          if (x > 0) C[d][c][i] = prefac*x*x;
        }
      }
    update_pml_arrays();
  }
}

void structure_chunk::update_pml_arrays() {
  FOR_DIRECTIONS(d)
    FOR_COMPONENTS(c) 
      if (C[d][c] != NULL) {
	bool all_zeros = true;
	for (int i=0;i<v.ntot();i++)
	  if (v.owns(v.iloc(c, i)))
	    if (C[d][c][i] != 0.0)
	      all_zeros = false;
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
	    if (Cdecay[d][c][d2][i] == 0.0)
	      abort("In update_Cdecay: Cdecay == 0\n");
	  }
	}
}

structure_chunk::structure_chunk(const structure_chunk *o) : gv(o->gv) {
  if (o->pb) pb = new polarizability(o->pb);
  else pb = NULL;
  a = o->a;
  v = o->v;
  the_proc = o->the_proc;
  the_is_mine = my_rank() == n_proc();
  if (is_mine()) {
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
      for (int i=0;i<v.ntot();i++) {
        const vec here = v.loc(c,i);
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

void structure_chunk::set_epsilon(material_function &epsilon, double minvol,
                            bool use_anisotropic_averaging) {
  if (!is_mine()) return;

  epsilon.set_volume(v.pad().surroundings());

  if (!eps)
       eps = new double[v.ntot()];
  for (int i=0;i<v.ntot();i++)
       eps[i] = epsilon.eps(v.loc(v.eps_component(),i));
  
  if (!use_anisotropic_averaging) {
      for (int i=0;i<v.ntot();i++) {
        FOR_ELECTRIC_COMPONENTS(c)
          if (v.has_field(c)) {
            const vec here = v.loc(c,i);
            LOOP_OVER_DIRECTIONS(v.dim,da)
              if (da != component_direction(c)) {
                vec dxa = zero_vec(v.dim);
                dxa.set_direction(da,0.5/a);
                bool have_other_direction = false;
                LOOP_OVER_DIRECTIONS(v.dim,db)
                  if (db != component_direction(c) && db != da) {
                    vec dxb = zero_vec(v.dim);
                    dxb.set_direction(db,0.5/a);
                    inveps[c][component_direction(c)][i] =
                      4.0/(epsilon.eps(here + dxa + dxb) +
                           epsilon.eps(here + dxa - dxb) +
                           epsilon.eps(here - dxa + dxb) +
                           epsilon.eps(here - dxa - dxb));
                    have_other_direction = true;
                  }
                if (!have_other_direction)
                  inveps[c][component_direction(c)][i] =
                    2.0/(epsilon.eps(here + dxa) + epsilon.eps(here - dxa));
                break;
              }
          }
      }
  } else {
    if (minvol == 0.0) minvol = v.dV(zero_ivec(v.dim)).full_volume()/100.0;
    FOR_ELECTRIC_COMPONENTS(c)
      if (v.has_field(c))
        LOOP_OVER_DIRECTIONS(v.dim,d) {
          if (!inveps[c][d]) inveps[c][d] = new double[v.ntot()];
          if (!inveps[c][d]) abort("Memory allocation error.\n");
          for (int i=0;i<v.ntot();i++)
            inveps[c][d][i] = anisoaverage(c, d, epsilon, v.dV(v.iloc(c,i)), minvol);
        }
  }
}

structure_chunk::structure_chunk(const volume &thev, material_function &epsilon,
                     const geometric_volume &vol_limit, int pr)
  : gv(thev.surroundings() & vol_limit) {
  pml_fmin = 0.2;
  pb = NULL;
  v = thev;
  a = thev.a;
  the_proc = pr;
  the_is_mine = n_proc() == my_rank();
  eps = NULL;
  FOR_COMPONENTS(c) kerr[c] = NULL;
  FOR_COMPONENTS(c) FOR_DIRECTIONS(d)
    if (is_mine() && v.has_field(c) && is_electric(c) &&
        d == component_direction(c)) {
      inveps[c][d] = new double[v.ntot()];
      // Initialize eps to 1;
      for (int i=0;i<v.ntot();i++) inveps[c][d][i] = 1;
    } else {
      inveps[c][d] = NULL;
    }
  set_epsilon(epsilon, 0.0, false);
  if (epsilon.has_kerr()) set_kerr(epsilon);
  // Allocate the conductivity arrays:
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

} // namespace meep
