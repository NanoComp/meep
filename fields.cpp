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

#include "dactyl.h"
#include "dactyl_internals.h"

fields::fields(const mat *ma, int tm) {
  verbosity = 0;
  v = ma->v;
  user_volume = ma->user_volume;
  outdir = ma->outdir;
  m = tm;
  phasein_time = 0;
  bands = NULL;
  S = ma->S;
  for (int d=0;d<5;d++) k[d] = 0.0;
  is_real = 0;
  a = v.a;
  inva = 1.0/a;
  t = 0;
  fluxes = NULL;
  // Time stuff:
  was_working_on = working_on = Other;
  for (int i=0;i<=Other;i++) times_spent[i] = 0.0;
  last_time = 0;
  am_now_working_on(Other);

  num_chunks = ma->num_chunks;
  chunks = new (fields_chunk *)[num_chunks];
  for (int i=0;i<num_chunks;i++)
    chunks[i] = new fields_chunk(ma->chunks[i], outdir, m);
  for (int ft=0;ft<2;ft++) {
    comm_sizes[ft] = new int[num_chunks*num_chunks];
    for (int i=0;i<num_chunks*num_chunks;i++) comm_sizes[ft][i] = 0;
    comm_blocks[ft] = new (double *)[num_chunks*num_chunks];
    for (int i=0;i<num_chunks*num_chunks;i++)
      comm_blocks[ft][i] = 0;
  }
  for (int b=0;b<2;b++) for (int d=0;d<5;d++) boundaries[b][d] = None;
  connect_chunks();
}

fields::~fields() {
  for (int i=0;i<num_chunks;i++) delete chunks[i];
  delete[] chunks;
  for (int ft=0;ft<2;ft++) {
    for (int i=0;i<num_chunks*num_chunks;i++)
      delete[] comm_blocks[ft][i];
    delete[] comm_blocks[ft];
    delete[] comm_sizes[ft];
  }
  delete fluxes;
  delete bands;
}
void fields::use_real_fields() {
  for (int d=0;d<5;d++)
    if (boundaries[High][d] == Periodic)
      abort("Can't use real fields_chunk with bloch boundary conditions!\n");
  is_real = 1;
  for (int i=0;i<num_chunks;i++) chunks[i]->use_real_fields();
}

bool fields::have_component(component c) {
  if (v.dim != D2) return v.has_field(c);
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      return chunks[i]->f[c][0] != NULL;
  return false;
}

fields_chunk::~fields_chunk() {
  delete ma;
  is_real = 0; // So that we can make sure to delete everything...
  DOCMP {
    for (int i=0;i<10;i++) delete[] f[i][cmp];
    for (int i=0;i<10;i++) delete[] f_backup[i][cmp];
    for (int i=0;i<10;i++) delete[] f_pml[i][cmp];
    for (int i=0;i<10;i++) delete[] f_backup_pml[i][cmp];
    for (int ft=0;ft<2;ft++)
      for (int io=0;io<2;io++)
        delete[] connections[ft][io][cmp];
  }
  for (int ft=0;ft<2;ft++) delete[] connection_phases[ft];
  delete h_sources;
  delete e_sources;
  delete pol;
  delete olpol;
  delete fluxes;
}

fields_chunk::fields_chunk(const mat_chunk *the_ma, const char *od, int tm) {
  ma = new mat_chunk(the_ma);
  verbosity = 0;
  v = ma->v;
  outdir = od;
  m = tm;
  new_ma = NULL;
  bands = NULL;
  is_real = 0;
  a = ma->a;
  inva = 1.0/a;
  fluxes = NULL;
  pol = polarization::set_up_polarizations(ma, is_real);
  olpol = polarization::set_up_polarizations(ma, is_real);
  h_sources = e_sources = NULL;
  DOCMP {
    for (int i=0;i<10;i++) f[i][cmp] = NULL;
    for (int i=0;i<10;i++) f_backup[i][cmp] = NULL;
    for (int i=0;i<10;i++) f_pml[i][cmp] = NULL;
    for (int i=0;i<10;i++) f_backup_pml[i][cmp] = NULL;

    for (int i=0;i<10;i++) if (v.dim != D2 && v.has_field((component)i))
      f[i][cmp] = new double[v.ntot()];
    for (int i=0;i<10;i++) if (f[i][cmp]) {
      f_pml[i][cmp] = new double[v.ntot()];
      if (f_pml[i][cmp] == NULL) abort("Out of memory!\n");
    }
  }
  DOCMP {
    for (int c=0;c<10;c++)
      if (f[c][cmp])
        for (int i=0;i<v.ntot();i++)
          f[c][cmp][i] = 0.0;
    // Now for pml extra fields_chunk...
    for (int c=0;c<10;c++)
      if (f[c][cmp])
        for (int i=0;i<v.ntot();i++)
          f_pml[c][cmp][i] = 0.0;
  }
  num_connections[E_stuff][Incoming] = num_connections[E_stuff][Outgoing] = 0;
  num_connections[H_stuff][Incoming] = num_connections[H_stuff][Outgoing] = 0;
  connection_phases[E_stuff] = connection_phases[H_stuff] = 0;
  for (int f=0;f<2;f++)
    for (int io=0;io<2;io++)
      for (int cmp=0;cmp<2;cmp++)
        connections[f][io][cmp] = 0;
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
  FOR_DIRECTIONS(d) have_pml_in_direction[d] = false;
  FOR_COMPONENTS(c1)
    if (f[c1][0]) {
      const direction dc1 = component_direction(c1);
      // Figure out which field components contribute.
      component c_p=Ex, c_m=Ex;
      bool var_have_p = false, var_have_m = false;
      FOR_COMPONENTS(c2)
        if ((is_electric(c1) && is_magnetic(c2)) ||
            (is_magnetic(c1) && is_electric(c2))) {
          const direction dc2 = component_direction(c2);
          if (dc1 != dc2 && v.has_field(c2) && v.has_field(c1) &&
              has_direction(v.dim,cross(dc1,dc2))) {
            direction d_deriv = cross(dc1,dc2);
            if (ma->C[d_deriv][c1]) have_pml_in_direction[d_deriv] = true;
            if (cross_negative(dc1, dc2)) {
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
}

void fields_chunk::alloc_f(component c) {
  DOCMP {
    f[c][cmp] = new double[v.ntot()];
    f_pml[c][cmp] = new double[v.ntot()];
    for (int i=0;i<v.ntot();i++) {
      f[c][cmp][i] = 0.0;
      f_pml[c][cmp][i] = 0.0;
    }
  }
  figure_out_step_plan();
}

void fields_chunk::use_real_fields() {
  is_real = 1;
  if (is_mine() && pol) pol->use_real_fields();
  if (is_mine() && olpol) olpol->use_real_fields();
}

int fields::phase_in_material(const mat *newma, double time) {
  if (newma->num_chunks != num_chunks)
    abort("Can only phase in similar sets of chunks...\n");
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->phase_in_material(newma->chunks[i]);
  phasein_time = (int) (time*a/c);
  printf("I'm going to take %d time steps to phase in the material.\n",
         phasein_time);
  return phasein_time;
}

void fields_chunk::phase_in_material(const mat_chunk *newma) {
  new_ma = newma;
}

int fields::is_phasing() {
  return phasein_time > 0;
}
