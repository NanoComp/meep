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

fields::fields(const mat *ma, int tm=0) {
  verbosity = 0;
  v = ma->v;
  outdir = ma->outdir;
  m = tm;
  phasein_time = 0;
  bands = NULL;
  k = -1;
  is_real = 0;
  a = v.a;
  inva = 1.0/a;
  t = 0;

  num_chunks = ma->num_chunks;
  chunks = new (fields_chunk *)[num_chunks];
  for (int i=0;i<num_chunks;i++)
    chunks[i] = new fields_chunk(ma->chunks[i], m);
  connect_chunks();
}
void fields::use_bloch(double kz) { // FIXME
  k = kz;
  const double cosknz = cos(k*2*pi*inva*v.nz());
  const double sinknz = sin(k*2*pi*inva*v.nz());
  eiknz = complex<double>(cosknz, sinknz);
  connect_chunks();
}

vec fields::lattice_vector() {
  if (v.dim == dcyl) {
    return vec(0,v.nz()*inva);
  } else if (v.dim == d1) {
    return vec(v.nz()*inva);
  } else {
    printf("Don't support lattice_vector in these dimensions.\n");
    exit(1);
  }
}

fields::~fields() {
  for (int i=0;i<num_chunks;i++) delete chunks[i];
  delete[] chunks;
  delete bands;
}
void fields::use_real_fields() {
  if (k >= 0.0) {
    printf("Can't use real fields_chunk with bloch boundary conditions!\n");
    exit(1);
  }
  is_real = 1;
  for (int i=0;i<num_chunks;i++) chunks[i]->use_real_fields();
}

fields_chunk::~fields_chunk() {
  delete ma;
  is_real = 0; // So that we can make sure to delete everything...
  DOCMP {
    for (int i=0;i<10;i++) delete[] f[i][cmp];
    for (int i=0;i<10;i++) delete[] f_backup[i][cmp];
    for (int i=0;i<10;i++) delete[] f_pml[i][cmp];
    for (int i=0;i<10;i++) delete[] f_backup_pml[i][cmp];
    delete[] e_connection_sinks[cmp];
    delete[] e_connection_sources[cmp];
    delete[] h_connection_sinks[cmp];
    delete[] h_connection_sources[cmp];
  }
  delete[] e_phases;
  delete[] h_phases;
  delete h_sources;
  delete e_sources;
  delete pol;
  delete olpol;
}

fields_chunk::fields_chunk(const mat_chunk *the_ma, int tm) {
  ma = new mat_chunk(the_ma);
  verbosity = 0;
  v = ma->v;
  outdir = ma->outdir;
  m = tm;
  new_ma = NULL;
  bands = NULL;
  is_real = 0;
  a = ma->a;
  inva = 1.0/a;
  pol = polarization::set_up_polarizations(ma, is_real);
  olpol = polarization::set_up_polarizations(ma, is_real);
  h_sources = e_sources = NULL;
  DOCMP {
    for (int i=0;i<10;i++) f[i][cmp] = NULL;
    for (int i=0;i<10;i++) f_backup[i][cmp] = NULL;
    for (int i=0;i<10;i++) f_pml[i][cmp] = NULL;
    for (int i=0;i<10;i++) f_backup_pml[i][cmp] = NULL;

    for (int i=0;i<10;i++) if (v.has_field((component)i))
      f[i][cmp] = new double[v.ntot()];
    for (int i=0;i<10;i++) if (v.has_field((component)i)) {
      f_pml[i][cmp] = new double[v.ntot()];
      if (f_pml[i][cmp] == NULL) {
        printf("Out of memory!\n");
        exit(1);
      }
    }
  }
  DOCMP {
    for (int c=0;c<10;c++)
      if (v.has_field((component)c))
        for (int i=0;i<v.ntot();i++)
          f[c][cmp][i] = 0.0;
    // Now for pml extra fields_chunk...
    for (int c=0;c<10;c++)
      if (v.has_field((component)c))
        for (int i=0;i<v.ntot();i++)
          f_pml[c][cmp][i] = 0.0;
  }
  num_e_connections = num_h_connections = 0;
  for (int cmp=0;cmp<2;cmp++) {
    e_connection_sources[cmp] = e_connection_sinks[cmp] = NULL;
    h_connection_sources[cmp] = h_connection_sinks[cmp] = NULL;
    h_phases = e_phases = NULL;
  }
}

void fields_chunk::use_real_fields() {
  is_real = 1;
  if (pol) pol->use_real_fields();
  if (olpol) olpol->use_real_fields();
}

int fields::phase_in_material(const mat *newma, double time) {
  if (newma->num_chunks != num_chunks) {
    printf("Can only phase in similar sets of chunks...\n");
    exit(1);
  }
  for (int i=0;i<num_chunks;i++)
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
