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
}
void fields::use_bloch(double kz) { // FIXME
  for (int i=0;i<num_chunks;i++) chunks[i]->use_bloch(kz);
}
fields::~fields() {
  for (int i=0;i<num_chunks;i++) delete chunks[i];
  delete[] chunks;
}
void fields::use_real_fields() {
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
  delete bands;
  delete pol;
  delete olpol;
}

void fields_chunk::use_bloch(double tk) {
  if (is_real) {
    printf("Can't do bloch boundaries with real fields_chunk!\n");
    // Note that I *could* implement bloch boundary conditions, at least
    // for gamma point and zone edge situations.
    exit(1);
  }
  k = tk;
  cosknz = cos(k*2*pi*inva*v.nz());
  sinknz = sin(k*2*pi*inva*v.nz());
  eiknz = complex<double>(cosknz, sinknz);
  complex<double> emiknz = complex<double>(cosknz, -sinknz);
  if (v.dim == d1) {
    delete[] e_phases;
    delete[] h_phases;
    num_h_connections = 0;
    num_e_connections = 1;
    e_phases = new complex<double>[num_e_connections];
    e_phases[0] = eiknz;
    DOCMP {
      delete[] e_connection_sinks[cmp];
      delete[] e_connection_sources[cmp];
      e_connection_sinks[cmp] = new (double *)[num_e_connections];
      e_connection_sources[cmp] = new (double *)[num_e_connections];
      e_connection_sources[cmp][0] = &f[Ex][cmp][0];
      e_connection_sinks[cmp][0] = &f[Ex][cmp][v.nz()];
    }
  } else if (v.dim == dcyl) {
    num_e_connections = 3*(v.nr()+1);
    num_h_connections = 3*(v.nr()+1);
    delete[] e_phases;
    delete[] h_phases;
    e_phases = new complex<double>[num_e_connections];
    h_phases = new complex<double>[num_h_connections];
    DOCMP {
      delete[] e_connection_sinks[cmp];
      delete[] e_connection_sources[cmp];
      delete[] h_connection_sinks[cmp];
      delete[] h_connection_sources[cmp];
      e_connection_sinks[cmp] = new (double *)[num_e_connections];
      h_connection_sinks[cmp] = new (double *)[num_h_connections];
      e_connection_sources[cmp] = new (double *)[num_e_connections];
      h_connection_sources[cmp] = new (double *)[num_h_connections];
      int econ = 0;
      int hcon = 0;
      for (int r=0;r<=v.nr();r++) {
        int right = v.nz() + r*(v.nz()+1);
        int left = 0 + r*(v.nz()+1);
        e_connection_sources[cmp][econ] = &f[Ep][cmp][right];
        e_connection_sinks[cmp][econ] = &f[Ep][cmp][left];
        e_phases[econ] = emiknz;
        econ++;
        e_connection_sources[cmp][econ] = &f[Er][cmp][right];
        e_connection_sinks[cmp][econ] = &f[Er][cmp][left];
        e_phases[econ] = emiknz;
        econ++;
        e_connection_sources[cmp][econ] = &f[Ez][cmp][left];
        e_connection_sinks[cmp][econ] = &f[Ez][cmp][right];
        e_phases[econ] = eiknz;
        econ++;

        h_connection_sources[cmp][hcon] = &f[Hz][cmp][right];
        h_connection_sinks[cmp][hcon] = &f[Hz][cmp][left];
        h_phases[hcon] = emiknz;
        hcon++;
        h_connection_sources[cmp][hcon] = &f[Hr][cmp][left];
        h_connection_sinks[cmp][hcon] = &f[Hr][cmp][right];
        h_phases[hcon] = eiknz;
        hcon++;
        h_connection_sources[cmp][hcon] = &f[Hp][cmp][left];
        h_connection_sinks[cmp][hcon] = &f[Hp][cmp][right];
        h_phases[hcon] = eiknz;
        hcon++;
        // FIXME: Need to add PML connections here!
      }
    }
  } else {
    printf("Unsupported dimension?!\n");
    exit(1);
  }
}

double zero = 0.0;

fields_chunk::fields_chunk(const mat_chunk *the_ma, int tm) {
  ma = new mat_chunk(the_ma);
  verbosity = 0;
  v = ma->v;
  outdir = ma->outdir;
  m = tm;
  phasein_time = 0;
  new_ma = NULL;
  bands = NULL;
  k = -1;
  is_real = 0;
  a = ma->a;
  inva = 1.0/a;
  preferred_fmax = 2.5; // Some sort of reasonable maximum
                        // frequency... (assuming a has a value on the
                        // order of your frequency).
  t = 0;
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
  if (v.dim == d1) {
    // Set up by default with metallic boundary conditions.
    num_h_connections = 0;
    num_e_connections = 1;
    e_phases = new complex<double>[num_e_connections];
    e_phases[0] = 0.0;
    h_phases = NULL;
    DOCMP {
      e_connection_sinks[cmp] = new (double *)[num_e_connections];
      e_connection_sources[cmp] = new (double *)[num_e_connections];
      e_connection_sources[cmp][0] = &zero;
      e_connection_sinks[cmp][0] = &f[Ex][cmp][v.nz()];
      h_connection_sinks[cmp] = NULL;
      h_connection_sources[cmp] = NULL;
    }
  } else if (v.dim == dcyl) {
    // Set up by default with metallic boundary conditions.
    num_e_connections = 2*(v.nr()+1);
    num_h_connections =   (v.nr()+1);
    if (f_pml[Ep][0]) num_e_connections += v.nr()+1;
    if (f_pml[Er][0]) num_e_connections += v.nr()+1;
    if (f_pml[Hz][0]) num_h_connections += v.nr()+1;
    e_phases = new complex<double>[num_e_connections];
    h_phases = new complex<double>[num_h_connections];
    DOCMP {
      e_connection_sinks[cmp] = new (double *)[num_e_connections];
      h_connection_sinks[cmp] = new (double *)[num_h_connections];
      e_connection_sources[cmp] = new (double *)[num_e_connections];
      h_connection_sources[cmp] = new (double *)[num_h_connections];
      int econ = 0;
      int hcon = 0;
      for (int i=0;i<num_e_connections;i++) e_phases[i] = 1.0;
      for (int i=0;i<num_h_connections;i++) h_phases[i] = 1.0;
      for (int r=0;r<=v.nr();r++) {
        int i = v.nz() + r*(v.nz()+1);
        e_connection_sources[cmp][econ] = &zero;
        e_connection_sinks[cmp][econ] = &f[Ep][cmp][i];
        econ++;
        e_connection_sources[cmp][econ] = &zero;
        e_connection_sinks[cmp][econ] = &f[Er][cmp][i];
        econ++;
        h_connection_sources[cmp][hcon] = &zero;
        h_connection_sinks[cmp][hcon] = &f[Hz][cmp][i];
        hcon++;
        // I have to add PML connections here too:
        if (f_pml[Ep][cmp]) {
          e_connection_sources[cmp][econ] = &zero;
          e_connection_sinks[cmp][econ] = &f_pml[Ep][cmp][i];
          econ++;
        }
        if (f_pml[Er][cmp]) {
          e_connection_sources[cmp][econ] = &zero;
          e_connection_sinks[cmp][econ] = &f_pml[Er][cmp][i];
          econ++;
        }
        if (f_pml[Hz][0]) {
          h_connection_sources[cmp][hcon] = &zero;
          h_connection_sinks[cmp][hcon] = &f_pml[Hz][cmp][i];
          hcon++;
        }
      }
    }
  }
}

void fields_chunk::use_real_fields() {
  if (k >= 0.0) {
    printf("Can't use real fields_chunk with bloch boundary conditions!\n");
    exit(1);
  }
  is_real = 1;
  if (pol) pol->use_real_fields();
  if (olpol) olpol->use_real_fields();
}

int fields_chunk::phase_in_material(const mat_chunk *newma, double time) {
  new_ma = newma;
  phasein_time = (int) (time*a/c);
  if (phasein_time == 0) phasein_time = 1;
  printf("I'm going to take %d time steps to phase in the material.\n", phasein_time);
  return phasein_time;
}

int fields_chunk::is_phasing() {
  return phasein_time > 0;
}
