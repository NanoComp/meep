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

#include <stdlib.h>
#include <complex>

#include "dactyl.h"
#include "dactyl_internals.h"

void fields::disconnect_chunks() {
  for (int i=0;i<num_chunks;i++) {
    DOCMP {
      delete[] chunks[i]->e_connection_sources[cmp];
      delete[] chunks[i]->e_connection_sinks[cmp];  
      delete[] chunks[i]->h_connection_sources[cmp];
      delete[] chunks[i]->h_connection_sinks[cmp];
      chunks[i]->e_connection_sources[cmp] = 0;
      chunks[i]->e_connection_sinks[cmp] = 0;
      chunks[i]->h_connection_sources[cmp] = 0;
      chunks[i]->h_connection_sinks[cmp] = 0;
    }
    delete[] chunks[i]->e_phases;
    delete[] chunks[i]->h_phases;
    chunks[i]->e_phases = 0;
    chunks[i]->h_phases = 0;
    chunks[i]->num_h_connections = 0;
    chunks[i]->num_e_connections = 0;
  }
}

void fields::connect_chunks() {
  disconnect_chunks();
  connect_adjacent_chunks();
  connect_periodic_chunks();
  connect_metallic_chunks();
}

void fields::connect_adjacent_chunks() {
  for (int i=0;i<num_chunks;i++) {
    // First count the border elements...
    int num_e_connections = 0, num_h_connections = 0;
    const volume vi = chunks[i]->v;
    for (int c=0;c<10;c++)
      if (vi.has_field((component)c))
        for (int n=0;n<vi.ntot();n++) {
          const vec here = vi.loc((component)c, n);
          if (!vi.owns(here))
            // Here we are, looking at a border element...
            for (int j=0;j<num_chunks;j++)
              if (j != i)
                if (chunks[j]->v.owns(here)) {
                  if (is_electric((component)c)) num_e_connections++;
                  else num_h_connections++;
                }
        }
    // Now allocate the connection arrays...
    chunks[i]->alloc_extra_connections(num_e_connections, num_h_connections);
    // Next start setting up the connections...
    int which_h = 0, which_e = 0;
    for (int c=0;c<10;c++)
      if (vi.has_field((component)c))
        for (int n=0;n<vi.ntot();n++) {
          const vec here = vi.loc((component)c, n);
          if (!vi.owns(here)) // This means we're looking at a border element...
            for (int j=0;j<num_chunks;j++)
              if (j != i)
                if (chunks[j]->v.owns(here)) {
                  // the following is deprecated, but ok here.
                  int m = chunks[j]->v.index((component)c, here);
                  if (is_electric((component)c)) {
                    DOCMP {
                      chunks[i]->e_connection_sinks[cmp][which_e] =
                        chunks[i]->f[c][cmp] + n;
                      chunks[i]->e_connection_sources[cmp][which_e] =
                        chunks[j]->f[c][cmp] + m;
                    }
                    chunks[i]->e_phases[which_e] = 1.0;
                    which_e++;
                  } else {
                    DOCMP {
                      chunks[i]->h_connection_sinks[cmp][which_h] =
                        chunks[i]->f[c][cmp] + n;
                      chunks[i]->h_connection_sources[cmp][which_h] =
                        chunks[j]->f[c][cmp] + m;
                    }
                    chunks[i]->h_phases[which_h] = 1.0;
                    which_h++;
                  }
                }
        }
  }
}

void fields::connect_periodic_chunks() {
  if (k == -1) return; // It isn't periodic!  :)
  for (int i=0;i<num_chunks;i++) {
    // First count the border elements in z direction...
    int num_e_connections = 0, num_h_connections = 0;
    const volume vi = chunks[i]->v;
    for (int c=0;c<10;c++)
      if (vi.has_field((component)c))
        for (int n=0;n<vi.ntot();n++) {
          const vec here = vi.loc((component)c, n);
          if (here.z() == v.origin.z() && !vi.owns(here)) {
            // It is on the z == 0 border...
            const vec there = here + lattice_vector();
            for (int j=0;j<num_chunks;j++)
              if (chunks[j]->v.owns(there)) {
                if (is_electric((component)c)) num_e_connections++;
                else num_h_connections++;
              }
          } else if (here.z() > v.origin.z() + lattice_vector().z() &&
                     !vi.owns(here)) {
            // This is on the high z border...
            const vec there = here - lattice_vector();
            for (int j=0;j<num_chunks;j++)
              if (chunks[j]->v.owns(there)) {
                if (is_electric((component)c)) num_e_connections++;
                else num_h_connections++;
              }
          }
        }
    // Now allocate the connection arrays...
    chunks[i]->alloc_extra_connections(num_e_connections, num_h_connections);
    // Next start setting up the connections...
    int which_h = 0, which_e = 0;
    for (int c=0;c<10;c++)
      if (vi.has_field((component)c))
        for (int n=0;n<vi.ntot();n++) {
          const vec here = vi.loc((component)c, n);
          vec there;
          int isborder = false;
          complex<double> phase;
          if (here.z() == 0.0) {
            // It is on the z == 0 border...
            isborder = true;
            phase = conj(eiknz);
            there = here + lattice_vector();
          } else if (here.z() >  + lattice_vector().z() && !vi.owns(here)) {
            // This is on the high z border...
            isborder = true;
            phase = eiknz;
            there = here - lattice_vector();
          }
          if (isborder && !vi.owns(here)) {
            for (int j=0;j<num_chunks;j++)
              if (chunks[j]->v.owns(there)) {
                int m = chunks[j]->v.index((component)c, there);
                //printf("Connected %4lg %4lg with %4lg %4lg -- %s\n",
                //       here.r(), here.z(), there.r(), there.z(),
                //       component_name((component)c));
                if (is_electric((component)c)) {
                  DOCMP {
                    chunks[i]->e_connection_sinks[cmp][which_e] =
                      chunks[i]->f[c][cmp] + n;
                    chunks[i]->e_connection_sources[cmp][which_e] =
                      chunks[j]->f[c][cmp] + m;
                  }
                  chunks[i]->e_phases[which_e] = phase;
                  which_e++;
                } else {
                  DOCMP {
                    chunks[i]->h_connection_sinks[cmp][which_h] =
                      chunks[i]->f[c][cmp] + n;
                    chunks[i]->h_connection_sources[cmp][which_h] =
                      chunks[j]->f[c][cmp] + m;
                  }
                  chunks[i]->h_phases[which_h] = phase;
                  which_h++;
                }
              }
          }
        }
  }
}

void fields::connect_metallic_chunks() {
  if (v.dim == dcyl) connect_metallic_bigr_chunks();
  if (k != -1.0) return;
  connect_metallic_bigz_chunks();
}

static double zero = 0.0;

void fields::connect_metallic_bigr_chunks() {
  // Note: Only call this function if you are in cylindrical coords!
  for (int i=0;i<num_chunks;i++) {
    // First count the radial border elements...
    int num_e_connections = 0, num_h_connections = 0;
    const volume vi = chunks[i]->v;
    for (int c=0;c<10;c++)
      if (vi.has_field((component)c))
        for (int n=0;n<vi.ntot();n++) {
          const vec here = vi.loc((component)c, n);
          if (here.r() > v.origin.r() + v.nr()*inva && !vi.owns(here)) {
            // It is on the big r border...
            if (is_electric((component)c)) num_e_connections++;
            else num_h_connections++;
          }
        }
    // Now allocate the connection arrays...
    chunks[i]->alloc_extra_connections(num_e_connections, num_h_connections);
    // Next start setting up the connections...
    int which_h = 0, which_e = 0;
    for (int c=0;c<10;c++)
      if (vi.has_field((component)c))
        for (int n=0;n<vi.ntot();n++) {
          const vec here = vi.loc((component)c, n);
          if (here.r() > v.origin.r() + v.nr()*inva && !vi.owns(here)) {
            // It is on the big r border...
            //printf("Connected %4lg %4lg with zero -- %s\n",
            //       here.r(), here.z(), component_name((component)c));
            if (is_electric((component)c)) {
              DOCMP {
                chunks[i]->e_connection_sinks[cmp][which_e] =
                  chunks[i]->f[c][cmp] + n;
                chunks[i]->e_connection_sources[cmp][which_e] = &zero;
              }
              chunks[i]->e_phases[which_e] = 1.0;
              which_e++;
            } else {
              DOCMP {
                chunks[i]->h_connection_sinks[cmp][which_h] =
                  chunks[i]->f[c][cmp] + n;
                chunks[i]->h_connection_sources[cmp][which_h] = &zero;
              }
              chunks[i]->h_phases[which_h] = 1.0;
              which_h++;
            }
          }
        }
  }
}

void fields::connect_metallic_bigz_chunks() {
  // Note: Only call this function if you are in cylindrical coords!
  for (int i=0;i<num_chunks;i++) {
    // First count the radial border elements...
    int num_e_connections = 0, num_h_connections = 0;
    const volume vi = chunks[i]->v;
    for (int c=0;c<10;c++)
      if (vi.has_field((component)c))
        for (int n=0;n<vi.ntot();n++) {
          const vec here = vi.loc((component)c, n);
          if (here.z() == v.origin.z() + (v.nz()+0.5)*inva && !vi.owns(here)) {
            // It is on the big z border...
            if (is_electric((component)c)) num_e_connections++;
            else num_h_connections++;
          }
        }
    // Now allocate the connection arrays...
    chunks[i]->alloc_extra_connections(num_e_connections, num_h_connections);
    // Next start setting up the connections...
    int which_h = 0, which_e = 0;
    for (int c=0;c<10;c++)
      if (vi.has_field((component)c))
        for (int n=0;n<vi.ntot();n++) {
          const vec here = vi.loc((component)c, n);
          if (here.z() == v.origin.z() + (v.nz()+0.5)*inva && !vi.owns(here)) {
            // It is on the big z border...
            if (is_electric((component)c)) {
              DOCMP {
                chunks[i]->e_connection_sinks[cmp][which_e] =
                  chunks[i]->f[c][cmp] + n;
                chunks[i]->e_connection_sources[cmp][which_e] = &zero;
              }
              chunks[i]->e_phases[which_e] = 1.0;
              which_e++;
            } else {
              DOCMP {
                chunks[i]->h_connection_sinks[cmp][which_h] =
                  chunks[i]->f[c][cmp] + n;
                chunks[i]->h_connection_sources[cmp][which_h] = &zero;
              }
              chunks[i]->h_phases[which_h] = 1.0;
              which_h++;
            }
          }
        }
  }
}

void fields_chunk::alloc_extra_connections(int nume, int numh) {
  const int toth = num_h_connections + numh;
  const int tote = num_e_connections + nume;
  complex<double> *eph = new complex<double>[tote];
  complex<double> *hph = new complex<double>[toth];
  if (!hph) {
    printf("Out of memory!\n");
    exit(1);
  }
  for (int x=0;x<num_e_connections;x++) eph[nume+x] = e_phases[x];
  for (int x=0;x<num_h_connections;x++) hph[numh+x] = h_phases[x];
  delete[] e_phases;
  delete[] h_phases;
  e_phases = eph;
  h_phases = hph;
  DOCMP {
    double **sink_e = new double *[tote];
    double **source_e = new double *[tote];
    double **sink_h = new double *[toth];
    double **source_h = new double *[toth];
    if (!source_h) {
      printf("Out of memory!\n");
      exit(1);
    }
    for (int x=0;x<num_e_connections;x++) {
      sink_e[nume+x] = e_connection_sinks[cmp][x];
      source_e[nume+x] = e_connection_sources[cmp][x];
    }
    for (int x=0;x<num_h_connections;x++) {
      sink_h[numh+x] = h_connection_sinks[cmp][x];
      source_h[numh+x] = h_connection_sources[cmp][x];
    }
    delete[] e_connection_sources[cmp];
    delete[] h_connection_sources[cmp];
    delete[] e_connection_sinks[cmp];
    delete[] h_connection_sinks[cmp];
    e_connection_sources[cmp] = source_e;
    h_connection_sources[cmp] = source_h;
    e_connection_sinks[cmp] = sink_e;
    h_connection_sinks[cmp] = sink_h;
  }
  num_e_connections += nume;
  num_h_connections += numh;
}
