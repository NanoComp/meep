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
      for (int f=0;f<2;f++)
        for (int io=0;io<2;io++) {
          delete[] chunks[i]->connections[f][io][cmp];
          chunks[i]->connections[f][io][cmp] = 0;
        }
    }
    delete[] chunks[i]->connection_phases[E_stuff];
    delete[] chunks[i]->connection_phases[H_stuff];
    chunks[i]->connection_phases[E_stuff] = 0;
    chunks[i]->connection_phases[H_stuff] = 0;
    for (int f=0;f<2;f++)
      for (int io=0;io<2;io++)
        chunks[i]->num_connections[f][io] = 0;
  }
  for (int ft=0;ft<2;ft++)
    for (int i=0;i<num_chunks*num_chunks;i++) {
      DOCMP {
        delete[] comm_blocks[ft][cmp][i];
        comm_blocks[ft][cmp][i] = 0;
      }
      comm_sizes[ft][i] = 0;
    }
}

void fields::connect_chunks() {
  disconnect_chunks();
  connect_adjacent_chunks();
  connect_periodic_chunks();
  connect_metallic_chunks();
}

void fields::connect_adjacent_chunks() {
  int *nc[2][2];
  for (int f=0;f<2;f++)
    for (int io=0;io<2;io++) {
      nc[f][io] = new int[num_chunks];
      for (int i=0;i<num_chunks;i++) nc[f][io][i] = 0;
    }
  int *new_comm_sizes[2];
  for (int ft=0;ft<2;ft++) new_comm_sizes[ft] = new int[num_chunks*num_chunks];
  for (int i=0;i<num_chunks;i++) {
    // First count the border elements...
    const volume vi = chunks[i]->v;
    for (int ft=0;ft<2;ft++) 
      for (int j=0;j<num_chunks;j++)
        new_comm_sizes[ft][j] = comm_sizes[ft][j+i*num_chunks];
    for (int c=0;c<10;c++)
      if (vi.has_field((component)c))
        for (int n=0;n<vi.ntot();n++) {
          const vec here = vi.loc((component)c, n);
          if (!vi.owns(here)) // This means we're looking at a border element...
            for (int j=0;j<num_chunks;j++)
              if (j != i)
                if (chunks[j]->v.owns(here)) {
                  nc[type((component)c)][Incoming][i]++;
                  nc[type((component)c)][Outgoing][j]++;
                  new_comm_sizes[type((component)c)][j]++;
                }
        }
    // Allocating comm blocks as we go...
    for (int ft=0;ft<2;ft++)
      for (int j=0;j<num_chunks;j++) {
        DOCMP {
          delete[] comm_blocks[ft][cmp][j+i*num_chunks];
          comm_blocks[ft][cmp][j+i*num_chunks] = new double[new_comm_sizes[ft][j]];
        }
        comm_sizes[ft][j+i*num_chunks] = new_comm_sizes[ft][j];
      }
  }
  for (int ft=0;ft<2;ft++) delete[] new_comm_sizes[ft];
  // Now allocate the connection arrays...
  for (int i=0;i<num_chunks;i++) {
    for (int f=0;f<2;f++)
      for (int io=0;io<2;io++)
        chunks[i]->alloc_extra_connections((field_type)f,(in_or_out)io,nc[f][io][i]);
  }
  for (int f=0;f<2;f++)
    for (int io=0;io<2;io++) delete[] nc[f][io];
  // Next start setting up the connections...
  int *which[2][2];
  for (int f=0;f<2;f++)
    for (int io=0;io<2;io++) {
      which[f][io] = new int[num_chunks];
      for (int i=0;i<num_chunks;i++) which[f][io][i] = 0;
    }
  for (int i=0;i<num_chunks;i++) {
    const volume vi = chunks[i]->v;
    for (int c=0;c<10;c++)
      if (vi.has_field((component)c))
        for (int n=0;n<vi.ntot();n++) {
          const vec here = vi.loc((component)c, n);
          if (!vi.owns(here)) // This means we're looking at a border element...
            for (int j=0;j<num_chunks;j++)
              if (j != i)
                if (chunks[j]->v.owns(here)) {
                  // index is deprecated, but ok in this case:
                  int m = chunks[j]->v.index((component)c, here);
                  int ft = type((component) c);
                  DOCMP {
                    chunks[i]->connections[ft][Incoming][cmp][which[ft][Incoming][i]]
                      = chunks[i]->f[c][cmp] + n;
                    chunks[j]->connections[ft][Outgoing][cmp][which[ft][Outgoing][j]]
                      = chunks[j]->f[c][cmp] + m;
                  }
                  chunks[i]->connection_phases[ft][which[ft][Incoming][i]] = 1.0;
                  which[ft][Incoming][i]++;
                  which[ft][Outgoing][j]++;
                }
        }
  }
  for (int f=0;f<2;f++)
    for (int io=0;io<2;io++) delete[] which[f][io];
}

void fields::connect_periodic_chunks() {
  if (k == -1) return; // It isn't periodic!  :)
  int *nc[2][2];
  for (int f=0;f<2;f++)
    for (int io=0;io<2;io++) {
      nc[f][io] = new int[num_chunks];
      for (int i=0;i<num_chunks;i++) nc[f][io][i] = 0;
    }
  int *new_comm_sizes[2];
  for (int ft=0;ft<2;ft++) new_comm_sizes[ft] = new int[num_chunks*num_chunks];
  for (int i=0;i<num_chunks;i++) {
    const volume vi = chunks[i]->v;
    // First count the border elements in z direction...
    for (int ft=0;ft<2;ft++)
      for (int j=0;j<num_chunks;j++)
        new_comm_sizes[ft][j] = comm_sizes[ft][j+i*num_chunks];
    for (int c=0;c<10;c++)
      if (vi.has_field((component)c))
        for (int n=0;n<vi.ntot();n++) {
          const vec here = vi.loc((component)c, n);
          if (here.z() == v.origin.z() && !vi.owns(here)) {
            // It is on the z == 0 border...
            const vec there = here + lattice_vector();
            for (int j=0;j<num_chunks;j++)
              if (chunks[j]->v.owns(there)) {
                nc[type((component)c)][Incoming][i]++;
                nc[type((component)c)][Outgoing][j]++;
                new_comm_sizes[type((component)c)][j]++;
              }
          } else if (here.z() > v.origin.z() + lattice_vector().z() &&
                     !vi.owns(here)) {
            // This is on the high z border...
            const vec there = here - lattice_vector();
            for (int j=0;j<num_chunks;j++)
              if (chunks[j]->v.owns(there)) {
                nc[type((component)c)][Incoming][i]++;
                nc[type((component)c)][Outgoing][j]++;
                new_comm_sizes[type((component)c)][j]++;
              }
          }
        }
    // Allocating comm blocks as we go...
    for (int ft=0;ft<2;ft++)
      for (int j=0;j<num_chunks;j++) {
        DOCMP {
          delete[] comm_blocks[ft][cmp][j+i*num_chunks];
          comm_blocks[ft][cmp][j+i*num_chunks] = new double[new_comm_sizes[ft][j]];
        }
        comm_sizes[ft][j+i*num_chunks] = new_comm_sizes[ft][j];
      }
  }
  for (int ft=0;ft<2;ft++) delete[] new_comm_sizes[ft];
  // Now allocate the connection arrays...
  for (int i=0;i<num_chunks;i++) {
    for (int f=0;f<2;f++)
      for (int io=0;io<2;io++)
        chunks[i]->alloc_extra_connections((field_type)f,(in_or_out)io,nc[f][io][i]);
  }
  for (int f=0;f<2;f++)
    for (int io=0;io<2;io++) delete[] nc[f][io];
  // Next start setting up the connections...
  int *which[2][2];
  for (int f=0;f<2;f++)
    for (int io=0;io<2;io++)
      which[f][io] = new int[num_chunks];
  for (int f=0;f<2;f++)
    for (int io=0;io<2;io++)
      for (int i=0;i<num_chunks;i++) which[f][io][i] = 0;
  for (int i=0;i<num_chunks;i++) {
    const volume vi = chunks[i]->v;
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
                // index is deprecated, but ok in this case:
                int m = chunks[j]->v.index((component)c, there);
                int ft = type((component) c);
                DOCMP {
                  chunks[i]->connections[ft][Incoming][cmp]
                    [which[ft][Incoming][i]] = chunks[i]->f[c][cmp] + n;
                  chunks[j]->connections[ft][Outgoing][cmp]
                    [which[ft][Outgoing][j]] = chunks[j]->f[c][cmp] + m;
                }
                chunks[i]->connection_phases[ft][which[ft][Incoming][i]] = 1.0;
                which[ft][Incoming][i]++;
                which[ft][Outgoing][j]++;
              }
          }
        }
  }
  for (int f=0;f<2;f++)
    for (int io=0;io<2;io++) delete[] which[f][io];
}

void fields::connect_metallic_chunks() {
  if (v.dim == dcyl) connect_metallic_bigr_chunks();
  if (k != -1.0) return;
  connect_metallic_bigz_chunks();
}

static double zero = 0.0;

void fields::connect_metallic_bigr_chunks() {
  int *nc[2][2];
  for (int f=0;f<2;f++)
    for (int io=0;io<2;io++) {
      nc[f][io] = new int[num_chunks];
      for (int i=0;i<num_chunks;i++) nc[f][io][i] = 0;
    }
  // Note: Only call this function if you are in cylindrical coords!
  for (int i=0;i<num_chunks;i++) {
    const volume vi = chunks[i]->v;
    // First count the radial border elements...
    int new_comm_size[2];
    for (int ft=0;ft<2;ft++) new_comm_size[ft] = comm_sizes[ft][i+i*num_chunks];
    for (int c=0;c<10;c++)
      if (vi.has_field((component)c))
        for (int n=0;n<vi.ntot();n++) {
          const vec here = vi.loc((component)c, n);
          if (here.r() > v.origin.r() + v.nr()*inva && !vi.owns(here)) {
            // It is on the big r border...
            nc[type((component)c)][Incoming][i]++;
            nc[type((component)c)][Outgoing][i]++;
            new_comm_size[type((component)c)]++;
          }
        }
    // Allocating comm blocks as we go...
    for (int ft=0;ft<2;ft++) {
      DOCMP {
        delete[] comm_blocks[ft][cmp][i+i*num_chunks];
        comm_blocks[ft][cmp][i+i*num_chunks] = new double[new_comm_size[ft]];
      }
      comm_sizes[ft][i+i*num_chunks] = new_comm_size[ft];
    }
  }
  // Now allocate the connection arrays...
  for (int i=0;i<num_chunks;i++) {
    for (int f=0;f<2;f++)
      for (int io=0;io<2;io++)
        chunks[i]->alloc_extra_connections((field_type)f,(in_or_out)io,nc[f][io][i]);
  }
  for (int f=0;f<2;f++)
    for (int io=0;io<2;io++) delete[] nc[f][io];
  // Next start setting up the connections...
  for (int i=0;i<num_chunks;i++) {
    const volume vi = chunks[i]->v;
    int *which[2][2];
    for (int f=0;f<2;f++)
      for (int io=0;io<2;io++)
        which[f][io] = new int[num_chunks];
    for (int f=0;f<2;f++)
      for (int io=0;io<2;io++)
        for (int i=0;i<num_chunks;i++) which[f][io][i] = 0;
    for (int c=0;c<10;c++)
      if (vi.has_field((component)c))
        for (int n=0;n<vi.ntot();n++) {
          const vec here = vi.loc((component)c, n);
          if (here.r() > v.origin.r() + v.nr()*inva && !vi.owns(here)) {
            int ft = type((component) c);
            DOCMP {
              chunks[i]->connections[ft][Incoming][cmp][which[ft][Incoming][i]]
                = chunks[i]->f[c][cmp] + n;
              chunks[i]->connections[ft][Outgoing][cmp][which[ft][Outgoing][i]]
                = &zero;
            }
            chunks[i]->connection_phases[ft][which[ft][Incoming][i]] = 1.0;
            which[ft][Incoming][i]++;
            which[ft][Outgoing][i]++;
          }
        }
    for (int f=0;f<2;f++)
      for (int io=0;io<2;io++) delete[] which[f][io];
  }
}

void fields::connect_metallic_bigz_chunks() {
  int *nc[2][2];
  for (int f=0;f<2;f++)
    for (int io=0;io<2;io++) {
      nc[f][io] = new int[num_chunks];
      for (int i=0;i<num_chunks;i++) nc[f][io][i] = 0;
    }
  for (int i=0;i<num_chunks;i++) {
    const volume vi = chunks[i]->v;
    // First count the border elements...
    int new_comm_size[2];
    for (int ft=0;ft<2;ft++) new_comm_size[ft] = comm_sizes[ft][i+i*num_chunks];
    for (int c=0;c<10;c++)
      if (vi.has_field((component)c))
        for (int n=0;n<vi.ntot();n++) {
          const vec here = vi.loc((component)c, n);
          if (here.z() == v.origin.z() + (v.nz()+0.5)*inva && !vi.owns(here)) {
            // It is on the big z border...
            nc[type((component)c)][Incoming][i]++;
            nc[type((component)c)][Outgoing][i]++;
            new_comm_size[type((component)c)]++;
          }
        }
    // Allocating comm blocks as we go...
    for (int ft=0;ft<2;ft++) {
      DOCMP {
        delete[] comm_blocks[ft][cmp][i+i*num_chunks];
        comm_blocks[ft][cmp][i+i*num_chunks] = new double[new_comm_size[ft]];
      }
      comm_sizes[ft][i+i*num_chunks] = new_comm_size[ft];
    }
  }
  // Now allocate the connection arrays...
  for (int i=0;i<num_chunks;i++) {
    for (int f=0;f<2;f++)
      for (int io=0;io<2;io++)
        chunks[i]->alloc_extra_connections((field_type)f,(in_or_out)io,nc[f][io][i]);
  }
  for (int f=0;f<2;f++)
    for (int io=0;io<2;io++) delete[] nc[f][io];
  // Next start setting up the connections...
  for (int i=0;i<num_chunks;i++) {
    const volume vi = chunks[i]->v;
    int *which[2][2];
    for (int f=0;f<2;f++)
      for (int io=0;io<2;io++)
        which[f][io] = new int[num_chunks];
    for (int f=0;f<2;f++)
      for (int io=0;io<2;io++)
        for (int i=0;i<num_chunks;i++) which[f][io][i] = 0;
    for (int c=0;c<10;c++)
      if (vi.has_field((component)c))
        for (int n=0;n<vi.ntot();n++) {
          const vec here = vi.loc((component)c, n);
          if (here.z() == v.origin.z() + (v.nz()+0.5)*inva && !vi.owns(here)) {
            // It is on the big z border...
            int ft = type((component) c);
            DOCMP {
              chunks[i]->connections[ft][Incoming][cmp][which[ft][Incoming][i]]
                = chunks[i]->f[c][cmp] + n;
              chunks[i]->connections[ft][Outgoing][cmp][which[ft][Outgoing][i]]
                = &zero;
            }
            chunks[i]->connection_phases[ft][which[ft][Incoming][i]] = 1.0;
            which[ft][Incoming][i]++;
            which[ft][Outgoing][i]++;
          }
        }
    for (int f=0;f<2;f++)
      for (int io=0;io<2;io++) delete[] which[f][io];
  }
}

void fields_chunk::alloc_extra_connections(field_type f, in_or_out io, int num) {
  if (num == 0) return; // No need to go to any bother...
  const int tot = num_connections[f][io] + num;
  if (io == Incoming) {
    complex<double> *ph = new complex<double>[tot];
    if (!ph) {
      printf("Out of memory!\n");
      exit(1);
    }
    for (int x=0;x<num_connections[f][io];x++)
      ph[num+x] = connection_phases[f][x];
    delete[] connection_phases[f];
    connection_phases[f] = ph;
  }
  DOCMP {
    double **conn = new double *[tot];
    if (!conn) {
      printf("Out of memory!\n");
      exit(1);
    }
    for (int x=0;x<num_connections[f][io];x++)
      conn[num+x] = connections[f][io][cmp][x];
    delete[] connections[f][io][cmp];
    connections[f][io][cmp] = conn;
  }
  num_connections[f][io] += num;
}
