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
      delete[] comm_blocks[ft][i];
      comm_blocks[ft][i] = 0;
      comm_sizes[ft][i] = 0;
    }
}

void fields::connect_chunks() {
  disconnect_chunks();
  connect_the_chunks();
}

static double zero = 0.0;

inline int fields::is_metal(const vec &here, const volume &vi) {
  if (!v.owns(here)) return 0;
  return
    // Check if it is on the big r border...
    (v.dim == dcyl && here.r() > v.origin.r() + (v.nr()-0.2)*inva) ||
    // Check if it is on the big z border...
    (k == -1 && here.z() > v.origin.z() + (v.nz()-0.2)*inva);
}

void fields::connect_the_chunks() {
  int *nc[2][2];
  for (int f=0;f<2;f++)
    for (int io=0;io<2;io++) {
      nc[f][io] = new int[num_chunks];
      for (int i=0;i<num_chunks;i++) nc[f][io][i] = 0;
    }
  int *new_comm_sizes[2];
  for (int ft=0;ft<2;ft++) new_comm_sizes[ft] = new int[num_chunks];
  for (int i=0;i<num_chunks;i++) {
    // First count the border elements...
    const volume vi = chunks[i]->v;
    for (int j=0;j<num_chunks;j++)
      for (int ft=0;ft<2;ft++)
        new_comm_sizes[ft][j] = comm_sizes[ft][j+i*num_chunks];
    for (int j=0;j<num_chunks;j++)
      for (int c=0;c<10;c++)
        if (vi.has_field((component)c))
          for (int n=0;n<vi.ntot();n++) {
            const vec here = vi.loc((component)c, n);
            if (!vi.owns(here)) {
              // We're looking at a border element...
              if (j != i && chunks[j]->v.owns(here)) {
                // Adjacent...
                nc[type((component)c)][Incoming][i]++;
                nc[type((component)c)][Outgoing][j]++;
                new_comm_sizes[type((component)c)][j]++;
              } else if (k != -1 && here.z() == v.origin.z() &&
                         chunks[j]->v.owns(here + lattice_vector())) {
                // Periodic...
                nc[type((component)c)][Incoming][i]++;
                nc[type((component)c)][Outgoing][j]++;
                new_comm_sizes[type((component)c)][j]++;
              } else if (k != -1 && here.z() > v.origin.z() + lattice_vector().z() &&
                         chunks[j]->v.owns(here - lattice_vector())) {
                // Periodic...
                nc[type((component)c)][Incoming][i]++;
                nc[type((component)c)][Outgoing][j]++;
                new_comm_sizes[type((component)c)][j]++;
              }
            } else if (j == i && is_metal(here, vi)) {
              // Try metallic...
              nc[type((component)c)][Incoming][i]++;
              nc[type((component)c)][Outgoing][i]++;
              new_comm_sizes[type((component)c)][j]++;
            }
        }
    // Allocating comm blocks as we go...
    for (int ft=0;ft<2;ft++)
      for (int j=0;j<num_chunks;j++) {
        delete[] comm_blocks[ft][j+i*num_chunks];
        comm_blocks[ft][j+i*num_chunks] = new double[2*new_comm_sizes[ft][j]];
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
  int *wh[2][2];
  for (int f=0;f<2;f++)
    for (int io=0;io<2;io++) {
      wh[f][io] = new int[num_chunks];
      for (int i=0;i<num_chunks;i++) wh[f][io][i] = 0;
    }
  for (int i=0;i<num_chunks;i++) {
    const volume vi = chunks[i]->v;
    for (int j=0;j<num_chunks;j++)
      for (int c=0;c<10;c++)
        if (vi.has_field((component)c))
          for (int n=0;n<vi.ntot();n++) {
            const field_type ft = type((component) c);
            const vec here = vi.loc((component)c, n);
            if (!vi.owns(here)) {
              // This means we're looking at a border element...
              if (j != i && chunks[j]->v.owns(here)) {
                // index is deprecated, but ok in this case:
                const int m = chunks[j]->v.index((component)c, here);
                DOCMP {
                  chunks[i]->connections[ft][Incoming][cmp][wh[ft][Incoming][i]]
                    = chunks[i]->f[c][cmp] + n;
                  chunks[j]->connections[ft][Outgoing][cmp][wh[ft][Outgoing][j]]
                    = chunks[j]->f[c][cmp] + m;
                }
                chunks[i]->connection_phases[ft][wh[ft][Incoming][i]] = 1.0;
                wh[ft][Incoming][i]++;
                wh[ft][Outgoing][j]++;
              } else if (k != -1 && here.z() == v.origin.z() &&
                         chunks[j]->v.owns(here + lattice_vector())) {
                // Periodic...
                // index is deprecated, but ok in this case:
                const int m = chunks[j]->v.index((component)c, here + lattice_vector());
                DOCMP {
                  chunks[i]->connections[ft][Incoming][cmp][wh[ft][Incoming][i]]
                    = chunks[i]->f[c][cmp] + n;
                  chunks[j]->connections[ft][Outgoing][cmp][wh[ft][Outgoing][j]]
                    = chunks[j]->f[c][cmp] + m;
                }
                chunks[i]->connection_phases[ft][wh[ft][Incoming][i]] = eiknz;
                wh[ft][Incoming][i]++;
                wh[ft][Outgoing][j]++;
              } else if (k != -1 && here.z() > v.origin.z() + lattice_vector().z() &&
                         chunks[j]->v.owns(here - lattice_vector())) {
                // Periodic...
                // index is deprecated, but ok in this case:
                const int m = chunks[j]->v.index((component)c, here - lattice_vector());
                DOCMP {
                  chunks[i]->connections[ft][Incoming][cmp][wh[ft][Incoming][i]]
                    = chunks[i]->f[c][cmp] + n;
                  chunks[j]->connections[ft][Outgoing][cmp][wh[ft][Outgoing][j]]
                    = chunks[j]->f[c][cmp] + m;
                }
                chunks[i]->connection_phases[ft][wh[ft][Incoming][i]] = conj(eiknz);
                wh[ft][Incoming][i]++;
                wh[ft][Outgoing][j]++;
              }
            } else if (j == i && is_metal(here, vi)) {
              // Try metallic...
              DOCMP {
                chunks[i]->connections[ft][Incoming][cmp][wh[ft][Incoming][i]]
                  = chunks[i]->f[c][cmp] + n;
                chunks[j]->connections[ft][Outgoing][cmp][wh[ft][Outgoing][j]]
                  = &zero;
              }
              chunks[i]->connection_phases[ft][wh[ft][Incoming][i]] = 1.0;
              wh[ft][Incoming][i]++;
              wh[ft][Outgoing][j]++;
            }
          }
  }
  for (int f=0;f<2;f++)
    for (int io=0;io<2;io++) delete[] wh[f][io];
}

void fields_chunk::alloc_extra_connections(field_type f, in_or_out io, int num) {
  if (num == 0) return; // No need to go to any bother...
  const int tot = num_connections[f][io] + num;
  if (io == Incoming) {
    complex<double> *ph = new complex<double>[tot];
    if (!ph) abort("Out of memory!\n");
    for (int x=0;x<num_connections[f][io];x++)
      ph[num+x] = connection_phases[f][x];
    delete[] connection_phases[f];
    connection_phases[f] = ph;
  }
  DOCMP {
    double **conn = new double *[tot];
    if (!conn) abort("Out of memory!\n");
    for (int x=0;x<num_connections[f][io];x++)
      conn[num+x] = connections[f][io][cmp][x];
    delete[] connections[f][io][cmp];
    connections[f][io][cmp] = conn;
  }
  num_connections[f][io] += num;
}
