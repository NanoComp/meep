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

void fields::set_boundary(boundary_side b,direction d,
                          boundary_condition cond, bool autoconnect,
                          complex<double> kcomponent) {
  boundaries[b][d] = cond;
  if (autoconnect) connect_chunks();
}

void fields::use_metal_everywhere() {
  for (int b=0;b<2;b++)
    for (int d=0;d<5;d++)
      if (v.has_boundary((boundary_side)b, (direction)d))
        set_boundary((boundary_side)b,(direction)d,Metallic, false);
  connect_chunks();
}

void fields::use_bloch(direction d, complex<double> kk, bool autoconnect) {
  k[d] = kk;
  for (int b=0;b<2;b++) boundaries[b][d] = Periodic;
  const complex<double> I = complex<double>(0.0,1.0);
  eikna[d] = exp(I*kk*((2*pi)*inva*v.num_direction(d)));
  coskna[d] = real(eikna[d]);
  sinkna[d] = imag(eikna[d]);
  if (autoconnect) connect_chunks();
}

void fields::use_bloch(const vec &k, bool autoconnect) {
  // Note that I allow a 1D k input when in cylindrical, since in that case
  // it is unambiguous.
  if (k.dim != v.dim && !(k.dim == D1 && v.dim == Dcyl))
    abort("Aaaack, k has wrong dimensions!\n");
  for (int dd=0;dd<5;dd++) {
    const direction d = (direction) dd;
    if (v.has_boundary(Low,d) && d != R)
      use_bloch(d, k.in_direction(d), false);
  }
  if (autoconnect) connect_chunks();
}

ivec fields::ilattice_vector(direction d) const {
  switch (v.dim) {
  case Dcyl: case D1: return ivec(0,2*v.nz()); // Only Z direction here...
  case D2:
    switch (d) {
    case X: return ivec2d(v.nx()*2,0);
    case Y: return ivec2d(0,v.ny()*2);
    }
  case D3:
    switch (d) {
    case X: return ivec(v.nx()*2,0,0);
    case Y: return ivec(0,v.ny()*2,0);
    case Z: return ivec(0,0,v.nz()*2);
    }
  }
}

vec fields::lattice_vector(direction d) const {
  if (v.dim == Dcyl) {
    return vec(0,v.nz()*inva); // Only Z direction here...
  } else if (v.dim == D1) {
    return vec(v.nz()*inva); // Only Z direction here...
  } else if (v.dim == D2) {
    switch (d) {
    case X: return vec2d(v.nx()*inva,0);
    case Y: return vec2d(0,v.ny()*inva);
    }
  } else if (v.dim == D3) {
    switch (d) {
    case X: return vec(v.nx()*inva,0,0);
    case Y: return vec(0,v.ny()*inva,0);
    case Z: return vec(0,0,v.nz()*inva);
    }
  }
}

void fields::disconnect_chunks() {
  for (int i=0;i<num_chunks;i++) {
    DOCMP {
      for (int f=0;f<2;f++)
        for (int io=0;io<2;io++) {
          delete[] chunks[i]->connections[f][io];
          chunks[i]->connections[f][io] = NULL;
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
  am_now_working_on(Connecting);
  disconnect_chunks();
  find_metals();
  connect_the_chunks();
  finished_working();
}

static double zero = 0.0;

inline int fields::is_metal(const ivec &here) {
  if (!user_volume.owns(here)) return 0;
  LOOP_OVER_DIRECTIONS(v.dim, d) {
    if (user_volume.has_boundary(High, d) &&
        here.in_direction(d) == user_volume.big_corner().in_direction(d)) {
      if (boundaries[High][d] == Metallic) return true;
    }
    if (user_volume.has_boundary(Low, d)) {
      if (boundaries[Low][d] == Magnetic &&
          here.in_direction(d) ==
          user_volume.little_corner().in_direction(d)+1)
        return true;
    }
  }
  return false;
}

bool fields::locate_point_in_user_volume(ivec *there, complex<double> *phase) const {
  // Check if a translational symmetry is needed to bring the point in...
  if (!user_volume.owns(*there))
    FOR_DIRECTIONS(d) {
      if (boundaries[Low][d] == Periodic &&
          there->in_direction(d) <= user_volume.origin.in_direction(d)) {
        while (there->in_direction(d) <=
               user_volume.origin.in_direction(d)) {
          *there += ilattice_vector(d);
          *phase *= conj(eikna[d]);
        }
      } else if (boundaries[High][d] == Periodic &&
                 there->in_direction(d)-ilattice_vector(d).in_direction(d)
                 > user_volume.origin.in_direction(d)) {
        while (there->in_direction(d)-ilattice_vector(d).in_direction(d)
               > user_volume.origin.in_direction(d)) {
          *there -= ilattice_vector(d);
          *phase *= eikna[d];
        }
      }
    }
}

bool fields::locate_component_point(component *c, ivec *there,
                                    complex<double> *phase) const {
  // returns true if this point and component exist in the user_volume.  If
  // that is the case, on return *c and *there store the component and
  // location of where the point actually is, and *phase determines holds
  // the phase needed to get the true field.  If the point is not located,
  // *c and *there will hold undefined values.

  // Check if nothing tricky is needed...
  *phase = 1.0;
  if (!locate_point_in_user_volume(there, phase)) return false;
  // Check if a rotation or inversion brings the point in...
  if (user_volume.owns(*there))
    for (int sn=0;sn<S.multiplicity();sn++) {
      const ivec here=S.transform(*there,sn);
      if (v.owns(here)) {
        *there = here;
        *phase *= S.phase_shift(*c,sn);
        *c = direction_component(*c,
                                 S.transform(component_direction(*c),sn).d);
        return true;
      }
    }
  return false;
}

void fields_chunk::zero_metal(field_type ft) {
  for (int i=0;i<num_zeroes[ft];i++) *(zeroes[ft][i]) = 0.0;
}

void fields::find_metals() {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine()) {
      const volume vi = chunks[i]->v;
      delete[] chunks[i]->zeroes[E_stuff];
      delete[] chunks[i]->zeroes[H_stuff];
      // First electric components...
      chunks[i]->num_zeroes[E_stuff] = 0;
      DOCMP FOR_ELECTRIC_COMPONENTS(c)
        if (chunks[i]->f[c][cmp])
          for (int n=0;n<vi.ntot();n++)
            if (vi.owns(vi.iloc(c,n)) && is_metal(vi.iloc(c,n)))
              chunks[i]->num_zeroes[E_stuff]++;
      chunks[i]->zeroes[E_stuff] =
        new (double *)[chunks[i]->num_zeroes[E_stuff]];
      int num = 0;
      DOCMP FOR_ELECTRIC_COMPONENTS(c)
        if (chunks[i]->f[c][cmp])
          for (int n=0;n<vi.ntot();n++)
            if (vi.owns(vi.iloc(c,n)) && is_metal(vi.iloc(c,n)))
              chunks[i]->zeroes[E_stuff][num++] = &(chunks[i]->f[c][cmp][n]);
      // Now magnetic components
      chunks[i]->num_zeroes[H_stuff] = 0;
      DOCMP FOR_MAGNETIC_COMPONENTS(c)
        if (chunks[i]->f[c][cmp])
          for (int n=0;n<vi.ntot();n++)
            if (vi.owns(vi.iloc(c,n)) && is_metal(vi.iloc(c,n)))
              chunks[i]->num_zeroes[H_stuff]++;
      chunks[i]->zeroes[H_stuff] =
        new (double *)[chunks[i]->num_zeroes[H_stuff]];
      num = 0;
      DOCMP FOR_MAGNETIC_COMPONENTS(c)
        if (chunks[i]->f[c][cmp])
          for (int n=0;n<vi.ntot();n++)
            if (vi.owns(vi.iloc(c,n)) && is_metal(vi.iloc(c,n)))
              chunks[i]->zeroes[H_stuff][num++] = &(chunks[i]->f[c][cmp][n]);
    }
}

void fields::connect_the_chunks() {
  int *nc[2][2];
  for (int f=0;f<2;f++)
    for (int io=0;io<2;io++) {
      nc[f][io] = new int[num_chunks];
      for (int i=0;i<num_chunks;i++) nc[f][io][i] = 0;
    }
  for (int i=0;i<num_chunks;i++) {
    // First count the border elements...
    const volume vi = chunks[i]->v;
    for (int j=0;j<num_chunks;j++)
      for (int ft=0;ft<2;ft++) {
        comm_sizes[ft][j+i*num_chunks] = 0;
        comm_num_complex[ft][j+i*num_chunks] = 0;
        comm_num_negate[ft][j+i*num_chunks] = 0;
      }
    for (int j=0;j<num_chunks;j++)
      FOR_COMPONENTS(corig)
        if (have_component(corig))
          for (int n=0;n<vi.ntot();n++) {
            component c = corig;
            ivec here = vi.iloc(c, n);
            if (!vi.owns(here)) {
              // We're looking at a border element...
              complex<double> thephase = 1.0;
              if (locate_component_point(&c,&here,&thephase))
                if (chunks[j]->v.owns(here) && !is_metal(here)) {
                  // Adjacent, periodic or rotational...
                  const int nn = is_real?1:2;
                  const int pair = j+i*num_chunks;
                  if (imag(thephase) != 0.0 || fabs(real(thephase)) != 1.0)
                    comm_num_complex[type(corig)][pair] += nn;
                  else if (thephase == -1.0)
                    comm_num_negate[type(corig)][pair] += nn;
                  nc[type(corig)][Incoming][i] += nn;
                  nc[type(c)][Outgoing][j] += nn;
                  comm_sizes[type(c)][pair] += nn;
                }
            }
        }
    // Allocating comm blocks as we go...
    for (int ft=0;ft<2;ft++)
      for (int j=0;j<num_chunks;j++) {
        delete[] comm_blocks[ft][j+i*num_chunks];
        comm_blocks[ft][j+i*num_chunks] =
          new double[comm_sizes[ft][j+i*num_chunks]];
      }
  }
  // Now allocate the connection arrays...
  for (int i=0;i<num_chunks;i++)
    for (int f=0;f<2;f++)
      for (int io=0;io<2;io++)
        chunks[i]->alloc_extra_connections((field_type)f,
                                           (in_or_out)io,nc[f][io][i]);
  for (int f=0;f<2;f++)
    for (int io=0;io<2;io++) delete[] nc[f][io];
  // Next start setting up the connections...
  int *wh[2][2];
  for (int f=0;f<2;f++)
    for (int io=0;io<2;io++) {
      wh[f][io] = new int[num_chunks];
      for (int i=0;i<num_chunks;i++) wh[f][io][i] = 0;
    }
  // First look for connections with complex phases...
  for (int i=0;i<num_chunks;i++) {
    const volume vi = chunks[i]->v;
    for (int j=0;j<num_chunks;j++) {
      FOR_COMPONENTS(corig)
        if (have_component(corig))
          for (int n=0;n<vi.ntot();n++) {
            component c = corig;
#define FT (type(c))
            ivec here = vi.iloc(c, n);
            if (!vi.owns(here)) {
              // We're looking at a border element...
              complex<double> thephase = 1.0;
              if (locate_component_point(&c,&here,&thephase))
                if (chunks[j]->v.owns(here) && !is_metal(here) &&
                    (imag(thephase) != 0.0 || fabs(real(thephase))!=1.0)) {
                  // Periodic, probably...
                  // index is deprecated, but ok in this case:
                  const int m = chunks[j]->v.index(c, here);
                  chunks[i]->connection_phases[FT][wh[FT][Incoming][i]/2] = thephase;
                  DOCMP {
                    chunks[i]->connections[FT][Incoming][wh[FT][Incoming][i]++]
                      = chunks[i]->f[corig][cmp] + n;
                    chunks[j]->connections[FT][Outgoing][wh[FT][Outgoing][j]++]
                      = chunks[j]->f[c][cmp] + m;
                  }
                }
            }
          }
      DOCMP FOR_COMPONENTS(corig)
        if (have_component(corig))
          for (int n=0;n<vi.ntot();n++) {
            component c = corig;
            ivec here = vi.iloc(c, n);
            if (!vi.owns(here)) {
              // We're looking at a border element...
              complex<double> thephase = 1.0;
              if (locate_component_point(&c,&here,&thephase))
                if (chunks[j]->v.owns(here) && !is_metal(here) &&
                    thephase == -1.0) {
                  // Adjacent, periodic or rotational...
                  // index is deprecated, but ok in this case:
                  const int m = chunks[j]->v.index(c, here);
                  chunks[i]->connections[FT][Incoming][wh[FT][Incoming][i]++]
                    = chunks[i]->f[corig][cmp] + n;
                  chunks[j]->connections[FT][Outgoing][wh[FT][Outgoing][j]++]
                    = chunks[j]->f[c][cmp] + m;
                }
            }
          }
      DOCMP FOR_COMPONENTS(corig)
        if (have_component(corig))
          for (int n=0;n<vi.ntot();n++) {
            component c = corig;
            ivec here = vi.iloc(c, n);
            if (!vi.owns(here)) {
              // We're looking at a border element...
              complex<double> thephase = 1.0;
              if (locate_component_point(&c,&here,&thephase))
                if (chunks[j]->v.owns(here) && !is_metal(here) &&
                    thephase == 1.0) {
                  // Adjacent, periodic or rotational...
                  // index is deprecated, but ok in this case:
                  const int m = chunks[j]->v.index(c, here);
                  chunks[i]->connections[FT][Incoming][wh[FT][Incoming][i]++]
                    = chunks[i]->f[corig][cmp] + n;
                  chunks[j]->connections[FT][Outgoing][wh[FT][Outgoing][j]++]
                    = chunks[j]->f[c][cmp] + m;
                }
            }
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
    delete[] connection_phases[f];
    connection_phases[f] = new complex<double>[tot]; // This is larger than necesary...
  }
  double **conn = new (double *)[tot];
  if (!conn) abort("Out of memory!\n");
  delete[] connections[f][io];
  connections[f][io] = conn;
  num_connections[f][io] += num;
}
