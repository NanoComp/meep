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

#include "meep.h"
#include "meep_internals.h"

#define UNUSED(x) (void) x // silence compiler warnings

namespace meep {

void fields::set_boundary(boundary_side b,direction d,
                          boundary_condition cond,
                          complex<double> kcomponent) {
  UNUSED(kcomponent);
  boundaries[b][d] = cond;
  chunk_connections_valid = false;
}

void fields::use_bloch(direction d, complex<double> kk) {
  k[d] = kk;
  for (int b=0;b<2;b++) boundaries[b][d] = Periodic;
  const complex<double> I = complex<double>(0.0,1.0);
  eikna[d] = exp(I*kk*((2*pi/a)*v.num_direction(d)));
  coskna[d] = real(eikna[d]);
  sinkna[d] = imag(eikna[d]);
  chunk_connections_valid = false;
}

void fields::use_bloch(const vec &k) {
  // Note that I allow a 1D k input when in cylindrical, since in that case
  // it is unambiguous.
  if (k.dim != v.dim && !(k.dim == D1 && v.dim == Dcyl))
    abort("Aaaack, k has wrong dimensions!\n");
  for (int dd=0;dd<5;dd++) {
    const direction d = (direction) dd;
    if (v.has_boundary(Low,d) && d != R)
      use_bloch(d, k.in_direction(d));
  }
  chunk_connections_valid = false;
}

ivec fields::ilattice_vector(direction d) const {
  switch (user_volume.dim) {
  case Dcyl: case D1: return ivec(0,2*user_volume.nz()); // Only Z direction here...
  case D2:
    switch (d) {
    case X: return ivec2d(user_volume.nx()*2,0);
    case Y: return ivec2d(0,user_volume.ny()*2);
    case Z: case R: case P: case NO_DIRECTION: break;
    }
  case D3:
    switch (d) {
    case X: return ivec(user_volume.nx()*2,0,0);
    case Y: return ivec(0,user_volume.ny()*2,0);
    case Z: return ivec(0,0,user_volume.nz()*2);
    case R: case P: case NO_DIRECTION: break;
    }
  }
  abort("Aaack in ilattice_vector.\n");
  return ivec(0);
}

vec fields::lattice_vector(direction d) const {
  return v[ilattice_vector(d)];
}

void fields::disconnect_chunks() {
  chunk_connections_valid = false;
  for (int i=0;i<num_chunks;i++) {
    DOCMP {
      FOR_FIELD_TYPES(f)
        for (int io=0;io<2;io++) {
          delete[] chunks[i]->connections[f][io];
          chunks[i]->connections[f][io] = NULL;
        }
    }
    delete[] chunks[i]->connection_phases[E_stuff];
    delete[] chunks[i]->connection_phases[H_stuff];
    chunks[i]->connection_phases[E_stuff] = 0;
    chunks[i]->connection_phases[H_stuff] = 0;
    FOR_FIELD_TYPES(f)
      for (int io=0;io<2;io++)
        chunks[i]->num_connections[f][io] = 0;
  }
  FOR_FIELD_TYPES(ft)
    for (int i=0;i<num_chunks*num_chunks;i++) {
      delete[] comm_blocks[ft][i];
      comm_blocks[ft][i] = 0;
      comm_sizes[ft][i] = 0;
    }
}

void fields::connect_chunks() {
  if (!chunk_connections_valid) {
    am_now_working_on(Connecting);
    disconnect_chunks();
    find_metals();
    connect_the_chunks();
    finished_working();
    chunk_connections_valid = true;
  }
}

inline bool fields::on_metal_boundary(const ivec &here) {
  LOOP_OVER_DIRECTIONS(v.dim, d) {
    if (user_volume.has_boundary(High, d) &&
        here.in_direction(d) == user_volume.big_corner().in_direction(d)) {
      if (boundaries[High][d] == Metallic) return true;
    }
    if (boundaries[Low][d] == Magnetic &&
        here.in_direction(d) ==
        user_volume.little_corner().in_direction(d)+1)
      return true;
    if (boundaries[Low][d] == Metallic &&
        here.in_direction(d) ==
        user_volume.little_corner().in_direction(d))
      return true;
  }
  return false;
}

bool fields::locate_point_in_user_volume(ivec *there, complex<double> *phase) const {
  // Check if a translational symmetry is needed to bring the point in...
  if (!user_volume.owns(*there))
    FOR_DIRECTIONS(d) {
      if (boundaries[High][d] == Periodic &&
          there->in_direction(d) <= user_volume.little_corner().in_direction(d)) {
        while (there->in_direction(d) <=
               user_volume.little_corner().in_direction(d)) {
          *there += ilattice_vector(d);
          *phase *= conj(eikna[d]);
        }
      } else if (boundaries[High][d] == Periodic &&
                 there->in_direction(d)-ilattice_vector(d).in_direction(d)
                 > user_volume.little_corner().in_direction(d)) {
        while (there->in_direction(d)-ilattice_vector(d).in_direction(d)
               > user_volume.little_corner().in_direction(d)) {
          *there -= ilattice_vector(d);
          *phase *= eikna[d];
        }
      }
    }
  return user_volume.owns(*there);
}

void fields::locate_volume_source_in_user_volume(const vec p1, const vec p2, vec newp1[8], vec newp2[8],
                                                  complex<double> kphase[8], int &ncopies) const {
  // For periodic boundary conditions, 
  // this function locates up to 8 translated copies of the initial volume specified by (p1,p2)
  // First bring center of volume inside
  ncopies = 1;
  newp1[0] = p1;
  newp2[0] = p2;
  kphase[0] = 1;
  vec cen = (newp1[0] + newp2[0]) * 0.5;
  LOOP_OVER_DIRECTIONS(v.dim, d) 
    if (boundaries[High][d] == Periodic)  {
      while (cen.in_direction(d) < v.boundary_location(Low, d)) {
        newp1[0] += lattice_vector(d);
        newp2[0] += lattice_vector(d);
        kphase[0] *= conj(eikna[d]);
        cen = (newp1[0] + newp2[0]) * 0.5;
      }
      while (cen.in_direction(d) > v.boundary_location(High, d)) {
        newp1[0] -= lattice_vector(d);
        newp2[0] -= lattice_vector(d);
        kphase[0] *= eikna[d];
        cen = (newp1[0] + newp2[0]) * 0.5;
      }
    }
  
  // if volume extends outside user_volume in any direction, we need to duplicate already existing copies
  LOOP_OVER_DIRECTIONS(v.dim, d) 
    if (boundaries[High][d] == Periodic) {
      if (newp1[0].in_direction(d) < v.boundary_location(Low, d) ||
          newp2[0].in_direction(d) < v.boundary_location(Low, d)) {
        for (int j=0; j<ncopies; j++) {
          newp1[ncopies+j] = newp1[j] + lattice_vector(d);
          newp2[ncopies+j] = newp2[j] + lattice_vector(d);
          kphase[ncopies+j] = kphase[j] * conj(eikna[d]);
        }
        ncopies *= 2;
      }
      else if (newp1[0].in_direction(d) > v.boundary_location(High, d) ||
               newp2[0].in_direction(d) > v.boundary_location(High, d)) {
        for (int j=0; j<ncopies; j++) {
          newp1[ncopies+j] = newp1[j] - lattice_vector(d);
          newp2[ncopies+j] = newp2[j] - lattice_vector(d);
          kphase[ncopies+j] = kphase[j] * eikna[d];
        }
        ncopies *= 2;
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
      FOR_FIELD_TYPES(ft) {
        delete[] chunks[i]->zeroes[ft];
        // First electric components...
        chunks[i]->num_zeroes[ft] = 0;
        DOCMP FOR_COMPONENTS(c)
          if (type(c) == ft && chunks[i]->f[c][cmp])
	    LOOP_OVER_VOL_OWNED(vi, c, n) 
	      if (IVEC_LOOP_AT_BOUNDARY) { // todo: just loop over boundaries
		IVEC_LOOP_ILOC(vi, here);
		if (on_metal_boundary(here))
		  chunks[i]->num_zeroes[ft]++;
	      }
        typedef double *double_ptr;
        chunks[i]->zeroes[ft] = new double_ptr[chunks[i]->num_zeroes[ft]];
        int num = 0;
        DOCMP FOR_COMPONENTS(c)
          if (type(c) == ft && chunks[i]->f[c][cmp])
	    LOOP_OVER_VOL_OWNED(vi, c, n)
	      if (IVEC_LOOP_AT_BOUNDARY) { // todo: just loop over boundaries
		IVEC_LOOP_ILOC(vi, here);
		if (on_metal_boundary(here))
		  chunks[i]->zeroes[ft][num++] = &(chunks[i]->f[c][cmp][n]);
	      }
      }
    }
}

void fields::connect_the_chunks() {
  int *nc[NUM_FIELD_TYPES][2];
  FOR_FIELD_TYPES(f)
    for (int io=0;io<2;io++) {
      nc[f][io] = new int[num_chunks];
      for (int i=0;i<num_chunks;i++) nc[f][io][i] = 0;
    }
  for (int i=0;i<num_chunks;i++) {
    // First count the border elements...
    const volume vi = chunks[i]->v;
    for (int j=0;j<num_chunks;j++)
      FOR_FIELD_TYPES(ft) {
        comm_sizes[ft][j+i*num_chunks] = 0;
        comm_num_complex[ft][j+i*num_chunks] = 0;
        comm_num_negate[ft][j+i*num_chunks] = 0;
      }
    for (int j=0;j<num_chunks;j++) {
      FOR_COMPONENTS(corig)
        if (have_component(corig))
	  LOOP_OVER_VOL_NOTOWNED(vi, corig, n) {
	    IVEC_LOOP_ILOC(vi, here);
            component c = corig;
	    // We're looking at a border element...
	    complex<double> thephase = 1.0;
	    if (locate_component_point(&c,&here,&thephase))
	      if (chunks[j]->v.owns(here) && !on_metal_boundary(here)) {
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
      FOR_ELECTRIC_COMPONENTS(corig)
        if (have_component(corig))
	  LOOP_OVER_VOL_NOTOWNED(vi, corig, n) {
	    IVEC_LOOP_ILOC(vi, here);
            component c = corig;
	    // We're looking at a border element...
	    complex<double> thephase = 1.0;
	    if (locate_component_point(&c,&here,&thephase))
	      if (chunks[j]->v.owns(here) && !on_metal_boundary(here)) {
		// Adjacent, periodic or rotational...
		int common_pols = 0;
		for (polarization *pi = chunks[i]->pol; pi; pi = pi->next)
		  for (polarization *pj = chunks[j]->pol; pj; pj = pj->next)
		    if (pi->pb->get_identifier() == pj->pb->get_identifier())
		      common_pols += 1;
		const int nn = (is_real?1:2) * common_pols;
		const int pair = j+i*num_chunks;
		if (imag(thephase) != 0.0 || fabs(real(thephase)) != 1.0)
		  comm_num_complex[P_stuff][pair] += 2*nn;
		else if (thephase == -1.0)
		  comm_num_negate[P_stuff][pair] += 2*nn;
		nc[P_stuff][Incoming][i] += 2*nn;
		nc[P_stuff][Outgoing][j] += 2*nn;
		comm_sizes[P_stuff][pair] += 2*nn;
		// Note above that the factor of two in 2*nn comes from
		// the fact that we have two polarization arrays, pol and
		// olpol.
	      }
        }
    }
    // Allocating comm blocks as we go...
    FOR_FIELD_TYPES(ft)
      for (int j=0;j<num_chunks;j++) {
        delete[] comm_blocks[ft][j+i*num_chunks];
        comm_blocks[ft][j+i*num_chunks] =
          new double[comm_sizes[ft][j+i*num_chunks]];
      }
  }
  // Now allocate the connection arrays...
  for (int i=0;i<num_chunks;i++)
    FOR_FIELD_TYPES(f)
      for (int io=0;io<2;io++)
        chunks[i]->alloc_extra_connections((field_type)f,
                                           (in_or_out)io,nc[f][io][i]);
  FOR_FIELD_TYPES(f)
    for (int io=0;io<2;io++) delete[] nc[f][io];
  // Next start setting up the connections...
  int *wh[NUM_FIELD_TYPES][2];
  FOR_FIELD_TYPES(f)
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
	  LOOP_OVER_VOL_NOTOWNED(vi, corig, n) {
	    IVEC_LOOP_ILOC(vi, here);
            component c = corig;
#define FT (type(c))
	    // We're looking at a border element...
	    complex<double> thephase = 1.0;
	    if (locate_component_point(&c,&here,&thephase))
	      if (chunks[j]->v.owns(here) && !on_metal_boundary(here) &&
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
      DOCMP FOR_COMPONENTS(corig)
        if (have_component(corig))
	  LOOP_OVER_VOL_NOTOWNED(vi, corig, n) {
	    IVEC_LOOP_ILOC(vi, here);
            component c = corig;
	    // We're looking at a border element...
	    complex<double> thephase = 1.0;
	    if (locate_component_point(&c,&here,&thephase))
	      if (chunks[j]->v.owns(here) && !on_metal_boundary(here) &&
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
      DOCMP FOR_COMPONENTS(corig)
        if (have_component(corig))
	  LOOP_OVER_VOL_NOTOWNED(vi, corig, n) {
	    IVEC_LOOP_ILOC(vi, here);
            component c = corig;
	    // We're looking at a border element...
	    complex<double> thephase = 1.0;
	    if (locate_component_point(&c,&here,&thephase))
	      if (chunks[j]->v.owns(here) && !on_metal_boundary(here) &&
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

      // Now connect up the polarizations...
#undef FT
#define FT P_stuff
      // First the complex connections:
      FOR_ELECTRIC_COMPONENTS(corig)
        if (have_component(corig))
	  LOOP_OVER_VOL_NOTOWNED(vi, corig, n) {
	    IVEC_LOOP_ILOC(vi, here);
            component c = corig;
	    // We're looking at a border element...
	    complex<double> thephase = 1.0;
	    if (locate_component_point(&c,&here,&thephase))
	      if (chunks[j]->v.owns(here) && !on_metal_boundary(here))
		if (imag(thephase) != 0.0 || fabs(real(thephase)) != 1.0) {
		  for (polarization *pi = chunks[i]->pol; pi; pi = pi->next)
		    for (polarization *pj = chunks[j]->pol; pj; pj = pj->next)
		      if (pi->pb->get_identifier() == pj->pb->get_identifier()) {
			const int m = chunks[j]->v.index(c, here);
			chunks[i]->connection_phases[FT][wh[FT][Incoming][i]/2] = thephase;
			DOCMP {
			  chunks[i]->connections[FT][Incoming][wh[FT][Incoming][i]++]
			    = pi->P[corig][cmp] + n;
			  chunks[j]->connections[FT][Outgoing][wh[FT][Outgoing][j]++]
			    = pj->P[c][cmp] + m;
			}
		      }
		  // Now do same thing with olpol (could this be uglier?)
		  for (polarization *pi = chunks[i]->olpol; pi; pi = pi->next)
		    for (polarization *pj = chunks[j]->olpol; pj; pj = pj->next)
		      if (pi->pb->get_identifier() == pj->pb->get_identifier()) {
			const int m = chunks[j]->v.index(c, here);
			chunks[i]->connection_phases[FT][wh[FT][Incoming][i]/2] = thephase;
			DOCMP {
			  chunks[i]->connections[FT][Incoming][wh[FT][Incoming][i]++]
			    = pi->P[corig][cmp] + n;
			  chunks[j]->connections[FT][Outgoing][wh[FT][Outgoing][j]++]
			    = pj->P[c][cmp] + m;
			}
		      }
		}
          }
      // Now the negative connections:
      FOR_ELECTRIC_COMPONENTS(corig)
        if (have_component(corig))
	  LOOP_OVER_VOL_NOTOWNED(vi, corig, n) {
	    IVEC_LOOP_ILOC(vi, here);
            component c = corig;
	    // We're looking at a border element...
	    complex<double> thephase = 1.0;
	    if (locate_component_point(&c,&here,&thephase))
	      if (chunks[j]->v.owns(here) && !on_metal_boundary(here))
		if (thephase == -1.0) {
		  for (polarization *pi = chunks[i]->pol; pi; pi = pi->next)
		    for (polarization *pj = chunks[j]->pol; pj; pj = pj->next)
		      if (pi->pb->get_identifier() == pj->pb->get_identifier()) {
			const int m = chunks[j]->v.index(c, here);
			DOCMP {
			  chunks[i]->connections[FT][Incoming][wh[FT][Incoming][i]++]
			    = pi->P[corig][cmp] + n;
			  chunks[j]->connections[FT][Outgoing][wh[FT][Outgoing][j]++]
			    = pj->P[c][cmp] + m;
			}
		      }
		  // Now do same thing with olpol (could this be uglier?)
		  for (polarization *pi = chunks[i]->olpol; pi; pi = pi->next)
		    for (polarization *pj = chunks[j]->olpol; pj; pj = pj->next)
		      if (pi->pb->get_identifier() == pj->pb->get_identifier()) {
			const int m = chunks[j]->v.index(c, here);
			DOCMP {
			  chunks[i]->connections[FT][Incoming][wh[FT][Incoming][i]++]
			    = pi->P[corig][cmp] + n;
			  chunks[j]->connections[FT][Outgoing][wh[FT][Outgoing][j]++]
			    = pj->P[c][cmp] + m;
			}
		      }
		}
          }
      // Finally the plain old copied connections:
      FOR_ELECTRIC_COMPONENTS(corig)
        if (have_component(corig))
	  LOOP_OVER_VOL_NOTOWNED(vi, corig, n) {
	    IVEC_LOOP_ILOC(vi, here);
            component c = corig;
	    // We're looking at a border element...
	    complex<double> thephase = 1.0;
	    if (locate_component_point(&c,&here,&thephase))
	      if (chunks[j]->v.owns(here) && !on_metal_boundary(here))
		if (thephase == 1.0) {
		  for (polarization *pi = chunks[i]->pol; pi; pi = pi->next)
		    for (polarization *pj = chunks[j]->pol; pj; pj = pj->next)
		      if (pi->pb->get_identifier() == pj->pb->get_identifier()) {
			const int m = chunks[j]->v.index(c, here);
			DOCMP {
			  chunks[i]->connections[FT][Incoming][wh[FT][Incoming][i]++]
			    = pi->P[corig][cmp] + n;
			  chunks[j]->connections[FT][Outgoing][wh[FT][Outgoing][j]++]
			    = pj->P[c][cmp] + m;
			}
		      }
		  // Now do same thing with olpol (could this be uglier?)
		  for (polarization *pi = chunks[i]->olpol; pi; pi = pi->next)
		    for (polarization *pj = chunks[j]->olpol; pj; pj = pj->next)
		      if (pi->pb->get_identifier() == pj->pb->get_identifier()) {
			const int m = chunks[j]->v.index(c, here);
			DOCMP {
			  chunks[i]->connections[FT][Incoming][wh[FT][Incoming][i]++]
			    = pi->P[corig][cmp] + n;
			  chunks[j]->connections[FT][Outgoing][wh[FT][Outgoing][j]++]
			    = pj->P[c][cmp] + m;
			}
		      }
		}
          }
      // Finished connecting up the polarizations...
    }
  }
  FOR_FIELD_TYPES(f)
    for (int io=0;io<2;io++) delete[] wh[f][io];
}

void fields_chunk::alloc_extra_connections(field_type f, in_or_out io, int num) {
  if (num == 0) return; // No need to go to any bother...
  const int tot = num_connections[f][io] + num;
  if (io == Incoming) {
    delete[] connection_phases[f];
    connection_phases[f] = new complex<double>[tot]; // This is larger than necesary...
  }
  typedef double *double_ptr;
  double **conn = new double_ptr[tot];
  if (!conn) abort("Out of memory!\n");
  delete[] connections[f][io];
  connections[f][io] = conn;
  num_connections[f][io] += num;
}

} // namespace meep
