/* Copyright (C) 2005-2009 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* This file implements dispersive materials for Meep via a polarization P = \chi(\omega) W,
   where W is e.g. E or H.  Each subclass of the susceptibility class should implement a different
   type of \chi(\omega).  The subclass knows how to timestep P given W at the current (and possibly
   previous) timestep, and any additional internal data that needs to be allocated along with P. 

   Each \chi(\omega) is spatially multiplied by a (scalar) sigma array.  The meep::fields class is
   responsible for allocating P and sigma and passing them to susceptibility::update_P. */


#include <string.h>
#include "meep.hpp"
#include "meep_internals.hpp"

namespace meep {

int susceptibility::cur_id = 0;

susceptibility *susceptibility::clone() const {
  susceptibility *sus = new susceptibility(*this);
  sus->next = 0;
  sus->ntot = ntot;
  sus->id = id;
  FOR_COMPONENTS(c) FOR_DIRECTIONS(d) {
    if (sigma[c][d]) {
      sus->sigma[c][d] = new realnum[ntot];
      memcpy(sus->sigma[c][d], sigma[c][d], sizeof(realnum) * ntot);
    }
    else sus->sigma[c][d] = NULL;
    sus->trivial_sigma[c][d] = trivial_sigma[c][d];
  }
  return sus;
}

/* Return whether or not we need to allocate P[c].  (We don't need to
   allocate P[c] if we can be sure it will be zero.)

   We are a bit wasteful because if sigma is nontrivial in *any* chunk,
   we allocate the corresponding P on *every* owned chunk.  This greatly
   simplifies communication in boundaries.cpp, because we can be sure that
   one chunk has a P then any chunk it borders has the same P, so we don't
   have to worry about communicating with something that doesn't exist.
   TODO: reduce memory usage (bookkeeping seem much harder, though).
*/
bool susceptibility::needs_P(component c, realnum *W[NUM_FIELD_COMPONENTS][2])
  const {
  if (!is_electric(c) && !is_magnetic(c)) return false;
  FOR_DIRECTIONS(d)
    if (!trivial_sigma[c][d] && W[direction_component(c, d)][0]) return true;
  return false;
}

/* return whether we need the notowned parts of the W field --
   by default, this is only the case if sigma has offdiagonal components
   coupling P to W.   (See needs_P: again, this true if the notowned
   W is needed in *any* chunk.) */
bool susceptibility::needs_W_notowned(component c,
				realnum *W[NUM_FIELD_COMPONENTS][2]) const {
  FOR_DIRECTIONS(d) if (d != component_direction(c)) {
    component cP = direction_component(c, d);
    if (needs_P(cP, W) && !trivial_sigma[cP][component_direction(c)])
      return true;
  }
  return false;
}

// for Lorentzian susc. the internal data is just a backup of P from
// the previous timestep.
int lorentzian_susceptibility::num_internal_data(
			 realnum *P[NUM_FIELD_COMPONENTS][2],
			 const volume &v) const {
  int num = 0;
  FOR_COMPONENTS(c) DOCMP2 if (P[c][cmp]) num += v.ntot();
  return num;
}

void lorentzian_susceptibility::update_P
       (realnum *P[NUM_FIELD_COMPONENTS][2],
	realnum *W[NUM_FIELD_COMPONENTS][2],
	realnum *W_prev[NUM_FIELD_COMPONENTS][2], 
	double dt, const volume &v, realnum *P_internal_data) const {
  const double omega2pi = 2*pi*omega_0, g2pi = gamma*2*pi;
  const double omega0dtsqr = omega2pi * omega2pi * dt * dt;
  const double gamma1inv = 1 / (1 + g2pi*dt/2), gamma1 = (1 - g2pi*dt/2);
  (void) W_prev;
  
  realnum *P_prev;
  P_prev = P_internal_data;
  FOR_COMPONENTS(c) DOCMP2 if (P[c][cmp]) {
    // FIXME: handle offdiagonal sigma
    const realnum *w = W[c][cmp], *s = sigma[c][component_direction(c)];
    if (w && s) {
      realnum *p = P[c][cmp], *pp = P_prev;
      for (int i = 0; i < v.ntot(); ++i) {
	realnum pcur = p[i];
	p[i] = gamma1inv * (pcur * (2 - omega0dtsqr) 
			    - gamma1 * pp[i] + w[i] * s[i] * omega0dtsqr);
	pp[i] = pcur;
      }
    }
    P_prev += v.ntot();
  }
}

} // namespace meep
