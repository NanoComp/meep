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

/* This file implements dispersive materials for Meep via a
   polarization P = \chi(\omega) W, where W is e.g. E or H.  Each
   subclass of the susceptibility class should implement a different
   type of \chi(\omega).  The subclass knows how to timestep P given W
   at the current (and possibly previous) timestep, and any additional
   internal data that needs to be allocated along with P.

   Each \chi(\omega) is spatially multiplied by a (scalar) sigma
   array.  The meep::fields class is responsible for allocating P and
   sigma and passing them to susceptibility::update_P. */


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

/* Return whether or not we need to allocate P[c][cmp].  (We don't need to
   allocate P[c] if we can be sure it will be zero.)

   We are a bit wasteful because if sigma is nontrivial in *any* chunk,
   we allocate the corresponding P on *every* owned chunk.  This greatly
   simplifies communication in boundaries.cpp, because we can be sure that
   one chunk has a P then any chunk it borders has the same P, so we don't
   have to worry about communicating with something that doesn't exist.
   TODO: reduce memory usage (bookkeeping seem much harder, though).
*/
bool susceptibility::needs_P(component c, int cmp,
			     realnum *W[NUM_FIELD_COMPONENTS][2])
  const {
  if (!is_electric(c) && !is_magnetic(c)) return false;
  FOR_DIRECTIONS(d)
    if (!trivial_sigma[c][d] && W[direction_component(c, d)][cmp]) return true;
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
    if (needs_P(cP, 0, W) && !trivial_sigma[cP][component_direction(c)])
      return true;
  }
  return false;
}

// for Lorentzian susc. the internal data is just a backup of P from
// the previous timestep.
int lorentzian_susceptibility::num_internal_data(
			 realnum *W[NUM_FIELD_COMPONENTS][2],
			 const grid_volume &gv) const {
  int num = 0;
  FOR_COMPONENTS(c) DOCMP2 if (needs_P(c, cmp, W)) num += 2 * gv.ntot();
  return num;
}

/* Return true if the discretized Lorentzian ODE is intrinsically unstable,
   i.e. if it corresponds to a filter with a pole z outside the unit circle.
   Note that the pole satisfies the quadratic equation:
            (z + 1/z - 2)/dt^2 + g*(z - 1/z)/(2*dt) + w^2 = 0
   where w = 2*pi*omega_0 and g = 2*pi*gamma.   It is just a little
   algebra from this to get the condition for a root with |z| > 1. */
static bool lorentzian_unstable(double omega_0, double gamma, double dt) {
  double w = 2*pi*omega_0, g = 2*pi*gamma;
  double g2 = g*dt/2, w2 = (w*dt)*(w*dt);
  double b = (1 - w2/2) / (1 + g2), c = (1 - g2) / (1 + g2);
  return b*b > c && 2*b*b - c + 2*fabs(b)*sqrt(b*b - c) > 1;
}

#define SWAP(t,a,b) { t SWAP_temp = a; a = b; b = SWAP_temp; }

  // stable averaging of offdiagonal components
#define OFFDIAG(u,g,sx,s) (0.25 * ((g[i]+g[i-sx])*u[i]		\
		   	         + (g[i+s]+g[(i+s)-sx])*u[i+s]))

void lorentzian_susceptibility::update_P
       (realnum *W[NUM_FIELD_COMPONENTS][2],
	realnum *W_prev[NUM_FIELD_COMPONENTS][2], 
	double dt, const grid_volume &gv, realnum *P_internal_data) const {
  const double omega2pi = 2*pi*omega_0, g2pi = gamma*2*pi;
  const double omega0dtsqr = omega2pi * omega2pi * dt * dt;
  const double gamma1inv = 1 / (1 + g2pi*dt/2), gamma1 = (1 - g2pi*dt/2);
  const double omega0dtsqr_denom = no_omega_0_denominator ? 0 : omega0dtsqr;
  (void) W_prev; // unused;

  if (!no_omega_0_denominator && gamma >= 0
      && lorentzian_unstable(omega_0, gamma, dt))
    abort("Lorentzian pole at too high a frequency %g for stability with dt = %g: reduce the Courant factor, increase the resolution, or use a different dielectric model\n", omega_0, dt);

  realnum *P = P_internal_data;
  realnum *P_prev = P_internal_data + num_internal_data(W, gv) / 2;
  FOR_COMPONENTS(c) DOCMP2 if (needs_P(c, cmp, W)) {
    const realnum *w = W[c][cmp], *s = sigma[c][component_direction(c)];
    if (w && s) {
      realnum *p = P, *pp = P_prev;

      // directions/strides for offdiagonal terms, similar to update_eh
      const direction d = component_direction(c);
      const int is = gv.stride(d) * (is_magnetic(c) ? -1 : +1);
      direction d1 = cycle_direction(gv.dim, d, 1);
      component c1 = direction_component(c, d1);
      int is1 = gv.stride(d1) * (is_magnetic(c) ? -1 : +1);
      const realnum *w1 = W[c1][cmp];
      const realnum *s1 = w1 ? sigma[c][d1] : NULL;
      direction d2 = cycle_direction(gv.dim, d, 2);
      component c2 = direction_component(c, d2);
      int is2 = gv.stride(d2) * (is_magnetic(c) ? -1 : +1);
      const realnum *w2 = W[c2][cmp];
      const realnum *s2 = w2 ? sigma[c][d2] : NULL;

      if (s2 && !s1) { // make s1 the non-NULL one if possible
	SWAP(direction, d1, d2);
	SWAP(component, c1, c2);
	SWAP(int, is1, is2);
	SWAP(const realnum *, w1, w2);
	SWAP(const realnum *, s1, s2);
      }
      if (s1 && s2) { // 3x3 anisotropic
	LOOP_OVER_VOL_OWNED(gv, c, i) {
	  realnum pcur = p[i];
	  p[i] = gamma1inv * (pcur * (2 - omega0dtsqr_denom) 
			      - gamma1 * pp[i] 
			      + omega0dtsqr * (s[i] * w[i]
					       + OFFDIAG(s1,w1,is1,is)
					       + OFFDIAG(s2,w2,is2,is)));
	  pp[i] = pcur;
	}
      }
      else if (s1) { // 2x2 anisotropic
	LOOP_OVER_VOL_OWNED(gv, c, i) {
	  realnum pcur = p[i];
	  p[i] = gamma1inv * (pcur * (2 - omega0dtsqr_denom) 
			      - gamma1 * pp[i] 
			      + omega0dtsqr * (s[i] * w[i]
					       + OFFDIAG(s1,w1,is1,is)));
	  pp[i] = pcur;
	}
      }
      else { // isotropic
	LOOP_OVER_VOL_OWNED(gv, c, i) {
	  realnum pcur = p[i];
	  p[i] = gamma1inv * (pcur * (2 - omega0dtsqr_denom) 
			      - gamma1 * pp[i] 
			      + omega0dtsqr * (s[i] * w[i]));
	  pp[i] = pcur;
	}
      }
    }
    P += gv.ntot();
    P_prev += gv.ntot();
  }
}

void lorentzian_susceptibility::subtract_P(field_type ft,
					   realnum *f[NUM_FIELD_COMPONENTS][2],
					   realnum *f_minus_p[NUM_FIELD_COMPONENTS][2], 
					   const grid_volume &gv,
					   realnum *P_internal_data) const {
  field_type ft2 = ft == E_stuff ? D_stuff : B_stuff; // for sources etc.
  realnum *P = P_internal_data;
  int ntot = gv.ntot();
  /* note: don't use FOR_FT_COMPONENTS, use FOR_COMPONENTS to ensure
     that we step through the P array in exactly the same way as for
     update_P */
  FOR_COMPONENTS(ec) DOCMP2 if (needs_P(ec, cmp, f)) {
    component dc = field_type_component(ft2, ec);
    if (type(ec) == ft && f_minus_p[dc][cmp]) {
      realnum *fmp = f_minus_p[dc][cmp];
      for (int i = 0; i < ntot; ++i) fmp[i] -= P[i];
    }
    P += ntot;
  }
}

int lorentzian_susceptibility::num_cinternal_notowned_needed(component c,
				   realnum *W[NUM_FIELD_COMPONENTS][2]) const {
  return needs_P(c, 0, W) ? 1 : 0;
}

int lorentzian_susceptibility::cinternal_notowned_offset(
				        int inotowned, component c0, int cmp0, 
					int n, 
					realnum *W[NUM_FIELD_COMPONENTS][2],
					const grid_volume &gv) const {
  (void) inotowned; // always = 0
  int offset = n;
  FOR_COMPONENTS(c) DOCMP2 if (needs_P(c, cmp, W)) {
    if (c0 == c && cmp0 == cmp) return offset;
    offset += gv.ntot();
  }
  abort("bug: notowned_offset called for unallocated component");
  return offset;
}


void noisy_lorentzian_susceptibility::update_P
       (realnum *W[NUM_FIELD_COMPONENTS][2],
	realnum *W_prev[NUM_FIELD_COMPONENTS][2], 
	double dt, const grid_volume &gv, realnum *P_internal_data) const {
  lorentzian_susceptibility::update_P(W, W_prev, dt, gv, P_internal_data);

  const double g2pi = gamma*2*pi;
  const double amp = noise_amp * sqrt(g2pi) * dt*dt / (1 + g2pi*dt/2);
  /* for uniform random numbers in [-amp,amp] below, multiply amp by sqrt(3) */

  realnum *P = P_internal_data;
  FOR_COMPONENTS(c) DOCMP2 if (needs_P(c, cmp, W)) {
    const realnum *s = sigma[c][component_direction(c)];
    if (s) {
      realnum *p = P;
      LOOP_OVER_VOL_OWNED(gv, c, i)
	p[i] += gaussian_random(0, amp * sqrt(s[i]));
      // for uniform random numbers, use uniform_random(-amp * sqrt(s[i]), +amp * sqrt(s[i]))
    }
    P += gv.ntot();
  }
}

} // namespace meep
