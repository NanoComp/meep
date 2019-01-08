/* Copyright (C) 2005-2019 Massachusetts Institute of Technology.
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

#include <stdlib.h>
#include <string.h>
#include "meep.hpp"
#include "meep_internals.hpp"

using namespace std;

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

void susceptibility::delete_internal_data(void *data) const {
  free(data);
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
  FOR_DIRECTIONS(d) {
    if (!trivial_sigma[c][d] && W[direction_component(c, d)][cmp]) return true;
  }
  return false;
}

/* return whether we need the notowned parts of the W field --
   by default, this is only the case if sigma has offdiagonal components
   coupling P to W.   (See needs_P: again, this true if the notowned
   W is needed in *any* chunk.) */
bool susceptibility::needs_W_notowned(component c,
				realnum *W[NUM_FIELD_COMPONENTS][2]) const {
  FOR_DIRECTIONS(d) { if (d != component_direction(c)) {
    component cP = direction_component(c, d);
    if (needs_P(cP, 0, W) && !trivial_sigma[cP][component_direction(c)])
      return true;
  }
  }
  return false;
}

typedef struct {
  size_t sz_data;
  size_t ntot;
  realnum *P[NUM_FIELD_COMPONENTS][2];
  realnum *P_prev[NUM_FIELD_COMPONENTS][2];
  realnum data[1];
} lorentzian_data;

// for Lorentzian susc. the internal data is just a backup of P from
// the previous timestep.
void *lorentzian_susceptibility::new_internal_data(
			 realnum *W[NUM_FIELD_COMPONENTS][2],
			 const grid_volume &gv) const {
  int num = 0;
  FOR_COMPONENTS(c) DOCMP2 { if (needs_P(c, cmp, W)) num += 2 * gv.ntot(); }
  size_t sz = sizeof(lorentzian_data) + sizeof(realnum) * (num - 1);
  lorentzian_data *d = (lorentzian_data *) malloc(sz);
  d->sz_data = sz;
  return (void*) d;
}

void lorentzian_susceptibility::init_internal_data(
			  realnum *W[NUM_FIELD_COMPONENTS][2],
			  double dt, const grid_volume &gv, void *data) const {
  (void) dt; // unused
  lorentzian_data *d = (lorentzian_data *) data;
  size_t sz_data = d->sz_data;
  memset(d, 0, sz_data);
  d->sz_data = sz_data;
  size_t ntot = d->ntot = gv.ntot();
  realnum *P = d->data;
  realnum *P_prev = d->data + ntot;
  FOR_COMPONENTS(c) DOCMP2 { if (needs_P(c, cmp, W)) {
    d->P[c][cmp] = P;
    d->P_prev[c][cmp] = P_prev;
    P += 2*ntot;
    P_prev += 2*ntot;
  }
  }
}

void *lorentzian_susceptibility::copy_internal_data(void *data) const {
  lorentzian_data *d = (lorentzian_data *) data;
  if (!d) return 0;
  lorentzian_data *dnew = (lorentzian_data *) malloc(d->sz_data);
  memcpy(dnew, d, d->sz_data);
  size_t ntot = d->ntot;
  realnum *P = dnew->data;
  realnum *P_prev = dnew->data + ntot;
  FOR_COMPONENTS(c) DOCMP2 { if (d->P[c][cmp]) {
    dnew->P[c][cmp] = P;
    dnew->P_prev[c][cmp] = P_prev;
    P += 2*ntot;
    P_prev += 2*ntot;
  }
  }
  return (void*) dnew;
}

/* Return true if the discretized Lorentzian ODE is intrinsically unstable,
   i.e. if it corresponds to a filter with a pole z outside the unit circle.
   Note that the pole satisfies the quadratic equation:
            (z + 1/z - 2)/dt^2 + g*(z - 1/z)/(2*dt) + w^2 = 0
   where w = 2*pi*omega_0 and g = 2*pi*gamma.   It is just a little
   algebra from this to get the condition for a root with |z| > 1.

   FIXME: this test seems to be too conservative (issue #12) */
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
	double dt, const grid_volume &gv, void *P_internal_data) const {
  lorentzian_data *d = (lorentzian_data *) P_internal_data;
  const double omega2pi = 2*pi*omega_0, g2pi = gamma*2*pi;
  const double omega0dtsqr = omega2pi * omega2pi * dt * dt;
  const double gamma1inv = 1 / (1 + g2pi*dt/2), gamma1 = (1 - g2pi*dt/2);
  const double omega0dtsqr_denom = no_omega_0_denominator ? 0 : omega0dtsqr;
  (void) W_prev; // unused;

  // TODO: add back lorentzian_unstable(omega_0, gamma, dt) if we can improve the stability test

  FOR_COMPONENTS(c) DOCMP2 { if (d->P[c][cmp]) {
    const realnum *w = W[c][cmp], *s = sigma[c][component_direction(c)];
    if (w && s) {
      realnum *p = d->P[c][cmp], *pp = d->P_prev[c][cmp];

      // directions/strides for offdiagonal terms, similar to update_eh
      const direction d = component_direction(c);
      const ptrdiff_t is = gv.stride(d) * (is_magnetic(c) ? -1 : +1);
      direction d1 = cycle_direction(gv.dim, d, 1);
      component c1 = direction_component(c, d1);
      ptrdiff_t is1 = gv.stride(d1) * (is_magnetic(c) ? -1 : +1);
      const realnum *w1 = W[c1][cmp];
      const realnum *s1 = w1 ? sigma[c][d1] : NULL;
      direction d2 = cycle_direction(gv.dim, d, 2);
      component c2 = direction_component(c, d2);
      ptrdiff_t is2 = gv.stride(d2) * (is_magnetic(c) ? -1 : +1);
      const realnum *w2 = W[c2][cmp];
      const realnum *s2 = w2 ? sigma[c][d2] : NULL;

      if (s2 && !s1) { // make s1 the non-NULL one if possible
	SWAP(direction, d1, d2);
	SWAP(component, c1, c2);
	SWAP(ptrdiff_t, is1, is2);
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
  }
  }
}

void lorentzian_susceptibility::subtract_P(field_type ft,
					   realnum *f_minus_p[NUM_FIELD_COMPONENTS][2],
					   void *P_internal_data) const {
  lorentzian_data *d = (lorentzian_data *) P_internal_data;
  field_type ft2 = ft == E_stuff ? D_stuff : B_stuff; // for sources etc.
  size_t ntot = d->ntot;
  FOR_FT_COMPONENTS(ft, ec) DOCMP2 { if (d->P[ec][cmp]) {
    component dc = field_type_component(ft2, ec);
    if (f_minus_p[dc][cmp]) {
      realnum *p = d->P[ec][cmp];
      realnum *fmp = f_minus_p[dc][cmp];
      for (size_t i = 0; i < ntot; ++i) fmp[i] -= p[i];
    }
  }
  }
}

int lorentzian_susceptibility::num_cinternal_notowned_needed(component c,
				   void *P_internal_data) const {
  lorentzian_data *d = (lorentzian_data *) P_internal_data;
  return d->P[c][0] ? 1 : 0;
}

realnum *lorentzian_susceptibility::cinternal_notowned_ptr(
				        int inotowned, component c, int cmp,
					int n,
					void *P_internal_data) const {
  lorentzian_data *d = (lorentzian_data *) P_internal_data;
  (void) inotowned; // always = 0
  if (!d || !d->P[c][cmp])
    return NULL;
  return d->P[c][cmp] + n;
}

void lorentzian_susceptibility::dump_params(h5file *h5f, size_t *start) {
  size_t num_params = 5;
  size_t params_dims[1] = {num_params};
  double params_data[] = {4, (double)get_id(), omega_0, gamma, (double)no_omega_0_denominator};
  h5f->write_chunk(1, start, params_dims, params_data);
  *start += num_params;
}

void noisy_lorentzian_susceptibility::update_P
       (realnum *W[NUM_FIELD_COMPONENTS][2],
	realnum *W_prev[NUM_FIELD_COMPONENTS][2],
	double dt, const grid_volume &gv, void *P_internal_data) const {
  lorentzian_susceptibility::update_P(W, W_prev, dt, gv, P_internal_data);
  lorentzian_data *d = (lorentzian_data *) P_internal_data;

  const double g2pi = gamma*2*pi;
  const double w2pi = omega_0*2*pi;
  const double amp = w2pi * noise_amp * sqrt(g2pi) * dt*dt / (1 + g2pi*dt/2);
  /* for uniform random numbers in [-amp,amp] below, multiply amp by sqrt(3) */

  FOR_COMPONENTS(c) DOCMP2 { if (d->P[c][cmp]) {
    const realnum *s = sigma[c][component_direction(c)];
    if (s) {
      realnum *p = d->P[c][cmp];
      LOOP_OVER_VOL_OWNED(gv, c, i) {
        p[i] += gaussian_random(0, amp * sqrt(s[i]));
      }
      // for uniform random numbers, use uniform_random(-1,1) * amp * sqrt(s[i])
      // for gaussian random numbers, use gaussian_random(0, amp * sqrt(s[i]))
    }
  }
  }
}

void noisy_lorentzian_susceptibility::dump_params(h5file *h5f, size_t *start) {
  size_t num_params = 6;
  size_t params_dims[1] = {num_params};
  double params_data[] = {5, (double)get_id(), noise_amp, omega_0, gamma, (double)no_omega_0_denominator};
  h5f->write_chunk(1, start, params_dims, params_data);
  *start += num_params;
}
} // namespace meep
