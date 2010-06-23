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

/* this file implements multilevel atomic materials for Meep */

#include <stdlib.h>
#include <string.h>
#include "meep.hpp"
#include "meep_internals.hpp"

namespace meep {

typedef realnum *realnumP;

typedef struct {
  size_t sz_data;
  int ntot;
  realnum *GammaInv; // inv(1 + Gamma * dt / 2)
  realnumP *P[NUM_FIELD_COMPONENTS][2]; // P[c][cmp][transition][i]
  realnumP *P_prev[NUM_FIELD_COMPONENTS][2];
  realnum *N; // ntot x L array of centered grid populations N[i*L + level]
  realnum *Ntmp; // temporary length L array of levels, used in updating
  realnum data[1];
} multilevel_data;

multilevel_susceptibility::multilevel_susceptibility(int theL, int theT,
			    const realnum *theGamma,
			    const realnum *theN0,
			    const realnum *thealpha,
			    const realnum *theomega,
			    const realnum *thegamma) {
  L = theL;
  T = theT;
  Gamma = new realnum[L*L];
  memcpy(Gamma, theGamma, sizeof(realnum) * L*L);
  N0 = new realnum[L];
  memcpy(N0, theN0, sizeof(realnum) * L);
  alpha = new realnum[L*T];
  memcpy(alpha, thealpha, sizeof(realnum) * L*T);
  omega = new realnum[T];
  memcpy(omega, theomega, sizeof(realnum) * T);
  gamma = new realnum[T];
  memcpy(gamma, thegamma, sizeof(realnum) * T);
}

void multilevel_susceptibility::subtract_P(field_type ft,
			  realnum *f_minus_p[NUM_FIELD_COMPONENTS][2], 
					   void *P_internal_data) const {
  multilevel_data *d = (multilevel_data *) P_internal_data;
  field_type ft2 = ft == E_stuff ? D_stuff : B_stuff; // for sources etc.
  int ntot = d->ntot;
  for (int t = 0; t < T; ++t) { 
    FOR_FT_COMPONENTS(ft, ec) DOCMP2 if (d->P[ec][cmp][t]) {
      component dc = field_type_component(ft2, ec);
      if (f_minus_p[dc][cmp]) {
	realnum *p = d->P[ec][cmp][t];
	realnum *fmp = f_minus_p[dc][cmp];
	for (int i = 0; i < ntot; ++i) fmp[i] -= p[i];
      }
    }
  }
}

/* U <- 1/U.  U must be real symmetric & positive-definite
   Returns 1 on success, 0 if failure (e.g. matrix singular) */
int sqmatrix_invert(realnum *S, int p)
{

  if (!lapackglue_potrf('U', p, S, p)) return 0;
  if (!lapackglue_potri('U', p, S, p)) return 0;

  int i,j;
  /* Now, copy the upper half onto the lower half of U */
  for (i = 0; i < p; ++i)
    for (j = i + 1; j < p; ++j)
       S[j * p + i] =  S[i * p + j];

  return 1;
}

void *multilevel_susceptibility::new_internal_data(
				    realnum *W[NUM_FIELD_COMPONENTS][2],
				    const grid_volume &gv) const {
  int num = 0; // number of P components
  FOR_COMPONENTS(c) DOCMP2 if (needs_P(c, cmp, W)) num += 2 * gv.ntot();
  size_t sz = sizeof(multilevel_data)
    + sizeof(realnum) * (L*L + L + gv.ntot()*L + num*T - 1);
  multilevel_data *d = (multilevel_data *) malloc(sz);
  memset(d, 0, sz);
  d->sz_data = sz;
  return (void*) d;
}

void multilevel_susceptibility::init_internal_data(
			  realnum *W[NUM_FIELD_COMPONENTS][2],
			  double dt, const grid_volume &gv, void *data) const {
  (void) dt;
  multilevel_data *d = (multilevel_data *) data;
  size_t sz_data = d->sz_data;
  memset(d, 0, sz_data);
  d->sz_data = sz_data;
  int ntot = d->ntot = gv.ntot();
  if (!sqmatrix_invert(Gamma, L)) abort("Gamma matrix singular");
  d->GammaInv = Gamma;
  realnum *P = d->data + L*L;
  realnum *P_prev = d->data + L*L + ntot;
  FOR_COMPONENTS(c) DOCMP2 if (needs_P(c, cmp, W)) {
    d->P[c][cmp] = new realnumP[T];
    d->P_prev[c][cmp] = new realnumP[T];
    for (int t = 0; t < T; ++t) {
      d->P[c][cmp][t] = P;
      d->P_prev[c][cmp][t] = P_prev;
      P += 2*ntot;
      P_prev += 2*ntot;
    }
  }
  d->Ntmp = P_prev;
  d->N = P_prev + L; // the last L*ntot block of the data

  // initial populations
  for (int i = 0; i < ntot; ++i)
    for (int l = 0; l < L; ++l)
      d->N[i*L + l] = N0[l];
}

void multilevel_susceptibility::delete_internal_data(void *data) const {
  if (data) {
    multilevel_data *d = (multilevel_data *) data;
    FOR_COMPONENTS(c) DOCMP2 {
      delete[] d->P[c][cmp];
      delete[] d->P_prev[c][cmp];
    }
    free(data);
  }
}

void *multilevel_susceptibility::copy_internal_data(void *data) const {
  multilevel_data *d = (multilevel_data *) data;
  if (!d) return 0;
  multilevel_data *dnew = (multilevel_data *) malloc(d->sz_data);
  memcpy(dnew, d, d->sz_data);
  int ntot = d->ntot;
  dnew->GammaInv = dnew->data;
  realnum *P = dnew->data + L*L;
  realnum *P_prev = dnew->data + L*L + ntot;
  FOR_COMPONENTS(c) DOCMP2 if (d->P[c][cmp]) {
    dnew->P[c][cmp] = new realnumP[T];
    dnew->P_prev[c][cmp] = new realnumP[T];
    for (int t = 0; t < T; ++t) {
      dnew->P[c][cmp][t] = P;
      dnew->P_prev[c][cmp][t] = P_prev;
      P += 2*ntot;
      P_prev += 2*ntot;
    }
  }
  dnew->Ntmp = P_prev;
  dnew->N = P_prev + L;
  return (void*) dnew;
}

int multilevel_susceptibility::num_cinternal_notowned_needed(component c,
				   void *P_internal_data) const {
  multilevel_data *d = (multilevel_data *) P_internal_data;
  return d->P[c][0] ? T : 0;
}

realnum *multilevel_susceptibility::cinternal_notowned_ptr(
				        int inotowned, component c, int cmp, 
					int n, 
					void *P_internal_data) const {
  multilevel_data *d = (multilevel_data *) P_internal_data;
  if (!d->P[c][cmp] || inotowned < 0 || inotowned >= T) // never true
    return NULL;
  return d->P[c][cmp][inotowned] + n;
}

void multilevel_susceptibility::update_P
       (realnum *W[NUM_FIELD_COMPONENTS][2],
	realnum *W_prev[NUM_FIELD_COMPONENTS][2], 
	double dt, const grid_volume &gv, void *P_internal_data) const {
  multilevel_data *d = (multilevel_data *) P_internal_data;
  double dtinv2 = 0.5 / dt;

  // field directions and offsets for E * dP dot product.
  component cdot[3] = {Dielectric,Dielectric,Dielectric};
  int o1[3], o2[3];
  int idot = 0;
  FOR_COMPONENTS(c) if (d->P[c][0]) {
    if (idot == 3) abort("bug in meep: too many polarization components");
    gv.yee2cent_offsets(c, o1[idot], o2[idot]);
    cdot[idot++] = c;
  }

  // update N from W and P
  realnum *GammaInv = d->GammaInv;
  realnum *Ntmp = d->Ntmp;
  LOOP_OVER_VOL_OWNED(gv, Centered, i) {
    realnum *N = d->N + i*L; // N at current point, to update
    
    // Ntmp = (I - Gamma * dt/2) * N
    for (int l1 = 0; l1 < L; ++l1) {
      Ntmp[l1] = (1.0 - Gamma[l1*L + l1]*dtinv2) * N[l1]; // diagonal term
      for (int l2 = 0; l2 < l1; ++l2) Ntmp[l1] -= Gamma[l1*L+l2]*dtinv2 * N[l2];
      for (int l2 = l1+1; l2 < L; ++l2) Ntmp[l1] -= Gamma[l1*L+l2]*dtinv2 * N[l2];
    }

    // compute E*8 at point i
    double E8[3][2];
    for (idot = 0; idot < 3 && cdot[idot] != Dielectric; ++idot) {
      realnum *w = W[cdot[idot]][0], *wp = W_prev[cdot[idot]][0];
      E8[idot][0] = w[i]+w[i+o1[idot]]+w[i+o2[idot]]+w[i+o1[idot]+o2[idot]]
	+ wp[i]+wp[i+o1[idot]]+wp[i+o2[idot]]+wp[i+o1[idot]+o2[idot]];
      if (W[cdot[idot]][1]) {
	w = W[cdot[idot]][1]; wp = W_prev[cdot[idot]][1];
	E8[idot][1] = w[i]+w[i+o1[idot]]+w[i+o2[idot]]+w[i+o1[idot]+o2[idot]]
	  + wp[i]+wp[i+o1[idot]]+wp[i+o2[idot]]+wp[i+o1[idot]+o2[idot]];
      }
      else
	E8[idot][1] = 0;
    }

    // Ntmp = Ntmp + alpha * E * dP
    for (int t = 0; t < T; ++t) {
      // compute 32 * E * dP at point i
      double EdP32 = 0;
      for (idot = 0; idot < 3 && cdot[idot] != Dielectric; ++idot) {
	realnum *p = d->P[cdot[idot]][0][t], *pp = d->P_prev[cdot[idot]][0][t];
	realnum dP = p[i]+p[i+o1[idot]]+p[i+o2[idot]]+p[i+o1[idot]+o2[idot]]
	  - (pp[i]+pp[i+o1[idot]]+pp[i+o2[idot]]+pp[i+o1[idot]+o2[idot]]);
	EdP32 += dP * E8[idot][0];
	if (d->P[cdot[idot]][1][t]) {
	  p = d->P[cdot[idot]][1][t]; pp = d->P_prev[cdot[idot]][1][t];
	  dP = p[i]+p[i+o1[idot]]+p[i+o2[idot]]+p[i+o1[idot]+o2[idot]]
	    + (pp[i]+pp[i+o1[idot]]+pp[i+o2[idot]]+pp[i+o1[idot]+o2[idot]]);
	  EdP32 += dP * E8[idot][1];
	}
      }

      for (int l = 0; l < L; ++l) Ntmp[l] += 0.03125 * alpha[l*T + t] * EdP32;
    }

    // N = GammaInv * Ntmp
    for (int l1 = 0; l1 < L; ++l1) {
      N[l1] = 0;
      for (int l2 = 0; l2 < L; ++l2) N[l1] += GammaInv[l1*L+l2] * Ntmp[l2];
    }
  }

  // each P is updated as a damped harmonic oscillator
  for (int t = 0; t < T; ++t) {
    const double omega2pi = 2*pi*omega[t], g2pi = gamma[t]*2*pi;
    const double omega0dtsqr = omega2pi * omega2pi * dt * dt;
    const double gamma1inv = 1 / (1 + g2pi*dt/2), gamma1 = (1 - g2pi*dt/2);

    // figure out which levels this transition couples
    int lp = -1, lm = -1;
    for (int l = 0; l < L; ++l) {
      if (alpha[l*T + t] > 0) lp = l;
      if (alpha[l*T + t] < 0) lm = l;
    }
    if (lp < 0 || lm < 0) abort("invalid alpha array for transition %d", t);

    FOR_COMPONENTS(c) DOCMP2 if (d->P[c][cmp][t]) {
      const realnum *w = W[c][cmp], *s = sigma[c][component_direction(c)];
      if (w && s) {
	realnum *p = d->P[c][cmp][t], *pp = d->P_prev[c][cmp][t];

	int o1, o2;
	gv.cent2yee_offsets(c, o1, o2);
	o1 *= L; o2 *= L;
	const realnum *N = d->N;
	
	// directions/strides for offdiagonal terms, similar to update_eh
	const direction d = component_direction(c);
	direction d1 = cycle_direction(gv.dim, d, 1);
	component c1 = direction_component(c, d1);
	const realnum *w1 = W[c1][cmp];
	const realnum *s1 = w1 ? sigma[c][d1] : NULL;
	direction d2 = cycle_direction(gv.dim, d, 2);
	component c2 = direction_component(c, d2);
	const realnum *w2 = W[c2][cmp];
	const realnum *s2 = w2 ? sigma[c][d2] : NULL;
	
	if (s1 || s2) {
	  abort("nondiagonal saturable gain is not yet supported");
	}
	else { // isotropic
	  LOOP_OVER_VOL_OWNED(gv, c, i) {
	    realnum pcur = p[i];
	    const realnum *Ni = N + i*L;
	    // dNi is population inversion for this transition
	    double dNi = -0.25 * (Ni[lp]+Ni[lp+o1]+Ni[lp+o2]+Ni[lp+o1+o2]
				  -Ni[lm]-Ni[lm+o1]-Ni[lm+o2]-Ni[lm+o1+o2]);
	    p[i] = gamma1inv * (pcur * (2 - omega0dtsqr) 
				- gamma1 * pp[i] 
				+ omega0dtsqr * (s[i] * w[i])) * dNi;
	    pp[i] = pcur;
	  }
	}
      }
    }
  }
}

}
