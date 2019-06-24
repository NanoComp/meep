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

#include <stdlib.h>
#include <string.h>
#include "meep.hpp"
#include "meep_internals.hpp"
#include "material_data.hpp"

using namespace std;

namespace meep {

// --------------------------------------------------------------- //
// S to Z transform routines
// --------------------------------------------------------------- //
/* The following routines are use to convert the IIR transfer function
 * in the s domain to the z domain.
 */

double factorial(int i) {
  double factorial = 1;
  for (i; i > 0; i--) factorial *= i;
  return factorial;
}

double comb(int n, int k) {
  return factorial(n) / (factorial(k) * factorial(n-k));
}


// Uses Tustins method to discretize the transfer function. S = (2/T)*(1-z^-1)/(1+z^-1)
// https://github.com/scipy/scipy/blob/v1.3.0/scipy/signal/filter_design.py#L1965-L2043
void tustins_method(
    double* vec_numS, 
    int N, 
    double* vec_denS, 
    int D, 
    double* vec_numZ,
    int N_z, 
    double* vec_denZ,
    int D_z, 
    double T
    ) {

  int M = D-1; //Denominator should always have highest rank
  int Np = M;
  int Dp = M;
  
  std::vector<double> bprime(Np + 1, 0.0);
  std::vector<double> aprime(Dp + 1, 0.0);
    
  int j, i, k, l;
  double val;
  for (j=0;j<=Np;j++) {
    val = 0.0;
    for(i=0;i<N;i++) {
      for (k=0;k<=i;k++) {
        for (l=0;l<=M-i;l++) {
          if (k+l == j) {
            val += (comb(i,k) * comb(M-i,l) * vec_numS[N-i-1] * pow(2.0/T,i) * pow(-1,k));
          }
        }
      }
    }
    bprime[j] = val;
  }
  for (j=0;j<=Dp;j++) {
    val = 0.0;
    for(i=0;i<D;i++) {
      for (k=0;k<=i;k++) {
        for (l=0;l<=M-i;l++) {
          if (k+l == j) {
            val += (comb(i,k) * comb(M-i,l) * vec_denS[D-i-1] * pow(2.0/T,i) * pow(-1,k));
          }
        }
      }
    }
    aprime[j] = val;
  }

  // Normalize and copy elements over
  double norm = aprime[0];
  for (i=0; i<=Np; i++) vec_numZ[i] = bprime[i]/norm;
  for (i=0; i<=Dp; i++) vec_denZ[i] = aprime[i]/norm;        
}

// --------------------------------------------------------------- //
// IIR transfer function susceptibility
// --------------------------------------------------------------- //

typedef struct {
  size_t sz_data;
  size_t ntot;
  realnum *P;
  realnum *W;
  realnum data[1];
} iir_data;

// The internal data is just a backup of P AND W from the previous timestep(s).
void *iir_susceptibility::new_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2], const grid_volume &gv) const {
  int num = 0;
  FOR_COMPONENTS(c) DOCMP2 {
    if (needs_P(c, cmp, W)) num += 2 * gv.ntot() * 2 * Dz; // for each component, we need 2 directions, ntot index points, and 2*Dz filter taps (1 for P, 1 for W)
  }
  size_t sz = sizeof(iir_susceptibility) + sizeof(realnum) * (num - 1);
  iir_data *d = (iir_data *)malloc(sz);
  d->sz_data = sz;
  return (void *)d;
}

void iir_susceptibility::init_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2], double dt, const grid_volume &gv, void *data) const {
  (void)dt; // unused
  iir_data *d = (iir_data *)data;
  size_t sz_data = d->sz_data;
  memset(d, 0, sz_data);
  d->sz_data = sz_data;
  size_t ntot = d->ntot = gv.ntot() * Dz;
  realnum *P = d->data;
  realnum *P_prev = d->data + ntot;
  realnum *W_prev = P_prev + ntot;
  FOR_COMPONENTS(c) DOCMP2 {
    if (needs_P(c, cmp, W)) {
      for (int dd = X; dd < R; dd++) {
      d->P[c][cmp] = P;
      d->P_prev[c][cmp] = P_prev;
      d->W_prev[c][cmp] = W_prev;
      P += 2 * ntot;
      P_prev += 2 * ntot;
      }
    }
  }
}

void iir_susceptibility::update_P(realnum *W[NUM_FIELD_COMPONENTS][2],
                                         realnum *W_prev[NUM_FIELD_COMPONENTS][2], double dt,
                                         const grid_volume &gv, void *P_internal_data) const {
  iir_data *d = (iir_data *)P_internal_data;
  (void)W_prev; // unused;


  FOR_COMPONENTS(c) DOCMP2 {
    if (d->P[c][cmp]) {
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
            // s[i] != 0 check is a bit of a hack to work around
            // some instabilities that occur near the boundaries
            // of materials; see PR #666
            if (s[i] != 0) {
              realnum pcur = p[i];
              p[i] = gamma1inv * (pcur * (2 - omega0dtsqr_denom) - gamma1 * pp[i] +
                                  omega0dtsqr * (s[i] * w[i] + OFFDIAG(s1, w1, is1, is) +
                                                 OFFDIAG(s2, w2, is2, is)));
              pp[i] = pcur;
            }
          }
        } else if (s1) { // 2x2 anisotropic
          LOOP_OVER_VOL_OWNED(gv, c, i) {
            if (s[i] != 0) { // see above
              realnum pcur = p[i];
              p[i] = gamma1inv * (pcur * (2 - omega0dtsqr_denom) - gamma1 * pp[i] +
                                  omega0dtsqr * (s[i] * w[i] + OFFDIAG(s1, w1, is1, is)));
              pp[i] = pcur;
            }
          }
        } else { // isotropic
          LOOP_OVER_VOL_OWNED(gv, c, i) {
            realnum pcur = p[i];
            // Assign current polarization by looping through first 
            // the numerator coefficients and then the denominator
            // coefficients.
            data->p[i] = numZ[0] * W[c2][cmp][i];
            for (int num_c=0; num_c<N_z-1; num_c++) data->p[i] += numZ[num_c+1] * data->W_prev[num_c];
            for (int num_c=0; num_c<N_z-1; num_c++) data->p[i] += numZ[num_c+1] * data->W_prev[num_c];

            // shift the W field
            // shift the P field
            p[i] = gamma1inv *
                   (pcur * (2 - omega0dtsqr_denom) - gamma1 * pp[i] + omega0dtsqr * (s[i] * w[i]));
            pp[i] = pcur;
          }
        }
      }
    }
  }
}

// constructor must take the s domain coefficients and transform
// them to the z domain and store them in the class.
iir_susceptibility::iir_susceptibility(std::vector<double> vec_numS, std::vector<double> vec_denS, double dt) {
    if (vec_numS.size() > vec_denS.size()) abort("error: the numerator rank cannot be higher than the denominator rank.\n");
    
    // copy over the S domain numerator and denominator
    N = vec_numS.size();
    D = vec_denS.size();
    N_z = D;
    D_z = D; // both numerator and denominator of Z domain are same rank.
    numS = new double[N];
    std::copy(vec_numS.begin(), vec_numS.end(), numS);
    denS = new double[D];
    std::copy(vec_denS.begin(), vec_denS.end(), denS);
    numZ = new double[N_z]; 
    denZ = new double[D_z];

    T = dt;

    // calculate the Z domain numerator and denominator
    tustins_method(numS, N, denS, D, numZ, N_z, denZ,D_z, T);
}

// be sure to clean up the numerator and denominator arrays
iir_susceptibility::~iir_susceptibility() {
  delete[] numZ;
  delete[] denZ;
  delete[] numS;
  delete[] denS;
}

}