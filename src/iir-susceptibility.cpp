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

// The internal data needs to contain the current polarizations,
// the previous polarizations, and the previous fields.
typedef struct {
  size_t sz_data;
  realnum *P[NUM_FIELD_COMPONENTS][2]; // current polarization value used in step equations
  size_t P_tot;
  realnum *P_prev[NUM_FIELD_COMPONENTS][2]; // previous polarization values used to update
  size_t P_prev_tot;
  realnum *W_prev[NUM_FIELD_COMPONENTS][2]; // previous field values used to update
  size_t W_prev_tot;
  realnum data[1];
} iir_data;

void *iir_susceptibility::new_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2], const grid_volume &gv) const {
  // Count to see how large we need our data structure to be
  int num = 0;
  FOR_COMPONENTS(c) DOCMP2 {
    // for each component, we need 2 directions, ntot index points, 
    // Dz denominator filter taps, and Nz-1 numerator filter taps
    if (needs_P(c, cmp, W)) num += 2 * gv.ntot() * D_z * (N_z-1); 
  }
  size_t sz = sizeof(iir_data) + sizeof(realnum) * (num);
  iir_data *d = (iir_data *)malloc(sz);
  d->sz_data = sz;
  return (void *)d;
}

void iir_susceptibility::init_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2], double dt, const grid_volume &gv, void *data) const {
  (void)dt; // unused
  iir_data *d = (iir_data *)data;

  // initialize everything to 0
  size_t sz_data = d->sz_data;
  memset(d, 0, sz_data);
  d->sz_data = sz_data;

  // Take that contiguous block of data passed by the user and
  // distribute it to the appropriate array locations
  d->P_tot = gv.ntot();
  d->P_prev_tot = gv.ntot() * D_z;
  d->W_prev_tot = gv.ntot() * (N_z -1);
  realnum *P = d->data;
  realnum *P_prev = d->data + d->P_tot;
  realnum *W_prev = P_prev + d->P_prev_tot;

  FOR_COMPONENTS(c) DOCMP2 {
    if (needs_P(c, cmp, W)) {
      for (int dd = X; dd < R; dd++) {
      d->P[c][cmp] = P;
      d->P_prev[c][cmp] = P_prev;
      d->W_prev[c][cmp] = W_prev;
      P += 2 * (d->P_tot + d->P_prev_tot + d->W_prev_tot);
      P_prev += 2 * (d->P_tot + d->P_prev_tot + d->W_prev_tot);
      W_prev += 2 * (d->P_tot + d->P_prev_tot + d->W_prev_tot);
      }
    }
  }
}

void *iir_susceptibility::copy_internal_data(void *data) const {
  iir_data *d = (iir_data *)data;
  if (!d) return 0;
  iir_data *dnew = (iir_data *)malloc(d->sz_data);
  memcpy(dnew, d, d->sz_data);

  dnew->P_tot = d->P_tot;
  dnew->P_prev_tot = d->P_prev_tot;
  dnew->W_prev_tot = d->W_prev_tot;
  realnum *P = d->data;
  realnum *P_prev = d->data + d->P_tot;
  realnum *W_prev = P_prev + d->P_prev_tot;
  FOR_COMPONENTS(c) DOCMP2 {
    if (d->P[c][cmp]) {
      dnew->P[c][cmp] = P;
      dnew->P_prev[c][cmp] = P_prev;
      dnew->W_prev[c][cmp] = W_prev;
      P += 2 * (d->P_tot + d->P_prev_tot + d->W_prev_tot);
      P_prev += 2 * (d->P_tot + d->P_prev_tot + d->W_prev_tot);
      W_prev += 2 * (d->P_tot + d->P_prev_tot + d->W_prev_tot);
    }
  }
  return (void *)dnew;
}

void iir_susceptibility::update_P(realnum *W[NUM_FIELD_COMPONENTS][2],
                                         realnum *W_prev[NUM_FIELD_COMPONENTS][2], double dt,
                                         const grid_volume &gv, void *P_internal_data) const {
  iir_data *d = (iir_data *)P_internal_data;
  (void)W_prev; // unused;
  (void)dt;

  FOR_COMPONENTS(c) DOCMP2 {
    if (d->P[c][cmp]) {
      const realnum *w = W[c][cmp], *s = sigma[c][component_direction(c)];
      if (w && s) {
        if (s != 0) {
          realnum *p = d->P[c][cmp], *pp = d->P_prev[c][cmp], *wp = d->W_prev[c][cmp];

          LOOP_OVER_VOL_OWNED(gv, c, i) {
              // Initialize current polarization by weighting the
              // current field value with the first numerator tap ...
              p[i] = numZ[0] * W[c][cmp][i];
              // ... then continue by looping through the rest of 
              // the numerator taps...
              for (int num_c=0; num_c<N_z-2; num_c++) p[i] += numZ[num_c+1] * wp[num_c][i];
              //... and now loop through the denominator taps,
              // which weight the previous polarizations 
              // (subtract since they must be moved to the RHS)
              for (int den_c=0; den_c<D_z-1; den_c++) p[i] -= denZ[den_c] * wp[den_c][i];
              // shift the cached W field
              for (int num_c=N_z-2; num_c>0; num_c--) wp[num_c][i] = wp[num_c-1][i];
              wp[0][i] = W[c][cmp][i];
              // shift the P field
              for (int den_c=D_z-1; den_c>0; den_c--) pp[den_c][i] = pp[den_c-1][i];
              pp[0][i] = p[i];
          }
        }
      }
    }
  }
}

void iir_susceptibility::subtract_P(field_type ft,
                                           realnum *f_minus_p[NUM_FIELD_COMPONENTS][2],
                                           void *P_internal_data) const {
  iir_data *d = (iir_data *)P_internal_data;
  field_type ft2 = ft == E_stuff ? D_stuff : B_stuff; // for sources etc.
  size_t ntot = d->ntot;
  FOR_FT_COMPONENTS(ft, ec) DOCMP2 {
    if (d->P[ec][cmp]) {
      component dc = field_type_component(ft2, ec);
      if (f_minus_p[dc][cmp]) {
        realnum *p = d->P[ec][cmp];
        realnum *fmp = f_minus_p[dc][cmp];
        for (size_t i = 0; i < ntot; ++i)
          fmp[i] -= p[i];
      }
    }
  }
}

// constructor must take the s domain coefficients and transform
// them to the z domain and store them in the class.
iir_susceptibility::iir_susceptibility(std::vector<double> vec_numS, std::vector<double> vec_denS, double dt) {
    if (vec_numS.size() > vec_denS.size()) abort("error: for iir-filter susceptibilities, the numerator rank cannot be higher than the denominator rank.\n");
    
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

int iir_susceptibility::num_cinternal_notowned_needed(component c,
                                                             void *P_internal_data) const {
  iir_data *d = (iir_data *)P_internal_data;
  return d->P[c][0][0] ? 1 : 0;
}

realnum *iir_susceptibility::cinternal_notowned_ptr(int inotowned, component c, int cmp,
                                                           int n, void *P_internal_data) const {
  iir_data *d = (iir_data *)P_internal_data;
  (void)inotowned; // always = 0
  if (!d || !d->P[c][cmp]) return NULL;
  return d->P[c][cmp] + n;
}

std::complex<double> iir_susceptibility::chi1(double freq, double sigma) {
  std::complex<double> num = 0;
  std::complex<double> den = 0;
  for (int i=0;i<N;i++) num += std::pow(std::complex<double>(0,freq),i) * numS[i];
  for (int i=0;i<D;i++) den += std::pow(std::complex<double>(0,freq),i) * denS[i];
  return sigma * num / den;
}

void iir_susceptibility::dump_params(h5file *h5f, size_t *start) {
  size_t num_params = N_z + D_z;
  size_t params_dims[1] = {num_params};
  double params_data[num_params];
  for (int i=0; i<N_z; i++) params_data[i] = numZ[i];
  for (int i=0; i<D_z; i++) params_data[i+N_z] = denZ[i];
  h5f->write_chunk(1, start, params_dims, params_data);
  *start += num_params;
}

// be sure to clean up the numerator and denominator arrays
iir_susceptibility::~iir_susceptibility() {
  delete[] numZ;
  delete[] denZ;
  delete[] numS;
  delete[] denS;
}

}