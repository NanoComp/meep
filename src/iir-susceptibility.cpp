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

 // generalized way to convert a polynomial in the s domain to the z 
// domain. Based on the paper:
// Scott, Dan M. "A simplified method for the bilinear sz transformation." IEEE Transactions on Education 37.3 (1994): 289-292.
// https://digitalcommons.unl.edu/cgi/viewcontent.cgi?referer=https://www.google.com/&httpsredir=1&article=1085&context=imsefacpub
void polynomial_transform(double *r, int N, double alpha, double beta, double delta, double gamma) {
    int k, j, n=N-1;
    double sum, beta_powers[n], gamma_powers[n];
    master_printf("entered PT -------------- \n");
    beta_powers[0] = gamma_powers[0] = 1;
    for (j=1;j<=n;j++)
    {
        beta_powers[j] = beta*beta_powers[j- 1];
        gamma_powers[j] = gamma*gamma_powers[j- 1];
    }

    for (k=0; k<n; k++)
    {
        sum = 0.0;
        for (j=0; j<=n-k; j++) sum += beta_powers[j] * gamma_powers[n-k-j] * r[j];
        for (j=0; j<n-k; j++) r[j] = ((n-k-j)*delta*r[j]+(j+1)*alpha*r[j+1])/(k+1);
        r[n-k] = sum;
    }
}
// Uses backward difference method to discretize the transfer function. S = (1-z^-1)/T
void backward_difference(
    realnum* vec_numS, 
    int N, 
    realnum* vec_denS, 
    int D, 
    realnum* vec_numZ,
    int N_z, 
    realnum* vec_denZ,
    int D_z, 
    realnum T) {
  master_printf("entered BD -------------- \n");
  double alpha = 1.0, beta = -1.0, delta = T, gamma = 0;

  // Process numerator
  for (int i=0; i< N; i++) vec_numZ[i] = vec_numS[i];
  polynomial_transform(vec_numZ,N_z,alpha,beta,delta,gamma);
  // Process denominator
  for (int i=0; i< D; i++) vec_denZ[i] = vec_denS[i];
  polynomial_transform(vec_denZ,D_z,alpha,beta,delta,gamma);
  // Normalize
  double norm = vec_denZ[0];
  for (int i=0; i<=N_z; i++) vec_numZ[i] = vec_numZ[i]/norm;
  for (int i=0; i<=D_z; i++) vec_denZ[i] = vec_denZ[i]/norm; 
  }
realnum factorial(int i) {
  realnum factorial = 1;
  for (int k=i; k > 0; k--) factorial *= k;
  return factorial;
}

realnum comb(int n, int k) {
  return factorial(n) / (factorial(k) * factorial(n-k));
}

// Uses Tustins method to discretize the transfer function. S = (2/T)*(1-z^-1)/(1+z^-1)
// https://github.com/scipy/scipy/blob/v1.3.0/scipy/signal/filter_design.py#L1965-L2043
void tustins_method(
    realnum* vec_numS, 
    int N, 
    realnum* vec_denS, 
    int D, 
    realnum* vec_numZ,
    int N_z, 
    realnum* vec_denZ,
    int D_z, 
    realnum T
    ) {
  (void)N_z;
  (void)D_z;

  int M = D-1; //Denominator should always have highest rank
  int Np = M;
  int Dp = M;
  
  std::vector<realnum> bprime(Np + 1, 0.0);
  std::vector<realnum> aprime(Dp + 1, 0.0);
    
  int j, i, k, l;
  realnum val;
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
  realnum norm = aprime[0];
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
  size_t ntot;
  realnum *P[NUM_FIELD_COMPONENTS][2]; // current polarization value used in step equations
  realnum *P_prev[NUM_FIELD_COMPONENTS][2]; // previous polarization values used to update
  realnum *W_prev[NUM_FIELD_COMPONENTS][2]; // previous field values used to update
  realnum data[1];

} iir_data;

void *iir_susceptibility::new_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2], const grid_volume &gv) const {
  size_t ntot = gv.ntot();
  // Count to see how large we need our data structure to be
  int num = 0;
  FOR_COMPONENTS(c) DOCMP2 {
    if (needs_P(c, cmp, W)) num += ntot + ntot * (D_z-1) + ntot * (N_z-1); 
  }
  size_t sz = sizeof(iir_data) + sizeof(realnum) * (num);
  iir_data *d = (iir_data *)malloc(sz);
  d->sz_data = sz;
  d->ntot = ntot;
  
  return (void *)d;
}

void iir_susceptibility::init_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2], double dt, const grid_volume &gv, void *data) const {
  (void)dt; // unused
  (void)gv;
  iir_data *d = (iir_data *)data;

  // initialize everything to 0
  size_t sz_data = d->sz_data;
  size_t ntot = d->ntot;
  memset(d, 0, sz_data);
  d->sz_data = sz_data;
  d->ntot = ntot;

  // Take that contiguous block of data passed by the user and
  // distribute it to the appropriate array locations
  realnum *P = d->data;

  FOR_COMPONENTS(c) DOCMP2 {
    if (needs_P(c, cmp, W)) {
      d->P[c][cmp] = P;
      P += ntot;
      d->P_prev[c][cmp] = P;
      P += ntot * (D_z-1);
      d->W_prev[c][cmp] = P;
      P += ntot * (N_z-1);
    }
  }
}

void *iir_susceptibility::copy_internal_data(void *data) const {
  iir_data *d = (iir_data *)data;
  if (!d) return 0;
  iir_data *dnew = (iir_data *)malloc(d->sz_data);
  memcpy(dnew, d, d->sz_data);

  dnew->sz_data = d->sz_data;
  dnew->ntot = d->ntot;
  realnum *P = dnew->data;

  FOR_COMPONENTS(c) DOCMP2 {
    if (d->P[c][cmp]) {
      dnew->P[c][cmp] = P;
      P += ntot;
      dnew->P_prev[c][cmp] = P;
      P += ntot * (D_z-1);
      dnew->W_prev[c][cmp] = P;
      P += ntot * (N_z-1);
    }
  }
  return (void *)dnew;
}

void iir_susceptibility::update_P(realnum *W[NUM_FIELD_COMPONENTS][2],
                                         realnum *W_prev[NUM_FIELD_COMPONENTS][2], double dt,
                                         const grid_volume &gv, void *P_internal_data) const {
  iir_data *d = (iir_data *)P_internal_data;
  size_t ntot = d->ntot;
  (void)dt;
  (void)W_prev; // unused;

  FOR_COMPONENTS(c) DOCMP2 {
    if (d->P[c][cmp]) {
      const realnum *w = W[c][cmp], *s = sigma[c][component_direction(c)];
      if (w && s) {
        realnum *p = d->P[c][cmp], *pp = d->P_prev[c][cmp], *wp = d->W_prev[c][cmp];
        LOOP_OVER_VOL_OWNED(gv, c, i) {
          //p[i] = 0;
          if (s[i] != 0) {
            for (int den_c=0; den_c<D_z-1; den_c++) {
              p[i] -= denZ[den_c+1] * pp[den_c * ntot + i];
            }
            // Initialize current polarization by weighting the
            // current field value with the first numerator tap ...
            p[i] = numZ[0] * w[i] * s[i];
            // ... then continue by looping through the rest of 
            // the numerator taps...
            for (int num_c=0; num_c<N_z-1; num_c++) {
              p[i] += numZ[num_c+1] * wp[num_c * ntot + i] * s[i];
            }
            // ..the entire numerator term is multiplied by the sigma weight...
            //p[i] *= s[i];
            //... and now loop through the denominator taps,
            // which weight the previous polarizations 
            // (subtract since they must be moved to the RHS)
          }
          // shift the cached W field
          for (int num_c=N_z-2; num_c>0; num_c--) wp[num_c * ntot + i] = wp[(num_c-1) * ntot + i];
          wp[0 * ntot+i] = w[i];
          // shift the P field
          for (int den_c=D_z-2; den_c>0; den_c--) pp[den_c * ntot + i] = pp[(den_c-1) * ntot + i];
          pp[0 * ntot+i] = p[i];
          //master_printf("W[%i]: %3.1e\n",i,w[i]);
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
iir_susceptibility::iir_susceptibility(std::vector<double> vec_numS, std::vector<double> vec_denS, double dt):
  N (vec_numS.size()), D (vec_denS.size()), N_z (vec_denS.size()), D_z (vec_denS.size()), T (dt),
  numS (new realnum[N]), denS (new realnum[D]),
  numZ (new realnum[N_z]), denZ (new realnum[D_z])
 {
   //if (vec_numS.size() > vec_denS.size()) abort("error: for iir-filter susceptibilities, the numerator rank cannot be higher than the denominator rank.\n");
    
    // copy over the S domain numerator and denominator vectors to arrays
    for (int i=0; i<N; i++) numS[i] = vec_numS[i];
    for (int i=0; i<D; i++) denS[i] = vec_denS[i];

    // calculate the Z domain numerator and denominator
    for (int i=0; i<N; i++) master_printf("N_s[%i]: %3.5e\n",i,numS[i]);
    for (int i=0; i<D; i++) master_printf("D_s[%i]: %3.5e\n",i,denS[i]);
    master_printf("+++++++++++++++\n");
    master_printf("dt: %3.5e\n",T);
    tustins_method(numS, N, denS, D, numZ, N_z, denZ, D_z, T);
    for (int i=0; i<N_z; i++) master_printf("N_z[%i]: %3.5e\n",i,numZ[i]);
    for (int i=0; i<D_z; i++) master_printf("D_z[%i]: %3.5e\n",i,denZ[i]);
}

iir_susceptibility::iir_susceptibility(const iir_susceptibility &from) : susceptibility(from) {
  T = from.T;
  N = from.N;
  D = from.D;
  N_z = from.N_z;
  D_z = from.D_z;
  
  numS = new realnum[N];
  memcpy(numS, from.numS, sizeof(realnum) * N);
  denS = new realnum[D];
  memcpy(denS, from.denS, sizeof(realnum) * D); 
  numZ = new realnum[N_z];
  memcpy(numZ, from.numZ, sizeof(realnum) * N_z);
  denZ = new realnum[D_z];
  memcpy(denZ, from.denZ, sizeof(realnum) * D_z);
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
  for (int i=0;i<N;i++) num += std::pow(std::complex<double>(0,freq),N-i) * numS[i];
  for (int i=0;i<D;i++) den += std::pow(std::complex<double>(0,freq),D-i) * denS[i];

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