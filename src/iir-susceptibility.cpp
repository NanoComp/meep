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
    int c, k, j, n=N-1;
    double sum, beta_powers[n], gamma_powers[n];

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
    /*
    // Flip array
    for (c = 0; c < N/2; c++) {
      j = r[c];
      r[c] = r[n];
      r[n] = j;
      n--;
  }*/
}

// Uses Tustins method to discretize the transfer function. S = (2/T)*(1-z^-1)/(1+z^-1)
void tustins_method(double *num, int num_N, double *den, int den_N, double dt) {
  double c = 2.0 / dt;
  double alpha = 2.0, beta = -2.0, delta = dt, gamma = dt;
  int k;

  if (den_N < num_N) abort("error: the numerator cannot be higher rank than the denominator");
  if (num_N < den_N) {
    double * temp = new double[den_N];
    for (k=0; k < den_N; k++) {
      temp[k] = (k <=num_N) ? num[k] : 0.0;
    }
    delete num;
    num = temp;
    num_N = den_N;
  }

  // Process numerator
  polynomial_transform(num,num_N,alpha,beta,delta,gamma);

  // Process denominator
  polynomial_transform(den,den_N,alpha,beta,delta,gamma);

  // Normalize
  double norm = den[0];
  for (k=0; k < num_N; k++) num[k] = num[k] / norm;
  for (k=0; k < den_N; k++) den[k] = den[k] / norm;
  for (k=0;k<num_N; k++) master_printf("num: %20.17le\n",num[k]);

}

// Uses backward difference method to discretize the transfer function. S = (1-z^-1)/T
void backward_difference(double *num, int num_N, double *den, int den_N, double dt) {
  double alpha = 1.0, beta = -1.0, delta = dt, gamma = 0;

  // Process numerator
  polynomial_transform(num,num_N,alpha,beta,delta,gamma);

  // Process denominator
  polynomial_transform(den,den_N,alpha,beta,delta,gamma);

  // Normalize
  /*
  double norm = den[0];
  int k;
  for (k=0; k < num_N; k++) num[k] = num[k] / norm;
  for (k=0; k < den_N; k++) den[k] = den[k] / norm;
  */
}

// --------------------------------------------------------------- //
// IIR transfer function susceptibility
// --------------------------------------------------------------- //

/*
typedef struct {
  size_t sz_data;
  size_t ntot;
  realnum *P;
  realnum *P_prev;
  realnum *W_prev;
  realnum data[1];
} iir_data;

// constructor must take the s domain coefficients and transform
// them to the z domain and store them in the class.
iir_susceptibility::iir_susceptibility(double *num, int num_N, double *den, int den_N, double dt){
    numz_N = num_N;
    denz_N = den_N;
    tustins_method(num,num_N,den,den_N, dt);
    numz = num;
    denz = den;
}

*/

}