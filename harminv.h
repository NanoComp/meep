/* Copyright (C) 2000 Massachusetts Institute of Technology.
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

#ifndef HARMINV_H
#define HARMINV_H

/* Require C99 complex number support; this is just too painful
   without it.  Alternatively, use the complex<double> STL class in
   C++. */

#include <complex>

#include "config.h"

/**************************************************************************/

typedef complex<double> cmplx;
#  define I cmplx(0,1)
#  define creal(c) real(c)
#  define cimag(c) imag(c)
#  define cabs(c) abs(c)
#  define carg(c) arg(c)
#  define cexp(c) exp(c)
#  define csqrt(c) sqrt(c)

typedef struct harminv_data_struct {
     const cmplx *c;
     int n, K, J, nfreqs;
     double fmin, fmax;
     cmplx *z;
     cmplx *U0, *U1;
     cmplx *B, *u;  /* eigen-solutions of U1*B = u*U0*B */
} *harminv_data;

/**************************************************************************/

extern harminv_data harminv_data_create(int n,
					const cmplx *signal,
					double fmin, double fmax, int nf);
extern void harminv_data_destroy(harminv_data d);

extern void harminv_solve(harminv_data d);
extern void harminv_solve_again(harminv_data d);

extern int harminv_get_num_freqs(harminv_data d);
extern double harminv_get_freq(harminv_data d, int k);
extern double harminv_get_decay(harminv_data d, int k);

extern double *harminv_compute_frequency_errors(harminv_data d);
extern cmplx *harminv_compute_amplitudes(harminv_data d);

/**************************************************************************/

#endif /* HARMINV_H */
