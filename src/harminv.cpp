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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "harminv.h"
#include "check.h"
#include "config.h"

/**************************************************************************/

/* The harminv routines are designed to perform "harmonic inversion."
   That is, given a signal (a set of samples as a function of time),
   they decompose the signal into a finite number of (possibly
   exponentially-decaying) sinusoids.

   Essentially, because we assume that the signal has this form, we
   can determine the frequencies, decay constants, and amplitudes of
   the sinusoids much more accurately than we could via taking the FFT
   and looking at the peaks, for the same number of samples.

   We use an "low-storage filter diagonalization" algorithm for finding
   the sinusoids near a given frequency interval, described in:

   V. A. Mandelshtam and H. S. Taylor, "Harmonic inversion of time
   signals," J. Chem. Phys., vol. 107, no. 17, p. 6756-6769 (Nov. 1
   1997).  See also erratum, ibid, vol. 109, no. 10, p. 4128 (Sep. 8
   1998).

   with a refinement (for generate_U below) described in:

   Rongqing Chen and Hua Guo, "Efficient calculation of matrix
   elements in low storate filter diagonalization," J. Chem. Phys.,
   vol. 111, no. 2, p. 464-471(Jul. 8 1999). 

   The seminal work (though less practical than the M&T algorithm) on
   this class of methods was done by:

   Michael R. Wall and Daniel Neuhauser, "Extraction, through
   filter-diagonalization, of general quantum eigenvalues or classical
   normal mode frequencies from a small number of residues or a
   short-time segment of a signal. I. Theory and application to a
   quantum-dynamics model," J. Chem. Phys., 102, no. 20, p. 8011-8022
   (May 22 1995). */

/**************************************************************************/

#define TWOPI 6.2831853071795864769252867665590057683943388

/**************************************************************************/

/* Crays have float == double, and don't have the Z* functions in
   LAPACK/BLAS...we have to use C*.  Sigh. */
#if defined(CRAY) || defined(_UNICOS) || defined(_CRAYMPP)
#  define BLAS_FUNC(x,X) F77_FUNC(c##x,C##X)
#else  /* ! CRAY */
#  define BLAS_FUNC(x,X) F77_FUNC(z##x,Z##X)
#endif /* ! CRAY */

#define ZGEEV BLAS_FUNC(geev,GEEV)
#define ZGEMM BLAS_FUNC(gemm,GEMM)
#define ZCOPY BLAS_FUNC(copy,COPY)
#define ZAXPY BLAS_FUNC(axpy,AXPY)
#define ZGEMV BLAS_FUNC(gemv,GEMV)
#define ZSCAL BLAS_FUNC(scal,SCAL)

#ifdef __cplusplus
extern "C" {
#endif

/* We have to pass strings in special ways on Crays, even
   for passing a single character as with LAPACK.  Sigh. */
#if defined(CRAY) || defined(_UNICOS) || defined(_CRAYMPP)
#  include <fortran.h>
#  define FCHARP _fcd
#  define F_(s) _cptofcd(s,1)  /* second argument is the string length */
#else  /* ! CRAY */
#  define FCHARP char*
#  define F_(s) (s)
#endif

extern void ZGEEV(FCHARP,FCHARP, int*, cmplx*,int*, cmplx*, cmplx*,int*,
		  cmplx*,int*, cmplx*,int*, double*, int*);
extern void ZGEMM(FCHARP,FCHARP, int*,int*,int*, cmplx*,
		  cmplx*,int*, cmplx*,int*, cmplx*, cmplx*,int*);
extern void ZCOPY(int*, cmplx*,int*, cmplx*,int*);
extern void ZAXPY(int*, cmplx*, cmplx*,int*, cmplx*,int*);
extern void ZGEMV(FCHARP, int*,int*, cmplx*, cmplx*,int*, cmplx*,int*,
		  cmplx*, cmplx*,int*);
extern void ZSCAL(int*, cmplx*, cmplx*,int*);

#ifdef __cplusplus
} /* extern "C" */
#endif

/**************************************************************************/

/* compute c^n, where n is an integer: */
static cmplx cpow_i(cmplx c, int n)
{
     if (n < 0)
	  return (1.0 / cpow_i(c, -n));
     else {
	  cmplx result = 1;
	  while (n > 1) {
	       if (n % 2 == 1)
		    result *= c;
	       c *= c;
	       n /= 2;
	  }
	  if (n > 0)
	       result *= c;
	  return result;
     }
}

#define SMALL (1e-12)
#define SMALL2 (SMALL*SMALL)
#define C_CLOSE(c1,c2) (creal(((c1) - (c2)) * conj((c1) - (c2))) < SMALL2)

/**************************************************************************/

/* Initialize the JxJ2 matrix U = U_p(z,z2), as described in M&T.
   Also, if U1 != NULL, then set U1 = U_{p+1}(z,z2).  If z == z2, it
   must be the case that no two elements of z are the same and that
   J2 == J1; in this case the matrix U will be symmetric.  Note that c
   must be an array whose size is at least 2*K+p-1 elements.  */
static void generate_U(cmplx *U, cmplx *U1,
		       int p, 
		       const cmplx *c,
		       int K,
		       int J, int J2, const cmplx *z, const cmplx *z2)
{
     int M = K - 1;
     int i, j, m;
     /* temp. arrays for 1/z, z^(-m), z^(-M), the G function of C&G,
	and the diagonal elements D[i] = U(z[i],z[i]): */
     cmplx *z_inv, *z_m, *z_M, *G, *G_M, *D; 
     cmplx *z2_inv, *z2_m, *z2_M, *G2, *G2_M; 

     CHECK(U && c && z && z2, "invalid arguments to generate_U");
     CHECK(z != z2 || J == J2, "invalid sizes passed to generate_U");
     
     /* Now, compute U according to eqs. 25-27 of Chen & Guo, but
        using the notation of eq. 25 of M&T.  This operation has
        complexity O(N*J + J*J).  At the same time, we can compute the
        matrix U1 as well by eqs. 29-30 of C&G, saving an extra pass
        over the input array. */

     /* first, set up some temporary arrays for caching things like
	z^m and 1/z, so we don't need to recompute them all the time. */
     CHK_MALLOC(z_inv, cmplx, J);
     CHK_MALLOC(z_m, cmplx, J);
     CHK_MALLOC(z_M, cmplx, J);
     CHK_MALLOC(G, cmplx, J);
     CHK_MALLOC(G_M, cmplx, J);
     CHK_MALLOC(D, cmplx, J);
     for (i = 0; i < J; ++i) {
	  z_inv[i] = 1.0/z[i];
	  z_m[i] = 1;
	  z_M[i] = 1.0 / cpow_i(z[i], M);
	  D[i] = G[i] = G_M[i] = 0;
     }
     if (z2 != z) {
	  CHK_MALLOC(z2_inv, cmplx, J2);
	  CHK_MALLOC(z2_m, cmplx, J2);
	  CHK_MALLOC(z2_M, cmplx, J2);
	  CHK_MALLOC(G2, cmplx, J2);
	  CHK_MALLOC(G2_M, cmplx, J2);
	  for (i = 0; i < J2; ++i) {
	       z2_inv[i] = 1.0/z2[i];
	       z2_m[i] = 1;
	       z2_M[i] = 1.0 / cpow_i(z2[i], M);
	       G2[i] = G2_M[i] = 0;
	  }
     }
     else {
	  z2_inv = z2_m = z2_M = G2 = G2_M = NULL;
     }

     /* First, loop over the signal array (c), building up the
	spectral functions G and G_M (corresponding to G_p and
	G_{p+M+1} in C&G), as well as the diagonal matrix entries: */
     for (m = 0; m <= M; ++m) {
	  cmplx c1 = c[m + p], c2;
	  double d = m + 1; /* M - fabs(M - m) + 1 */
	  double D2 = M - m; /* M - fabs(M - (m + M + 1)) + 1 */

	  if (m < M)
	       c2 = c[m + p + M + 1];
	  else
	       c2 = 0;

	  for (i = 0; i < J; ++i) {
	       cmplx x1 = z_m[i] * c1;
	       cmplx x2 = z_m[i] * c2;
	       G[i] += x1;
	       G_M[i] += x2;
	       D[i] += x1 * d + x2 * D2 * z_M[i] * z_inv[i];
	       z_m[i] *= z_inv[i];
	  }
	  if (z2 != z)
	       for (i = 0; i < J2; ++i) {
		    G2[i] += z2_m[i] * c1;
		    G2_M[i] += z2_m[i] * c2;
		    z2_m[i] *= z2_inv[i];
	       }
     }

     /* Compute U (or just the upper part if U is symmetric), via the
        formula from C&G; compute U1 at the same time as in C&G. */
     if (z2 != z) {
	  for (i = 0; i < J; ++i)
	       for (j = 0; j < J2; ++j) {
		    if (C_CLOSE(z[i], z2[j]))
			 U[i*J2 + j] = D[i];
		    else
			 U[i*J2 + j] = (z[i] * G2[j] - z2[j] * G[i] +
				       z2_M[j] * G_M[i] - z_M[i] * G2_M[j])
			      / (z[i] - z2[j]);
	       }

	  if (U1)
	       for (i = 0; i < J; ++i)
		    for (j = 0; j < J2; ++j) {
			 if (C_CLOSE(z[i], z2[j]))
			      U1[i*J2 + j] = z[i] * (D[i] - G[i]) +
				   z_M[i] * G_M[i];
			 else
			      U1[i*J2 + j] = (z[i] * z2[j] * (G2[j] - G[i])
					     + z2_M[j] * z[i] * G_M[i] 
					     - z_M[i] * z2[j] * G2_M[j])
				   / (z[i] - z2[j]);
		    }
     }
     else {  /* z == z2 */
	  for (i = 0; i < J; ++i) {
	       U[i*J + i] = D[i];
	       for (j = i + 1; j < J; ++j) {
		    U[i*J + j] = (z[i] * G[j] - z[j] * G[i] +
				  z_M[j] * G_M[i] - z_M[i] * G_M[j])
			 / (z[i] - z[j]);
	       }
	  }

	  if (U1)
	       for (i = 0; i < J; ++i) {
		    U1[i*J + i] = z[i] * (D[i] - G[i]) + z_M[i] * G_M[i];
		    for (j = i + 1; j < J; ++j) {
			 U1[i*J + j] = (z[i] * z[j] * (G[j] - G[i])
					+ z_M[j] * z[i] * G_M[i] 
					- z_M[i] * z[j] * G_M[j])
			      / (z[i] - z[j]);
		    }
	       }
     }


     /* finally, copy the upper to the lower triangle if U is symmetric: */
     if (z == z2) {
	  for (i = 0; i < J; ++i)
	       for (j = i + 1; j < J; ++j)
		    U[j*J + i] = U[i*J + j];
	  if (U1)
	       for (i = 0; i < J; ++i)
		    for (j = i + 1; j < J; ++j)
			 U1[j*J + i] = U1[i*J + j];
     }

     free(G2_M);
     free(G2);
     free(z2_M);
     free(z2_m);
     free(z2_inv);

     free(D);
     free(G_M);
     free(G);
     free(z_M);
     free(z_m);
     free(z_inv);
}

/**************************************************************************/

static void init_z(harminv_data d, int J, cmplx *z)
{
     d->J = J;
     d->z = z;
     CHK_MALLOC(d->U0, cmplx, J*J);
     CHK_MALLOC(d->U1, cmplx, J*J);
     generate_U(d->U0, d->U1, 0, d->c, d->K, d->J, d->J, d->z, d->z);
}

/**************************************************************************/

harminv_data harminv_data_create(int n,
				 const cmplx *signal,
				 double fmin, double fmax, int nf)
{
     int i;
     harminv_data d;

     CHECK(nf == 0 || nf > 1, "# frequencies must be zero or > 1");
     CHECK(n > 0, "invalid number of data points");
     CHECK(signal, "invalid NULL signal array");
     CHECK(fmin < fmax, "should have fmin < fmax");

     if (!nf) {
	  /* use "reasonable choice" suggested by M&T: */
	  nf = (int) (n*(fmax-fmin)/2 + 0.5);
	  if (nf < 2)
	       nf = 2;
     }

     CHK_MALLOC(d, struct harminv_data_struct, 1);
     d->c = signal;
     d->n = n;
     d->K = (n - 1) / 2;
     d->fmin = fmin;
     d->fmax = fmax;
     
     CHK_MALLOC(d->z, cmplx, nf);
     for (i = 0; i < nf; ++i)
	  d->z[i] = cexp(-I * TWOPI * (fmin + i * ((fmax - fmin) / (nf - 1))));

     init_z(d, nf, d->z);

     d->nfreqs = 0;
     d->B = d->u = NULL;  /* we haven't computed eigen-solutions yet */

     return d;
}

/**************************************************************************/

void harminv_data_destroy(harminv_data d)
{
     if (d) {
	  free(d->u); free(d->B);
	  free(d->U1); free(d->U0);
	  free(d->z);
	  free(d);
     }
}

/**************************************************************************/

/* Compute the symmetric dot product of x and y, both vectors of
   length n.  If they are column-vectors, this is: transpose(x) * y. */
static cmplx symmetric_dot(int n, cmplx *x, cmplx *y)
{
     cmplx dot = 0;
     int i;
     for (i = 0; i < n; ++i)
	  dot += x[i] * y[i];
     return dot;
}

/* Solve for the eigenvalues (v) and eigenvectors (rows of V) of the
   complex-symmetric n x n matrix A.  The eigenvectors are normalized
   to 1 according to the symmetric dot product (i.e. no complex
   conjugation). */
static void solve_eigenvects(int n, cmplx *A, cmplx *V, cmplx *v)
{
     int lwork, info;
     cmplx *work;
     double *rwork;

     /* Unfortunately, LAPACK doesn't have a special solver for the
	complex-symmetric eigenproblem.  For now, just use the general
	non-symmetric solver, and realize that the left eigenvectors
	are the complex-conjugates of the right eigenvectors. */

#if 0  /* LAPACK seems to be buggy here, returning ridiculous sizes at times */
     cmplx wsize;
     lwork = -1; /* compute optimal workspace size */
     ZGEEV(F_("N"), F_("V"), &n, A, &n, v, V, &n, V, &n, &wsize, &lwork, rwork, &info);
     if (info == 0)
	  lwork = floor(creal(wsize) + 0.5);
     else
	  lwork = 2*n;
     CHECK(lwork > 0, "zgeev is not returning a positive work size!");
#else
     lwork = 4*n; /* minimum is 2*n; we'll be generous. */
#endif

     CHK_MALLOC(rwork, double, 2*n);
     CHK_MALLOC(work, cmplx, lwork);

     ZGEEV(F_("N"), F_("V"), &n, A, &n, v, V, &n, V, &n, work, &lwork, rwork, &info);

     free(work);
     free(rwork);

     CHECK(info >= 0, "invalid argument to ZGEEV");
     CHECK(info <= 0, "failed convergence in ZGEEV");

     /* Finally, we need to fix the normalization of the eigenvectors,
        since LAPACK normalizes them under the ordinary dot product,
        i.e. with complex conjugation.  (In principle, do we also need
        to re-orthogonalize, for the case of degenerate eigenvalues?) */
     {
	  int i, one = 1;
	  for (i = 0; i < n; ++i) {
	       cmplx norm = 1.0 / csqrt(symmetric_dot(n, V+i*n, V+i*n));
	       ZSCAL(&n, &norm, V+i*n, &one);
	  }
     }
}

/**************************************************************************/

/* how conservative do we need to be for this? */
#define SINGULAR_THRESHOLD 1e-5

/* Solve the eigenvalue problem U1 b = u U0 b, where b is the eigenvector
   and u is the eigenvalue.  u = exp(iwt - at) then contains both the
   frequency and the decay constant. */
void harminv_solve(harminv_data d)
{
     int J, i, one=1;
     cmplx zone = 1.0, zzero = 0.0;
     cmplx *V0, *v0, *H1, *V1; /* for eigensolutions of U0 and U1 */
     double max_v0 = 0.0;
     
     J = d->J;
     CHK_MALLOC(V0, cmplx, J*J);
     CHK_MALLOC(v0, cmplx, J);

     /* Unfortunately, U0 is very likely to be singular, so we must
	first extract the non-singular eigenvectors and only work in
	that sub-space.  See the Wall & Neuhauser paper.  */
     solve_eigenvects(J, d->U0, V0, v0);
     
     /* find maximum |eigenvalue| */
     for (i = 0; i < J; ++i) {
	  double v = cabs(v0[i]);
      //      printf("v[%d] is %lg\n", i, v);
	  if (v > max_v0)
	       max_v0 = v;
     }

     /* we must remove the singular components of U0, those
	that are less than some threshold times the maximum eigenvalue.
        Also, we need to scale the eigenvectors by 1/sqrt(eigenval). */
     d->nfreqs = J;
     for (i = 0; i < J; ++i) {
	  if (cabs(v0[i]) < SINGULAR_THRESHOLD * max_v0) {
	       v0[i] = 0; /* tag as a "hole" */
	       d->nfreqs -= 1;
	  }
	  else { /* not singular */
	       cmplx s;
	       int j;
	       /* move the eigenvector to the first "hole" left by
		  deleting singular eigenvalues: */
	       for (j = 0; j < i && v0[j] != 0.0; ++j)
		    ;
	       if (j < i) {
		    ZCOPY(&J, V0 + i*J, &one, V0 + j*J, &one);
		    v0[j] = v0[i];
		    v0[i] = 0; /* tag as a "hole" */
	       }
	       s = 1.0 / csqrt(v0[j]);
	       ZSCAL(&J, &s, V0 + j*J, &one);
	  }
     }

     CHK_MALLOC(d->B, cmplx, d->nfreqs * J);
     CHK_MALLOC(d->u, cmplx, d->nfreqs);
     CHK_MALLOC(V1, cmplx, d->nfreqs * d->nfreqs);
     CHK_MALLOC(H1, cmplx, d->nfreqs * d->nfreqs);

     /* compute H1 = V0 * U1 * V0': */

     /* B = V0 * U1: */
     ZGEMM(F_("N"), F_("N"), &J, &d->nfreqs, &J,
	   &zone, d->U1, &J, V0, &J, &zzero, d->B, &J);
     /* H1 = B * transpose(V0) */
     ZGEMM(F_("T"), F_("N"), &d->nfreqs, &d->nfreqs, &J,
	   &zone, V0, &J, d->B, &J, &zzero, H1, &d->nfreqs);

     /* Finally, we can find the eigenvalues and eigenvectors: */
     solve_eigenvects(d->nfreqs, H1, V1, d->u);
     /* B = V1 * V0: */
     ZGEMM(F_("N"), F_("N"), &J, &d->nfreqs, &d->nfreqs,
	   &zone, V0, &J, V1, &d->nfreqs, &zzero, d->B, &J);

     free(H1);
     free(V1);
     free(v0);
     free(V0);
}

/**************************************************************************/

/* After solving once, solve again using the solutions from last
   time as the input to the spectra estimator this time. */
void harminv_solve_again(harminv_data d)
{
     int i;
     CHECK(d->B && d->u, "haven't computed eigensolutions yet");

     free(d->B);
     free(d->U1); free(d->U0);
     free(d->z);

     /* Spectral grid needs to be on the unit circle or system is unstable: */
     for (i = 0; i < d->nfreqs; ++i)
	  d->u[i] /= cabs(d->u[i]);

     init_z(d, d->nfreqs, d->u);

     d->nfreqs = 0;
     d->B = d->u = NULL;

     harminv_solve(d);
}

/**************************************************************************/

#define NORMSQR(c) (creal(c) * creal(c) + cimag(c) * cimag(c))

/* Returns an array (of size harminv_get_num_freqs(d)) of estimates
   for the |error| in the solution frequencies.  Solutions with
   errors not << 1 are likely to be spurious. */
double *harminv_compute_frequency_errors(harminv_data d)
{
     int i, J2, one = 1;
     cmplx *U2, *R, *r;
     double *freq_err;

     CHECK(d->B && d->u, "haven't computed eigensolutions yet");

     CHK_MALLOC(freq_err, double, d->nfreqs);
     
     J2 = d->J*d->J;
     CHK_MALLOC(U2, cmplx, J2);
     generate_U(U2, NULL, 2, d->c, d->K, d->J, d->J, d->z, d->z);
     CHK_MALLOC(R, cmplx, J2);
     CHK_MALLOC(r, cmplx, d->J);

     /* For each eigenstate, compute an estimate of the error, as
	suggested in W&N, eq. (2.19).  (Note that this equation is
	properly normalized, so that the error / u^2 should be
	dimensionless.) */

     for (i = 0; i < d->nfreqs; ++i) {
	  cmplx u2m = -(d->u[i] * d->u[i]);
	  cmplx zone = 1.0, zzero = 0.0;

	  /* Compute R = (U2 - u^2 U0): */
	  ZCOPY(&J2, U2, &one, R, &one);
	  ZAXPY(&J2, &u2m, d->U0, &one, R, &one);

	  /* Compute residual r = R*b.  ("T" transpose argument to zgemv
	     is for row/column-major conversion.  In some sense it doesn't
	     matter since R is symmetric, but this should be more efficient
	     since it can then read R row-wise (contiguously).)*/
	  ZGEMV(F_("T"), &d->J, &d->J, 
		&zone, R, &d->J, d->B + i * d->J, &one,
		&zzero, r, &one);

	  freq_err[i] =
	       2 * cabs(symmetric_dot(d->J, d->B + i * d->J, r)) / cabs(u2m);
     }

     free(r);
     free(R);
     free(U2);
     return freq_err;
}

/**************************************************************************/

/* Return an array (of size harminv_get_num_freqs(d)) of complex
   amplitudes of each sinusoid in the solution. */
cmplx *harminv_compute_amplitudes(harminv_data d)
{
     int k, j;
     cmplx *Uu;
     cmplx *a; /* the amplitudes of the eigenfrequencies */

     CHECK(d->B, "haven't computed eigensolutions yet");
     CHK_MALLOC(a, cmplx, d->nfreqs);

     CHK_MALLOC(Uu, cmplx, d->J * d->nfreqs);
     generate_U(Uu, NULL, 0, d->c, d->K, d->J, d->nfreqs, d->z, d->u);

     /* compute the amplitudes via eq. 27 of M&T: */
     for (k = 0; k < d->nfreqs; ++k) {
	  a[k] = 0;
	  for (j = 0; j < d->J; ++j)
	       a[k] += d->B[k * d->J + j] * Uu[j * d->nfreqs + k];
	  a[k] /= d->K;
	  a[k] *= a[k];
     }

     free(Uu);
     return a;
}

/**************************************************************************/

int harminv_get_num_freqs(harminv_data d)
{
     return d->nfreqs;
}

double harminv_get_freq(harminv_data d, int k)
{
     CHECK(d->u, "haven't computed eigensolutions yet");
     CHECK(k >= 0 && k < d->nfreqs,
	   "argument out of range in harminv_get_freq");
     return(-carg(d->u[k]) / TWOPI);
}

double harminv_get_decay(harminv_data d, int k)
{
     CHECK(d->u, "haven't computed eigensolutions yet");
     CHECK(k >= 0 && k < d->nfreqs,
	   "argument out of range in harminv_get_freq");
     return(-log(cabs(d->u[k])));
}

/**************************************************************************/

