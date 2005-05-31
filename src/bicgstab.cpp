#include <math.h>

#include "meep/mympi.hpp"
#include "bicgstab.hpp"

#include "config.h"

namespace meep {

#ifdef HAVE_CBLAS_DDOT
   extern "C" double cblas_ddot(int, const double*, int, const double*, int);
#  define dot(n, x, y) sum_to_all(cblas_ddot(n, x, 1, y, 1))
#else
static double dot(int n, const double *x, const double *y)
{
  double sum = 0;
  for (int i = 0; i < n; ++i) sum += x[i] * y[i];
  return sum_to_all(sum);
}
#endif

static double norm2(int n, const double *x) { return sqrt(dot(n, x, x)); }

/* Stablilized Biconjugate Gradient algorithm for Ax = b,
   based on book _Templates for the Solution of Linear Systems_ */

int bicgstab(int n, double *x,
	     bicgstab_op A, void *Adata, const double *b,
	     double tol,
	     int *iters,
	     double *work)
{
  if (!work) return 5*n; // required workspace
  double *r = work, *p = work+n, *v = work+2*n, *rtilde = work+3*n, *t = work + 4*n;

  double *s = r, *phat = p, *shat = s; // superfluous aliases
  double rho, prev_rho = 0.0, alpha, omega;
  const double breaktol = 1e-30;
  double bnrm;
  
  bnrm = norm2(n, b);
  if (bnrm == 0.0) bnrm = 1.0;
  
  // rtilde = r = b - Ax
  A(x, r, Adata);
  for (int i = 0; i < n; ++i) rtilde[i] = r[i] = b[i] - r[i];
  if (norm2(n, r) < tol)
       return 0;

  int iter = 0;
  do {
    ++iter;

    rho = dot(n, rtilde, r);
    if (fabs(rho) < breaktol) return -1;
    if (prev_rho == 0.0)
      for (int i = 0; i < n; ++i) p[i] = r[i];
    else {
      double beta = (rho / prev_rho) * (alpha / omega);
      for (int i = 0; i < n; ++i) p[i] = r[i] + beta * (p[i] - omega * v[i]);
    }
    prev_rho = rho;
    
    // OMITTED: precondition phat = (1/M) p

    A(phat, v, Adata);

    alpha = rho / dot(n, rtilde, v);
    for (int i = 0; i < n; ++i) s[i] = r[i] - alpha * v[i];

    if (norm2(n, s) < tol) { // CONVERGENCE
      for (int i = 0; i < n; ++i) x[i] += alpha * phat[i];
      break;
    }
    
    // OMITTED: precondition shat = (1/M) s

    A(shat, t, Adata);
    omega = dot(n, t, s) / dot(n, t, t);
    for (int i = 0; i < n; ++i) x[i] += alpha * phat[i] + omega * shat[i];
    for (int i = 0; i < n; ++i) r[i] = s[i] - omega * t[i];

    if (norm2(n, r) < tol * bnrm) // CONVERGENCE TEST #2 from book
	 break;
    if (fabs(omega) < breaktol) return -2;
  } while (iter < *iters);
  *iters = iter;
  return 0;
}

} // namespace meep
