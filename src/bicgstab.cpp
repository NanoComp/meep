/* Copyright (C) 2005-2026 Massachusetts Institute of Technology
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

#include <math.h>
#include <string.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "meep/mympi.hpp"
#include "bicgstab.hpp"

#include "config.h"

/* bicgstab() implements an iterative solver for non-symmetric linear
   operators, using the algorithm described in:

      Gerard L. G. Sleijpen and Diederik R. Fokkema, "BiCGSTAB(L) for
      linear equations involving unsymmetric matrices with complex
      spectrum," Electronic Trans. on Numerical Analysis 1, 11-32
      (1993).

   and also:

      Gerard L.G. Sleijpen, Henk A. van der Vorst, and Diederik
      R. Fokkema, " BiCGstab(L) and other Hybrid Bi-CG Methods,"
      Numerical Algorithms 7, 75-109 (1994).

   This is a generalization of the stabilized biconjugate-gradient
   (BiCGSTAB) algorithm proposed by van der Vorst (and described
   in the book _Templates for the Solution of Linear Systems_ by
   Barrett et al.)  BiCGSTAB(1) is equivalent to BiCGSTAB, and
   BiCGSTAB(2) is a slightly more efficient version of the BiCGSTAB2
   algorithm by Gutknecht, while BiCGSTAB(L>2) is a further
   generalization.

   The reason that we use this generalization of BiCGSTAB is that the
   BiCGSTAB(1) algorithm was observed by Sleijpen and Fokkema to have
   poor (or even failing) convergence when the linear operator has
   near-pure imaginary eigenvalues.  This is precisely the case for
   our problem (the eigenvalues of the timestep operator are i*omega),
   and we observed precisely such stagnation of convergence.  The
   BiCGSTAB(2) algorithm was reported to fix most such convergence
   problems, and indeed L > 1 seems to converge well for us. */

/* Other variations to explore:

   G. L. G. Sleijpen and H. A. van der Vorst, "Reliable updated
   residuals in hybrid Bi-CG methods," Computing 56 (2), 141-163
   (1996).

   G. L. G. Sleijpen and H. A. van der Vorst, "Maintaining convergence
   properties of BiCGstab methods in finite precision arithmetic,"
   Numerical Algorithms 10, 203-223 (1995).

   See also code on Sleijpen's web page:
                 http://www.math.uu.nl/people/sleijpen/

*/

using namespace std;

namespace meep {

/* The BLAS-1 operations below dominate the cost of solve_cw unless they are
   vectorized and threaded.  For the reductions, we do not use an OpenMP
   "reduction" clause, because the order in which it combines the per-thread
   partial sums varies from run to run, which makes the results (and, in single
   precision, even the iteration counts) non-reproducible.  Instead, blocked_sum
   returns sum_{i<n} term(i), summing fixed-size blocks in parallel (with a fixed
   number of independent accumulators, which allows SIMD vectorization) and then
   adding the block sums in order, so the result is bitwise identical regardless
   of the number of threads.  Terms are accumulated in double precision even
   when realnum is float. */
template <typename F> static double blocked_sum(size_t n, F term) {
  constexpr size_t B = 8192, W = 32; // block size and number of accumulators
  const size_t nb = (n + B - 1) / B;
  std::vector<double> partial(nb);
#pragma omp parallel for schedule(static)
  for (size_t k = 0; k < nb; ++k) {
    const size_t lo = k * B, hi = std::min(n, lo + B);
    double acc[W] = {0};
    size_t i = lo;
    for (; i + W <= hi; i += W) {
#pragma omp simd
      for (size_t j = 0; j < W; ++j)
        acc[j] += term(i + j);
    }
    double s = 0;
    for (; i < hi; ++i)
      s += term(i);
    for (size_t j = 0; j < W; ++j)
      s += acc[j];
    partial[k] = s;
  }
  double sum = 0;
  for (size_t k = 0; k < nb; ++k)
    sum += partial[k];
  return sum;
}

static double dot(size_t n, const realnum *x, const realnum *y) {
  return sum_to_all(blocked_sum(n, [=](size_t i) { return double(x[i]) * double(y[i]); }));
}

static double norm2(size_t n, const realnum *x) {
  const double sum = dot(n, x, x);
  if (!std::isinf(sum)) return sqrt(sum); // (NaN propagates)

  // overflow in sum of squares: fall back to the slower but robust scaled computation
  double xmax = 0;
#pragma omp parallel for reduction(max : xmax) // (max is exact, hence order-independent)
  for (size_t i = 0; i < n; ++i)
    xmax = std::max(xmax, fabs(double(x[i])));
  xmax = max_to_all(xmax);
  if (xmax == 0) return 0;
  const double scale = 1.0 / xmax;
  return xmax * sqrt(sum_to_all(blocked_sum(n, [=](size_t i) {
           const double xs = scale * x[i];
           return xs * xs;
         })));
}

// x += a*y
static void xpay(size_t n, realnum *x, double a, const realnum *y) {
#pragma omp parallel for simd
  for (size_t m = 0; m < n; ++m)
    x[m] += a * y[m];
}

// y = x - b*y
static void xmby(size_t n, const realnum *x, double b, realnum *y) {
#pragma omp parallel for simd
  for (size_t m = 0; m < n; ++m)
    y[m] = x[m] - b * y[m];
}

#define MEEP_MIN_OUTPUT_TIME 4.0 // output no more often than this many seconds

/* BiCGSTAB(L) algorithm for the n-by-n problem Ax = b */
ptrdiff_t bicgstabL(const int L, const size_t n, realnum *x, bicgstab_op A, void *Adata,
                    const realnum *b, const double tol, int *iters, realnum *work,
                    const bool quiet) {
  if (!work) return (2 * L + 3) * n; // required workspace

  // (std::vector rather than new[] so that nothing leaks if A throws, e.g. via meep::abort)
  std::vector<realnum *> r(L + 1), u(L + 1);
  for (int i = 0; i <= L; ++i) {
    r[i] = work + i * n;
    u[i] = work + (L + 1 + i) * n;
  }

  double bnrm = norm2(n, b);
  if (bnrm == 0.0) bnrm = 1.0;

  int iter = 0;
  double last_output_wall_time = wall_time();

  std::vector<double> gamma(L + 1), gamma_p(L + 1), gamma_pp(L + 1), tau(L * L), sigma(L + 1);

  int ierr = 0; // error code to return, if any
  const double breaktol = 1e-30;

  /* In finite precision (especially single precision), the residual eventually stops
     decreasing; if tol is below that floor, BiCGSTAB(L) does not merely stagnate but
     diverges exponentially, eventually overflowing to Inf/NaN.  So we keep a copy of the
     iterate xbest with the smallest (recursively updated) residual, stop if the residual
     grows by more than divtol relative to that, and return xbest on any failure.  The
     residual of a converging BiCGSTAB(L) can temporarily spike far above its running minimum
     (by up to ~4e4 in solve_cw tests), so divtol must be large; diverging iterations grow
     exponentially, so a large divtol costs only a modst number of extra iterations. */
  const double divtol = 1e8;
  realnum *xbest = work + (2 * L + 3) * n;
  double best_resid = HUGE_VAL;
  bool have_best = false;

  // rtilde = r[0] = b - Ax
  realnum *rtilde = work + (2 * L + 2) * n;
  A(x, r[0], Adata);
  {
    realnum *r0 = r[0];
#pragma omp parallel for simd
    for (size_t m = 0; m < n; ++m)
      rtilde[m] = r0[m] = b[m] - r0[m];
  }

  { /* Sleipjen normalizes rtilde in his code; it seems to help slightly */
    const double s = 1.0 / norm2(n, rtilde);
#pragma omp parallel for simd
    for (size_t m = 0; m < n; ++m)
      rtilde[m] *= s;
  }

  memset(u[0], 0, sizeof(realnum) * n); // u[0] = 0

  double rho = 1.0, alpha = 0, omega = 1;

  double resid;
  while (true) {
    resid = norm2(n, r[0]);
    if (!std::isfinite(resid) || resid > divtol * best_resid) {
      ierr = 2; // diverged
      break;
    }
    if (resid < best_resid) {
      best_resid = resid;
      have_best = true;
      memcpy(xbest, x, sizeof(realnum) * n);
    }
    if (resid <= tol * bnrm) break; // converged
    if (iter == *iters) {
      ierr = 1; // maximum number of iterations reached
      break;
    }

    ++iter;
    if (!quiet && wall_time() > last_output_wall_time + MEEP_MIN_OUTPUT_TIME) {
      master_printf("residual[%d] = %g\n", iter, resid / bnrm);
      last_output_wall_time = wall_time();
    }

    rho = -omega * rho;
    for (int j = 0; j < L; ++j) {
      if (fabs(rho) < breaktol) {
        ierr = -1; // breakdown
        goto finish;
      }
      double rho1 = dot(n, r[j], rtilde);
      double beta = alpha * rho1 / rho;
      rho = rho1;
      for (int i = 0; i <= j; ++i)
        xmby(n, r[i], beta, u[i]); // u[i] = r[i] - beta * u[i]
      A(u[j], u[j + 1], Adata);
      alpha = rho / dot(n, u[j + 1], rtilde);
      for (int i = 0; i <= j; ++i)
        xpay(n, r[i], -alpha, u[i + 1]);
      A(r[j], r[j + 1], Adata);
      xpay(n, x, alpha, u[0]);
    }

    for (int j = 1; j <= L; ++j) {
      for (int i = 1; i < j; ++i) {
        int ij = (j - 1) * L + (i - 1);
        tau[ij] = dot(n, r[j], r[i]) / sigma[i];
        xpay(n, r[j], -tau[ij], r[i]);
      }
      sigma[j] = dot(n, r[j], r[j]);
      gamma_p[j] = dot(n, r[0], r[j]) / sigma[j];
    }

    omega = gamma[L] = gamma_p[L];
    for (int j = L - 1; j >= 1; --j) {
      gamma[j] = gamma_p[j];
      for (int i = j + 1; i <= L; ++i)
        gamma[j] -= tau[(i - 1) * L + (j - 1)] * gamma[i];
    }
    for (int j = 1; j < L; ++j) {
      gamma_pp[j] = gamma[j + 1];
      for (int i = j + 1; i < L; ++i)
        gamma_pp[j] += tau[(i - 1) * L + (j - 1)] * gamma[i + 1];
    }

    /* Update x, r[0], and u[0] in a single pass over memory (rather than 3L separate
       xpay passes, since these memory-bound loops are limited by bandwidth):
          x    += gamma[1] r[0] + sum_{j=1}^{L-1} gamma_pp[j] r[j]
          r[0] -= sum_{j=1}^{L} gamma_p[j] r[j]
          u[0] -= sum_{j=1}^{L} gamma[j] u[j]
       where the x update uses r[0] from before its update. */
    {
      realnum *const *rp = r.data();
      realnum *const *up = u.data();
      const double *g = gamma.data(), *gp = gamma_p.data(), *gpp = gamma_pp.data();
#pragma omp parallel for
      for (size_t m = 0; m < n; ++m) {
        const double r0 = rp[0][m];
        double xm = x[m] + g[1] * r0;
        double rm = r0 - gp[L] * rp[L][m];
        double um = up[0][m] - g[L] * up[L][m];
        for (int j = 1; j < L; ++j) {
          const double rj = rp[j][m];
          xm += gpp[j] * rj;
          rm -= gp[j] * rj;
          um -= g[j] * up[j][m];
        }
        x[m] = xm;
        rp[0][m] = rm;
        up[0][m] = um;
      }
    }
  }

finish:
  if (ierr != 0 && have_best) { // on failure, return the best iterate found
    memcpy(x, xbest, sizeof(realnum) * n);
    resid = best_resid;
  }
  if (!quiet) {
    if (ierr == 2) master_printf("residual diverged: returning best iterate\n");
    master_printf("final residual = %g\n", resid / bnrm);
  }

  *iters = iter;
  return ierr;
}

} // namespace meep
