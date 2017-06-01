/*
 * Copyright (c) 2005 Steven G. Johnson
 *
 * Portions (see comments) based on HIntLib (also distributed under
 * the GNU GPL, v2 or later), copyright (c) 2002-2005 Rudolf Schuerer.
 *     (http://www.cosy.sbg.ac.at/~rschuer/hintlib/)
 *
 * Portions (see comments) based on GNU GSL (also distributed under
 * the GNU GPL, v2 or later), copyright (c) 1996-2000 Brian Gough.
 *     (http://www.gnu.org/software/gsl/)
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
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>

/* Adaptive multidimensional integration on hypercubes (or, really,
   hyper-rectangles) using cubature rules.

   A cubature rule takes a function and a hypercube and evaluates
   the function at a small number of points, returning an estimate
   of the integral as well as an estimate of the error, and also
   a suggested dimension of the hypercube to subdivide.

   Given such a rule, the adaptive integration is simple:

   1) Evaluate the cubature rule on the hypercube(s).
      Stop if converged.

   2) Pick the hypercube with the largest estimated error,
      and divide it in two along the suggested dimension.

   3) Goto (1).

*/

typedef double (*integrand) (unsigned ndim, const double *x, void *);

/* Integrate the function f from xmin[dim] to xmax[dim], with at most
   maxEval function evaluations (0 for no limit), until the given
   absolute or relative error is achieved.  val returns the integral,
   and err returns the estimate for the absolute error in val.  The
   return value of the function is 0 on success and non-zero if there
   was an error. */
int adapt_integrate(integrand f, void *fdata,
		    unsigned dim, const double *xmin, const double *xmax, 
		    unsigned maxEval, 
		    double reqAbsError, double reqRelError, 
		    double *val, double *err);

/***************************************************************************/
/* Basic datatypes */

typedef struct {
     double val, err;
} esterr;

static double relError(esterr ee)
{
     return (ee.val == 0.0 ? HUGE_VAL : fabs(ee.err / ee.val));
}

typedef struct {
     unsigned dim;
     double *data;	/* length 2*dim = center followed by half-widths */
     double vol;	/* cache volume = product of widths */
} hypercube;

static double compute_vol(const hypercube *h)
{
     unsigned i;
     double vol = 1;
     for (i = 0; i < h->dim; ++i)
	  vol *= 2 * h->data[i + h->dim];
     return vol;
}

static hypercube make_hypercube(unsigned dim, const double *center, const double *halfwidth)
{
     unsigned i;
     hypercube h;
     h.dim = dim;
     h.data = (double *) malloc(sizeof(double) * dim * 2);
     for (i = 0; i < dim; ++i) {
	  h.data[i] = center[i];
	  h.data[i + dim] = halfwidth[i];
     }
     h.vol = compute_vol(&h);
     return h;
}

static hypercube make_hypercube_range(unsigned dim, const double *xmin, const double *xmax)
{
     hypercube h = make_hypercube(dim, xmin, xmax);
     unsigned i;
     for (i = 0; i < dim; ++i) {
	  h.data[i] = 0.5 * (xmin[i] + xmax[i]);
	  h.data[i + dim] = 0.5 * (xmax[i] - xmin[i]);
     }
     h.vol = compute_vol(&h);
     return h;
}

static void destroy_hypercube(hypercube *h)
{
     free(h->data);
     h->dim = 0;
}

typedef struct {
     hypercube h;
     esterr ee;
     unsigned splitDim;
} region;

static region make_region(const hypercube *h)
{
     region R;
     R.h = make_hypercube(h->dim, h->data, h->data + h->dim);
     R.splitDim = 0;
     return R;
}

static void destroy_region(region *R)
{
     destroy_hypercube(&R->h);
}

static void cut_region(region *R, region *R2)
{
     unsigned d = R->splitDim, dim = R->h.dim;
     *R2 = *R;
     R->h.data[d + dim] *= 0.5;
     R->h.vol *= 0.5;
     R2->h = make_hypercube(dim, R->h.data, R->h.data + dim);
     R->h.data[d] -= R->h.data[d + dim];
     R2->h.data[d] += R->h.data[d + dim];
}

typedef struct rule_s {
     unsigned dim;              /* the dimensionality */
     unsigned num_points;       /* number of evaluation points */
     unsigned (*evalError)(struct rule_s *r, integrand f, void *fdata,
			   const hypercube *h, esterr *ee);
     void (*destroy)(struct rule_s *r);
} rule;

static void destroy_rule(rule *r)
{
     if (r->destroy) r->destroy(r);
     free(r);
}

static region eval_region(region R, integrand f, void *fdata, rule *r)
{
     R.splitDim = r->evalError(r, f, fdata, &R.h, &R.ee);
     return R;
}

/***************************************************************************/
/* Functions to loop over points in a hypercube. */

/* Based on orbitrule.cpp in HIntLib-0.0.10 */

/* ls0 returns the least-significant 0 bit of n (e.g. it returns
   0 if the LSB is 0, it returns 1 if the 2 LSBs are 01, etcetera). */

#if (defined(__GNUC__) || defined(__ICC)) && (defined(__i386__) || defined (__x86_64__))
/* use x86 bit-scan instruction, based on count_trailing_zeros()
   macro in GNU GMP's longlong.h. */
static unsigned ls0(unsigned n)
{
     unsigned count;
     n = ~n;
     __asm__("bsfl %1,%0": "=r"(count):"rm"(n));
     return count;
}
#else
static unsigned ls0(unsigned n)
{
     const unsigned bits[256] = {
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 7,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 8,
     };
     unsigned bit = 0;
     while ((n & 0xff) == 0xff) {
	  n >>= 8;
	  bit += 8;
     }
     return bit + bits[n & 0xff];
}
#endif

/**
 *  Evaluate the integral on all 2^n points (+/-r,...+/-r)
 *
 *  A Gray-code ordering is used to minimize the number of coordinate updates
 *  in p.
 */
static double evalR_Rfs(integrand f, void *fdata, unsigned dim, double *p, const double *c, const double *r)
{
     double sum = 0;
     unsigned i;
     unsigned signs = 0; /* 0/1 bit = +/- for corresponding element of r[] */

     /* We start with the point where r is ADDed in every coordinate
        (this implies signs=0). */
     for (i = 0; i < dim; ++i)
	  p[i] = c[i] + r[i];

     /* Loop through the points in Gray-code ordering */
     for (i = 0;; ++i) {
	  unsigned mask, d;

	  sum += f(dim, p, fdata);

	  d = ls0(i);	/* which coordinate to flip */
	  if (d >= dim)
	       break;

	  /* flip the d-th bit and add/subtract r[d] */
	  mask = 1U << d;
	  signs ^= mask;
	  p[d] = (signs & mask) ? c[d] - r[d] : c[d] + r[d];
     }
     return sum;
}

static double evalRR0_0fs(integrand f, void *fdata, unsigned dim, double *p, const double *c, const double *r)
{
     unsigned i, j;
     double sum = 0;

     for (i = 0; i < dim - 1; ++i) {
	  p[i] = c[i] - r[i];
	  for (j = i + 1; j < dim; ++j) {
	       p[j] = c[j] - r[j];
	       sum += f(dim, p, fdata);
	       p[i] = c[i] + r[i];
	       sum += f(dim, p, fdata);
	       p[j] = c[j] + r[j];
	       sum += f(dim, p, fdata);
	       p[i] = c[i] - r[i];
	       sum += f(dim, p, fdata);

	       p[j] = c[j];	/* Done with j -> Restore p[j] */
	  }
	  p[i] = c[i];		/* Done with i -> Restore p[i] */
     }
     return sum;
}

static unsigned evalR0_0fs4d(integrand f, void *fdata, unsigned dim, double *p, const double *c, double *sum0_, const double *r1, double *sum1_, const double *r2, double *sum2_)
{
     double maxdiff = 0;
     unsigned i, dimDiffMax = 0;
     double sum0, sum1 = 0, sum2 = 0; /* copies for aliasing, performance */

     double ratio = r1[0] / r2[0];

     ratio *= ratio;
     sum0 = f(dim, p, fdata);

     for (i = 0; i < dim; i++) {
	  double f1a, f1b, f2a, f2b, diff;

	  p[i] = c[i] - r1[i];
	  sum1 += (f1a = f(dim, p, fdata));
	  p[i] = c[i] + r1[i];
	  sum1 += (f1b = f(dim, p, fdata));
	  p[i] = c[i] - r2[i];
	  sum2 += (f2a = f(dim, p, fdata));
	  p[i] = c[i] + r2[i];
	  sum2 += (f2b = f(dim, p, fdata));
	  p[i] = c[i];

	  diff = fabs(f1a + f1b - 2 * sum0 - ratio * (f2a + f2b - 2 * sum0));

	  if (diff > maxdiff) {
	       maxdiff = diff;
	       dimDiffMax = i;
	  }
     }

     *sum0_ += sum0;
     *sum1_ += sum1;
     *sum2_ += sum2;

     return dimDiffMax;
}

#define num0_0(dim) (1U)
#define numR0_0fs(dim) (2 * (dim))
#define numRR0_0fs(dim) (2 * (dim) * (dim-1))
#define numR_Rfs(dim) (1U << (dim))

/***************************************************************************/
/* Based on rule75genzmalik.cpp in HIntLib-0.0.10: An embedded
   cubature rule of degree 7 (embedded rule degree 5) due to A. C. Genz
   and A. A. Malik.  See:

         A. C. Genz and A. A. Malik, "An imbedded [sic] family of fully
         symmetric numerical integration rules," SIAM
         J. Numer. Anal. 20 (3), 580-588 (1983).
*/

typedef struct {
     rule parent;

     /* temporary arrays of length dim */
     double *widthLambda, *widthLambda2, *p;

     /* dimension-dependent constants */
     double weight1, weight3, weight5;
     double weightE1, weightE3;
} rule75genzmalik;

#define real(x) ((double)(x))
#define to_int(n) ((int)(n))

static int isqr(int x)
{
     return x * x;
}

static void destroy_rule75genzmalik(rule *r_)
{
     rule75genzmalik *r = (rule75genzmalik *) r_;
     free(r->p);
}

static unsigned rule75genzmalik_evalError(rule *r_, integrand f, void *fdata, const hypercube *h, esterr *ee)
{
     /* lambda2 = sqrt(9/70), lambda4 = sqrt(9/10), lambda5 = sqrt(9/19) */
     const double lambda2 = 0.3585685828003180919906451539079374954541;
     const double lambda4 = 0.9486832980505137995996680633298155601160;
     const double lambda5 = 0.6882472016116852977216287342936235251269;
     const double weight2 = 980. / 6561.;
     const double weight4 = 200. / 19683.;
     const double weightE2 = 245. / 486.;
     const double weightE4 = 25. / 729.;

     rule75genzmalik *r = (rule75genzmalik *) r_;
     unsigned i, dimDiffMax, dim = r_->dim;
     double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4, sum5, result, res5th;
     const double *center = h->data;
     const double *halfwidth = h->data + dim;

     for (i = 0; i < dim; ++i)
	  r->p[i] = center[i];

     for (i = 0; i < dim; ++i)
	  r->widthLambda2[i] = halfwidth[i] * lambda2;
     for (i = 0; i < dim; ++i)
	  r->widthLambda[i] = halfwidth[i] * lambda4;

     /* Evaluate function in the center, in f(lambda2,0,...,0) and
        f(lambda3=lambda4, 0,...,0).  Estimate dimension with largest error */
     dimDiffMax = evalR0_0fs4d(f, fdata, dim, r->p, center, &sum1, r->widthLambda2, &sum2, r->widthLambda, &sum3);

     /* Calculate sum4 for f(lambda4, lambda4, 0, ...,0) */
     sum4 = evalRR0_0fs(f, fdata, dim, r->p, center, r->widthLambda);

     /* Calculate sum5 for f(lambda5, lambda5, ..., lambda5) */
     for (i = 0; i < dim; ++i)
	  r->widthLambda[i] = halfwidth[i] * lambda5;
     sum5 = evalR_Rfs(f, fdata, dim, r->p, center, r->widthLambda);

     /* Calculate fifth and seventh order results */

     result = h->vol * (r->weight1 * sum1 + weight2 * sum2 + r->weight3 * sum3 + weight4 * sum4 + r->weight5 * sum5);
     res5th = h->vol * (r->weightE1 * sum1 + weightE2 * sum2 + r->weightE3 * sum3 + weightE4 * sum4);

     ee->val = result;
     ee->err = fabs(res5th - result);

     return dimDiffMax;
}

static rule *make_rule75genzmalik(unsigned dim)
{
     rule75genzmalik *r;

     if (dim < 2) return 0; /* this rule does not support 1d integrals */

     /* Because of the use of a bit-field in evalR_Rfs, we are limited
	to be < 32 dimensions (or however many bits are in unsigned).
	This is not a practical limitation...long before you reach
	32 dimensions, the Genz-Malik cubature becomes excruciatingly
	slow and is superseded by other methods (e.g. Monte-Carlo). */
     if (dim >= sizeof(unsigned) * 8) return 0;

     r = (rule75genzmalik *) malloc(sizeof(rule75genzmalik));
     r->parent.dim = dim;

     r->weight1 = (real(12824 - 9120 * to_int(dim) + 400 * isqr(to_int(dim)))
		   / real(19683));
     r->weight3 = real(1820 - 400 * to_int(dim)) / real(19683);
     r->weight5 = real(6859) / real(19683) / real(1U << dim);
     r->weightE1 = (real(729 - 950 * to_int(dim) + 50 * isqr(to_int(dim)))
		    / real(729));
     r->weightE3 = real(265 - 100 * to_int(dim)) / real(1458);

     r->p = (double *) malloc(sizeof(double) * dim * 3);
     r->widthLambda = r->p + dim;
     r->widthLambda2 = r->p + 2 * dim;

     r->parent.num_points = num0_0(dim) + 2 * numR0_0fs(dim)
	  + numRR0_0fs(dim) + numR_Rfs(dim);

     r->parent.evalError = rule75genzmalik_evalError;
     r->parent.destroy = destroy_rule75genzmalik;

     return (rule *) r;
}

/***************************************************************************/
/* 1d 15-point Gaussian quadrature rule, based on qk15.c and qk.c in
   GNU GSL (which in turn is based on QUADPACK). */

static unsigned rule15gauss_evalError(rule *r, integrand f, void *fdata,
				      const hypercube *h, esterr *ee)
{
     /* Gauss quadrature weights and kronrod quadrature abscissae and
	weights as evaluated with 80 decimal digit arithmetic by
	L. W. Fullerton, Bell Labs, Nov. 1981. */
     const unsigned n = 8;
     const double xgk[8] = {  /* abscissae of the 15-point kronrod rule */
	  0.991455371120812639206854697526329,
	  0.949107912342758524526189684047851,
	  0.864864423359769072789712788640926,
	  0.741531185599394439863864773280788,
	  0.586087235467691130294144838258730,
	  0.405845151377397166906606412076961,
	  0.207784955007898467600689403773245,
	  0.000000000000000000000000000000000
	  /* xgk[1], xgk[3], ... abscissae of the 7-point gauss rule. 
	     xgk[0], xgk[2], ... to optimally extend the 7-point gauss rule */
     };
     static const double wg[4] = {  /* weights of the 7-point gauss rule */
	  0.129484966168869693270611432679082,
	  0.279705391489276667901467771423780,
	  0.381830050505118944950369775488975,
	  0.417959183673469387755102040816327
     };
     static const double wgk[8] = { /* weights of the 15-point kronrod rule */
	  0.022935322010529224963732008058970,
	  0.063092092629978553290700663189204,
	  0.104790010322250183839876322541518,
	  0.140653259715525918745189590510238,
	  0.169004726639267902826583426598550,
	  0.190350578064785409913256402421014,
	  0.204432940075298892414161999234649,
	  0.209482141084727828012999174891714
     };

     const double center = h->data[0];
     const double halfwidth = h->data[1];
     double fv1[7], fv2[7];
     const double f_center = f(1, &center, fdata);
     double result_gauss = f_center * wg[n/2 - 1];
     double result_kronrod = f_center * wgk[n - 1];
     double result_abs = fabs(result_kronrod);
     double result_asc, mean, err;
     unsigned j;

     for (j = 0; j < (n - 1) / 2; ++j) {
	  int j2 = 2*j + 1;
	  double x, f1, f2, fsum, w = halfwidth * xgk[j2];
	  x = center - w; fv1[j2] = f1 = f(1, &x, fdata);
	  x = center + w; fv2[j2] = f2 = f(1, &x, fdata);
	  fsum = f1 + f2;
	  result_gauss += wg[j] * fsum;
	  result_kronrod += wgk[j2] * fsum;
	  result_abs += wgk[j2] * (fabs(f1) + fabs(f2));
     }

     for (j = 0; j < n/2; ++j) {
	  int j2 = 2*j;
	  double x, f1, f2, w = halfwidth * xgk[j2];
	  x = center - w; fv1[j2] = f1 = f(1, &x, fdata);
          x = center + w; fv2[j2] = f2 = f(1, &x, fdata);
          result_kronrod += wgk[j2] * (f1 + f2);
          result_abs += wgk[j2] * (fabs(f1) + fabs(f2));
     }

     ee->val = result_kronrod * halfwidth;

     /* compute error estimate: */
     mean = result_kronrod * 0.5;
     result_asc = wgk[n - 1] * fabs(f_center - mean);
     for (j = 0; j < n - 1; ++j)
	  result_asc += wgk[j] * (fabs(fv1[j]-mean) + fabs(fv2[j]-mean));
     err = fabs(result_kronrod - result_gauss) * halfwidth;
     result_abs *= halfwidth;
     result_asc *= halfwidth;
     if (result_asc != 0 && err != 0) {
	  double scale = pow((200 * err / result_asc), 1.5);
	  if (scale < 1)
	       err = result_asc * scale;
	  else
	       err = result_asc;
     }
     if (result_abs > DBL_MIN / (50 * DBL_EPSILON)) {
	  double min_err = 50 * DBL_EPSILON * result_abs;
	  if (min_err > err)
	       err = min_err;
     }
     ee->err = err;
     
     return 0; /* no choice but to divide 0th dimension */
}

static rule *make_rule15gauss(unsigned dim)
{
     rule *r;
     if (dim != 1) return 0; /* this rule is only for 1d integrals */
     r = (rule *) malloc(sizeof(rule));
     r->dim = dim;
     r->num_points = 15;
     r->evalError = rule15gauss_evalError;
     r->destroy = 0;
     return r;
}

/***************************************************************************/
/* binary heap implementation (ala _Introduction to Algorithms_ by
   Cormen, Leiserson, and Rivest), for use as a priority queue of
   regions to integrate. */

typedef region heap_item;
#define KEY(hi) ((hi).ee.err)

typedef struct {
     unsigned n, nalloc;
     heap_item *items;
     esterr ee;
} heap;

static void heap_resize(heap *h, unsigned nalloc)
{
     h->nalloc = nalloc;
     h->items = (heap_item *) realloc(h->items, sizeof(heap_item) * nalloc);
}

static heap heap_alloc(unsigned nalloc)
{
     heap h;
     h.n = 0;
     h.nalloc = 0;
     h.items = 0;
     h.ee.val = h.ee.err = 0;
     heap_resize(&h, nalloc);
     return h;
}

/* note that heap_free does not deallocate anything referenced by the items */
static void heap_free(heap *h)
{
     h->n = 0;
     heap_resize(h, 0);
}

static void heap_push(heap *h, heap_item hi)
{
     int insert;

     h->ee.val += hi.ee.val;
     h->ee.err += hi.ee.err;
     insert = h->n;
     if (++(h->n) > h->nalloc)
	  heap_resize(h, h->n * 2);

     while (insert) {
	  int parent = (insert - 1) / 2;
	  if (KEY(hi) <= KEY(h->items[parent]))
	       break;
	  h->items[insert] = h->items[parent];
	  insert = parent;
     }
     h->items[insert] = hi;
}

static heap_item heap_pop(heap *h)
{
     heap_item ret;
     int i, n, child;

     if (!(h->n)) {
	  fprintf(stderr, "attempted to pop an empty heap\n");
	  exit(EXIT_FAILURE);
     }

     ret = h->items[0];
     h->items[i = 0] = h->items[n = --(h->n)];
     while ((child = i * 2 + 1) < n) {
	  int largest;
	  heap_item swap;

	  if (KEY(h->items[child]) <= KEY(h->items[i]))
	       largest = i;
	  else
	       largest = child;
	  if (++child < n && KEY(h->items[largest]) < KEY(h->items[child]))
	       largest = child;
	  if (largest == i)
	       break;
	  swap = h->items[i];
	  h->items[i] = h->items[largest];
	  h->items[i = largest] = swap;
     }

     h->ee.val -= ret.ee.val;
     h->ee.err -= ret.ee.err;
     return ret;
}

/***************************************************************************/

/* adaptive integration, analogous to adaptintegrator.cpp in HIntLib */

static int ruleadapt_integrate(rule *r, integrand f, void *fdata, const hypercube *h, unsigned maxEval, double reqAbsError, double reqRelError, esterr *ee)
{
     unsigned maxIter;		/* maximum number of adaptive subdivisions */
     heap regions;
     unsigned i;
     int status = -1; /* = ERROR */

     if (maxEval) {
	  if (r->num_points > maxEval)
	       return status; /* ERROR */
	  maxIter = (maxEval - r->num_points) / (2 * r->num_points);
     }
     else
	  maxIter = UINT_MAX;

     regions = heap_alloc(1);

     heap_push(&regions, eval_region(make_region(h), f, fdata, r));
     /* another possibility is to specify some non-adaptive subdivisions: 
	if (initialRegions != 1)
	   partition(h, initialRegions, EQUIDISTANT, &regions, f,fdata, r); */

     for (i = 0; i < maxIter; ++i) {
	  region R, R2;
	  if (regions.ee.err <= reqAbsError 
	      || relError(regions.ee) <= reqRelError) {
	       status = 0; /* converged! */
	       break;
	  }
	  R = heap_pop(&regions); /* get worst region */
	  cut_region(&R, &R2);
	  heap_push(&regions, eval_region(R, f, fdata, r));
	  heap_push(&regions, eval_region(R2, f, fdata, r));
     }

     ee->val = ee->err = 0;  /* re-sum integral and errors */
     for (i = 0; i < regions.n; ++i) {
	  ee->val += regions.items[i].ee.val;
	  ee->err += regions.items[i].ee.err;
	  destroy_region(&regions.items[i]);
     }
     /* printf("regions.nalloc = %d\n", regions.nalloc); */
     heap_free(&regions);

     return status;
}

int adapt_integrate(integrand f, void *fdata, 
		    unsigned dim, const double *xmin, const double *xmax, 
		    unsigned maxEval, double reqAbsError, double reqRelError, 
		    double *val, double *err)
{
     rule *r;
     hypercube h;
     esterr ee;
     int status;
     
     if (dim == 0) { /* trivial integration */
	  *val = f(0, xmin, fdata);
	  *err = 0;
	  return 0;
     }
     r = dim == 1 ? make_rule15gauss(dim) : make_rule75genzmalik(dim);
     if (!r) { *val = 0; *err = HUGE_VAL; return -2; /* ERROR */ }
     h = make_hypercube_range(dim, xmin, xmax);
     status = ruleadapt_integrate(r, f, fdata, &h,
				  maxEval, reqAbsError, reqRelError,
				  &ee);
     *val = ee.val;
     *err = ee.err;
     destroy_hypercube(&h);
     destroy_rule(r);
     return status;
}

double adaptive_integration(integrand f, double *xmin, double *xmax,
			    int n, void *fdata,
			    double abstol, double reltol, int maxnfe,
			    double *esterr, int *errflag)
{
     double val;
     *errflag = adapt_integrate(f, fdata, n, xmin, xmax,
				maxnfe, abstol, reltol, &val, esterr);
     return val;
}
