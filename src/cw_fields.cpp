/* Copyright (C) 2004 Massachusetts Institute of Technology.
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

#include "meep_internals.hpp"
#include "bicgstab.hpp"

namespace meep {

static void fields_to_array(const fields &f, complex<double> *x)
{
  int ix = 0;
  for (int i=0;i<f.num_chunks;i++)
    if (f.chunks[i]->is_mine())
      FOR_COMPONENTS(c)
        if (f.chunks[i]->f[c][0] && (is_D(c) || is_magnetic(c))) {
	  double *fs[3][2];
	  // TODO: this is somewhat wasteful, because at least
	  // one of f/f_p_pml/f_m_pml is redundant.
	  DOCMP2 { fs[0][cmp] = f.chunks[i]->f[c][cmp];
 	           fs[1][cmp] = f.chunks[i]->f_p_pml[c][cmp];
 	           fs[2][cmp] = f.chunks[i]->f_m_pml[c][cmp]; }
	  for (int k = 0; k < 3; ++k)
	    if (fs[k][0]) {
	      LOOP_OVER_VOL_OWNED(f.chunks[i]->v, c, idx)
		x[ix++] = complex<double>(fs[k][0][idx], fs[k][1][idx]);
	    }
        }
}

static void array_to_fields(const complex<double> *x, fields &f)
{
  int ix = 0;
  for (int i=0;i<f.num_chunks;i++)
    if (f.chunks[i]->is_mine())
      FOR_COMPONENTS(c)
        if (f.chunks[i]->f[c][0] && (is_D(c) || is_magnetic(c))) {
	  double *fs[3][2];
	  DOCMP2 { fs[0][cmp] = f.chunks[i]->f[c][cmp];
 	           fs[1][cmp] = f.chunks[i]->f_p_pml[c][cmp];
 	           fs[2][cmp] = f.chunks[i]->f_m_pml[c][cmp]; }
	  for (int k = 0; k < 3; ++k)
	    if (fs[k][0]) {
	      LOOP_OVER_VOL_OWNED(f.chunks[i]->v, c, idx) {
	        fs[k][0][idx] = real(x[ix]);
		fs[k][1][idx] = imag(x[ix]);
		++ix;
	      }
	    }
	}
  f.force_consistency(H_stuff);
  f.force_consistency(D_stuff);
}

typedef struct {
  int n;
  fields *f;
  complex<double> iomega;
  int iters;
} fieldop_data;

static void fieldop(const double *xr, double *yr, void *data_)
{
  const complex<double> *x = reinterpret_cast<const complex<double>*>(xr);
  complex<double> *y = reinterpret_cast<complex<double>*>(yr);
  fieldop_data *data = (fieldop_data *) data_;
  array_to_fields(x, *data->f);
  data->f->step();
  fields_to_array(*data->f, y);
  int n = data->n;
  double dt_inv = 1.0 / data->f->dt;
  complex<double> iomega = data->iomega;
  for (int i = 0; i < n; ++i) y[i] = (y[i] - x[i]) * dt_inv + iomega * x[i];
  data->iters++;
}

/* Solve for the CW (constant frequency) field response at the given
   frequency to the sources (with amplitude given by the current sources
   at the current time).  The solver halts at a fractional convergence
   of tol, or when maxiters is reached, or when convergence fails;
   returns true if convergence succeeds and false if it fails.

   The parameter L determines the order of the iterative algorithm
   that is used.  L should always be positive and should normally be
   >= 2.  Larger values of L will often lead to faster convergence, at
   the expense of more memory and more work per iteration. */
bool fields::solve_cw(double tol, int maxiters, complex<double> frequency,
		      int L) {
  if (is_real) abort("solve_cw is incompatible with use_real_fields()");
  if (L < 1) abort("solve_cw called with L = %d < 1", L);

  int N = 0; // size of linear system (on this processor, at least)
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      FOR_COMPONENTS(c)
	if (chunks[i]->f[c][0] && (is_D(c) || is_magnetic(c))) {
	  N += 2 * chunks[i]->v.nowned(c)
	    * (1 + !!chunks[i]->f_p_pml[c][0] + !!chunks[i]->f_m_pml[c][0]);
	}

  int nwork = bicgstabL(L, N, 0, 0, 0, 0, tol, &maxiters, 0, true);
  double *work = new double[nwork + 2*N];
  complex<double> *x = reinterpret_cast<complex<double>*>(work + nwork);
  complex<double> *b = reinterpret_cast<complex<double>*>(work + nwork + N);

  int tsave = t; // save time (gets incremented by iterations)

  fields_to_array(*this, x); // initial guess = initial fields

  // get J amplitudes from current time step
  zero_fields(); // note that we've saved the fields in x above
  calc_sources(time());
  step_h_source();
  step_boundaries(H_stuff);
  step_d_source();
  step_boundaries(D_stuff);
  fields_to_array(*this, b);
  double mdt_inv = -1.0 / dt;
  for (int i = 0; i < N/2; ++i) b[i] *= mdt_inv;
  {
    double bmax = 0;
    for (int i = 0; i < N/2; ++i) {
      double babs = abs(b[i]);
      if (babs > bmax) bmax = babs;
    }
    if (max_to_all(bmax) == 0.0) abort("zero current amplitudes in solve_cw");
  }

  fieldop_data data;
  data.f = this;
  data.n = N / 2;
  data.iomega = ((1.0 - exp(complex<double>(0.,-1.) * (2*pi*frequency) * dt))
		 * (1.0 / dt));
  data.iters = 0;

  bool save_disable_sources = disable_sources;
  disable_sources = true;

  int ierr = bicgstabL(L, N, reinterpret_cast<double*>(x),
		       fieldop, &data, reinterpret_cast<double*>(b),
		       tol, &maxiters, work, quiet);

  if (!quiet) {
    master_printf("Finished solve_cw after %d steps and %d CG iters.\n", 
		  data.iters, maxiters);
    if (ierr)
      master_printf(" -- CONVERGENCE FAILURE (%d) in solve_cw!\n", ierr);
  }

  array_to_fields(x, *this);
  
  delete[] work;
  t = tsave;
  disable_sources = save_disable_sources;

  return !ierr;
}

/* as solve_cw, but infers frequency from sources */
bool fields::solve_cw(double tol, int maxiters, int L) {
  complex<double> freq = 0.0;
  for (src_time *s = sources; s; s = s->next) {
    complex<double> sf = s->frequency();
    if (sf != freq && freq != 0.0 && sf != 0.0)
      abort("must pass frequency to solve_cw if sources do not agree");
    if (sf != 0.0)
      freq = sf;
  }
  if (freq == 0.0)
    abort("must pass frequency to solve_cw if sources do not specify one");
  return solve_cw(tol, maxiters, freq, L);
}

} // namespace meep
