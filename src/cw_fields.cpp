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

#include "meep_internals.hpp"
#include "bicgstab.hpp"

#include <vector>

using namespace std;

namespace meep {

/* The unknowns of the frequency-domain linear system are D and B, but in PML regions
   the U, W, E/H, C (conductivity), and bfast auxiliary fields are also unknowns (in
   principle, we might be able to compute these extra fields in frequency domain via
   scaling by the appropriate s factors, rather than storing them).
   for_each_cw_unknown calls fn(gv, c, re, im) for the real and imaginary arrays of
   each unknown field on this process, in a fixed order; the system size, fields_to_array,
   and array_to_fields all use this one traversal so that they are guaranteed to agree. */
template <typename F> static void for_each_cw_unknown(const fields &f, F fn) {
  for (int i = 0; i < f.num_chunks; i++) {
    fields_chunk *fc = f.chunks[i];
    if (!fc->is_mine()) continue;
    FOR_COMPONENTS(c) {
      if (!is_D(c) && !is_B(c)) continue;
      const component c2 = field_type_component(is_D(c) ? E_stuff : H_stuff, c);
      auto visit = [&](realnum *const fld[2]) {
        if (fld[0] && fld[1]) fn(fc->gv, c, fld[0], fld[1]);
      };
      visit(fc->f[c]);
      visit(fc->f_u[c]);
      visit(fc->f_cond[c]);
      visit(fc->f_bfast[c]);
      visit(fc->f_w[c2]);
      if (fc->f_w[c2][0]) visit(fc->f[c2]);
    }
  }
}

// number of complex unknowns on this process
static size_t cw_system_size(const fields &f) {
  size_t n = 0;
  for_each_cw_unknown(
      f, [&](const grid_volume &gv, component c, realnum *, realnum *) { n += gv.nowned(c); });
  return n;
}

/* Calls body(idx, offset) for each point owned by gv for component c (in parallel), where
   idx is the index of the point in the field arrays and offset is its index in the
   (contiguous) sequence of owned points. Uses the stride-1 loop when possible so that the
   inner loop is contiguous and can be vectorized. */
#define CW_LOOP_OFFSET ((loop_i1 * loop_n2 + loop_i2) * loop_n3 + loop_i3)
template <typename F> static void loop_over_owned(const grid_volume &gv, component c, F body) {
  if (LOOPS_ARE_STRIDE1(gv)) {
    PS1LOOP_OVER_VOL_OWNED(gv, c, idx) { body(idx, CW_LOOP_OFFSET); }
  }
  else {
    PLOOP_OVER_VOL_OWNED(gv, c, idx) { body(idx, CW_LOOP_OFFSET); }
  }
}
#undef CW_LOOP_OFFSET

static void fields_to_array(const fields &f, complex<realnum> *x) {
  size_t ix = 0;
  for_each_cw_unknown(
      f, [&](const grid_volume &gv, component c, const realnum *fr, const realnum *fi) {
        complex<realnum> *xc = x + ix;
        loop_over_owned(
            gv, c, [=](ptrdiff_t idx, ptrdiff_t k) { xc[k] = complex<realnum>(fr[idx], fi[idx]); });
        ix += gv.nowned(c);
      });
}

static void array_to_fields(const complex<realnum> *x, fields &f) {
  size_t ix = 0;
  for_each_cw_unknown(f, [&](const grid_volume &gv, component c, realnum *fr, realnum *fi) {
    const complex<realnum> *xc = x + ix;
    loop_over_owned(gv, c, [=](ptrdiff_t idx, ptrdiff_t k) {
      fr[idx] = real(xc[k]);
      fi[idx] = imag(xc[k]);
    });
    ix += gv.nowned(c);
  });

  f.step_boundaries(D_stuff);
  f.update_eh(E_stuff, true);
  f.step_boundaries(E_stuff);

  // H is computed from B by f.step() before updating D, or at the end of solve_cw
}

typedef struct {
  size_t n;
  fields *f;
  complex<double> iomega;
} fieldop_data;

// y = A x, where A = (step - 1)/dt + i*omega, applied to complex vectors of length n.
static void fieldop(const realnum *xr, realnum *yr, void *data_) {
  fieldop_data *data = (fieldop_data *)data_;
  array_to_fields(reinterpret_cast<const complex<realnum> *>(xr), *data->f);
  data->f->step();
  fields_to_array(*data->f, reinterpret_cast<complex<realnum> *>(yr));
  const size_t n = data->n;
  const double dt_inv = 1.0 / data->f->dt;
  const double wr = real(data->iomega), wi = imag(data->iomega);
#pragma omp parallel for simd
  for (size_t i = 0; i < n; ++i) {
    const double a = xr[2 * i], b = xr[2 * i + 1];
    yr[2 * i] = (yr[2 * i] - a) * dt_inv + (wr * a - wi * b);
    yr[2 * i + 1] = (yr[2 * i + 1] - b) * dt_inv + (wr * b + wi * a);
  }
}

// Rayleigh-quotient estimate <x,Ax>/<x,x> for eigenfrequency given approximate eigenvector x
// (length n), overwriting x with Ax and b with x/|x|.
static complex<double> estimate_eigfreq(complex<realnum> *b, complex<realnum> *x, size_t n,
                                        fieldop_data *data) {
  memcpy(b, x, n * sizeof(complex<realnum>));
  fieldop(reinterpret_cast<realnum *>(b), reinterpret_cast<realnum *>(x), (void *)data);
  complex<double> bdotx(0, 0);
  double bnorm2 = 0;
  for (size_t i = 0; i < n; ++i) {
    complex<realnum> bi = b[i];
    bnorm2 += real(bi) * real(bi) + imag(bi) * imag(bi);
    complex<realnum> bx = conj(bi) * x[i];
    bdotx += complex<double>(real(bx), imag(bx));
  }
  bnorm2 = sum_to_all(bnorm2);
  bdotx = sum_to_all(bdotx);
  double bnorminv = 1 / sqrt(bnorm2);
  for (size_t i = 0; i < n; ++i) {
    b[i] *= bnorminv; // normalize b for subsequent shift-and-invert iterations
  }
  complex<double> iomega = data->iomega - bdotx / bnorm2; // unshifted eigenvalue
  // now, invert: iomega = (1 - exp(-i * (2 * pi * frequency) * dt)) / dt)
  // to get frequency = log(1 - iomega * dt) / (-2 pi i * dt)
  double dt = data->f->dt;
  return log(1.0 - iomega * dt) / complex<double>(0, -2 * pi * dt);
}

/* Solve for the CW (constant frequency) field response at the given
   frequency to the sources (with amplitude given by the current sources
   at the current time).  The solver halts at a fractional convergence
   of tol, or when maxiters is reached, or when convergence fails;
   returns true if convergence succeeds and false if it fails.

   The parameter L determines the order of the iterative algorithm
   that is used.  L should always be positive and should normally be
   >= 2.  Larger values of L will often lead to faster convergence, at
   the expense of more memory and more work per iteration.

   If the optional argument eigfreq is non-NULL, then the solver is used for a
   shift-and-invert power iteration to find the closest eigenfrequency and
   eigenvector to frequency: the solver is iterated up to eigiters times,
   or until the estimated eigenfreq stops changing by <= eigtol (relative).
   In this case, the return value is true only if the eigenfrequency
   converged (failures of individual solves are reported but tolerated,
   since shift-and-invert iteration tolerates inexact solves.) */
namespace {
/* While solve_cw is running, the fields are in a special "solve_cw" state (sources and
   DFT updates are disabled in step(), and step() increments t).  This guard restores the
   time step and normal time-domain operation when solve_cw exits, including when it exits
   via an exception (e.g. meep::abort, which raises a RuntimeError in Python), so that the
   fields are not left with sources and DFTs silently disabled. */
class solve_cw_guard {
public:
  solve_cw_guard(fields &f, complex<double> omega) : f(f), tsave(f.t), active(true) {
    f.set_solve_cw_omega(omega);
  }
  ~solve_cw_guard() { restore(); }
  void restore() {
    if (active) {
      f.t = tsave;
      f.unset_solve_cw_omega();
      active = false;
    }
  }
  fields &f;
  const int tsave; // time step at which solve_cw was called
private:
  bool active;
};
} // namespace

bool fields::solve_cw(double tol, int maxiters, complex<double> frequency, int L,
                      complex<double> *eigfreq, double eigtol, int eigiters) {
  if (is_real) meep::abort("solve_cw is incompatible with use_real_fields()");
  if (L < 1) meep::abort("solve_cw called with L = %d < 1", L);

  solve_cw_guard guard(*this, 2 * pi * frequency);

  step(); // step once to make sure everything is allocated
  /* restore the time, which step() incremented, so that the source amplitudes below are
     those of the time at which solve_cw was called, and the solution therefore corresponds
     to the fields at that time */
  t = guard.tsave;

  // size of linear system (on this processor, at least), in real numbers
  const size_t N = 2 * cw_system_size(*this);

  int iters = maxiters;
  const size_t nwork = (size_t)bicgstabL(L, N, 0, 0, 0, 0, tol, &iters, 0, true);
  std::vector<realnum> work(nwork + 2 * N);
  complex<realnum> *x = reinterpret_cast<complex<realnum> *>(work.data() + nwork);
  complex<realnum> *b = reinterpret_cast<complex<realnum> *>(work.data() + nwork + N);

  fields_to_array(*this, x); // initial guess = initial fields

  // get J amplitudes from current time step
  zero_fields(); // note that we've saved the fields in x above
  calc_sources(time());
  step_source(B_stuff, true);
  step_boundaries(B_stuff);
  update_eh(H_stuff);
  calc_sources(time() + 0.5 * dt);
  step_source(D_stuff, true);
  step_boundaries(D_stuff);
  update_eh(E_stuff);
  fields_to_array(*this, b);
  double mdt_inv = -1.0 / dt;
  for (size_t i = 0; i < N / 2; ++i)
    b[i] *= mdt_inv;
  {
    double bmax = 0;
    for (size_t i = 0; i < N / 2; ++i) {
      double babs = abs(b[i]);
      if (babs > bmax) bmax = babs;
    }
    am_now_working_on(MpiAllTime);
    if (max_to_all(bmax) == 0.0) meep::abort("zero current amplitudes in solve_cw");
    finished_working();
  }

  fieldop_data data;
  data.f = this;
  data.n = N / 2;
  data.iomega = ((1.0 - exp(complex<double>(0., -1.) * (2 * pi * frequency) * dt)) * (1.0 / dt));

  // solve A x = b, starting from the current x; returns error code from bicgstabL;
  auto solve = [&]() {
    iters = maxiters;
    return (int)bicgstabL(L, N, reinterpret_cast<realnum *>(x), fieldop, &data,
                          reinterpret_cast<realnum *>(b), tol, &iters, work.data(), verbosity == 0);
  };

  auto report_failure = [&](int ierr) {
    if (ierr && verbosity > 0)
      master_printf(" -- CONVERGENCE FAILURE (%d) in solve_cw!%s\n", ierr,
                    ierr == 2 ? " (residual diverged: tol may be below the floating-point "
                                "precision floor; using the best solution found)"
                              : "");
  };

  int ierr = solve();
  bool success = !ierr;

  if (verbosity > 0)
    master_printf("Finished solve_cw after %d CG iters (~ %d timesteps).\n", iters, iters * 2 * L);
  report_failure(ierr);

  // do additional shift-and-invert iterations to find eigenfrequency
  if (eigfreq) {
    /* The first solve only provides the starting vector for the shift-and-invert iteration,
       so its accuracy is not critical, and we continue even if it failed. */
    *eigfreq = estimate_eigfreq(b, x, data.n, &data);
    if (verbosity > 0) {
      master_printf("Initial eigen-frequency estimate = %g%+gi\n", real(*eigfreq), imag(*eigfreq));
    }
    bool converged = false;
    for (int eigiter = 0; eigiter < eigiters; ++eigiter) {
      ierr = solve();
      // (if the solve failed, x is its best iterate, so this estimate is still meaningful)
      complex<double> newfreq = estimate_eigfreq(b, x, data.n, &data);
      complex<double> dfreq = newfreq - *eigfreq;
      if (verbosity > 0)
        master_printf("Eigensolver step %d: %d CG iters, freq = %g%+gi (change = %g%+gi).\n",
                      eigiter + 1, iters, real(newfreq), imag(newfreq), real(dfreq), imag(dfreq));
      report_failure(ierr);
      *eigfreq = newfreq;
      /* If the solve failed, we continue anyway: shift-and-invert iteration tolerates inexact
         solves, and x is the best iterate found.  Failures are reported above, and the
         eigenfrequency convergence test below still determines success. */
      if (abs(dfreq) <= eigtol * abs(newfreq)) {
        converged = true;
        break;
      }
    }
    if (!converged && verbosity > 0)
      master_printf(" --- eigensolver did not converge within %d iterations\n", eigiters);
    success = converged;
    memcpy(x, b, N * sizeof(realnum));
  }

  array_to_fields(x, *this); // sets D, B, auxiliary fields, and E

  // compute H from B, as step() would do before updating D (without advancing the time)
  step_boundaries(B_stuff);
  update_eh(H_stuff, true);
  step_boundaries(H_stuff);

  guard.restore(); // re-enable DFT updates, which are skipped while doing_solve_cw
  update_dfts();

  return success;
}

/* as solve_cw, but infers frequency from sources */
bool fields::solve_cw(double tol, int maxiters, int L, complex<double> *eigfreq, double eigtol,
                      int eigiters) {
  complex<double> freq = 0.0;
  for (src_time *s = sources; s; s = s->next) {
    complex<double> sf = s->frequency();
    if (sf != freq && freq != 0.0 && sf != 0.0)
      meep::abort("must pass frequency to solve_cw if sources do not agree");
    if (sf != 0.0) freq = sf;
  }
  if (freq == 0.0) meep::abort("must pass frequency to solve_cw if sources do not specify one");
  return solve_cw(tol, maxiters, freq, L, eigfreq, eigtol, eigiters);
}

} // namespace meep
