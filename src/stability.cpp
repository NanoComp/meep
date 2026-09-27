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

#include <cmath>
#include <complex>
#include <vector>

#include "meep.hpp"
#include "stability.hpp"

namespace meep {

/* Diagnose unstable dispersive timesteps (issue #12).

   E, H and the polarizations are stepped together, so a Lorentzian whose
   own recurrence is stable can still break the Courant condition; checking
   each pole on its own is not enough.  This analyses the coupled step, for a
   uniform bulk medium only; it says nothing about interfaces, PML or sources.

   For a single plane wave, one timestep is a small matrix acting on E, H
   and the current and previous value of each polarization, built from the
   same coefficients as update_P.  The step is stable for that wave if no
   eigenvalue of the matrix has modulus above 1.

   The wavevector enters the matrix through a single number, which is largest
   at k*dx = pi on every axis.  For the media medium_is_modelled accepts, an
   eigenvalue can only leave the unit circle as that number grows, never come
   back in, so the largest wavevector decides whether any wavevector is
   unstable, and it is the only one evaluated.  The growth rate itself can
   peak at a smaller wavevector, so it is not reported. */

namespace {

// Row-major n x n complex matrix.
typedef std::vector<std::complex<double> > cmatrix;

// Matrix infinity norm (max row sum, operator norm).
double matrix_norm(const cmatrix &A, int n) {
  double best = 0;
  for (int i = 0; i < n; ++i) {
    double s = 0;
    for (int j = 0; j < n; ++j)
      s += std::abs(A[i * n + j]);
    if (s > best) best = s;
  }
  return best;
}

void matrix_square(cmatrix &A, int n) {
  cmatrix B(n * n, std::complex<double>(0, 0));
  for (int i = 0; i < n; ++i)
    for (int k = 0; k < n; ++k)
      for (int j = 0; j < n; ++j)
        B[i * n + j] += A[i * n + k] * A[k * n + j];
  A.swap(B);
}

/* Spectral radius via Gelfand's formula, |A^m|^(1/m) with m = 2^30, by
   repeated squaring with renormalization.  It needs no eigensolver, copes
   with defective eigenvalues, and never falls below rho, since
   rho(A)^m <= |A^m| for an operator norm.  The overshoot is a few times
   1e-6 at worst, for the undamped DC mode of a Drude term. */
double spectral_radius(cmatrix A, int n) {
  const int nsquare = 30;
  double log_rho = 0;
  for (int j = 1; j <= nsquare; ++j) {
    matrix_square(A, n);
    const double nrm = matrix_norm(A, n);
    if (nrm != nrm) return infinity; // NaN: report as unstable, never as safe
    if (!(nrm > 0)) return 0;
    for (int i = 0; i < n * n; ++i)
      A[i] /= nrm;
    log_rho += log(nrm) / ldexp(1.0, j);
  }
  return exp(log_rho);
}

} // namespace

double dispersive_spectral_radius(double eps_inf, double mu_inf,
                                  const std::vector<lorentzian_pole> &poles, double dt, double dx,
                                  int ndim) {
  const int N = int(poles.size()), n = 2 + 2 * N;
  const double q = 2 * sqrt(double(ndim)) * dt / dx;

  cmatrix M(n * n, std::complex<double>(0, 0));
  double csum = 0;
  for (int j = 0; j < N; ++j) { // as in lorentzian_susceptibility::update_P
    const double w = 2 * pi * poles[j].omega_0, g = 2 * pi * poles[j].gamma;
    const double w2dt2 = (w * dt) * (w * dt);
    const double gamma1inv = 1 / (1 + g * dt / 2), gamma1 = 1 - g * dt / 2;
    const double a = gamma1inv * (2 - (poles[j].no_omega_0_denominator ? 0 : w2dt2));
    const double b = -gamma1inv * gamma1;
    const double c = gamma1inv * w2dt2 * poles[j].sigma;
    csum += c;
    M[0 * n + (2 + 2 * j)] = (1 - a) / eps_inf;
    M[0 * n + (3 + 2 * j)] = -b / eps_inf;
    M[(2 + 2 * j) * n + 0] = c;
    M[(2 + 2 * j) * n + (2 + 2 * j)] = a;
    M[(2 + 2 * j) * n + (3 + 2 * j)] = b;
    M[(3 + 2 * j) * n + (2 + 2 * j)] = 1;
  }
  // E sees the H just updated from it, hence the q^2/(eps*mu) term.
  M[0 * n + 0] = (eps_inf - csum - q * q / mu_inf) / eps_inf;
  M[0 * n + 1] = std::complex<double>(0, -q) / eps_inf;
  M[1 * n + 0] = std::complex<double>(0, -q) / mu_inf;
  M[1 * n + 1] = 1;
  return spectral_radius(M, n);
}

} // namespace meep

namespace meep_geom {

static bool is_isotropic(vector3 v) { return v.x == v.y && v.y == v.z; }
static bool is_zero(vector3 v) { return v.x == 0 && v.y == 0 && v.z == 0; }
static bool is_zero(cvector3 v) {
  return v.x.re == 0 && v.y.re == 0 && v.z.re == 0 && v.x.im == 0 && v.y.im == 0 && v.z.im == 0;
}

/* An isotropic Lorentzian or Drude term with sigma, gamma >= 0.  Noise only
   adds a source, so a noisy Lorentzian qualifies. */
static bool is_modelled(const susceptibility &ss) {
  return !ss.is_file && !ss.saturated_gyrotropy && is_zero(ss.bias) && ss.transitions.empty() &&
         ss.initial_populations.empty() && is_zero(ss.sigma_offdiag) &&
         is_isotropic(ss.sigma_diag) && ss.sigma_diag.x >= 0 && ss.gamma >= 0;
}

static bool all_modelled(const susceptibility_list &suscs) {
  for (size_t i = 0; i < suscs.size(); ++i)
    if (!is_modelled(suscs[i])) return false;
  return true;
}

/* True if the matrix represents this medium on this grid: a Cartesian grid,
   isotropic eps_inf, mu_inf > 0 (anisotropy couples components for oblique
   k), and passive Lorentzian/Drude terms on E or on H only, with no
   conductivity or nonlinearity.  The out-of-plane wavevector of a 2D run
   (fields::beta) is not known yet at this point and is taken to be zero. */
static bool medium_is_modelled(const medium_struct *mm, const meep::grid_volume &gv) {
  const bool cartesian = gv.dim == meep::D1 || gv.dim == meep::D2 || gv.dim == meep::D3;
  const bool electric = !mm->E_susceptibilities.empty();
  const bool magnetic = !mm->H_susceptibilities.empty();
  return cartesian && electric != magnetic && all_modelled(mm->E_susceptibilities) &&
         all_modelled(mm->H_susceptibilities) && is_isotropic(mm->epsilon_diag) &&
         is_isotropic(mm->mu_diag) && is_zero(mm->epsilon_offdiag) && is_zero(mm->mu_offdiag) &&
         mm->epsilon_diag.x > 0 && mm->mu_diag.x > 0 && is_zero(mm->D_conductivity_diag) &&
         is_zero(mm->B_conductivity_diag) && is_zero(mm->E_chi2_diag) && is_zero(mm->E_chi3_diag) &&
         is_zero(mm->H_chi2_diag) && is_zero(mm->H_chi3_diag);
}

stability_status medium_stability(const medium_struct *mm, const meep::grid_volume &gv, double dt) {
  if (!medium_is_modelled(mm, gv)) return STABILITY_UNCHECKED;

  // Well above the estimator's error, and still e^10 of growth in 1e5 steps.
  const double margin = 1e-4;

  // A magnetic medium runs the same recursion with eps and mu swapped.
  const bool electric = !mm->E_susceptibilities.empty();
  const susceptibility_list &suscs = electric ? mm->E_susceptibilities : mm->H_susceptibilities;
  std::vector<meep::lorentzian_pole> poles;
  for (size_t i = 0; i < suscs.size(); ++i) {
    const meep::lorentzian_pole p = {suscs[i].frequency, suscs[i].gamma, suscs[i].sigma_diag.x,
                                     suscs[i].drude};
    if (p.sigma != 0) poles.push_back(p);
  }
  const double eps_inf = electric ? mm->epsilon_diag.x : mm->mu_diag.x;
  const double mu_inf = electric ? mm->mu_diag.x : mm->epsilon_diag.x;
  const double rho = meep::dispersive_spectral_radius(eps_inf, mu_inf, poles, dt, 1.0 / gv.a,
                                                      meep::number_of_directions(gv.dim));
  return rho > 1 + margin ? STABILITY_UNSTABLE : STABILITY_STABLE;
}

void check_medium_stability_once(const medium_struct *mm, const meep::grid_volume &gv, double dt,
                                 std::vector<const medium_struct *> &checked) {
  if (mm->E_susceptibilities.empty() && mm->H_susceptibilities.empty()) return;
  for (size_t i = 0; i < checked.size(); ++i)
    if (checked[i] == mm || medium_struct_equal(checked[i], mm)) return;
  checked.push_back(mm);

  switch (medium_stability(mm, gv, dt)) {
    case STABILITY_UNSTABLE:
      meep::master_printf_stderr("warning: dispersive material may be unstable for the current "
                                 "timestep; reduce Courant or increase resolution\n");
      break;
    case STABILITY_UNCHECKED:
      if (meep::verbosity > 1)
        meep::master_printf("dispersive material stability not checked: outside the "
                            "isotropic Cartesian Lorentzian/Drude model\n");
      break;
    case STABILITY_STABLE: break;
  }
}

} // namespace meep_geom
