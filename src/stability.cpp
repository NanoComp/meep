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
#include <vector>

#include "meep.hpp"
#include "stability.hpp"

namespace meep {

/* Diagnose unstable dispersive timesteps (issue #12).

   E, H and the polarizations are stepped together, so a Lorentzian whose
   own recurrence is stable can still break the Courant condition; checking
   each pole on its own is not enough.  This analyses the coupled step, for a
   uniform bulk medium only; it says nothing about interfaces, PML or sources.

   For passive poles an eigenvalue of the plane-wave step can only leave the
   unit circle through z = -1, and only as k grows, so it suffices to check
   z = -1 at k*dx = pi on every axis.  There the damping drops out of each
   pole's recurrence (1+g) z^2 - (2 - w) z + (1-g), where g = gamma dt/2 and
   w = (omega_0 dt)^2 (w is dropped from the z term for Drude), and the step
   is stable iff every Lorentzian with sigma > 0 has w < 4 and
       ndim (dt/dx)^2 < mu * eps_N,  eps_N = eps - sum_j sigma_j w_j / (4 - w_j),
   the Courant condition with the permittivity the update sees at z = -1.
   tests/dispersive_stability.cpp checks this against the full amplification
   matrix. */

namespace {

/* eps_N above, or -infinity if a Lorentzian's own recurrence already has a
   root at z = -1 or beyond. */
double nyquist_epsilon(double eps_inf, const std::vector<lorentzian_pole> &poles, double dt) {
  double eps = eps_inf;
  for (size_t j = 0; j < poles.size(); ++j) {
    if (poles[j].sigma == 0) continue;
    const double w2dt2 = pow(2 * pi * poles[j].omega_0 * dt, 2);
    const double denom = 4 - (poles[j].no_omega_0_denominator ? 0 : w2dt2);
    if (!(denom > 0)) return -infinity;
    eps -= poles[j].sigma * w2dt2 / denom;
  }
  return eps;
}

} // namespace

bool dispersive_unstable(double eps_inf, double mu_inf, const std::vector<lorentzian_pole> &poles,
                         double dt, double dx, int ndim) {
  const double S = dt / dx;
  return !(ndim * S * S < mu_inf * nyquist_epsilon(eps_inf, poles, dt));
}

double dispersive_max_dt(double eps_inf, double mu_inf, const std::vector<lorentzian_pole> &poles,
                         double dx, int ndim) {
  /* Stable below it and unstable above, since eps_N only falls as dt grows.
     The poles lower eps_N, so the non-dispersive limit bounds it. */
  double lo = 0, hi = dx * sqrt(mu_inf * eps_inf / ndim);
  for (int i = 0; i < 60; ++i) {
    const double mid = (lo + hi) / 2;
    (dispersive_unstable(eps_inf, mu_inf, poles, mid, dx, ndim) ? hi : lo) = mid;
  }
  return hi;
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

stability_status medium_stability(const medium_struct *mm, const meep::grid_volume &gv, double dt,
                                  double *max_courant) {
  if (!medium_is_modelled(mm, gv)) return STABILITY_UNCHECKED;

  // A magnetic medium runs the same recursion with eps and mu swapped.
  const bool electric = !mm->E_susceptibilities.empty();
  const susceptibility_list &suscs = electric ? mm->E_susceptibilities : mm->H_susceptibilities;
  std::vector<meep::lorentzian_pole> poles;
  for (size_t i = 0; i < suscs.size(); ++i) {
    const meep::lorentzian_pole p = {suscs[i].frequency, suscs[i].gamma, suscs[i].sigma_diag.x,
                                     suscs[i].drude};
    poles.push_back(p);
  }
  const double eps_inf = electric ? mm->epsilon_diag.x : mm->mu_diag.x;
  const double mu_inf = electric ? mm->mu_diag.x : mm->epsilon_diag.x;
  const double dx = 1.0 / gv.a;
  const int ndim = meep::number_of_directions(gv.dim);
  if (max_courant) *max_courant = meep::dispersive_max_dt(eps_inf, mu_inf, poles, dx, ndim) / dx;
  return meep::dispersive_unstable(eps_inf, mu_inf, poles, dt, dx, ndim) ? STABILITY_UNSTABLE
                                                                         : STABILITY_STABLE;
}

void check_medium_stability_once(const medium_struct *mm, const meep::grid_volume &gv, double dt,
                                 std::vector<const medium_struct *> &checked) {
  if (mm->E_susceptibilities.empty() && mm->H_susceptibilities.empty()) return;
  for (size_t i = 0; i < checked.size(); ++i)
    if (checked[i] == mm || medium_struct_equal(checked[i], mm)) return;
  checked.push_back(mm);

  double max_courant;
  switch (medium_stability(mm, gv, dt, &max_courant)) {
    case STABILITY_UNSTABLE:
      meep::master_printf_stderr("warning: dispersive material may be unstable for the current "
                                 "timestep; it needs Courant < %g at this resolution, or a "
                                 "finer resolution\n",
                                 max_courant);
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
