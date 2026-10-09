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

/* Setup-time stability diagnostic for Lorentzian/Drude media (issue #12).
   Not installed; declared here so that tests/dispersive_stability.cpp can
   reach it. */

#ifndef MEEP_STABILITY_H
#define MEEP_STABILITY_H

#include <vector>

#include "meepgeom.hpp"

namespace meep {

// One Lorentzian or Drude pole, as lorentzian_susceptibility sees it.
struct lorentzian_pole {
  double omega_0, gamma, sigma;
  bool no_omega_0_denominator; // true for a Drude term
};

/* Spectral radius of one timestep of a uniform, isotropic medium on a
   Cartesian grid with ndim axes, at the largest Yee spatial frequency
   (k*dx = pi on every axis).  Requires eps_inf > 0 and mu_inf > 0. */
double dispersive_spectral_radius(double eps_inf, double mu_inf,
                                  const std::vector<lorentzian_pole> &poles, double dt, double dx,
                                  int ndim);

} // namespace meep

namespace meep_geom {

enum stability_status { STABILITY_UNCHECKED, STABILITY_STABLE, STABILITY_UNSTABLE };

/* Verdict of the bulk von Neumann estimate for this medium on this grid, or
   STABILITY_UNCHECKED if the estimate does not model it (or it has no
   susceptibilities). */
stability_status medium_stability(const medium_struct *mm, const meep::grid_volume &gv, double dt);

/* Warn if medium_stability is STABILITY_UNSTABLE, skipping media already in
   `checked`. */
void check_medium_stability_once(const medium_struct *mm, const meep::grid_volume &gv, double dt,
                                 std::vector<const medium_struct *> &checked);

} // namespace meep_geom

#endif /* MEEP_STABILITY_H */
