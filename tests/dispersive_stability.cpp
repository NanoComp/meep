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

/* Checks on the setup-time stability diagnostic for dispersive media. */

#include <cmath>
#include <cstdio>
#include <vector>

#include <meep.hpp>
#include <stability.hpp>

using namespace meep;
using meep_geom::medium_stability;

static const double margin = 1e-4; // as in medium_stability

static int failures = 0;

static void check(bool ok, const char *what) {
  master_printf("%-60s %s\n", what, ok ? "ok" : "FAILED");
  if (!ok) ++failures;
}

static bool unstable(double eps, double mu, const std::vector<lorentzian_pole> &poles, double dt,
                     double dx, int ndim) {
  return dispersive_spectral_radius(eps, mu, poles, dt, dx, ndim) > 1 + margin;
}

/* For lossless poles the analysis reduces to the Courant condition with the
   discrete permittivity at the Nyquist frequency, z = -1:
       ndim (dt/dx)^2 <= mu * [eps - sum_j c_j/(2 + a_j)]
   Return the dx at which that holds with equality. */
static double lossless_critical_dx(double eps, double mu, const std::vector<lorentzian_pole> &poles,
                                   double dt, int ndim) {
  double eps_nyquist = eps;
  for (size_t j = 0; j < poles.size(); ++j) {
    const double w2dt2 = pow(2 * pi * poles[j].omega_0 * dt, 2);
    const double a = 2 - (poles[j].no_omega_0_denominator ? 0 : w2dt2);
    eps_nyquist -= w2dt2 * poles[j].sigma / (2 + a);
  }
  return dt * sqrt(ndim / (mu * eps_nyquist));
}

static double seed_resolution, seed_sigma;

/* Excite the grid-scale mode that first becomes unstable.  Hy sits half a
   cell off the origin, where a cosine of this period would sample zero. */
static std::complex<double> nyquist_seed(const vec &p) {
  return std::complex<double>(sin(pi * seed_resolution * p.z()), 0.0);
}

// Deterministic white noise, which excites every wavevector in the cell.
static std::complex<double> broadband_seed(const vec &p) {
  const double x = sin(12.9898 * p.z() * seed_resolution + 1.0) * 43758.5453;
  return std::complex<double>(x - floor(x) - 0.5, 0.0);
}

/* Run the actual 1D field update at Courant 0.5 with one susceptibility on the
   E or H side, and return the field energy relative to its initial value. */
static double energy_growth(double resolution, int steps, field_type ft,
                            const lorentzian_pole &pole, std::complex<double> seed(const vec &)) {
  grid_volume gv = volone(1.0, resolution);
  gv.center_origin();
  structure s(gv, [](const vec &) { return 1.0; });
  seed_resolution = resolution;
  seed_sigma = pole.sigma;
  s.add_susceptibility(
      [](const vec &) { return seed_sigma; }, ft,
      lorentzian_susceptibility(pole.omega_0, pole.gamma, pole.no_omega_0_denominator));
  fields f(&s);
  f.use_bloch(0.0);
  f.initialize_field(Hy, seed);
  const double e0 = f.field_energy();
  for (int i = 0; i < steps; ++i)
    f.step();
  return f.field_energy() / e0;
}

static meep_geom::medium_struct medium(const lorentzian_pole &p, bool magnetic) {
  meep_geom::medium_struct mm;
  meep_geom::susceptibility ss = {}; // vector3 members are not zeroed otherwise
  ss.frequency = p.omega_0;
  ss.gamma = p.gamma;
  ss.drude = p.no_omega_0_denominator;
  ss.sigma_diag.x = ss.sigma_diag.y = ss.sigma_diag.z = p.sigma;
  (magnetic ? mm.H_susceptibilities : mm.E_susceptibilities).push_back(ss);
  return mm;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  const double courant = 0.5;

  // Without dispersion the estimate reproduces dt <= dx*sqrt(eps*mu/ndim).
  const std::vector<lorentzian_pole> none;
  const double vdx = 0.05;
  for (int ndim = 1; ndim <= 3; ++ndim) {
    const double dt = vdx / sqrt(double(ndim));
    char what[64];
    snprintf(what, sizeof(what), "vacuum brackets the Courant limit in %dD", ndim);
    check(!unstable(1, 1, none, 0.999 * dt, vdx, ndim) &&
              unstable(1, 1, none, 1.001 * dt, vdx, ndim),
          what);
  }
  check(!unstable(2, 0.25, none, 0.999 * vdx / sqrt(2.0), vdx, 1) &&
            unstable(2, 0.25, none, 1.001 * vdx / sqrt(2.0), vdx, 1),
        "eps = 2, mu = 1/4 bracket dt = dx/sqrt(2)");

  // A strong Lorentzian, unstable on a coarse grid and stable on a fine one.
  const lorentzian_pole lorentz = {2.0, 0.0, 60.0, false};
  const double coarse = 24, fine = 400;

  /* Its own recurrence is stable (omega_0 dt = 0.26, well below 2), yet the
     coupled update blows up, in the prediction and in FDTD. */
  check(2 * pi * lorentz.omega_0 * courant / coarse < 2 &&
            unstable(1, 1, {lorentz}, courant / coarse, 1 / coarse, 1) &&
            energy_growth(coarse, 8, E_stuff, lorentz, nyquist_seed) > 1e4,
        "a stable pole still makes the coarse FDTD run blow up");

  /* Lossless poles have a closed-form limit; evaluating only the largest Yee
     wavevector has to land exactly on it. */
  const std::vector<lorentzian_pole> two_poles = {{2.0, 0.0, 60.0, false}, {1.0, 0.0, 5.0, true}};
  for (int ndim = 1; ndim <= 3; ndim += 2) {
    const double dt = courant / coarse;
    const double dx = lossless_critical_dx(1.5, 1, two_poles, dt, ndim);
    char what[64];
    snprintf(what, sizeof(what), "lossless poles bracket the Nyquist limit in %dD", ndim);
    check(!unstable(1.5, 1, two_poles, dt, 1.001 * dx, ndim) &&
              unstable(1.5, 1, two_poles, dt, 0.999 * dx, ndim),
          what);
  }

  /* A Drude metal with omega_p dt ~ 1 sits close to its limit in 3D: stable at
     resolution 24, unstable at 20.  Its undamped |z| = 1 mode is where the
     estimator is least accurate, so it must not be mistaken for growth. */
  const lorentzian_pole drude = {1.0, 0.0, 50.0, true};
  check(!unstable(1, 1, {drude}, courant / 24, 1.0 / 24, 3), "undamped Drude is stable at res 24");
  check(unstable(1, 1, {drude}, courant / 20, 1.0 / 20, 3), "undamped Drude is unstable at res 20");

  check(energy_growth(coarse, 8, H_stuff, lorentz, nyquist_seed) > 1e4,
        "FDTD grows for the magnetic Lorentzian too");

  /* Every wavevector in the cell, not just the Nyquist one, has to stay bounded
     just past the predicted threshold resolution, and blow up just short of it. */
  const lorentzian_pole lossy = {2.0, 0.2, 60.0, false};
  for (const lorentzian_pole &p : {lossy, drude}) {
    double lo = 10, hi = 1000;
    for (int i = 0; i < 50; ++i) {
      const double mid = sqrt(lo * hi);
      (unstable(1, 1, {p}, courant / mid, 1 / mid, 1) ? lo : hi) = mid;
    }
    check(energy_growth(ceil(1.03 * hi), 4000, E_stuff, p, broadband_seed) < 2 &&
              energy_growth(floor(0.97 * hi), 20, E_stuff, p, broadband_seed) > 1e4,
          p.no_omega_0_denominator ? "FDTD agrees with the Drude threshold on both sides"
                                   : "FDTD agrees with the Lorentzian threshold on both sides");
  }

  // The medium_struct path, as used at setup.
  const meep_geom::medium_struct em = medium(lorentz, false), hm = medium(lorentz, true);
  const grid_volume gc = vol3d(1.0, 1.0, 1.0, coarse), gf = vol3d(1.0, 1.0, 1.0, fine);
  const double dtc = courant / coarse;
  check(medium_stability(&em, gc, dtc) == meep_geom::STABILITY_UNSTABLE,
        "coarse electric medium is unstable");
  check(medium_stability(&hm, gc, dtc) == meep_geom::STABILITY_UNSTABLE,
        "coarse magnetic medium is unstable");
  check(medium_stability(&em, gf, courant / fine) == meep_geom::STABILITY_STABLE,
        "refined medium is stable");
  check(medium_stability(&em, volcyl(1.0, 1.0, coarse), dtc) == meep_geom::STABILITY_UNCHECKED,
        "cylindrical grids are not checked");

  // Diagonal anisotropy mixes components for oblique k, so it is declined.
  meep_geom::medium_struct aniso = em;
  aniso.E_susceptibilities[0].sigma_diag.z = 0;
  check(medium_stability(&aniso, gc, dtc) == meep_geom::STABILITY_UNCHECKED,
        "anisotropic sigma is not checked");

  // The same medium, or an equal copy of it, is checked (and warned about) once.
  std::vector<const meep_geom::medium_struct *> checked;
  const meep_geom::medium_struct copy = medium(lorentz, false);
  meep_geom::check_medium_stability_once(&em, gc, dtc, checked);
  meep_geom::check_medium_stability_once(&em, gc, dtc, checked);
  meep_geom::check_medium_stability_once(&copy, gc, dtc, checked);
  check(checked.size() == 1, "each distinct medium is checked once");

  master_printf("%s\n", failures ? "FAILED" : "all dispersive stability checks passed");
  return failures ? 1 : 0;
}
