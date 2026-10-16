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
#include <complex>
#include <cstdio>
#include <vector>

#include <meep.hpp>
#include <stability.hpp>

using namespace meep;
using meep_geom::medium_stability;

static int failures = 0;

static void check(bool ok, const char *what) {
  master_printf("%-60s %s\n", what, ok ? "ok" : "FAILED");
  if (!ok) ++failures;
}

static bool unstable(double eps, double mu, const std::vector<lorentzian_pole> &poles, double dt,
                     double dx, int ndim) {
  return dispersive_unstable(eps, mu, poles, dt, dx, ndim);
}

/* An independent check on the closed form in stability.cpp: the full
   amplification matrix of one timestep on [E, H, P_j, P_j_prev], for a plane
   wave at k*dx = pi on every axis, built from the coefficients of update_P. */

// Row-major n x n complex matrix.
typedef std::vector<std::complex<double> > cmatrix;

// Matrix infinity norm (max row sum, operator norm).
static double matrix_norm(const cmatrix &A, int n) {
  double best = 0;
  for (int i = 0; i < n; ++i) {
    double s = 0;
    for (int j = 0; j < n; ++j)
      s += std::abs(A[i * n + j]);
    if (s > best) best = s;
  }
  return best;
}

static void matrix_square(cmatrix &A, int n) {
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
   rho(A)^m <= |A^m| for an operator norm. */
static double spectral_radius(cmatrix A, int n) {
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

static double matrix_spectral_radius(double eps_inf, double mu_inf,
                                     const std::vector<lorentzian_pole> &poles, double dt,
                                     double dx, int ndim) {
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

// Undamped modes sit on the unit circle, where the estimate overshoots slightly.
static bool matrix_unstable(double eps, double mu, const std::vector<lorentzian_pole> &poles,
                            double dt, double dx, int ndim) {
  return matrix_spectral_radius(eps, mu, poles, dt, dx, ndim) > 1 + 1e-4;
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

  /* The full amplification matrix has to agree with the closed form on both
     sides of its threshold, with or without damping, and also where a pole
     above 1/(pi dt) is what sets the limit. */
  const lorentzian_pole drude = {1.0, 0.0, 50.0, true};
  const lorentzian_pole lossy = {2.0, 0.2, 60.0, false};
  struct {
    const char *name;
    double eps, mu;
    std::vector<lorentzian_pole> poles;
  } media[] = {
      {"lossless Lorentzian", 1, 1, {lorentz}},
      {"lossy Lorentzian", 1, 1, {lossy}},
      {"undamped Drude", 1, 1, {drude}},
      {"damped Drude", 1, 1, {{1.0, 0.5, 50.0, true}}},
      {"Lorentzian + Drude", 1.5, 1, {lorentz, {1.0, 0.0, 5.0, true}}},
      {"lossy Lorentzian, eps = 2, mu = 1/2", 2, 0.5, {lossy}},
      {"weak pole above 1/(pi dt)", 1, 1, {{12.0, 0.1, 0.01, false}}},
  };
  for (const auto &m : media)
    for (int ndim = 1; ndim <= 3; ++ndim) {
      const double dx = 1 / coarse, dt = dispersive_max_dt(m.eps, m.mu, m.poles, dx, ndim);
      char what[96];
      snprintf(what, sizeof(what), "matrix brackets the limit: %s, %dD", m.name, ndim);
      check(!matrix_unstable(m.eps, m.mu, m.poles, 0.999 * dt, dx, ndim) &&
                matrix_unstable(m.eps, m.mu, m.poles, 1.001 * dt, dx, ndim),
            what);
    }
  check(dispersive_max_dt(1, 1, {lorentz}, 1 / coarse, 3) ==
            dispersive_max_dt(1, 1, {lossy}, 1 / coarse, 3),
        "damping does not move the limit");

  /* A Drude metal with omega_p dt ~ 1 sits close to its limit in 3D: stable at
     resolution 24, unstable at 20.  Its undamped z = 1 mode sits on the unit
     circle, and must not be mistaken for growth. */
  check(!unstable(1, 1, {drude}, courant / 24, 1.0 / 24, 3), "undamped Drude is stable at res 24");
  check(unstable(1, 1, {drude}, courant / 20, 1.0 / 20, 3), "undamped Drude is unstable at res 20");

  check(energy_growth(coarse, 8, H_stuff, lorentz, nyquist_seed) > 1e4,
        "FDTD grows for the magnetic Lorentzian too");

  /* Every wavevector in the cell, not just the Nyquist one, has to stay bounded
     just past the predicted threshold resolution, and blow up just short of it. */
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
  double max_courant;
  medium_stability(&em, gc, dtc, &max_courant);
  check(
      max_courant < courant &&
          medium_stability(&em, gc, 0.999 * max_courant / coarse) == meep_geom::STABILITY_STABLE &&
          medium_stability(&em, gc, 1.001 * max_courant / coarse) == meep_geom::STABILITY_UNSTABLE,
      "the reported Courant number is the limit");
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
