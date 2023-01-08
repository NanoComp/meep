/* Copyright (C) 2005-2023 Massachusetts Institute of Technology
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

/* Nonlinear test program checking 2nd and 3rd harmonic generation */

#include <meep.hpp>
using namespace meep;
using std::max;

double the_value = 1.0;
double value(const vec &) { return the_value; }

void harmonics(double freq, double chi2, double chi3, double J, double &A2, double &A3) {
  const double dpml = 5.0;
  const double res = 20;
  const double sz = 100 + 2 * dpml;
  grid_volume gv = vol1d(sz, res);
  gv.center_origin();

  the_value = 1.0;
  structure s(gv, value, pml(dpml));
  the_value = chi2;
  s.set_chi2(value);
  the_value = chi3;
  s.set_chi3(value);

  fields f(&s);
  f.use_real_fields();

  gaussian_src_time src(freq, freq / 20);
  f.add_point_source(Ex, src, vec(-0.5 * sz + dpml), J);

  vec fpt(0.5 * sz - dpml - 0.5);
  dft_flux d1 =
      f.add_dft_flux(Z, volume(fpt), freq, freq, 1, true, true, 1 /* decimation_factor */);
  dft_flux d2 =
      f.add_dft_flux(Z, volume(fpt), 2 * freq, 2 * freq, 1, true, true, 1 /* decimation_factor */);
  dft_flux d3 =
      f.add_dft_flux(Z, volume(fpt), 3 * freq, 3 * freq, 1, true, true, 1 /* decimation_factor */);

  double emax = 0;

  while (f.time() < f.last_source_time()) {
    emax = max(emax, abs(f.get_field(Ex, fpt)));
    f.step();
  }
  do {
    double emaxcur = 0;
    double T = f.time() + 50;
    while (f.time() < T) {
      double e = abs(f.get_field(Ex, fpt));
      emax = max(emax, e);
      emaxcur = max(emaxcur, e);
      f.step();
    }
    if (emaxcur < 1e-6 * emax) break;
  } while (1);

  double *d1f = d1.flux();
  double *d2f = d2.flux();
  double *d3f = d3.flux();

  A2 = *d2f / *d1f;
  A3 = *d3f / *d1f;

  master_printf("harmonics(%g,%g,%g) = %g, %g\n", chi2, chi3, J, A2, A3);

  delete[] d1f;
  delete[] d2f;
  delete[] d3f;
}

int different(double a, double a0, double thresh, const char *msg) {
  if (fabs(a - a0) > thresh * fabs(a0)) {
    master_printf("error: %s\n --- %g vs. %g (%g error > %g)\n", msg, a, a0,
                  fabs(a - a0) / fabs(a0), thresh);
    return 1;
  }
  else
    return 0;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 0;
  const double freq = 1.0 / 3.0;

  double a2, a3, a2_2, a3_2;

  double thresh = sizeof(realnum) == sizeof(float) ? 1e-3 : 1e-5;
  harmonics(freq, 0.27e-4, 1e-4, 1.0, a2, a3);
  if (different(a2, 9.80330e-07, thresh, "2nd harmonic mismatches known val")) return 1;
  if (sizeof(realnum) == sizeof(float)) {
    if (different(a3, 9.99349e-07, thresh, "3rd harmonic mismatches known val")) return 1;
  }
  else {
    if (different(a3, 9.97747e-07, thresh, "3rd harmonic mismatches known val")) return 1;
  }

  harmonics(freq, 0.54e-4, 2e-4, 1.0, a2_2, a3_2);
  master_printf("doubling chi2, chi3 = %g x 2nd harmonic, %g x 3rd\n", a2_2 / a2, a3_2 / a3);
  if (different(a2_2 / a2, 4.0, 0.01, "incorrect chi2 scaling")) return 1;
  if (different(a3_2 / a3, 4.0, 0.01, "incorrect chi3 scaling")) return 1;

  harmonics(freq, 0.27e-4, 1e-4, 2.0, a2_2, a3_2);
  master_printf("doubling J = %g x 2nd harmonic, %g x 3rd\n", a2_2 / a2, a3_2 / a3);
  if (different(a2_2 / a2, 4.0, 0.01, "incorrect J scaling for 2nd harm.")) return 1;
  if (different(a3_2 / a3, 16.0, 0.01, "incorrect J scaling for 3rd harm.")) return 1;

  harmonics(freq, 0.27e-4, 0.0, 1.0, a2_2, a3_2);
  if (different(a2, a2_2, 1e-2, "chi3 has too big effect on 2nd harmonic")) return 1;
  if (a3_2 / a3 > 1e-4) {
    master_printf("error: too much 3rd harmonic without chi3\n");
    return 1;
  }

  harmonics(freq, 0.0, 1e-4, 1.0, a2_2, a3_2);
  if (sizeof(realnum) == sizeof(float)) {
    if (different(a3, a3_2, 0.0017, "chi2 has too big effect on 3rd harmonic")) return 1;
  }
  else {
    if (different(a3, a3_2, 0.001, "chi2 has too big effect on 3rd harmonic")) return 1;
  }

  if (a2_2 / a2 > 1e-5) {
    master_printf("error: too much 2nd harmonic without chi3\n");
    return 1;
  }

  return 0;
}
