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

/* Compute the transmission spectrum through a 4-layer 1d Bragg mirror,
   and compare to the result from the analytical transfer matrices.
   The transmission spectrum is computed via the dft_flux feature,
   which dynamically updates the DFTs of the fields on the flux plane
   as we go along. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <meep.hpp>
using namespace meep;
using std::complex;

const double nhi = 3, nlo = 1;
const double wlo = nhi / (nlo + nhi);
const int Nperiods = 4;
const double zsize = 10;

double eps_nlo(const vec &) { return nlo * nlo; }

double eps_bragg(const vec &pt) {
  double z = pt.z() - zsize * 0.5;

  if (fabs(z) * 2 > Nperiods)
    return nlo * nlo;
  else {
    double zi;
    double zf = modf(z, &zi);
    if (zf < 0) zf += 1;
    if (zf < wlo)
      return (nlo * nlo);
    else
      return (nhi * nhi);
  }
}

typedef complex<double> matrix2x2[2][2];

/* multiply m by transfer matrix from n1 to n2 */
inline void byT12(matrix2x2 m, double n1, double n2) {
  complex<double> m00, m01, m10, m11;
  double td, tod;
  double n12 = n1 / n2;

  td = 0.5 * (1 + n12);
  tod = 0.5 * (1 - n12);

  m00 = m[0][0];
  m01 = m[0][1];
  m10 = m[1][0];
  m11 = m[1][1];

  m[0][0] = m00 * td + m01 * tod;
  m[0][1] = m00 * tod + m01 * td;
  m[1][0] = m10 * td + m11 * tod;
  m[1][1] = m10 * tod + m11 * td;
}

/* multiply m by propagation matrix through dz of index n, frequency w */
inline void byP(matrix2x2 m, double n, double w, double dz) {
  complex<double> p, pc;

  p = exp(complex<double>(0, n * w * dz));
  pc = conj(p);

  m[0][0] *= p;
  m[0][1] *= pc;
  m[1][0] *= p;
  m[1][1] *= pc;
}

inline double abs2(complex<double> x) {
  double ax = abs(x);
  return ax * ax;
}

void bragg_transmission_analytic(double freq_min, double freq_max, int nfreq, double *T,
                                 double *R) {
  for (int i = 0; i < nfreq; ++i) {
    double omega = 2 * pi * (freq_min + i * (freq_max - freq_min) / (nfreq - 1));
    matrix2x2 Tm = {{1, 0}, {0, 1}};
    for (int j = 0; j < Nperiods; ++j) {
      byT12(Tm, nlo, nhi);
      byP(Tm, nhi, omega, 1 - wlo);
      byT12(Tm, nhi, nlo);
      byP(Tm, nlo, omega, wlo);
    }
    complex<double> refl = -Tm[1][0] / Tm[1][1];
    T[i] = abs2(Tm[0][0] + refl * Tm[0][1]);
    R[i] = abs2(refl);
  }
}

void bragg_transmission(double a, double freq_min, double freq_max, int nfreq, double *T, double *R,
                        bool use_hdf5) {
  const grid_volume gv = volone(zsize, a);

  structure *s = new structure(gv, eps_bragg, pml(0.5));
  fields f(s);
  f.use_real_fields();

  structure s0(gv, eps_nlo, pml(0.5));
  fields f0(&s0);
  f0.use_real_fields();

  vec srcpt(0.1), Tfluxpt(zsize - 0.1), Rfluxpt(0.1);

  gaussian_src_time src((freq_min + freq_max) * 0.5, 0.5 / fabs(freq_max - freq_min), 0,
                        5 / fabs(freq_max - freq_min));
  f.add_point_source(Ex, src, srcpt);
  f0.add_point_source(Ex, src, srcpt);

  dft_flux ft = f.add_dft_flux_plane(Tfluxpt, freq_min, freq_max, nfreq);
  dft_flux fr = f.add_dft_flux_plane(Rfluxpt, freq_min, freq_max, nfreq);
  dft_flux ft0 = f0.add_dft_flux_plane(Tfluxpt, freq_min, freq_max, nfreq);
  dft_flux fr0 = f0.add_dft_flux_plane(Rfluxpt, freq_min, freq_max, nfreq);

  while (f0.time() < nfreq / fabs(freq_max - freq_min) / 2)
    f0.step();

  /* we want to subtract the fields for the reflection... */
  if (use_hdf5) {
    /* simulate a case where the normalization is done
       by a separate run and saved to a file */
    fr0.save_hdf5(f, "flux", "reflection");
    fr.load_hdf5(f, "flux", "reflection");
    fr.scale_dfts(-1.0);

    // clean up after ourselves: delete the file
    h5file *ff = f.open_h5file("flux", h5file::READONLY);
    ff->remove();
    delete ff;
  }
  else
    fr -= fr0;

  while (f.time() < nfreq / fabs(freq_max - freq_min) / 2)
    f.step();

  double *flux = ft.flux();
  double *flux0 = ft0.flux();
  for (int i = 0; i < nfreq; ++i)
    T[i] = flux[i] / flux0[i];
  delete[] flux;
  flux = fr.flux();
  for (int i = 0; i < nfreq; ++i)
    R[i] = -flux[i] / flux0[i];
  delete[] flux;
  delete[] flux0;
  delete s; // tests whether okay to delete s before f
}

inline double max2(double a, double b) { return (a > b ? a : b); }
inline double min2(double a, double b) { return (a < b ? a : b); }
inline double max2a(double a, double b) { return max2(abs(a), abs(b)); }
inline double sqr(double x) { return x * x; }

/* The discretization errors tend to result in a *shift* of the spectral
   features more than a change in their amplitude.  Because these features
   are very sharp (e.g. at the gap edges), it is more appropriate to compute
   errors via the distance from a point to the curve, rather than just
   the difference of the abscissae.  That's what this function does. */
double distance_from_curve(int n, double dx, double ys[], double x, double y) {
  double d = meep::infinity;
  for (int i = 1; i < n; ++i) {
    double theta = atan2(ys[i] - ys[i - 1], dx);
    double L = sqrt(sqr(dx) + sqr(ys[i] - ys[i - 1]));
    double x0 = x - (i - 1) * dx;
    double y0 = y - ys[i - 1];
    double x0p = x0 * cos(theta) + y0 * sin(theta);
    double y0p = y0 * cos(theta) - x0 * sin(theta);
    if (x0p < 0)
      d = min2(sqrt(sqr(x0) + sqr(y0)), d);
    else if (x0p > L)
      d = min2(sqrt(sqr(x - i * dx) + sqr(y - ys[i])), d);
    else
      d = min2(abs(y0p), d);
  }
  return d;
}

void doit(bool use_hdf5) {
  const int nfreq = 100;
  const double freq_min = 0.1, freq_max = 0.5;

  double *T = new double[nfreq];
  double *R = new double[nfreq];
  bragg_transmission(40.0, freq_min, freq_max, nfreq, T, R, use_hdf5);
  double *T0 = new double[nfreq];
  double *R0 = new double[nfreq];
  bragg_transmission_analytic(freq_min, freq_max, nfreq, T0, R0);

  double dfreq = (freq_max - freq_min) / (nfreq - 1);

  if (0) { // output transmission & reflection spectra for debugging
    master_printf("transmission:, freq (c/a), T, R, T0, R0\n");
    for (int i = 0; i < nfreq; ++i)
      master_printf("transmission:, %g, %g, %g, %g, %g\n", freq_min + i * dfreq, T[i], R[i], T0[i],
                    R0[i]);
  }

  double maxerrT = 0, maxerrR = 0;
  for (int i = 0; i < nfreq; ++i) {
    double errT = distance_from_curve(nfreq, dfreq, T0, i * dfreq, T[i]);
    double errR = distance_from_curve(nfreq, dfreq, R0, i * dfreq, R[i]);
    if (errT > maxerrT) maxerrT = errT;
    if (errR > maxerrR) maxerrR = errR;
    if (errT * sqr(freq_min / (freq_min + i * dfreq)) > 0.01)
      meep::abort("large error %g at freq = %g: T = %g instead of %g\n", errT, freq_min + i * dfreq,
                  T[i], T0[i]);
    if (errR * sqr(freq_min / (freq_min + i * dfreq)) > 0.01)
      meep::abort("large error %g at freq = %g: R = %g instead of %g\n", errR, freq_min + i * dfreq,
                  R[i], R0[i]);
  }
  master_printf("Done (max. err in T = %e, in R = %e)\n", maxerrT, maxerrR);

  delete[] R0;
  delete[] T0;
  delete[] R;
  delete[] T;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 0;

#ifdef HAVE_HDF5
  doit(true);
#endif
  doit(false);

  return 0;
}
