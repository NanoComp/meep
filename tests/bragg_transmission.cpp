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

/* Compute the transmission spectrum through a 4-layer 1d Bragg mirror,
   and compare to the result from the analytical transfer matrices.
   The transmission spectrum is computed via the dft_flux feature,
   which dynamically updates the DFTs of the fields on the flux plane
   as we go along. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <meep.h>
using namespace meep;

const double nhi = 3, nlo = 1;
const double wlo = nhi / (nlo + nhi);
const int Nperiods = 4;
const double zsize = 10;

double one(const vec &) { return 1.0; }

double eps_bragg(const vec &v) {
  double z = v.z() - zsize * 0.5;
  
  if (fabs(z)*2 > Nperiods)
    return 1.0;
  else {
    double zi;
    double zf = modf(z, &zi);
    if (zf < 0) zf += 1;
    if (zf < wlo)
      return (nlo*nlo);
    else
      return (nhi*nhi);
  }
}

typedef complex<double> matrix2x2[2][2];

/* multiply m by transfer matrix from n1 to n2 */
inline void byT12(matrix2x2 m, double n1, double n2)
{
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
inline void byP(matrix2x2 m, double n, double w, double dz)
{
     complex<double> p, pc;

     p = exp(complex<double>(0, n * w * dz));
     pc = conj(p);

     m[0][0] *= p;
     m[0][1] *= pc;
     m[1][0] *= p;
     m[1][1] *= pc;
}

inline double abs2(complex<double> x) { double ax = abs(x); return ax*ax; }

void bragg_transmission_analytic(double freq_min, double freq_max, int nfreq,
				 double *T)
{
  for (int i = 0; i < nfreq; ++i) {
    double omega = 2*pi * (freq_min + i * (freq_max - freq_min) / (nfreq - 1));
    matrix2x2 Tm = { { 1, 0 }, { 0, 1 } };
    for (int j = 0; j < Nperiods; ++j) {
      byT12(Tm, nlo, nhi);
      byP(Tm, nhi, omega, 1 - wlo);
      byT12(Tm, nhi, nlo);
      byP(Tm, nlo, omega, wlo);
    }
    complex<double> refl = - Tm[1][0] / Tm[1][1];
    T[i] = abs2(Tm[0][0] + refl * Tm[0][1]);
  }
}

void bragg_transmission(double a, double freq_min, double freq_max, int nfreq,
			double *T) {
  const volume v = volone(zsize, a);

  structure s(v, eps_bragg);
  s.use_pml_everywhere(0.5);
  fields f(&s);
  f.use_real_fields();

  structure s0(v, one);
  s0.use_pml_everywhere(0.5);
  fields f0(&s0);
  f0.use_real_fields();

  vec srcpt(0.1), fluxpt(zsize - 0.1);

  gaussian_src_time src((freq_min + freq_max) * 0.5,
			0.5 / fabs(freq_max - freq_min),
			0, 5 / fabs(freq_max - freq_min));
  f.add_point_source(Ex, src, srcpt);
  f0.add_point_source(Ex, src, srcpt);

  dft_flux fp = f.add_dft_flux_plane(geometric_volume(fluxpt, fluxpt),
				     freq_min, freq_max, nfreq);
  dft_flux fp0 = f0.add_dft_flux_plane(geometric_volume(fluxpt, fluxpt),
				       freq_min, freq_max, nfreq);

  int dindex = 0;
  f.output_hdf5(Dielectric, f0.v.surroundings(), a, true, dindex, false, true, "f");
  while (f.time() < nfreq / fabs(freq_max - freq_min) / 2) {
    f.step();
    f0.step();

    if (0 && f.t % 20 == 0) { // output fields for debugging
      f.output_hdf5(Ex, f0.v.surroundings(), a, true, dindex, false, true, "f");
      f0.output_hdf5(Ex, f0.v.surroundings(), a, true, dindex, false, true, "f0");
      dindex++;
    }
  }

  double *flux = fp.flux();
  double *flux0 = fp0.flux();
  for (int i = 0; i < nfreq; ++i)
    T[i] = flux[i] / flux0[i];
  delete[] flux;
  delete[] flux0;
}

inline double max2(double a, double b) { return (a > b ? a : b); }
inline double min2(double a, double b) { return (a < b ? a : b); }
inline double max2a(double a, double b) { return max2(abs(a), abs(b)); }
inline double sqr(double x) { return x*x; }

/* The discretization errors tend to result in a *shift* of the spectral
   features more than a change in their amplitude.  Because these features
   are very sharp (e.g. at the gap edges), it is more appropriate to compute
   errors via the distance from a point to the curve, rather than just
   the difference of the abscissae.  That's what this function does. */
double distance_from_curve(int n, double dx, double ys[], double x, double y)
{
  double d = infinity;
  for (int i = 1; i < n; ++i) {
    double theta = atan2(ys[i] - ys[i-1], dx);
    double L = sqrt(sqr(dx) + sqr(ys[i]-ys[i-1]));
    double x0 = x - (i-1) * dx;
    double y0 = y - ys[i-1];
    double x0p = x0 * cos(theta) + y0 * sin(theta);
    double y0p = y0 * cos(theta) - x0 * sin(theta);
    if (x0p < 0)
      d = min2(sqrt(sqr(x0) + sqr(y0)), d);
    else if (x0p > L)
      d = min2(sqrt(sqr(x-i*dx) + sqr(y-ys[i])), d);
    else
      d = min2(abs(y0p), d);
  }
  return d;
}

int main(int argc, char **argv) {
  const int nfreq = 100;
  const double freq_min = 0.1, freq_max = 0.5;
  initialize mpi(argc, argv);

  double *T = new double[nfreq];
  bragg_transmission(40.0, freq_min, freq_max, nfreq, T);
  double *T0 = new double[nfreq];
  bragg_transmission_analytic(freq_min, freq_max, nfreq, T0);

  double dfreq = (freq_max - freq_min) / (nfreq - 1);

  for (int i = 0; i < nfreq; ++i) {
    double err = distance_from_curve(nfreq, dfreq, T0, i * dfreq, T[i])
      / T0[i];
    if (err * sqr(freq_min / (freq_min + i*dfreq)) > 0.01)
      abort("large rel. error %g at freq = %g: T = %g instead of %g\n",
	    err, freq_min + i*dfreq, T[i], T0[i]);
  }

  if (0) { // output transmissions for debugging
    master_printf("transmission:, freq (c/a), T\n");
    for (int i = 0; i < nfreq; ++i)
      master_printf("transmission:, %g, %g, %g\n",
		    freq_min + i * dfreq,
		    T[i], T0[i]);
  }

  delete[] T0;
  delete[] T;

  return 0;
}
