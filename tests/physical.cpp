/* Copyright (C) 2005-2023 Massachusetts Institute of Technology
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2, or (at your option)
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#include <stdio.h>
#include <stdlib.h>

#include <meep.hpp>
using namespace meep;
using std::complex;

double one(const vec &) { return 1.0; }

int radiating_2D(const double xmax) {
  const double a = 10.0;
  const double ymax = 3.0;

  grid_volume gv = voltwo(xmax, ymax, a);
  structure s(gv, one, pml(ymax / 3));

  fields f(&s);
  double w = 0.30;
  double dx = 2.0;
  continuous_src_time src(w);
  f.add_point_source(Ez, src, vec(xmax / 2 - dx, ymax / 2));

  vec p1(xmax / 2 + 0 * dx, ymax / 2);
  vec p2(xmax / 2 + 1 * dx, ymax / 2);

  // let the source reach steady state
#if 1
  f.solve_cw(sizeof(realnum) == sizeof(float) ? 1e-5 : 1e-6);
#else
  while (f.time() < 400)
    f.step();
#endif

  complex<double> amp1 = f.get_field(Ez, p1);
  complex<double> amp2 = f.get_field(Ez, p2);
  double ratio = pow(abs(amp1) / abs(amp2), 2.0);
  master_printf("Ratio is %g from (%g %g) and (%g %g)\n", ratio, real(amp1), imag(amp1), real(amp2),
                imag(amp2));
  if (ratio > 2.12 || ratio < 1.88)
    meep::abort(
        "Failed: amp1 = (%g, %g), amp2 = (%g, %g)\n abs(amp1/amp2)^2 = %g, too far from 2.0\n",
        real(amp1), imag(amp1), real(amp2), imag(amp2), ratio);
  return 1;
}

int radiating_3D(const double xmax) {
  const double a = 10.0;
  const double ymax = 3.0;

  grid_volume gv = vol3d(xmax, ymax, ymax, a);
  symmetry S = mirror(Y, gv) - mirror(Z, gv);
  structure s(gv, one, pml(ymax / 3));

  fields f(&s);
  double w = 0.30;
  double dx = 2.0;
  continuous_src_time src(w);
  f.add_point_source(Ez, src, vec(xmax / 2 - dx, ymax / 2, ymax / 2));

  vec p1(xmax / 2 + 0 * dx, ymax / 2, ymax / 2);
  vec p2(xmax / 2 + 1 * dx, ymax / 2, ymax / 2);

  // let the source reach steady state
#if 1
  f.solve_cw(1e-3);
#else
  while (f.time() < 400)
    f.step();
#endif

  complex<double> amp1 = f.get_field(Ez, p1);
  complex<double> amp2 = f.get_field(Ez, p2);
  double ratio = abs(amp1) / abs(amp2);
  master_printf("Ratio is %g from (%g %g) and (%g %g)\n", ratio, real(amp1), imag(amp1), real(amp2),
                imag(amp2));
  if (ratio > 2.12 || ratio < 1.88)
    meep::abort(
        "Failed: amp1 = (%g, %g), amp2 = (%g, %g)\n abs(amp1/amp2) = %g, too far from 2.0\n",
        real(amp1), imag(amp1), real(amp2), imag(amp2), ratio);
  return 1;
}

void attempt(const char *name, int allright) {
  if (allright)
    master_printf("Passed %s\n", name);
  else
    meep::abort("Failed %s!\n", name);
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  // verbosity = 0;
  master_printf("Trying out some physical tests...\n");

  attempt("radiating source should decay spatially as 1/sqrt(r) in 2D.", radiating_2D(8.0));
  attempt("radiating source should decay spatially as 1/r in 3D.", radiating_3D(7.0));
  return 0;
}
