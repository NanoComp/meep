/* Copyright (C) 2003 Massachusetts Institute of Technology  
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

#include "dactyl.h"

double one(const vec &) { return 1.0; }
static double width = 20.0;
double bump(const vec &v) { return (fabs(v.z()-50.0) > width)?1.0:12.0; }

int compare(double a, double b, double eps, const char *n) {
  if (fabs(a-b) > fabs(b)*eps) {
    printf("%s differs by\t%lg out of\t%lg\n", n, a-b, b);
    printf("This gives a fractional error of %lg\n", fabs(a-b)/fabs(b));
    return 0;
  } else {
    return 1;
  }
}

int flux_1d(const double zmax,
                  double eps(const vec &)) {
  const double a = 10.0;
  const double gridpts = a*zmax;
  const double ttot = 10.0 + 1e5/zmax;

  volume v = volone(zmax,a);
  mat ma(v, eps);
  ma.use_pml_left(zmax/6);
  ma.use_pml_right(zmax/6);

  fields f(&ma);
  f.use_real_fields();
  f.add_point_source(Ex, 0.7, 2.5, 0.0, 3.0, vec(zmax/2+0.3), 1.0);
  flux_plane *left = f.add_flux_plane(vec(zmax/3.0), vec(zmax/3.0));
  flux_plane *right = f.add_flux_plane(vec(zmax*2.0/3.0), vec(zmax*2.0/3.0));

  while (f.time() <= f.find_last_source()) f.step();

  volume mid = volone(zmax/3,a);
  mid.origin = vec(zmax/3);
  double flux_energy=0.0;
  double delta_energy = f.energy_in_box(mid);
  while (f.time() < ttot) {
    f.step();
    flux_energy += (c/a)*(right->flux() - left->flux());
  }
  delta_energy -= f.energy_in_box(mid);
  master_printf("  Energy change:  \t%lg\n  Integrated flux:\t%lg\n  Ratio:\t%lg\n",
                delta_energy, flux_energy, flux_energy/delta_energy);
  return compare(delta_energy, flux_energy, 0.15, "Flux");
}

void attempt(const char *name, int allright) {
  if (allright) master_printf("Passed %s\n", name);
  else abort("Failed %s!\n", name);
}

int main(int argc, char **argv) {
  initialize(argc, argv);
  master_printf("Trying out the fluxes...\n");

  width = 20.0;
  attempt("Flux 1D 100", flux_1d(100.0, bump));
  width = 10.0;
  attempt("Flux 1D 100", flux_1d(100.0, bump));
  width = 300.0;
  attempt("Flux 1D 100", flux_1d(100.0, bump));

  finished();
  exit(0);
}

