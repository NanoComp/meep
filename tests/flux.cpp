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

double cavity(const vec &v) {
  const double zz = fabs(v.z() - 7.5);
  if (zz > 4.0) return 1.0;
  if (zz < 1.0) return 1.0;
  double norm = zz;
  while (norm > 1.0) norm -= 1.0;
  if (norm > 0.3) return 1.0;
  return 12.0;
}

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

  volume v = volone(zmax,a);
  mat ma(v, eps);
  ma.use_pml_left(zmax/6);
  ma.use_pml_right(zmax/6);

  fields f(&ma);
  f.use_real_fields();
  f.add_point_source(Ex, 0.25, 3.5, 0.0, 8.0, vec(zmax/6+0.3), 1.0);
  flux_plane *left = f.add_flux_plane(vec(zmax/3.0), vec(zmax/3.0));
  flux_plane *right = f.add_flux_plane(vec(zmax*2.0/3.0), vec(zmax*2.0/3.0));

  const double ttot = min(10.0 + 1e5/zmax,f.find_last_source());

  f.step();
  volume mid = volone(zmax/3,a);
  mid.origin = vec(zmax/3);
  double flux_left=0.0, flux_right=0.0;
  double delta_energy = f.energy_in_box(mid);
  while (f.time() < ttot) {
    f.step();
    flux_left  +=  (c/a)*left->flux();
    flux_right +=  (c/a)*right->flux();
  }
  delta_energy -= f.energy_in_box(mid);
  const double del = flux_left;
  const double der = flux_right - delta_energy;
  master_printf("  Delta E:\t%lg\n  Flux left:\t%lg\n  Flux right:\t%lg\n  Ratio:\t%lg\n",
                delta_energy, del, der, del/der);
  return compare(del, der, 0.05, "Flux");
}

double oned_flux(const fields &f, const vec &loc) {
  monitor_point p;
  f.get_point(&p, loc);
  master_printf("Ex is %lg\n", real(p.get_component(Ex)));
  master_printf("Hy is %lg\n", real(p.get_component(Hy)));
  return (real(p.get_component(Ex))*real(p.get_component(Hy)) +
          imag(p.get_component(Ex))*imag(p.get_component(Hy)))
    *(1.0/pi/8/c);
}

int cavity_1d(const double boxwidth, const double timewait,
              double eps(const vec &)) {
  const double zmax = 15.0;
  const double a = 10.0;
  const double gridpts = a*zmax;

  volume v = volone(zmax,a);
  mat ma(v, eps);
  ma.use_pml_left(2.0);
  ma.use_pml_right(2.0);

  fields f(&ma);
  f.use_real_fields();
  f.add_point_source(Ex, 0.25, 4.5, 0.0, 8.0, vec(zmax/6+0.3), 1.0e2);
  flux_plane *left  = f.add_flux_plane(vec(zmax*.5-boxwidth),
                                       vec(zmax*.5-boxwidth));
  flux_plane *right = f.add_flux_plane(vec(zmax*.5+boxwidth),
                                       vec(zmax*.5+boxwidth));
  volume mid = volone(2*boxwidth,a);
  mid.origin = vec(zmax*.5-boxwidth-0.25/a);

  while (f.time() < f.find_last_source()) f.step();
  const double ttot = f.time() + timewait;
  double flux_left=0.0, flux_right=0.0;
  const double start_energy = f.energy_in_box(mid);
  master_printf("  Energy starts at\t%lg\n", start_energy);
  while (f.time() < ttot) {
    f.step();
    flux_left  +=  (c/a)*left->flux();
    flux_right +=  (c/a)*right->flux();
  }
  const double delta_energy = start_energy - f.energy_in_box(mid);
  const double defl = flux_right - flux_left;
  master_printf("  Delta E:         \t%lg\n  Integrated Flux:\t%lg\n",
                delta_energy, defl);
  master_printf("  Ratio:         \t%lg\n", delta_energy/defl);
  master_printf("  Fractional error:\t%lg\n",
                (delta_energy - defl)/start_energy);
  return compare(start_energy - delta_energy,
                 start_energy - defl, 0.001, "Flux");
}

void attempt(const char *name, int allright) {
  if (allright) master_printf("Passed %s\n", name);
  else abort("Failed %s!\n", name);
}

int main(int argc, char **argv) {
  initialize(argc, argv);
  master_printf("Trying out the fluxes...\n");

  attempt("Cavity 1D 3.01 73", cavity_1d(3.01, 73.0, cavity));
  attempt("Cavity 1D 5.0   1", cavity_1d(5.0, 1.0, cavity));
  attempt("Cavity 1D 3.85 55", cavity_1d(3.85, 55.0, cavity));

  width = 20.0;
  attempt("Flux 1D 20", flux_1d(100.0, bump));
  width = 10.0;
  attempt("Flux 1D 10", flux_1d(100.0, bump));
  width = 300.0;
  attempt("Flux 1D 300", flux_1d(100.0, bump));

  finished();
  exit(0);
}

