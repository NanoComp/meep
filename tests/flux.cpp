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

#include "meep.h"

double one(const vec &) { return 1.0; }
static double width = 20.0;
double bump(const vec &v) { return (fabs(v.z()-50.0) > width)?1.0:12.0; }

double cavity(const vec &v) {
  const double zz = fabs(v.z() - 7.5) + 0.3001;
  if (zz > 5.0) return 1.0;
  if (zz < 2.0) return 1.0;
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

static inline double min(double a, double b) { return (a<b)?a:b; }

int flux_1d(const double zmax,
            double eps(const vec &)) {
  const double a = 10.0;
  const double gridpts = a*zmax;

  volume v = volone(zmax,a);
  mat ma(v, eps);
  ma.use_pml_everywhere(zmax/6);

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
  double delta_energy = f.energy_in_box(mid.surroundings());
  master_printf("Initial energy is %lg\n", f.energy_in_box(mid.surroundings()));
  master_printf("Initial electric energy is %lg\n",
                f.electric_energy_in_box(mid.surroundings()));
  while (f.time() < ttot) {
    f.step();
    flux_left  +=  -(c/a)*left->flux();
    flux_right +=  -(c/a)*right->flux();
  }
  delta_energy -= f.energy_in_box(mid.surroundings());
  master_printf("Final energy is %lg\n", f.energy_in_box(mid.surroundings()));
  master_printf("Final electric energy is %lg\n",
                f.electric_energy_in_box(mid.surroundings()));
  const double del = flux_left;
  const double der = flux_right - delta_energy;
  master_printf("  Delta E:\t%lg\n  Flux left:\t%lg\n  Flux right:\t%lg\n  Ratio:\t%lg\n",
                delta_energy, del, der, del/der);
  return compare(del, der, 0.06, "Flux");
}

int split_1d(double eps(const vec &), int splitting) {
  const double boxwidth = 5.0, timewait = 1.0;
  const double zmax = 15.0, a = 10.0, gridpts = a*zmax;

  volume v = volone(zmax,a);
  mat ma1(v, eps, 1);
  mat ma(v, eps, splitting);
  ma1.use_pml_everywhere(2.0);
  ma.use_pml_everywhere(2.0);

  fields f1(&ma1);
  fields f(&ma);
  f1.use_real_fields();
  f.use_real_fields();
  f1.add_point_source(Ex, 0.25, 4.5, 0.0, 8.0, vec(zmax/2+0.3), 1.0e2);
  f.add_point_source(Ex, 0.25, 4.5, 0.0, 8.0, vec(zmax/2+0.3), 1.0e2);
  flux_plane *left1  = f1.add_flux_plane(vec(zmax*.5-boxwidth),
                                         vec(zmax*.5-boxwidth));
  flux_plane *left  = f.add_flux_plane(vec(zmax*.5-boxwidth),
                                       vec(zmax*.5-boxwidth));
  volume mid = volone(2*boxwidth,a);
  mid.origin = vec(zmax*.5-boxwidth-0.25/a);

  const double ttot = f.find_last_source() + timewait;
  while (f.time() < ttot) {
    f1.step();
    f.step();
    if (!compare((c/a)*left1->flux(), (c/a)*left->flux(), 0.0, "Flux"))
      return 0;
  }
  return 1;
}

int cavity_1d(const double boxwidth, const double timewait,
              double eps(const vec &)) {
  const double zmax = 15.0;
  const double a = 10.0;
  const double gridpts = a*zmax;

  volume v = volone(zmax,a);
  mat ma(v, eps);
  ma.use_pml_everywhere(2.0);

  fields f(&ma);
  f.use_real_fields();
  f.add_point_source(Ex, 0.25, 4.5, 0.0, 8.0, vec(zmax/2+0.3), 1.0e2);
  flux_plane *left  = f.add_flux_plane(vec(zmax*.5-boxwidth),
                                       vec(zmax*.5-boxwidth));
  flux_plane *right = f.add_flux_plane(vec(zmax*.5+boxwidth),
                                       vec(zmax*.5+boxwidth));
  volume mid = volone(2*boxwidth,a);
  mid.origin = vec(zmax*.5-boxwidth-0.25/a);

  while (f.time() < f.find_last_source()) f.step();
  const double ttot = f.time() + timewait;
  double flux_left=0.0, flux_right=0.0;
  const double start_energy = f.energy_in_box(mid.surroundings());
  master_printf("  Energy starts at\t%lg\n", start_energy);
  while (f.time() < ttot) {
    f.step();
    flux_left  +=  -(c/a)*left->flux();
    flux_right +=  -(c/a)*right->flux();
  }
  const double delta_energy = start_energy - f.energy_in_box(mid.surroundings());
  const double defl = flux_right - flux_left;
  master_printf("  Delta E:         \t%lg\n  Integrated Flux:\t%lg\n",
                delta_energy, defl);
  master_printf("  Ratio:         \t%lg\n", delta_energy/defl);
  master_printf("  Fractional error:\t%lg\n",
                (delta_energy - defl)/start_energy);
  return compare(start_energy - delta_energy,
                 start_energy - defl,
                 (timewait>50)?0.03:0.003, "Flux"); // Yuck, problem with flux.
}

void attempt(const char *name, int allright) {
  if (allright) master_printf("Passed %s\n", name);
  else abort("Failed %s!\n", name);
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  master_printf("Trying out the fluxes...\n");

  attempt("Split flux plane split by 7...", split_1d(cavity, 7));

  attempt("Cavity 1D 6.01 73", cavity_1d(6.01, 137.0, cavity));
  attempt("Cavity 1D 5.0   1", cavity_1d(5.0, 1.0, cavity));
  attempt("Cavity 1D 3.85 55", cavity_1d(3.85, 55.0, cavity));

  width = 20.0;
  attempt("Flux 1D 20", flux_1d(100.0, bump));
  width = 10.0;
  attempt("Flux 1D 10", flux_1d(100.0, bump));
  width = 300.0;
  attempt("Flux 1D 300", flux_1d(100.0, bump));

  exit(0);
}

