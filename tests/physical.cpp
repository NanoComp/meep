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

double slab(const vec &v) {
  abort("not yet\n");
}

int radiating_2D(const double xmax) {
  const double a = 10.0;
  const double gridpts = a*xmax;
  const double ymax = 3.0;

  volume v = voltwo(xmax,ymax,a);
  mat ma(v, one);
  ma.use_pml_everywhere(ymax/3);

  fields f(&ma);
  double w = 0.30;
  double dx = 2.0;
  f.add_point_source(Ez, w, 3.0, 0.0, 2.0, vec2d(xmax/2 - dx, ymax/2), 1.0, 1); //continuous
  const double t1 = f.find_last_source();

  // let the source reach steady state
  while (f.time() < t1)
    f.step();

  monitor_point p1, p2;
  f.get_point(&p1, vec2d(xmax/2, ymax/2));
  f.get_point(&p2, vec2d(xmax/2 + dx, ymax/2));

  complex<double> amp1 = p1.get_component(Ez);
  complex<double> amp2 = p2.get_component(Ez);

  double ratio = pow(abs(amp1)/abs(amp2), 2.0) ;
  
  if (ratio > 2.02 || ratio < 1.98)
    abort("Failed: amp1 = (%lg, %lg), amp2 = (%lg, %lg)\n abs(amp1/amp2)^2 = %lg, too far from 2.0\n",
	  real(amp1), imag(amp1), real(amp2), imag(amp2), ratio);
  return 1;
}

int radiating_3D() {
  const double a = 10.0;
  const double ymax = 3.0;
  const double w = 0.30;
  const double xmax = 12.0;
  const double dx = 2.0;

  volume v = vol3d(xmax,ymax,ymax,a);
  symmetry S = mirror(X,v) + mirror(Y,v) + mirror(Z,v)*(-1.0);
  mat ma(v, one, 0, S);
  ma.use_pml_everywhere(ymax/3);

  fields f(&ma);
  f.add_point_source(Ez, w, 3.0, 0.0, 2.0,
                     vec(xmax/2, ymax/2, ymax/2), 1.0, 1); //continuous
  const double t1 = f.find_last_source();

  // let the source reach steady state
  double next_print_time = 1.0;
  while (f.time() < t1*2 + 2*xmax + 2*ymax) {
    f.step();
    if (f.time() > next_print_time) {
      monitor_point p1, p2;
      f.get_point(&p1, vec(xmax/2 + dx, ymax/2, ymax/2));
      f.get_point(&p2, vec(xmax/2 + 2*dx, ymax/2, ymax/2));
      complex<double> amp1 = p1.get_component(Ez);
      complex<double> amp2 = p2.get_component(Ez);
      const double ratio = abs(amp1)/abs(amp2);
      //printf("At time %lg ratio is %lg from %lg and %lg\n",
      //       f.time(), ratio, abs(amp1), abs(amp2));
      next_print_time += 1.0;
    }
  }

  monitor_point p1, p2;
  f.get_point(&p1, vec(xmax/2 + dx, ymax/2, ymax/2));
  f.get_point(&p2, vec(xmax/2 + 2*dx, ymax/2, ymax/2));

  complex<double> amp1 = p1.get_component(Ez);
  complex<double> amp2 = p2.get_component(Ez);

  const double ratio = abs(amp1)/abs(amp2);
  printf("Ratio is %lg from (%lg %lg) and (%lg %lg)\n",
         ratio, real(amp1), imag(amp1), real(amp2), imag(amp2));
  if (ratio > 2.05 || ratio < 1.9)
    abort("Failed: amp1 = (%lg, %lg), amp2 = (%lg, %lg)\n abs(amp1/amp2)^2 = %lg, too far from 2.0\n",
	  real(amp1), imag(amp1), real(amp2), imag(amp2), ratio);
  return 1;
}

void attempt(const char *name, int allright) {
  if (allright) master_printf("Passed %s\n", name);
  else abort("Failed %s!\n", name);
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  master_printf("Trying out some physical tests...\n");

  attempt("radiating source should decay spatially as 1/sqrt(r) in 2D.", radiating_2D(8.0));
  attempt("radiating source should decay spatially as 1/r in 3D.", radiating_3D());
  return 0;
}

