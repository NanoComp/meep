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

#include "config.h"

double one(const vec &) { return 1.0; }

double rods(const vec &r) {
  vec p = r;
  while (p.x() < -0.5)
    p.set_direction(X, p.x() + 1.0);
  while (p.x() > 0.5)
    p.set_direction(X, p.x() - 1.0);
  while (p.y() < -0.5)
    p.set_direction(Y, p.y() + 1.0);
  while (p.y() > 0.5)
    p.set_direction(Y, p.y() - 1.0);
  if (p.x() * p.x() + p.y() * p.y() < 0.3) return 12.0;
  return 1.0;
}

void compare(double b, double a, const char *n) {
  double thresh = sizeof(realnum) == sizeof(float) ? 1e-4 : 1e-5;
  if (fabs(a - b) > fabs(b) * thresh || b != b) {
    meep::abort("Failed %s (%g instead of %g, relerr %0.2g)\n", n, a, b, fabs(a - b) / fabs(b));
  }
  else { master_printf("Passed %s\n", n); }
}

static double dpml = 1.0;

double using_pml_ez(const grid_volume &gv, double eps(const vec &)) {
  const double ttot = 30.0;
  structure s(gv, eps, pml(dpml));
  fields f(&s);
  f.add_point_source(Ez, 0.2, 3.0, 0.0, 2.0, gv.center(), complex<double>(0, -2 * pi * 0.2));
  while (f.round_time() < ttot)
    f.step();
  monitor_point p;
  f.get_point(&p, gv.center());
  return real(p.get_component(Ez));
}

double x_periodic_y_pml(const grid_volume &gv, double eps(const vec &)) {
  const double ttot = 30.0;
  structure s(gv, eps, pml(dpml, Y));
  fields f(&s);
  f.add_point_source(Ez, 0.2, 3.0, 0.0, 2.0, gv.center(), complex<double>(0, -2 * pi * 0.2));
  f.use_bloch(X, 0.1);
  while (f.round_time() < ttot)
    f.step();
  monitor_point p;
  f.get_point(&p, gv.center());
  return real(p.get_component(Ez));
}

double x_periodic(const grid_volume &gv, double eps(const vec &)) {
  const double ttot = 30.0;
  structure s(gv, eps);
  fields f(&s);
  f.add_point_source(Ez, 0.2, 3.0, 0.0, 2.0, gv.center(), complex<double>(0, -2 * pi * 0.2));
  f.use_bloch(X, 0.1);
  while (f.round_time() < ttot)
    f.step();
  monitor_point p;
  f.get_point(&p, gv.center());
  return real(p.get_component(Ez));
}

double periodic_ez(const grid_volume &gv, double eps(const vec &)) {
  const double ttot = 30.0;
  structure s(gv, eps);
  fields f(&s);
  f.add_point_source(Ez, 0.2, 3.0, 0.0, 2.0, gv.center(), complex<double>(0, -2 * pi * 0.2));
  vec k;
  switch (gv.dim) {
    case D1: k = vec(0.3); break;
    case D2: k = vec(0.3, 0.4); break;
    case D3: k = vec(0.3, 0.5, 0.8); break;
    case Dcyl: k = veccyl(0.3, 0.2); break;
  }
  f.use_bloch(k);
  while (f.round_time() < ttot)
    f.step();
  monitor_point p;
  f.get_point(&p, gv.center());
  return real(p.get_component(Ez));
}

double metallic_ez(const grid_volume &gv, double eps(const vec &)) {
  const double ttot = 10.0;
  structure s(gv, eps);
  fields f(&s);
  f.add_point_source(Ez, 0.2, 3.0, 0.0, 2.0, gv.center(), complex<double>(0, -2 * pi * 0.2));
  while (f.round_time() < ttot)
    f.step();
  monitor_point p;
  f.get_point(&p, gv.center());
  return real(p.get_component(Ez));
}

double sigma(const vec &) { return 7.63; }

double polariton_ex(const grid_volume &gv, double eps(const vec &)) {
  const double ttot = 10.0;
  structure s(gv, eps);
  s.add_susceptibility(sigma, E_stuff, lorentzian_susceptibility(0.3, 0.1));
  fields f(&s);
  f.add_point_source(Ex, 0.2, 3.0, 0.0, 2.0, gv.center(), complex<double>(0, -2 * pi * 0.2));
  while (f.round_time() < ttot)
    f.step();
  monitor_point p;
  f.get_point(&p, gv.center());
  return real(p.get_component(Ex));
}

double polariton_energy(const grid_volume &gv, double eps(const vec &)) {
  const double ttot = 10.0;
  structure s(gv, eps);
  s.add_susceptibility(sigma, E_stuff, lorentzian_susceptibility(0.3, 0.1));
  fields f(&s, 0);
  f.add_point_source(Ex, 0.2, 3.0, 0.0, 2.0, gv.center(), complex<double>(0, -2 * pi * 0.2));
  while (f.round_time() < ttot)
    f.step();
  return f.field_energy();
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 0;
  master_printf("Testing with some known results...\n");
  const double a = 10.0;

  compare(-0.0894851, polariton_ex(volone(1.0, a), one), "1D polariton");
  compare(0.0863443, polariton_energy(volone(1.0, a), one), "1D polariton energy");
  compare(5.20605, metallic_ez(voltwo(1.0, 1.0, a), one), "1x1 metallic 2D TM");
  compare(0.883776, using_pml_ez(voltwo(1.0 + 2 * dpml, 1.0 + 2 * dpml, a), one), "1x1 PML 2D TM");
  compare(0.110425, x_periodic(voltwo(1.0, 1.0, a), one), "1x1 X periodic 2D TM");
  compare(-4.78767, periodic_ez(voltwo(1.0, 3.0, a), rods), "1x1 fully periodic 2D TM rods");
  compare(1.12502, periodic_ez(voltwo(1.0, 3.0, a), one), "1x1 fully periodic 2D TM");
  compare(0.608815, x_periodic_y_pml(voltwo(1.0, 1.0 + 2 * dpml, a), one),
          "1x1 X periodic Y PML 2D TM");
  compare(-41.8057, metallic_ez(vol3d(1.0, 1.0, 1.0, a), one), "1x1x1 metallic 3D");
  compare(-100.758, x_periodic(vol3d(1.0, 1.0, 1.0, a), one), "1x1x1 X periodic 3D");
  compare(-101.398, x_periodic_y_pml(vol3d(1.0, 1.0 + 2 * dpml, 1.0, a), one),
          "1x1x1 X periodic Y PML 3D");
  compare(-103.844, periodic_ez(vol3d(1.0, 1.0, 1.0, a), rods), "1x1x1 fully periodic 3D rods");
  compare(-99.1618, periodic_ez(vol3d(1.0, 1.0, 1.0, a), one), "1x1x1 fully periodic 3D");

  return 0;
}
