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
#include <signal.h>

#include <meep.h>
using namespace meep;

double one(const vec &) { return 1.0; }

int compare(double a, double b, const char *n) {
  if (fabs(a-b) > fabs(b)*8e-15) {
    master_printf("%s differs by\t%lg out of\t%lg\n", n, a-b, b);
    master_printf("This gives a fractional error of %lg\n", fabs(a-b)/fabs(b));
    return 0;
  } else {
    return 1;
  }
}

int compare_point(fields &f1, fields &f2, const vec &p) {
  monitor_point m1, m_test;
  f1.get_point(&m_test, p);
  f2.get_point(&m1, p);
  for (int i=0;i<10;i++) {
    component c = (component) i;
    if (f1.v.has_field(c)) {
      complex<double> v1 = m_test.get_component(c), v2 = m1.get_component(c);
      if (abs(v1 - v2) > 0.0*2e-15*abs(v2)) {
        master_printf("%s differs:  %lg %lg out of %lg %lg\n",
               component_name(c), real(v2-v1), imag(v2-v1), real(v2), imag(v2));
        master_printf("This comes out to a fractional error of %lg\n",
               abs(v1 - v2)/abs(v2));
        master_printf("Right now I'm looking at %lg, time %lg\n", p.z(), f1.time());
        f1.output_real_imaginary_slices("multi");
        f2.output_real_imaginary_slices("single");
        return 0;
      }
    }
  }
  return 1;
}

int test_simple_periodic(double eps(const vec &), int splitting, const char *dirname) {
  double a = 10.0;
  double ttot = 170.0;
  
  volume v = volone(6.0,a);
  mat ma1(v, eps, 1);
  mat ma(v, eps, splitting);
  ma.set_output_directory(dirname);
  ma1.set_output_directory(dirname);

  master_printf("Trying splitting into %d chunks...\n", splitting);
  fields f(&ma);
  f.use_bloch(0.0);
  f.add_point_source(Hy, 0.7, 2.5, 0.0, 4.0, vec(0.5), 1.0);
  f.add_point_source(Ex, 0.8, 0.6, 0.0, 4.0, vec(0.401), 1.0);
  fields f1(&ma1);
  f1.use_bloch(0.0);
  f1.add_point_source(Hy, 0.7, 2.5, 0.0, 4.0, vec(0.5), 1.0);
  f1.add_point_source(Ex, 0.8, 0.6, 0.0, 4.0, vec(0.401), 1.0);
  if (!compare(f1.count_volume(Ex), f.count_volume(Ex), "volume")) return 0;
  double total_energy_check_time = 29.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.5  ))) return 0;
    if (!compare_point(f, f1, vec(0.46 ))) return 0;
    if (!compare_point(f, f1, vec(1.0  ))) return 0;
    if (!compare_point(f, f1, vec(0.01 ))) return 0;
    if (!compare_point(f, f1, vec(0.601))) return 0;
    if (f.time() >= total_energy_check_time) {
      if (!compare(f.total_energy(), f1.total_energy(),
                   "   total energy")) return 0;
      if (!compare(f.electric_energy_in_box(v.surroundings()),
                   f1.electric_energy_in_box(v.surroundings()),
                   "electric energy")) return 0;
      if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                   f1.magnetic_energy_in_box(v.surroundings()),
                   "magnetic energy")) return 0;
      
      total_energy_check_time += 5.0;
    }
  }
  return 1;
}

complex<double> checkers(const vec &v) {
  const double thez = v.z()+0.00001;
  int z = (int) (thez*5.0);
  int zz = (int) (thez*10.0);
  if (z & 1) return cos(thez);
  if (zz & 1) return 2.0;
  return 1.0;
}

int test_pattern(double eps(const vec &), int splitting,
                 const char *dirname) {
  double a = 10.0;
  volume v = volone(6.0,a);
  mat ma1(v, eps, 1);
  mat ma(v, eps, splitting);
  ma.set_output_directory(dirname);
  ma1.set_output_directory(dirname);

  master_printf("Trying test pattern with %d chunks...\n", splitting);
  fields f(&ma);
  f.use_bloch(0.0);
  fields f1(&ma1);
  f1.use_bloch(0.0);
  if (!compare(f1.count_volume(Ex), f.count_volume(Ex), "volume")) return 0;
  f1.initialize_field(Hy, checkers);
  f.initialize_field(Hy, checkers);

  f.step();
  f1.step();
  if (!compare_point(f, f1, vec(27.99))) return 0;
  if (!compare_point(f, f1, vec(42.01))) return 0;
  if (!compare_point(f, f1, vec(0.751))) return 0;
  if (!compare_point(f, f1, vec(0.01 ))) return 0;
  if (!compare_point(f, f1, vec(1.0  ))) return 0;
  if (!compare(f.total_energy(), f1.total_energy(),
               "   total energy")) return 0;
  if (!compare(f.electric_energy_in_box(v.surroundings()),
               f1.electric_energy_in_box(v.surroundings()),
               "electric energy")) return 0;
  if (!compare(f.magnetic_energy_in_box(v.surroundings()),
               f1.magnetic_energy_in_box(v.surroundings()),
               "magnetic energy")) return 0;
  return 1;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  const char *dirname = "one_dimensional-out";
  master_printf("Testing one dimension under different splittings...\n");

  for (int s=2;s<7;s++)
    if (!test_pattern(one, s, dirname)) abort("error in test_pattern\n");

  for (int s=2;s<7;s++)
    if (!test_simple_periodic(one, s, dirname)) abort("error in test_simple_periodic\n");
  exit(0);
}

