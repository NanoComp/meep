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

#include "dactyl.h"

double one(const vec &) { return 1.0; }

int compare(double a, double b, const char *n, double eps=4e-15) {
  if (fabs(a-b) > fabs(b)*eps) {
    printf("%s differs by\t%lg out of\t%lg\n", n, a-b, b);
    printf("This gives a fractional error of %lg\n", fabs(a-b)/fabs(b));
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
        printf("%s differs:  %lg %lg out of %lg %lg\n",
               component_name(c), real(v2-v1), imag(v2-v1), real(v2), imag(v2));
        printf("This comes out to a fractional error of %lg\n",
               abs(v1 - v2)/abs(v2));
        printf("Right now I'm looking at %lg %lg, time %lg\n", p.r(), p.z(), f1.time());
        f1.output_real_imaginary_slices("multi");
        f2.output_real_imaginary_slices("single");
        all_wait();
        return 0;
      }
    }
  }
  return 1;
}

int test_simple_periodic(double eps(const vec &), int splitting, const char *dirname) {
  double a = 10.0;
  double ttot = 30.0;
  
  volume v = volcyl(1.5,0.8,a);
  mat ma1(v, eps, 1);
  mat ma(v, eps, splitting);
  ma.set_output_directory(dirname);
  ma1.set_output_directory(dirname);
  for (int m=0;m<3;m++) {
    char m_str[10];
    snprintf(m_str, 10, "%d", m);
    master_printf("Trying with m = %d and a splitting into %d chunks...\n",
                  m, splitting);
    fields f(&ma, m);
    f.use_bloch(0.0);
    f.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, vec(0.5, 0.4), 1.0);
    f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.401, 0.301), 1.0);
    fields f1(&ma1, m);
    f1.use_bloch(0.0);
    f1.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, vec(0.5, 0.4), 1.0);
    f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.401, 0.301), 1.0);
    if (!compare(f1.count_volume(Ep), f.count_volume(Ep), "volume")) return 0;
    master_printf("Chunks are %lg by %lg\n",
                  f.chunks[0]->v.nr()/a, f.chunks[0]->v.nz()/a);
    double total_energy_check_time = 29.0;
    while (f.time() < ttot) {
      f.step();
      f1.step();
      if (!compare_point(f, f1, vec(0.5, 0.4))) return 0;
      if (!compare_point(f, f1, vec(0.46, 0.36))) return 0;
      if (!compare_point(f, f1, vec(1.0, 0.4))) return 0;
      if (!compare_point(f, f1, vec(0.01, 0.02))) return 0;
      if (!compare_point(f, f1, vec(0.601, 0.701))) return 0;
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
  }
  return 1;
}

int test_simple_metallic(double eps(const vec &), int splitting, const char *dirname) {
  double a = 10.0;
  double ttot = 30.0;
  
  volume v = volcyl(1.5,0.8,a);
  mat ma1(v, eps, 1);
  mat ma(v, eps, splitting);
  ma.set_output_directory(dirname);
  ma1.set_output_directory(dirname);
  for (int m=0;m<3;m++) {
    char m_str[10];
    snprintf(m_str, 10, "%d", m);
    master_printf("Metallic with m = %d and a splitting into %d chunks...\n",
                  m, splitting);
    fields f(&ma, m);
    f.use_metal_everywhere();
    f.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, vec(0.5, 0.4), 1.0);
    f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.401, 0.301), 1.0);
    fields f1(&ma1, m);
    f1.use_metal_everywhere();
    f1.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, vec(0.5, 0.4), 1.0);
    f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.401, 0.301), 1.0);
    if (!compare(f1.count_volume(Ep), f.count_volume(Ep), "volume")) return 0;
    master_printf("Chunks are %lg by %lg\n",
                  f.chunks[0]->v.nr()/a, f.chunks[0]->v.nz()/a);
    double total_energy_check_time = 29.0;
    while (f.time() < ttot) {
      f.step();
      f1.step();
      if (!compare_point(f, f1, vec(0.5, 0.4))) return 0;
      if (!compare_point(f, f1, vec(0.46, 0.36))) return 0;
      if (!compare_point(f, f1, vec(1.0, 0.4))) return 0;
      if (!compare_point(f, f1, vec(0.01, 0.02))) return 0;
      if (!compare_point(f, f1, vec(0.601, 0.701))) return 0;
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
  }
  return 1;
}

int test_pml(double eps(const vec &), int splitting, const char *dirname) {
  double a = 8;
  double ttot = 25.0;
  
  volume v = volcyl(3.5,10.0,a);
  mat ma1(v, eps, 1);
  mat ma(v, eps, splitting);
  ma.use_pml_everywhere(2.0);
  ma1.use_pml_everywhere(2.0);
  ma.set_output_directory(dirname);
  ma1.set_output_directory(dirname);
  for (int m=0;m<3;m++) {
    char m_str[10];
    snprintf(m_str, 10, "%d", m);
    master_printf("PML with m = %d and a splitting into %d chunks...\n",
                  m, splitting);
    fields f(&ma, m);
    f.use_metal_everywhere();
    f.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, vec(0.3, 7.0), 1.0);
    f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.3, 7.0), 1.0);
    fields f1(&ma1, m);
    f1.use_metal_everywhere();
    f1.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, vec(0.3, 7.0), 1.0);
    f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.3, 7.0), 1.0);
    f.eps_slices("multi");
    f1.eps_slices("single");
    if (!compare(f1.count_volume(Ep), f.count_volume(Ep), "volume", 3e-14)) return 0;
    master_printf("Chunks are %lg by %lg\n",
                  f.chunks[0]->v.nr()/a, f.chunks[0]->v.nz()/a);
    double total_energy_check_time = 10.0;
    while (f.time() < ttot) {
      f.step();
      f1.step();
      //f.output_real_imaginary_slices("multi");
      //f1.output_real_imaginary_slices("single");
      if (!compare_point(f, f1, vec(0.5, 7.0))) return 0;
      if (!compare_point(f, f1, vec(0.46, 0.36))) return 0;
      if (!compare_point(f, f1, vec(1.0, 0.4))) return 0;
      if (!compare_point(f, f1, vec(0.01, 0.02))) return 0;
      if (!compare_point(f, f1, vec(0.601, 0.701))) return 0;
      if (f.time() >= total_energy_check_time) {
        if (!compare(f.total_energy(), f1.total_energy(),
                     "pml total energy", 1e-13)) return 0;
        if (!compare(f.electric_energy_in_box(v.surroundings()),
                     f1.electric_energy_in_box(v.surroundings()),
                     "electric energy", 1e-13)) return 0;
        if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                     f1.magnetic_energy_in_box(v.surroundings()),
                     "magnetic energy", 1e-13)) return 0;

        total_energy_check_time += 10.0;
      }
    }
  }
  return 1;
}

complex<double> checkers(const vec &v) {
  const double ther = v.r() + 0.0001; // Just to avoid roundoff issues.
  const double thez = v.r() + 0.0001; // Just to avoid roundoff issues.
  int z = (int) (thez*5.0);
  int r = (int) (ther*5.0);
  int zz = (int) (thez*10.0);
  int rr = (int) (ther*10.0);
  if ((r & 1) ^ (z & 1)) return cos(thez*ther);
  if ((rr & 1) ^ (zz & 1)) return 1.0;
  return 0.0;
}

int test_pattern(double eps(const vec &), int splitting,
                 const char *dirname) {
  double a = 10.0;
  volume v = volcyl(1.5,0.8,a);
  mat ma1(v, eps, 1);
  mat ma(v, eps, splitting);
  ma.set_output_directory(dirname);
  ma1.set_output_directory(dirname);
  for (int m=0;m<1;m++) {
    char m_str[10];
    snprintf(m_str, 10, "%d", m);
    master_printf("Trying test pattern with m = %d and %d chunks...\n",
                  m, splitting);
    fields f(&ma, m);
    f.use_bloch(0.0);
    fields f1(&ma1, m);
    f1.use_bloch(0.0);
    if (!compare(f1.count_volume(Ep), f.count_volume(Ep), "volume")) return 0;
    master_printf("First chunk is %lg by %lg\n",
                  f.chunks[0]->v.nr()/a, f.chunks[0]->v.nz()/a);
    f1.initialize_field(Hp, checkers);
    f.initialize_field(Hp, checkers);

    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.751, 0.401))) return 0;
    if (!compare_point(f, f1, vec(0.01, 0.02))) return 0;
    if (!compare_point(f, f1, vec(1.0, 0.7))) return 0;
    if (!compare(f.total_energy(), f1.total_energy(),
                 "   total energy")) return 0;
    if (!compare(f.electric_energy_in_box(v.surroundings()),
                 f1.electric_energy_in_box(v.surroundings()),
                 "electric energy")) return 0;
    if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                 f1.magnetic_energy_in_box(v.surroundings()),
                 "magnetic energy")) return 0;
  }
  return 1;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  const char *dirname = "cylindrical-out";
  trash_output_directory(dirname);
  master_printf("Testing cylindrical coords under different splittings...\n");

  for (int s=2;s<8;s++)
    if (!test_pattern(one, s, dirname)) abort("error in test_pattern\n");
  if (!test_pattern(one, 8, dirname)) abort("error in crazy test_pattern\n");
  if (!test_pattern(one, 120, dirname)) abort("error in crazy test_pattern\n");
  
  for (int s=2;s<7;s++)
    if (!test_simple_periodic(one, s, dirname)) abort("error in test_simple_periodic\n");
  if (!test_simple_periodic(one, 8, dirname))
    abort("error in crazy test_simple_periodic\n");
  if (!test_simple_periodic(one, 120, dirname))
    abort("error in crazy test_simple_periodic\n");
  
  for (int s=2;s<7;s++)
    if (!test_simple_metallic(one, s, dirname)) abort("error in test_simple_metallic\n");
  if (!test_simple_metallic(one, 8, dirname))
    abort("error in crazy test_simple_metallic\n");
  if (!test_simple_metallic(one, 120, dirname))
    abort("error in crazy test_simple_metallic\n");
  
  for (int s=2;s<6;s++)
    if (!test_pml(one, s, dirname)) abort("error in test_pml\n");

  exit(0);
}

