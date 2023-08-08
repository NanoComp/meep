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
#include <signal.h>

#include <meep.hpp>
using namespace meep;
using std::complex;

double one(const vec &) { return 1.0; }

int compare(double a, double b, const char *n, double eps = 4e-15) {
  if (sizeof(realnum) == sizeof(float)) eps = sqrt(eps) * 10;
  if (fabs(a - b) > fabs(b) * eps && fabs(b) > 1e-14) {
    master_printf("%s differs by\t%g out of\t%g\n", n, a - b, b);
    master_printf("This gives a fractional error of %g\n", fabs(a - b) / fabs(b));
    return 0;
  }
  else { return 1; }
}

int compare_point(fields &f1, fields &f2, const vec &p, double eps = 4e-8) {
  if (sizeof(realnum) == sizeof(float)) eps = sqrt(eps);
  monitor_point m1, m_test;
  f1.get_point(&m_test, p);
  f2.get_point(&m1, p);
  for (int i = 0; i < 10; i++) {
    component c = (component)i;
    if (f1.gv.has_field(c)) {
      complex<double> v1 = m_test.get_component(c), v2 = m1.get_component(c);
      if (abs(v1 - v2) > eps * abs(v2) && abs(v2) > eps * 100) {
        master_printf("%s differs:  %g %g out of %g %g\n", component_name(c), real(v2 - v1),
                      imag(v2 - v1), real(v2), imag(v2));
        master_printf("This comes out to a fractional error of %g\n", abs(v1 - v2) / abs(v2));
        master_printf("Right now I'm looking at %g %g, time %g\n", p.r(), p.z(), f1.time());
        all_wait();
        return 0;
      }
    }
  }
  return 1;
}

int test_simple_periodic(double eps(const vec &), int splitting) {
  double a = 10.0;
  double ttot = 30.0;

  grid_volume gv = volcyl(1.5, 0.8, a);
  structure s1(gv, eps, no_pml(), identity(), 0, 0.4);
  structure s(gv, eps, no_pml(), identity(), splitting, 0.4);
  for (int m = 0; m < 3; m++) {
    char m_str[10];
    snprintf(m_str, 10, "%d", m);
    master_printf("Trying with m = %d and a splitting into %d chunks...\n", m, splitting);
    fields f(&s, m);
    f.use_bloch(0.0);
    f.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, veccyl(0.5, 0.4), 1.0);
    f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, veccyl(0.401, 0.301), 1.0);
    fields f1(&s1, m);
    f1.use_bloch(0.0);
    f1.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, veccyl(0.5, 0.4), 1.0);
    f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, veccyl(0.401, 0.301), 1.0);
    if (!compare(f1.count_volume(Ep), f.count_volume(Ep), "grid_volume")) return 0;
    master_printf("Chunks are %g by %g\n", f.chunks[0]->gv.nr() / a, f.chunks[0]->gv.nz() / a);
    double field_energy_check_time = 29.0;
    while (f.time() < ttot) {
      f.step();
      f1.step();
      if (!compare_point(f, f1, veccyl(0.5, 0.4))) return 0;
      if (!compare_point(f, f1, veccyl(0.46, 0.36))) return 0;
      if (!compare_point(f, f1, veccyl(1.0, 0.4))) return 0;
      if (!compare_point(f, f1, veccyl(0.01, 0.02), sizeof(realnum) == sizeof(float) ? 2e-6 : 4e-8))
        return 0;
      if (!compare_point(f, f1, veccyl(0.601, 0.701))) return 0;
      if (f.time() >= field_energy_check_time) {
        if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) return 0;
        if (!compare(f.electric_energy_in_box(gv.surroundings()),
                     f1.electric_energy_in_box(gv.surroundings()), "electric energy"))
          return 0;
        if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                     f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy"))
          return 0;

        field_energy_check_time += 5.0;
      }
    }
  }
  return 1;
}

int test_simple_metallic(double eps(const vec &), int splitting) {
  double a = 10.0;
  double ttot = 30.0;

  grid_volume gv = volcyl(1.5, 0.8, a);
  structure s1(gv, eps, no_pml(), identity(), 0, 0.4);
  structure s(gv, eps, no_pml(), identity(), splitting, 0.4);
  for (int m = 0; m < 3; m++) {
    char m_str[10];
    snprintf(m_str, 10, "%d", m);
    master_printf("Metallic with m = %d and a splitting into %d chunks...\n", m, splitting);
    fields f(&s, m);
    f.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, veccyl(0.5, 0.4), 1.0);
    f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, veccyl(0.401, 0.301), 1.0);
    fields f1(&s1, m);
    f1.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, veccyl(0.5, 0.4), 1.0);
    f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, veccyl(0.401, 0.301), 1.0);
    if (!compare(f1.count_volume(Ep), f.count_volume(Ep), "grid_volume")) return 0;
    master_printf("Chunks are %g by %g\n", f.chunks[0]->gv.nr() / a, f.chunks[0]->gv.nz() / a);
    double field_energy_check_time = 29.0;
    while (f.time() < ttot) {
      f.step();
      f1.step();
      if (!compare_point(f, f1, veccyl(0.5, 0.4))) return 0;
      if (!compare_point(f, f1, veccyl(0.46, 0.36))) return 0;
      if (!compare_point(f, f1, veccyl(1.0, 0.4))) return 0;
      if (!compare_point(f, f1, veccyl(0.01, 0.02), sizeof(realnum) == sizeof(float) ? 2e-6 : 4e-8))
        return 0;
      if (!compare_point(f, f1, veccyl(0.601, 0.701))) return 0;
      if (f.time() >= field_energy_check_time) {
        if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) return 0;
        if (!compare(f.electric_energy_in_box(gv.surroundings()),
                     f1.electric_energy_in_box(gv.surroundings()), "electric energy"))
          return 0;
        if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                     f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy"))
          return 0;

        field_energy_check_time += 5.0;
      }
    }
  }
  return 1;
}

static bool issmall(std::complex<double> x) { return abs(x) < 1e-16; }

int test_r_equals_zero(double eps(const vec &)) {
  double a = 10.0;
  double ttot = 3.0;
  grid_volume gv = volcyl(1.5, 0.8, a);
  structure s(gv, eps, no_pml(), identity(), 0, 0.4);
  for (int m = 0; m < 3; m++) {
    char m_str[10];
    snprintf(m_str, 10, "%d", m);
    master_printf("Checking at r == 0 with m = %d...\n", m);
    fields f(&s, m);
    f.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, veccyl(0.5, 0.4), 1.0);
    f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, veccyl(0.401, 0.301), 1.0);
    while (f.time() < ttot)
      f.step();
    monitor_point p;
    f.get_point(&p, veccyl(0.0, 0.5));
    if (!issmall(p.get_component(Ez)) && (m & 1)) {
      master_printf("Got non-zero Ez with m == %d\n", m);
      return 0;
    }
    if (!issmall(p.get_component(Hz)) && (m & 1)) {
      master_printf("Got non-zero Hz with m == %d\n", m);
      return 0;
    }
    if (!issmall(p.get_component(Er)) && !(m & 1)) {
      master_printf("Got non-zero Er with m == %d\n", m);
      return 0;
    }
    if (!issmall(p.get_component(Ep)) && !(m & 1)) {
      master_printf("Got non-zero Ep with m == %d\n", m);
      return 0;
    }
    if (!issmall(p.get_component(Hr)) && !(m & 1)) {
      master_printf("Got non-zero Hr with m == %d\n", m);
      return 0;
    }
    if (!issmall(p.get_component(Hp)) && !(m & 1)) {
      master_printf("Got non-zero Hp of %g %g with m == %d\n", real(p.get_component(Hp)),
                    imag(p.get_component(Hp)), m);
      return 0;
    }
  }
  return 1;
}

int test_pml(double eps(const vec &), int splitting) {
  double a = 8;
  double ttot = 25.0;

  grid_volume gv = volcyl(3.5, 10.0, a);
  structure s1(gv, eps, pml(2.0), identity(), 0, 0.4);
  structure s(gv, eps, pml(2.0), identity(), splitting, 0.4);
  for (int m = 0; m < 3; m++) {
    char m_str[10];
    snprintf(m_str, 10, "%d", m);
    master_printf("PML with m = %d and a splitting into %d chunks...\n", m, splitting);
    fields f(&s, m);
    f.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, veccyl(0.3, 7.0), 1.0);
    f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, veccyl(0.3, 7.0), 1.0);
    fields f1(&s1, m);
    f1.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, veccyl(0.3, 7.0), 1.0);
    f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, veccyl(0.3, 7.0), 1.0);
    if (!compare(f1.count_volume(Ep), f.count_volume(Ep), "grid_volume", 3e-14)) return 0;
    master_printf("Chunks are %g by %g\n", f.chunks[0]->gv.nr() / a, f.chunks[0]->gv.nz() / a);
    double field_energy_check_time = 10.0;
    while (f.time() < ttot) {
      f.step();
      f1.step();
      // f.output_real_imaginary_slices("multi");
      // f1.output_real_imaginary_slices("single");
      if (!compare_point(f, f1, veccyl(0.5, 7.0))) return 0;
      if (!compare_point(f, f1, veccyl(0.46, 0.36))) return 0;
      if (!compare_point(f, f1, veccyl(1.0, 0.4))) return 0;
      if (!compare_point(f, f1, veccyl(0.01, 0.02))) return 0;
      if (!compare_point(f, f1, veccyl(0.601, 0.701))) return 0;
      if (f.time() >= field_energy_check_time) {
        if (!compare(f.field_energy(), f1.field_energy(), "pml total energy", 1e-13)) return 0;
        if (!compare(f.electric_energy_in_box(gv.surroundings()),
                     f1.electric_energy_in_box(gv.surroundings()), "electric energy", 1e-13))
          return 0;
        if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                     f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy", 1e-13))
          return 0;

        field_energy_check_time += 10.0;
      }
    }
  }
  return 1;
}

complex<double> checkers(const vec &pt) {
  const double ther = pt.r() + 0.0001; // Just to avoid roundoff issues.
  const double thez = pt.r() + 0.0001; // Just to avoid roundoff issues.
  int z = (int)(thez * 5.0);
  int r = (int)(ther * 5.0);
  int zz = (int)(thez * 10.0);
  int rr = (int)(ther * 10.0);
  if ((r & 1) ^ (z & 1)) return cos(thez * ther);
  if ((rr & 1) ^ (zz & 1)) return 1.0;
  return 0.0;
}

int test_pattern(double eps(const vec &), int splitting) {
  double a = 10.0;
  grid_volume gv = volcyl(1.5, 0.8, a);
  structure s1(gv, eps);
  structure s(gv, eps, no_pml(), identity(), splitting);
  for (int m = 0; m < 1; m++) {
    char m_str[10];
    snprintf(m_str, 10, "%d", m);
    master_printf("Trying test pattern with m = %d and %d chunks...\n", m, splitting);
    fields f(&s, m);
    f.use_bloch(0.0);
    fields f1(&s1, m);
    f1.use_bloch(0.0);
    if (!compare(f1.count_volume(Ep), f.count_volume(Ep), "grid_volume")) return 0;
    master_printf("First chunk is %g by %g\n", f.chunks[0]->gv.nr() / a, f.chunks[0]->gv.nz() / a);
    f1.initialize_field(Hp, checkers);
    f.initialize_field(Hp, checkers);

    f.step();
    f1.step();
    if (!compare_point(f, f1, veccyl(0.751, 0.401))) return 0;
    if (!compare_point(f, f1, veccyl(0.01, 0.02))) return 0;
    if (!compare_point(f, f1, veccyl(1.0, 0.7))) return 0;
    if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) return 0;
    if (!compare(f.electric_energy_in_box(gv.surroundings()),
                 f1.electric_energy_in_box(gv.surroundings()), "electric energy"))
      return 0;
    if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                 f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy"))
      return 0;
  }
  return 1;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 0;
  master_printf("Testing cylindrical coords under different splittings...\n");

  if (!test_r_equals_zero(one)) meep::abort("error in test_r_equals_zero");

  for (int s = 2; s < 6; s++)
    if (!test_pattern(one, s)) meep::abort("error in test_pattern\n");
  // if (!test_pattern(one, 8)) meep::abort("error in crazy test_pattern\n");
  // if (!test_pattern(one, 120)) meep::abort("error in crazy test_pattern\n");

  for (int s = 2; s < 4; s++)
    if (!test_simple_periodic(one, s)) meep::abort("error in test_simple_periodic\n");
  // if (!test_simple_periodic(one, 8))
  //  meep::abort("error in crazy test_simple_periodic\n");
  // if (!test_simple_periodic(one, 120))
  //  meep::abort("error in crazy test_simple_periodic\n");

  for (int s = 2; s < 5; s++)
    if (!test_simple_metallic(one, s)) meep::abort("error in test_simple_metallic\n");
  // if (!test_simple_metallic(one, 8))
  //  meep::abort("error in crazy test_simple_metallic\n");
  // if (!test_simple_metallic(one, 120))
  //  meep::abort("error in crazy test_simple_metallic\n");

  for (int s = 2; s < 6; s++)
    if (!test_pml(one, s)) meep::abort("error in test_pml\n");

  return 0;
}
