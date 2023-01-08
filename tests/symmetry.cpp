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
vec the_center;
double rods_2d(const vec &pp) {
  vec p = pp - the_center;
  while (p.x() > 0.5)
    p -= vec(1.0, 0);
  while (p.x() < -0.5)
    p += vec(1.0, 0);
  while (p.y() > 0.5)
    p -= vec(0, 1.0);
  while (p.y() < -0.5)
    p += vec(0, 1.0);
  if (fabs(p.x()) < 0.314) return 12.0;
  if (fabs(p.y()) < 0.314) return 12.0;
  return 1.0;
}

#if MEEP_SINGLE
static double eps_compare = 1e-3;
static double thresh_compare = 1e-3;
#else
static double eps_compare = 1e-8;
static double thresh_compare = 1e-8;
#endif

static inline double max(double a, double b) { return a > b ? a : b; }

int compare(double a, double b, const char *n) {
  if (fabs(a - b) > fabs(b) * eps_compare && max(fabs(a), fabs(b)) > thresh_compare) {
    master_printf("%s = %g differs by %g from %g\n", n, a, a - b, b);
    master_printf("This gives a fractional error of %g\n", fabs(a - b) / fabs(b));
    return 0;
  }
  else { return 1; }
}

int compare_point(fields &f1, fields &f2, const vec &p) {
  monitor_point m1, m_test;
  f1.get_point(&m_test, p);
  f2.get_point(&m1, p);
  for (int i = 0; i < 10; i++) {
    component c = (component)i;
    if (f1.gv.has_field(c)) {
      complex<double> v1 = m_test.get_component(c), v2 = m1.get_component(c);
      if (!compare(real(v1), real(v2), "real part") ||
          !compare(imag(v1), imag(v2), "imaginary part")) {
        master_printf("%s differs by %g%+gi from %g%+gi\n", component_name(c), real(v2 - v1),
                      imag(v2 - v1), real(v2), imag(v2));
        master_printf("This comes out to a fractional error of %g\n", abs(v1 - v2) / abs(v2));
        master_printf("Right now I'm looking at ");
        LOOP_OVER_DIRECTIONS(p.dim, d) {
          master_printf("%s = %g, ", direction_name(d), p.in_direction(d));
        }
        master_printf("time %g\n", f1.time());
        return 0;
      }
    }
  }
  return 1;
}

void check_unequal_layout(const fields &f1, const fields &f2) {
  if (f1.equal_layout(f2) || !f1.equal_layout(f1) || !f2.equal_layout(f2))
    meep::abort("fields::equal_layout did not return expected result");
}

int test_cyl_metal_mirror(double eps(const vec &)) {
  master_printf("Testing Z mirror symmetry in Cylindrical...\n");
  double a = 8.0;
  double ttot = 3.0;

  const grid_volume gv = volcyl(1.0, 1.0, a);
  the_center = gv.center();
  const symmetry S = mirror(Z, gv);
  structure s(gv, eps, no_pml(), S);
  structure s1(gv, eps);

  fields f1(&s1);
  f1.add_point_source(Er, 0.7, 2.5, 0.0, 4.0, veccyl(0.5, 0.5));
  f1.add_point_source(Ep, 0.8, 0.6, 0.0, 4.0, veccyl(0.401, 0.5));
  fields f(&s);
  f.add_point_source(Er, 0.7, 2.5, 0.0, 4.0, veccyl(0.5, 0.5));
  f.add_point_source(Ep, 0.8, 0.6, 0.0, 4.0, veccyl(0.401, 0.5));
  check_unequal_layout(f, f1);
  double field_energy_check_time = 1.0;
  while (f.round_time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, veccyl(0.01, 0.5))) return 0;
    if (!compare_point(f, f1, veccyl(0.21, 0.5))) return 0;
    if (!compare_point(f, f1, veccyl(0.501, 0.5))) return 0;
    if (!compare_point(f, f1, veccyl(0.33, 0.46))) return 0;
    if (!compare_point(f, f1, veccyl(0.2, 0.2))) return 0;
    if (f.round_time() >= field_energy_check_time) {
      if (!compare(f.electric_energy_in_box(gv.surroundings()),
                   f1.electric_energy_in_box(gv.surroundings()), "electric energy"))
        return 0;
      if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                   f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy"))
        return 0;
      if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) return 0;
      field_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_cyl_metal_mirror_nonlinear(double eps(const vec &)) {
  master_printf("Testing Z mirror symmetry in Cylindrical...\n");
  double a = 16.0;
  double ttot = 3.0;

  const grid_volume gv = volcyl(1.0, 1.0, a);
  the_center = gv.center();
  const symmetry S = mirror(Z, gv);
  structure s(gv, eps, no_pml(), S);
  structure s1(gv, eps);
  s.set_chi3(one);
  s1.set_chi3(one);

  fields f1(&s1);
  f1.add_point_source(Er, 0.7, 2.5, 0.0, 4.0, veccyl(0.5, 0.5));
  // f1.add_point_source(Ep, 0.8, 0.6, 0.0, 4.0, veccyl(0.401,0.5));
  fields f(&s);
  f.add_point_source(Er, 0.7, 2.5, 0.0, 4.0, veccyl(0.5, 0.5));
  // f.add_point_source(Ep, 0.8, 0.6, 0.0, 4.0, veccyl(0.401,0.5));
  check_unequal_layout(f, f1);
  double field_energy_check_time = 1.0;
  while (f.round_time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, veccyl(0.01, 0.5))) return 0;
    if (!compare_point(f, f1, veccyl(0.21, 0.5))) return 0;
    if (!compare_point(f, f1, veccyl(0.501, 0.5))) return 0;
    if (!compare_point(f, f1, veccyl(0.33, 0.46))) return 0;
    if (!compare_point(f, f1, veccyl(0.2, 0.2))) return 0;
    if (f.round_time() >= field_energy_check_time) {
      if (!compare(f.electric_energy_in_box(gv.surroundings()),
                   f1.electric_energy_in_box(gv.surroundings()), "electric energy"))
        return 0;
      if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                   f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy"))
        return 0;
      if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) return 0;
      field_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_1d_periodic_mirror(double eps(const vec &)) {
  master_printf("Testing Z mirror symmetry in 1D...\n");
  double a = 16.0;
  double ttot = 3.0;

  const grid_volume gv = volone(1.0, a);
  the_center = gv.center();
  const symmetry S = mirror(Z, gv);
  structure s(gv, eps, no_pml(), S);
  structure s1(gv, eps);

  fields f1(&s1);
  f1.use_bloch(0.0);
  f1.add_point_source(Ex, 0.7, 2.5, 0.0, 4.0, vec(0.5));
  fields f(&s);
  f.use_bloch(0.0);
  f.add_point_source(Ex, 0.7, 2.5, 0.0, 4.0, vec(0.5));
  check_unequal_layout(f, f1);
  double field_energy_check_time = 1.0;
  while (f.round_time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.01))) return 0;
    if (!compare_point(f, f1, vec(0.33))) return 0;
    if (!compare_point(f, f1, vec(0.50))) return 0;
    if (f.round_time() >= field_energy_check_time) {
      if (!compare(f.electric_energy_in_box(gv.surroundings()),
                   f1.electric_energy_in_box(gv.surroundings()), "electric energy"))
        return 0;
      if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                   f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy"))
        return 0;
      if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) return 0;
      field_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_origin_shift(void) {
  master_printf("Testing origin shift in 2D...\n");
  double a = 8.0;
  double ttot = 3.0;

  const grid_volume gv = voltwo(1.0, 1.0, a);
  grid_volume vcentered = gv;
  vcentered.shift_origin(-gv.center());
  structure s(vcentered, one);
  structure s1(gv, one);

  fields f1(&s1);
  fields f(&s);
  f1.add_point_source(Ey, 0.7, 2.5, 0.0, 4.0, gv.center());
  f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, gv.center());
  f.add_point_source(Ey, 0.7, 2.5, 0.0, 4.0, vec(0.0, 0.0));
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.0, 0.0));
  check_unequal_layout(f, f1);
  while (f.round_time() < ttot) {
    f.step();
    f1.step();
    if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) {
      master_printf("Time is %g\n", f.time());
      return 0;
    }
  }
  return 1;
}

int test_metal_xmirror(double eps(const vec &)) {
  master_printf("Testing X mirror symmetry...\n");
  double a = 8.0;
  double ttot = 3.0;

  const grid_volume gv = voltwo(1.0, 1.0, a);
  the_center = gv.center();
  const symmetry S = mirror(X, gv);
  structure s(gv, eps, no_pml(), S);
  structure s1(gv, eps);

  fields f1(&s1);
  f1.add_point_source(Ey, 0.7, 2.5, 0.0, 4.0, vec(0.5, 0.5));
  f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.5, 0.401));
  fields f(&s);
  f.add_point_source(Ey, 0.7, 2.5, 0.0, 4.0, vec(0.5, 0.5));
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.5, 0.401));
  check_unequal_layout(f, f1);
  double field_energy_check_time = 1.0;
  while (f.round_time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.5, 0.01))) return 0;
    if (!compare_point(f, f1, vec(0.5, 0.21))) return 0;
    if (!compare_point(f, f1, vec(0.5, 0.501))) return 0;
    if (!compare_point(f, f1, vec(0.46, 0.33))) return 0;
    if (!compare_point(f, f1, vec(0.2, 0.2))) return 0;
    if (f.round_time() >= field_energy_check_time) {
      if (!compare(f.electric_energy_in_box(gv.surroundings()),
                   f1.electric_energy_in_box(gv.surroundings()), "electric energy"))
        return 0;
      if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                   f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy"))
        return 0;
      if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) return 0;
      field_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_3D_metal_xmirror(double eps(const vec &)) {
  double a = 8.0;
  double ttot = 3.0;

  const grid_volume gv = vol3d(1.0, 1.0, 1.0, a);
  const symmetry S = mirror(X, gv);
  structure s(gv, eps, no_pml(), S);
  structure s1(gv, eps);
  master_printf("Testing X mirror symmetry in 3D...\n");

  fields f1(&s1);
  f1.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.5, 0.51, 0.55));
  f1.add_point_source(Hx, 0.8, 0.6, 0.0, 4.0, vec(0.5, 0.401, 0.43));
  fields f(&s);
  f.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.5, 0.51, 0.55));
  f.add_point_source(Hx, 0.8, 0.6, 0.0, 4.0, vec(0.5, 0.401, 0.43));
  check_unequal_layout(f, f1);
  double field_energy_check_time = 1.0;
  while (f.round_time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.5, 0.01, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.5, 0.21, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.5, 0.501, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.46, 0.33, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.2, 0.2, 0.5))) return 0;
    if (f.round_time() >= field_energy_check_time) {
      if (!compare(f.electric_energy_in_box(gv.surroundings()),
                   f1.electric_energy_in_box(gv.surroundings()), "electric energy"))
        return 0;
      if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                   f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy"))
        return 0;
      if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) return 0;
      field_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_3D_metal_zmirror(double eps(const vec &)) {
  double a = 8.0;
  double ttot = 3.0;

  const grid_volume gv = vol3d(1.1, 0.6, 1.0, a);
  const symmetry S = mirror(Z, gv);
  structure s(gv, eps, no_pml(), S);
  structure s1(gv, eps);
  master_printf("Testing Z mirror symmetry in 3D...\n");

  fields f1(&s1);
  f1.add_point_source(Ex, 0.7, 2.5, 0.0, 4.0, vec(0.55, 0.51, 0.5));
  f1.add_point_source(Ey, 0.8, 0.6, 0.0, 4.0, vec(0.43, 0.401, 0.5));
  fields f(&s);
  f.add_point_source(Ex, 0.7, 2.5, 0.0, 4.0, vec(0.55, 0.51, 0.5));
  f.add_point_source(Ey, 0.8, 0.6, 0.0, 4.0, vec(0.43, 0.401, 0.5));
  check_unequal_layout(f, f1);
  double field_energy_check_time = 1.0;
  while (f.round_time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.5, 0.01, 0.75))) return 0;
    if (!compare_point(f, f1, vec(0.5, 0.21, 0.15))) return 0;
    if (!compare_point(f, f1, vec(0.5, 0.501, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.46, 0.33, 0.51))) return 0;
    if (!compare_point(f, f1, vec(0.2, 0.2, 0.05))) return 0;
    if (f.round_time() >= field_energy_check_time) {
      if (!compare(f.electric_energy_in_box(gv.surroundings()),
                   f1.electric_energy_in_box(gv.surroundings()), "electric energy"))
        return 0;
      if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                   f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy"))
        return 0;
      if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) return 0;
      field_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_3D_metal_odd_zmirror(double eps(const vec &)) {
  double a = 8.0;
  double ttot = 3.0;

  const grid_volume gv = vol3d(1.1, 0.6, 1.0, a);
  const symmetry S = -mirror(Z, gv);
  structure s(gv, eps, no_pml(), S);
  structure s1(gv, eps);
  master_printf("Testing odd Z mirror symmetry in 3D...\n");

  fields f1(&s1);
  f1.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.55, 0.51, 0.5));
  fields f(&s);
  f.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.55, 0.51, 0.5));
  check_unequal_layout(f, f1);
  double field_energy_check_time = 1.0;
  while (f.round_time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.5, 0.01, 0.75))) return 0;
    if (!compare_point(f, f1, vec(0.5, 0.21, 0.15))) return 0;
    if (!compare_point(f, f1, vec(0.5, 0.501, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.46, 0.33, 0.51))) return 0;
    if (!compare_point(f, f1, vec(0.2, 0.2, 0.05))) return 0;
    if (f.round_time() >= field_energy_check_time) {
      if (!compare(f.electric_energy_in_box(gv.surroundings()),
                   f1.electric_energy_in_box(gv.surroundings()), "electric energy"))
        return 0;
      if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                   f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy"))
        return 0;
      if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) return 0;
      field_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_3D_metal_rot4z(double eps(const vec &)) {
  double a = 8.0;
  double ttot = 3.0;

  const grid_volume gv = vol3d(1.0, 1.0, 1.0, a);
  const symmetry S = rotate4(Z, gv);
  structure s(gv, eps, no_pml(), S);
  structure s1(gv, eps);
  master_printf("Testing Z fourfold rotational symmetry in 3D...\n");

  fields f1(&s1);
  f1.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.5, 0.5, 0.52));
  f1.add_point_source(Hz, 0.8, 0.6, 0.0, 4.0, vec(0.5, 0.5, 0.43));
  fields f(&s);
  f.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.5, 0.5, 0.52));
  f.add_point_source(Hz, 0.8, 0.6, 0.0, 4.0, vec(0.5, 0.5, 0.43));
  check_unequal_layout(f, f1);
  double field_energy_check_time = 1.0;
  while (f.round_time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.5, 0.01, 0.75))) return 0;
    if (!compare_point(f, f1, vec(0.5, 0.21, 0.15))) return 0;
    if (!compare_point(f, f1, vec(0.5, 0.501, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.46, 0.33, 0.51))) return 0;
    if (!compare_point(f, f1, vec(0.2, 0.2, 0.05))) return 0;
    if (f.round_time() >= field_energy_check_time) {
      if (!compare(f.electric_energy_in_box(gv.surroundings()),
                   f1.electric_energy_in_box(gv.surroundings()), "electric energy"))
        return 0;
      if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                   f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy"))
        return 0;
      if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) return 0;
      field_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_3D_metal_rot4z_mirror(double eps(const vec &)) {
  double a = 8.0;
  double ttot = 3.0;

  const grid_volume gv = vol3d(1.0, 1.0, 1.0, a);
  const symmetry S = rotate4(Z, gv) + mirror(Z, gv);
  structure s(gv, eps, no_pml(), S);
  structure s1(gv, eps);
  master_printf("Testing Z fourfold rotational symmetry in 3D with horizontal mirror...\n");

  fields f1(&s1);
  f1.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec(0.5, 0.5, 0.5));
  fields f(&s);
  f.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec(0.5, 0.5, 0.5));
  check_unequal_layout(f, f1);
  double field_energy_check_time = 1.0;
  while (f.round_time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.5, 0.01, 0.75))) return 0;
    if (!compare_point(f, f1, vec(0.5, 0.21, 0.15))) return 0;
    if (!compare_point(f, f1, vec(0.5, 0.501, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.46, 0.33, 0.51))) return 0;
    if (!compare_point(f, f1, vec(0.2, 0.2, 0.05))) return 0;
    if (f.round_time() >= field_energy_check_time) {
      if (!compare(f.electric_energy_in_box(gv.surroundings()),
                   f1.electric_energy_in_box(gv.surroundings()), "electric energy"))
        return 0;
      if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                   f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy"))
        return 0;
      if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) return 0;
      field_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_3D_metal_3mirror(double eps(const vec &)) {
  double a = 8.0;
  double ttot = 3.0;

  const grid_volume gv = vol3d(1.0, 1.0, 1.0, a);
  const symmetry S = mirror(Z, gv) - mirror(Y, gv) - mirror(X, gv);
  structure s(gv, eps, no_pml(), S);
  structure s1(gv, eps);
  master_printf("Testing three mirror planes in 3D...\n");

  fields f1(&s1);
  f1.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec(0.5, 0.5, 0.5));
  fields f(&s);
  f.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec(0.5, 0.5, 0.5));
  check_unequal_layout(f, f1);
  double field_energy_check_time = 1.0;
  while (f.round_time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.5, 0.01, 0.75))) return 0;
    if (!compare_point(f, f1, vec(0.5, 0.21, 0.15))) return 0;
    if (!compare_point(f, f1, vec(0.5, 0.501, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.46, 0.33, 0.51))) return 0;
    if (!compare_point(f, f1, vec(0.2, 0.2, 0.05))) return 0;
    if (f.round_time() >= field_energy_check_time) {
      if (!compare(f.electric_energy_in_box(gv.surroundings()),
                   f1.electric_energy_in_box(gv.surroundings()), "electric energy"))
        return 0;
      if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                   f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy"))
        return 0;
      if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) return 0;
      field_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_metal_ymirror(double eps(const vec &)) {
  double a = 8.0;
  double ttot = 5.0;

  const grid_volume gv = voltwo(1.0, 1.0, a);
  the_center = gv.center();
  const symmetry S = mirror(Y, gv);
  structure s(gv, eps, no_pml(), S);
  structure s1(gv, eps);
  master_printf("Testing Y mirror symmetry...\n");

  fields f1(&s1);
  f1.add_point_source(Ex, 0.7, 2.5, 0.0, 4.0, vec(0.85, 0.5));
  f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.401, 0.5));
  fields f(&s);
  f.add_point_source(Ex, 0.7, 2.5, 0.0, 4.0, vec(0.85, 0.5));
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.401, 0.5));
  check_unequal_layout(f, f1);
  double field_energy_check_time = 1.0;
  while (f.round_time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.01, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.21, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.46, 0.33))) return 0;
    if (!compare_point(f, f1, vec(0.2, 0.2))) return 0;
    if (f.round_time() >= field_energy_check_time) {
      if (!compare(f.electric_energy_in_box(gv.surroundings()),
                   f1.electric_energy_in_box(gv.surroundings()), "electric energy"))
        return 0;
      if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                   f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy"))
        return 0;
      if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) return 0;
      field_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_yperiodic_ymirror(double eps(const vec &)) {
  double a = 8.0;
  double ttot = 5.0;

  const grid_volume gv = voltwo(1.0, 1.0, a);
  the_center = gv.center();
  const symmetry S = mirror(Y, gv);
  structure s(gv, eps, no_pml(), S);
  structure s1(gv, eps);
  master_printf("Testing Y periodic with mirror symmetry...\n");

  fields f1(&s1);
  f1.use_bloch(vec(0.1 * pi / 2, 0.0));
  // f1.add_point_source(Ex, 0.7, 2.5, 0.0, 4.0, vec(0.85 ,0.5));
  f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.401, 0.5));
  fields f(&s);
  f.use_bloch(vec(0.1 * pi / 2, 0.0));
  // f.add_point_source(Ex, 0.7, 2.5, 0.0, 4.0, vec(0.85 ,0.5));
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.401, 0.5));
  check_unequal_layout(f, f1);
  double field_energy_check_time = 1.0;
  while (f.round_time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.951, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.01, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.21, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.46, 0.33))) return 0;
    if (!compare_point(f, f1, vec(0.2, 0.2))) return 0;
    if (f.round_time() >= field_energy_check_time) {
      if (!compare(f.electric_energy_in_box(gv.surroundings()),
                   f1.electric_energy_in_box(gv.surroundings()), "electric energy")) {
        return 0;
      }
      if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                   f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy"))
        return 0;
      if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) return 0;
      field_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_metal_rot2y(double eps(const vec &)) {
  double a = 16.0;
  double ttot = 5.0;

  const grid_volume gv = voltwo(1.0, 1.0, a);
  the_center = gv.center();
  const symmetry S = rotate2(Y, gv);
  structure s(gv, eps, no_pml(), S);
  structure s1(gv, eps);
  master_printf("Testing Y twofold rotational symmetry...\n");

  fields f1(&s1);
  f1.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec(0.25, 0.875), 1.0);
  f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.25, 0.375), 1.0);
  f1.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec(0.75, 0.875), -1.0);
  f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.75, 0.375), -1.0);
  fields f(&s);
  f.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec(0.25, 0.875), 1.0);
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.25, 0.375), 1.0);
  f.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec(0.75, 0.875), -1.0);
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.75, 0.375), -1.0);
  check_unequal_layout(f, f1);
  double field_energy_check_time = 1.0;
  while (f.round_time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.01, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.21, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.46, 0.33))) return 0;
    if (!compare_point(f, f1, vec(0.2, 0.2))) return 0;
    if (f.round_time() >= field_energy_check_time) {
      if (!compare(f.electric_energy_in_box(gv.surroundings()),
                   f1.electric_energy_in_box(gv.surroundings()), "electric energy"))
        return 0;
      if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                   f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy"))
        return 0;
      if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) return 0;
      field_energy_check_time += 1.0;
    }
  }
  return 1;
}

int exact_metal_rot2y(double eps(const vec &)) {
  double a = 16.0;
  double ttot = 5.0;

  const grid_volume gv = voltwo(1.0, 1.5, a);
  the_center = gv.center();
  const symmetry S = rotate2(Y, gv);
  structure s(gv, eps, no_pml(), S);
  structure s1(gv, eps);
  master_printf("Testing exact Y twofold rotational symmetry...\n");

  fields f1(&s1);
  f1.add_point_source(Ey, 0.7, 2.5, 0.0, 4.0, vec(0.5, 0.875));
  f1.add_point_source(Hy, 0.8, 0.6, 0.0, 4.0, vec(0.5, 0.375));
  fields f(&s);
  f.add_point_source(Ey, 0.7, 2.5, 0.0, 4.0, vec(0.5, 0.875));
  f.add_point_source(Hy, 0.8, 0.6, 0.0, 4.0, vec(0.5, 0.375));
  check_unequal_layout(f, f1);
  double field_energy_check_time = 1.0;
  while (f.round_time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.01, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.21, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.46, 0.33))) return 0;
    if (!compare_point(f, f1, vec(0.2, 0.2))) return 0;
    if (f.round_time() >= field_energy_check_time) {
      if (!compare(f.electric_energy_in_box(gv.surroundings()),
                   f1.electric_energy_in_box(gv.surroundings()), "electric energy"))
        return 0;
      if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                   f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy"))
        return 0;
      if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) return 0;
      field_energy_check_time += 1.0;
    }
  }
  return 1;
}

int pml_twomirrors(double eps(const vec &)) {
  double a = 16.0;
  double ttot = 10.0;

  const grid_volume gv = voltwo(2.0, 2.0, a);
  the_center = gv.center();
  const symmetry S = mirror(X, gv) + mirror(Y, gv);

  structure s_mm(gv, eps, pml(0.5), S);
  structure s1(gv, eps, pml(0.5), identity());
  master_printf("Testing two mirrors with PML...\n");
  fields f_mm(&s_mm);
  fields f1(&s1);

  f_mm.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(1.0, 1.0), -1.5);
  f_mm.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.75, 0.75));
  f_mm.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.75, 1.25));
  f_mm.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(1.25, 0.75));
  f_mm.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(1.25, 1.25));

  f1.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(1.0, 1.0), -1.5);
  f1.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.75, 0.75));
  f1.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.75, 1.25));
  f1.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(1.25, 0.75));
  f1.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(1.25, 1.25));

  check_unequal_layout(f_mm, f1);
  double field_energy_check_time = 3.0;
  while (f_mm.round_time() < ttot) {
    f_mm.step();
    f1.step();
    if (!compare_point(f1, f_mm, vec(0.01, 0.5))) return 0;
    if (!compare_point(f1, f_mm, vec(0.21, 0.5))) return 0;
    if (!compare_point(f1, f_mm, vec(0.46, 0.33))) return 0;
    if (!compare_point(f1, f_mm, vec(0.2, 0.2))) return 0;
    if (f_mm.round_time() >= field_energy_check_time) {
      if (!compare(f_mm.electric_energy_in_box(gv.surroundings()),
                   f1.electric_energy_in_box(gv.surroundings()), "electric energy"))
        return 0;
      field_energy_check_time += 3.0;
    }
  }
  return 1;
}

int exact_metal_rot4z(double eps(const vec &)) {
  double a = 8.0;
  double ttot = 5.0;

  const grid_volume gv = voltwo(1.0, 1.0, a);
  the_center = gv.center();
  const symmetry S = rotate4(Z, gv);

  structure s(gv, eps, no_pml(), S);
  structure s1(gv, eps);
  master_printf("Testing Z fourfold rotational symmetry...\n");

  fields f1(&s1);
  f1.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.5, 0.5));
  f1.add_point_source(Hz, 0.8, 0.6, 0.0, 4.0, vec(0.5, 0.5));
  fields f(&s);
  f.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.5, 0.5));
  f.add_point_source(Hz, 0.8, 0.6, 0.0, 4.0, vec(0.5, 0.5));
  check_unequal_layout(f, f1);
  double field_energy_check_time = 1.0;
  while (f.round_time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.01, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.21, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.46, 0.33))) return 0;
    if (!compare_point(f, f1, vec(0.2, 0.2))) return 0;
    if (f.round_time() >= field_energy_check_time) {
      if (!compare(f.electric_energy_in_box(gv.surroundings()),
                   f1.electric_energy_in_box(gv.surroundings()), "electric energy"))
        return 0;
      if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                   f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy"))
        return 0;
      if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) return 0;
      field_energy_check_time += 1.0;
    }
  }
  return 1;
}

int exact_metal_rot4z_nonlinear(double eps(const vec &)) {
  double a = 8.0;
  double ttot = 5.0;

  const grid_volume gv = voltwo(1.0, 1.0, a);
  the_center = gv.center();
  const symmetry S = rotate4(Z, gv);

  structure s(gv, eps, no_pml(), S);
  structure s1(gv, eps);
  s.set_chi3(one);
  s1.set_chi3(one);
  master_printf("Testing nonlinear Z fourfold rotational symmetry...\n");

  fields f1(&s1);
  // f1.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.5,0.5));
  f1.add_point_source(Hz, 0.8, 0.6, 0.0, 4.0, vec(0.5, 0.5));
  fields f(&s);
  // f.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.5,0.5));
  f.add_point_source(Hz, 0.8, 0.6, 0.0, 4.0, vec(0.5, 0.5));
  check_unequal_layout(f, f1);
  double field_energy_check_time = 1.0;
  while (f.round_time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.01, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.21, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.46, 0.33))) return 0;
    if (!compare_point(f, f1, vec(0.2, 0.2))) return 0;
    if (f.round_time() >= field_energy_check_time) {
      if (!compare(f.electric_energy_in_box(gv.surroundings()),
                   f1.electric_energy_in_box(gv.surroundings()), "electric energy"))
        return 0;
      if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                   f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy"))
        return 0;
      if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) return 0;
      field_energy_check_time += 1.0;
    }
  }
  return 1;
}

int exact_pml_rot2x_tm(double eps(const vec &)) {
  double a = 8.0;
  double ttot = 30.0;

  const grid_volume gv = voltwo(3.0, 3.0, a);
  the_center = gv.center();
  const symmetry S = rotate2(X, gv);

  structure s(gv, eps, pml(1.0), S);
  structure s1(gv, eps, pml(1.0), identity());
  master_printf("Testing X twofold rotational symmetry with PML...\n");

  fields f1(&s1);
  f1.add_point_source(Hx, 0.7, 2.5, 0.0, 4.0, vec(1.3, 1.5));
  fields f(&s);
  f.add_point_source(Hx, 0.7, 2.5, 0.0, 4.0, vec(1.3, 1.5));
  check_unequal_layout(f, f1);
  double field_energy_check_time = 1.0;
  while (f.round_time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.01, 1.5))) return 0;
    if (!compare_point(f, f1, vec(1.21, 1.5))) return 0;
    if (!compare_point(f, f1, vec(1.46, 0.33))) return 0;
    if (!compare_point(f, f1, vec(1.2, 1.2))) return 0;
    if (f.round_time() >= field_energy_check_time) {
      if (!compare(f.electric_energy_in_box(gv.surroundings()),
                   f1.electric_energy_in_box(gv.surroundings()), "electric energy"))
        return 0;
      if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                   f1.magnetic_energy_in_box(gv.surroundings()), "magnetic energy"))
        return 0;
      if (!compare(f.field_energy(), f1.field_energy(), "   total energy")) return 0;
      field_energy_check_time += 1.0;
    }
  }
  return 1;
}

double sigma(const vec &) { return 7.63; }

double polariton_ex(const grid_volume &gv, double eps(const vec &)) {
  const double ttot = 10.0;
  master_printf("Testing polariton in %s...\n", dimension_name(gv.dim));
  the_center = gv.center();
  const symmetry S = mirror(Z, gv);
  structure s(gv, eps);
  structure sS(gv, eps, no_pml(), S);
  s.add_susceptibility(sigma, E_stuff, lorentzian_susceptibility(0.3, 0.1));
  sS.add_susceptibility(sigma, E_stuff, lorentzian_susceptibility(0.3, 0.1));
  fields f(&s);
  f.use_real_fields();
  f.add_point_source(Ex, 0.2, 3.0, 0.0, 2.0, gv.center());
  fields fS(&sS);
  fS.use_real_fields();
  fS.add_point_source(Ex, 0.2, 3.0, 0.0, 2.0, gv.center());
  f.use_bloch(zero_vec(gv.dim));
  fS.use_bloch(zero_vec(gv.dim));
  check_unequal_layout(f, fS);
  while (f.round_time() < ttot) {
    f.step();
    fS.step();
    if (!compare_point(fS, f, gv.center())) return 0;
    if (!compare_point(fS, f, zero_vec(gv.dim))) return 0;
    if (!compare_point(fS, f, gv.center() * 0.3)) return 0;
  }
  return 1;
}

double nonlinear_ex(const grid_volume &gv, double eps(const vec &)) {
  const double ttot = 10.0;
  master_printf("Testing nonlinear in %s...\n", dimension_name(gv.dim));
  the_center = gv.center();
  const symmetry S = mirror(Z, gv);
  structure s(gv, eps);
  structure sS(gv, eps, no_pml(), S);
  s.set_chi3(one);
  sS.set_chi3(one);
  fields f(&s);
  f.use_real_fields();
  f.add_point_source(Ex, 0.2, 3.0, 0.0, 2.0, gv.center());
  fields fS(&sS);
  fS.use_real_fields();
  fS.add_point_source(Ex, 0.2, 3.0, 0.0, 2.0, gv.center());
  f.use_bloch(zero_vec(gv.dim));
  fS.use_bloch(zero_vec(gv.dim));
  check_unequal_layout(f, fS);
  while (f.round_time() < ttot) {
    f.step();
    fS.step();
    if (!compare_point(fS, f, gv.center())) return 0;
    if (!compare_point(fS, f, zero_vec(gv.dim))) return 0;
    if (!compare_point(fS, f, gv.center() * 0.3)) return 0;
  }
  return 1;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 0;
  master_printf("Testing with various kinds of symmetry...\n");

  if (!test_1d_periodic_mirror(one)) meep::abort("error in test_1d_periodic_mirror vacuum\n");

  // disable for parallel runs due to a bug in splitting cylindrical cell
  // with z-mirror symmetry into multiple chunks
  if ((count_processors() == 1) && !test_cyl_metal_mirror(one))
    meep::abort("error in test_cyl_metal_mirror vacuum\n");

  if (!test_yperiodic_ymirror(one)) meep::abort("error in test_yperiodic_ymirror vacuum\n");
  if (!test_yperiodic_ymirror(rods_2d)) meep::abort("error in test_yperiodic_ymirror rods2d\n");

  if (!pml_twomirrors(one)) meep::abort("error in pml_twomirrors vacuum\n");

  if (!test_origin_shift()) meep::abort("error in test_origin_shift\n");

  if (!exact_pml_rot2x_tm(one)) meep::abort("error in exact_pml_rot2x_tm vacuum\n");

  if (!test_metal_xmirror(rods_2d)) meep::abort("error in test_metal_xmirror rods_2d\n");

  if (!test_metal_xmirror(one)) meep::abort("error in test_metal_xmirror vacuum\n");
  if (!test_metal_ymirror(one)) meep::abort("error in test_metal_ymirror vacuum\n");
  if (!test_metal_ymirror(rods_2d)) meep::abort("error in test_metal_ymirror rods_2d\n");

  if (!test_metal_rot2y(one)) meep::abort("error in test_metal_rot2y vacuum\n");
  if (!test_metal_rot2y(rods_2d)) meep::abort("error in test_metal_rot2y rods_2d\n");

  if (!exact_metal_rot2y(one)) meep::abort("error in exact_metal_rot2y vacuum\n");
  if (!exact_metal_rot2y(rods_2d)) meep::abort("error in exact_metal_rot2y rods_2d\n");

  if (!exact_metal_rot4z(one)) meep::abort("error in exact_metal_rot4z vacuum\n");
  if (!exact_metal_rot4z(rods_2d)) meep::abort("error in exact_metal_rot4z rods_2d\n");

  if (!test_3D_metal_xmirror(one)) meep::abort("error in test_3D_metal_xmirror vacuum\n");

  if (!test_3D_metal_zmirror(one)) meep::abort("error in test_3D_metal_zmirror vacuum\n");

  if (!test_3D_metal_odd_zmirror(one)) meep::abort("error in test_3D_metal_odd_zmirror vacuum\n");

  if (!test_3D_metal_rot4z(one)) {
    all_wait();
    meep::abort("error in test_3D_metal_rot4z vacuum\n");
  }

  if (!test_3D_metal_rot4z_mirror(one)) meep::abort("error in test_3D_metal_rot4z_mirror vacuum\n");

  if (!test_3D_metal_3mirror(one)) meep::abort("error in test_3D_metal_3mirror\n");

    /**************************************************************************/
    /* For the following tests, we increase the check tolerance slightly.
       Floating-point errors can cause these tests to have slightly different
       results with and without symmetry.

       Note also that symmetry is tricky with nonlinearity, since in
       general a nonlinear system does *not* conserve the irreducible
       representation of the symmetry group (i.e. symmetry doesn't work).
       The simulations here are chosen to preserve the symmetry, however. */
#if !MEEP_SINGLE
  thresh_compare = 1e-10;
#endif

  if (!nonlinear_ex(vol1d(1.0, 30.0), one)) meep::abort("error in 1D nonlinear vacuum\n");
  if (sizeof(realnum) == sizeof(double)) {
    if (!nonlinear_ex(vol3d(1.0, 1.2, 0.8, 10.0), one))
      meep::abort("error in 3D nonlinear vacuum\n");
  }

  // disable for parallel runs due to a bug in splitting cylindrical cell
  // with z-mirror symmetry into multiple chunks
  if ((count_processors() == 1) && (!test_cyl_metal_mirror_nonlinear(one)))
    meep::abort("error in test_cyl_metal_mirror nonlinear vacuum\n");

  if (!exact_metal_rot4z_nonlinear(one))
    meep::abort("error in exact_metal_rot4z nonlinear vacuum\n");
  if (!exact_metal_rot4z_nonlinear(rods_2d))
    meep::abort("error in exact_metal_rot4z nonlinear rods_2d\n");

  // I'm not sure why the polariton tests require increased tolerances...?
  if (!polariton_ex(vol1d(1.0, 30.0), one)) meep::abort("error in 1D polariton vacuum\n");
  if (!polariton_ex(vol3d(1.0, 1.2, 0.8, 10.0), one)) meep::abort("error in 3D polariton vacuum\n");

  return 0;
}
