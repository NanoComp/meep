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

const char *dirname = "symmetry-out";

double one(const vec &) { return 1.0; }
vec the_center;
double rods_2d(const vec &pp) {
  vec p = pp - the_center;
  while (p.x() > 0.5) p -= vec2d(1.0,0);
  while (p.x() <-0.5) p += vec2d(1.0,0);
  while (p.y() > 0.5) p -= vec2d(0,1.0);
  while (p.y() <-0.5) p += vec2d(0,1.0);
  if (fabs(p.x()) < 0.314) return 12.0;
  if (fabs(p.y()) < 0.314) return 12.0;
  return 1.0;
}

static const double eps_compare = 1.9e-14;
static const double thresh_compare = 1e-16;

int compare(double a, double b, const char *n) {
  if (fabs(a-b) > fabs(b)*eps_compare && fabs(b) > thresh_compare) {
    master_printf("%s differs by\t%g out of\t%g\n", n, a-b, b);
    master_printf("This gives a fractional error of %g\n", fabs(a-b)/fabs(b));
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
      if (abs(v1 - v2) > eps_compare*abs(v2) && abs(v2) > thresh_compare) {
        master_printf("%s differs:  %g %g out of %g %g\n",
               component_name(c), real(v2-v1), imag(v2-v1), real(v2), imag(v2));
        master_printf("This comes out to a fractional error of %g\n",
               abs(v1 - v2)/abs(v2));
        master_printf("Right now I'm looking at ");
        LOOP_OVER_DIRECTIONS(p.dim,d)
          master_printf("%s = %g, ", direction_name(d), p.in_direction(d));
        master_printf("time %g\n", f1.time());
        f1.output_real_imaginary_slices("multi");
        f2.output_real_imaginary_slices("single");
        f1.eps_slices("multi");
        f2.eps_slices("single");
        return 0;
      }
    }
  }
  return 1;
}

int test_cyl_metal_mirror(double eps(const vec &)) {
  master_printf("Testing Z mirror symmetry in Cylindrical...\n");
  double a = 10.0;
  double ttot = 3.0;

  const volume v = volcyl(1.0, 1.0, a);
  the_center = v.center();
  const symmetry S = mirror(Z,v);
  structure s(v, eps, 0, S);
  structure s1(v, eps, 0, identity());

  fields f1(&s1);
  f1.add_point_source(Er, 0.7, 2.5, 0.0, 4.0, vec(0.5,0.5));
  f1.add_point_source(Ep, 0.8, 0.6, 0.0, 4.0, vec(0.401,0.5));
  fields f(&s);
  f.add_point_source(Er, 0.7, 2.5, 0.0, 4.0, vec(0.5,0.5));
  f.add_point_source(Ep, 0.8, 0.6, 0.0, 4.0, vec(0.401,0.5));
  double total_energy_check_time = 1.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.01,  0.5  ))) return 0;
    if (!compare_point(f, f1, vec(0.21,  0.5  ))) return 0;
    if (!compare_point(f, f1, vec(0.501, 0.5  ))) return 0;
    if (!compare_point(f, f1, vec(0.33,  0.46 ))) return 0;
    if (!compare_point(f, f1, vec(0.2,   0.2  ))) return 0;
    if (f.time() >= total_energy_check_time) {
      if (!compare(f.electric_energy_in_box(v.surroundings()),
                   f1.electric_energy_in_box(v.surroundings()),
                   "electric energy")) return 0;
      if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                   f1.magnetic_energy_in_box(v.surroundings()),
                   "magnetic energy")) return 0;
      if (!compare(f.total_energy(), f1.total_energy(),
                   "   total energy")) return 0;
      total_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_cyl_metal_mirror_nonlinear(double eps(const vec &)) {
  master_printf("Testing Z mirror symmetry in Cylindrical...\n");
  double a = 10.0;
  double ttot = 3.0;

  const volume v = volcyl(1.0, 1.0, a);
  the_center = v.center();
  const symmetry S = mirror(Z,v);
  structure s(v, eps, 0, S);
  structure s1(v, eps, 0, identity());
  s.set_kerr(one);
  s1.set_kerr(one);

  fields f1(&s1);
  f1.add_point_source(Er, 0.7, 2.5, 0.0, 4.0, vec(0.5,0.5));
  f1.add_point_source(Ep, 0.8, 0.6, 0.0, 4.0, vec(0.401,0.5));
  fields f(&s);
  f.add_point_source(Er, 0.7, 2.5, 0.0, 4.0, vec(0.5,0.5));
  f.add_point_source(Ep, 0.8, 0.6, 0.0, 4.0, vec(0.401,0.5));
  double total_energy_check_time = 1.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.01,  0.5  ))) return 0;
    if (!compare_point(f, f1, vec(0.21,  0.5  ))) return 0;
    if (!compare_point(f, f1, vec(0.501, 0.5  ))) return 0;
    if (!compare_point(f, f1, vec(0.33,  0.46 ))) return 0;
    if (!compare_point(f, f1, vec(0.2,   0.2  ))) return 0;
    if (f.time() >= total_energy_check_time) {
      if (!compare(f.electric_energy_in_box(v.surroundings()),
                   f1.electric_energy_in_box(v.surroundings()),
                   "electric energy")) return 0;
      if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                   f1.magnetic_energy_in_box(v.surroundings()),
                   "magnetic energy")) return 0;
      if (!compare(f.total_energy(), f1.total_energy(),
                   "   total energy")) return 0;
      total_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_1d_periodic_mirror(double eps(const vec &)) {
  master_printf("Testing Z mirror symmetry in 1D...\n");
  double a = 10.0;
  double ttot = 3.0;

  const volume v = volone(1.0, a);
  the_center = v.center();
  const symmetry S = mirror(Z,v);
  structure s(v, eps, 0, S);
  structure s1(v, eps, 0, identity());

  fields f1(&s1);
  f1.use_bloch(0.0);
  f1.add_point_source(Ex, 0.7, 2.5, 0.0, 4.0, vec(0.5));
  fields f(&s);
  f.use_bloch(0.0);
  f.add_point_source(Ex, 0.7, 2.5, 0.0, 4.0, vec(0.5));
  double total_energy_check_time = 1.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.01))) return 0;
    if (!compare_point(f, f1, vec(0.33))) return 0;
    if (!compare_point(f, f1, vec(0.50))) return 0;
    if (f.time() >= total_energy_check_time) {
      if (!compare(f.electric_energy_in_box(v.surroundings()),
                   f1.electric_energy_in_box(v.surroundings()),
                   "electric energy")) return 0;
      if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                   f1.magnetic_energy_in_box(v.surroundings()),
                   "magnetic energy")) return 0;
      if (!compare(f.total_energy(), f1.total_energy(),
                   "   total energy")) return 0;
      total_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_origin_shift(const char *dirname) {
  master_printf("Testing origin shift in 2D...\n");
  double a = 10.0;
  double ttot = 3.0;

  const volume v = voltwo(1.0, 1.0, a);
  volume vcentered = v;
  vcentered.origin -= v.center();
  structure s(vcentered, one);
  structure s1(v, one);
  s.set_output_directory(dirname);
  s1.set_output_directory(dirname);

  fields f1(&s1);
  fields f(&s);
  f1.add_point_source(Ey, 0.7, 2.5, 0.0, 4.0, v.center());
  f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, v.center());
  f.add_point_source(Ey, 0.7, 2.5, 0.0, 4.0, vec2d(0.0,0.0));
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec2d(0.0,0.0));
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare(f.total_energy(), f1.total_energy(), "   total energy")) {
      master_printf("Time is %g\n", f.time());
      f1.output_real_imaginary_slices("unshifted");
      f.output_real_imaginary_slices("shifted");
      f1.eps_slices("unshifted");
      f.eps_slices("shifted");
      return 0;
    }
  }
  return 1;
}

int test_metal_xmirror(double eps(const vec &)) {
  master_printf("Testing X mirror symmetry...\n");
  double a = 10.0;
  double ttot = 3.0;

  const volume v = voltwo(1.0, 1.0, a);
  the_center = v.center();
  const symmetry S = mirror(X,v);
  structure s(v, eps, 0, S);
  structure s1(v, eps, 0, identity());

  fields f1(&s1);
  f1.add_point_source(Ey, 0.7, 2.5, 0.0, 4.0, vec2d(0.5,0.5));
  f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec2d(0.5,0.401));
  fields f(&s);
  f.add_point_source(Ey, 0.7, 2.5, 0.0, 4.0, vec2d(0.5,0.5));
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec2d(0.5,0.401));
  double total_energy_check_time = 1.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec2d(0.5  , 0.01))) return 0;
    if (!compare_point(f, f1, vec2d(0.5  , 0.21))) return 0;
    if (!compare_point(f, f1, vec2d(0.5  , 0.501))) return 0;
    if (!compare_point(f, f1, vec2d(0.46 , 0.33))) return 0;
    if (!compare_point(f, f1, vec2d(0.2  , 0.2 ))) return 0;
    if (f.time() >= total_energy_check_time) {
      if (!compare(f.electric_energy_in_box(v.surroundings()),
                   f1.electric_energy_in_box(v.surroundings()),
                   "electric energy")) return 0;
      if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                   f1.magnetic_energy_in_box(v.surroundings()),
                   "magnetic energy")) return 0;
      if (!compare(f.total_energy(), f1.total_energy(),
                   "   total energy")) return 0;
      total_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_3D_metal_xmirror(double eps(const vec &)) {
  double a = 10.0;
  double ttot = 3.0;

  const volume v = vol3d(1.0, 1.0, 1.0, a);
  const symmetry S = mirror(X,v);
  structure s(v, eps, 0, S);
  structure s1(v, eps, 0, identity());
  master_printf("Testing X mirror symmetry in 3D...\n");

  fields f1(&s1);
  f1.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.5,0.51,0.55));
  fields f(&s);
  f.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.5,0.51,0.55));
  double total_energy_check_time = 1.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.5  , 0.01 , 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.5  , 0.21 , 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.5  , 0.501, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.46 , 0.33 , 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.2  , 0.2  , 0.5))) return 0;
    if (f.time() >= total_energy_check_time) {
      if (!compare(f.electric_energy_in_box(v.surroundings()),
                   f1.electric_energy_in_box(v.surroundings()),
                   "electric energy")) return 0;
      if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                   f1.magnetic_energy_in_box(v.surroundings()),
                   "magnetic energy")) return 0;
      if (!compare(f.total_energy(), f1.total_energy(),
                   "   total energy")) return 0;
      total_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_3D_metal_zmirror(double eps(const vec &)) {
  double a = 10.0;
  double ttot = 3.0;

  const volume v = vol3d(1.1, 0.6, 1.0, a);
  const symmetry S = mirror(Z,v);
  structure s(v, eps, 0, S);
  structure s1(v, eps, 0, identity());
  master_printf("Testing Z mirror symmetry in 3D...\n");

  fields f1(&s1);
  f1.add_point_source(Ex, 0.7, 2.5, 0.0, 4.0, vec(0.55,0.51,0.5));
  fields f(&s);
  f.add_point_source(Ex, 0.7, 2.5, 0.0, 4.0, vec(0.55,0.51,0.5));
  double total_energy_check_time = 1.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.5  , 0.01 , 0.75))) return 0;
    if (!compare_point(f, f1, vec(0.5  , 0.21 , 0.15))) return 0;
    if (!compare_point(f, f1, vec(0.5  , 0.501, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.46 , 0.33 , 0.51))) return 0;
    if (!compare_point(f, f1, vec(0.2  , 0.2  , 0.05))) return 0;
    if (f.time() >= total_energy_check_time) {
      if (!compare(f.electric_energy_in_box(v.surroundings()),
                   f1.electric_energy_in_box(v.surroundings()),
                   "electric energy")) return 0;
      if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                   f1.magnetic_energy_in_box(v.surroundings()),
                   "magnetic energy")) return 0;
      if (!compare(f.total_energy(), f1.total_energy(),
                   "   total energy")) return 0;
      total_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_3D_metal_odd_zmirror(double eps(const vec &)) {
  double a = 10.0;
  double ttot = 3.0;

  const volume v = vol3d(1.1, 0.6, 1.0, a);
  const symmetry S = mirror(Z,v)*(-1.0);
  structure s(v, eps, 0, S);
  structure s1(v, eps, 0, identity());
  master_printf("Testing odd Z mirror symmetry in 3D...\n");

  fields f1(&s1);
  f1.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.55,0.51,0.5));
  fields f(&s);
  f.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.55,0.51,0.5));
  double total_energy_check_time = 1.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.5  , 0.01 , 0.75))) return 0;
    if (!compare_point(f, f1, vec(0.5  , 0.21 , 0.15))) return 0;
    if (!compare_point(f, f1, vec(0.5  , 0.501, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.46 , 0.33 , 0.51))) return 0;
    if (!compare_point(f, f1, vec(0.2  , 0.2  , 0.05))) return 0;
    if (f.time() >= total_energy_check_time) {
      if (!compare(f.electric_energy_in_box(v.surroundings()),
                   f1.electric_energy_in_box(v.surroundings()),
                   "electric energy")) return 0;
      if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                   f1.magnetic_energy_in_box(v.surroundings()),
                   "magnetic energy")) return 0;
      if (!compare(f.total_energy(), f1.total_energy(),
                   "   total energy")) return 0;
      total_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_3D_metal_rot4z(double eps(const vec &)) {
  double a = 10.0;
  double ttot = 3.0;

  const volume v = vol3d(1.0, 1.0, 1.0, a);
  const symmetry S = rotate4(Z,v);
  structure s(v, eps, 0, S);
  structure s1(v, eps, 0, identity());
  master_printf("Testing Z fourfold rotational symmetry in 3D...\n");

  fields f1(&s1);
  f1.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.5,0.5,0.52));
  fields f(&s);
  f.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.5,0.5,0.52));
  double total_energy_check_time = 1.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.5  , 0.01 , 0.75))) return 0;
    if (!compare_point(f, f1, vec(0.5  , 0.21 , 0.15))) return 0;
    if (!compare_point(f, f1, vec(0.5  , 0.501, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.46 , 0.33 , 0.51))) return 0;
    if (!compare_point(f, f1, vec(0.2  , 0.2  , 0.05))) return 0;
    if (f.time() >= total_energy_check_time) {
      if (!compare(f.electric_energy_in_box(v.surroundings()),
                   f1.electric_energy_in_box(v.surroundings()),
                   "electric energy")) return 0;
      if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                   f1.magnetic_energy_in_box(v.surroundings()),
                   "magnetic energy")) return 0;
      if (!compare(f.total_energy(), f1.total_energy(),
                   "   total energy")) return 0;
      total_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_3D_metal_rot4z_mirror(double eps(const vec &)) {
  double a = 10.0;
  double ttot = 3.0;

  const volume v = vol3d(1.0, 1.0, 1.0, a);
  const symmetry S = rotate4(Z,v) + mirror(Z,v);
  structure s(v, eps, 0, S);
  structure s1(v, eps, 0, identity());
  master_printf("Testing Z fourfold rotational symmetry in 3D with horizontal mirror...\n");

  fields f1(&s1);
  f1.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec(0.5,0.5,0.5));
  fields f(&s);
  f.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec(0.5,0.5,0.5));
  double total_energy_check_time = 1.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.5  , 0.01 , 0.75))) return 0;
    if (!compare_point(f, f1, vec(0.5  , 0.21 , 0.15))) return 0;
    if (!compare_point(f, f1, vec(0.5  , 0.501, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.46 , 0.33 , 0.51))) return 0;
    if (!compare_point(f, f1, vec(0.2  , 0.2  , 0.05))) return 0;
    if (f.time() >= total_energy_check_time) {
      if (!compare(f.electric_energy_in_box(v.surroundings()),
                   f1.electric_energy_in_box(v.surroundings()),
                   "electric energy")) return 0;
      if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                   f1.magnetic_energy_in_box(v.surroundings()),
                   "magnetic energy")) return 0;
      if (!compare(f.total_energy(), f1.total_energy(),
                   "   total energy")) return 0;
      total_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_metal_ymirror(double eps(const vec &)) {
  double a = 10.0;
  double ttot = 5.0;

  const volume v = voltwo(1.0, 1.0, a);
  the_center = v.center();
  const symmetry S = mirror(Y,v);
  structure s(v, eps, 0, S);
  structure s1(v, eps, 0, identity());
  master_printf("Testing Y mirror symmetry...\n");

  fields f1(&s1);
  f1.add_point_source(Ex, 0.7, 2.5, 0.0, 4.0, vec2d(0.85 ,0.5));
  f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec2d(0.401,0.5));
  fields f(&s);
  f.add_point_source(Ex, 0.7, 2.5, 0.0, 4.0, vec2d(0.85 ,0.5));
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec2d(0.401,0.5));
  double total_energy_check_time = 1.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec2d(0.01 ,  0.5))) return 0;
    if (!compare_point(f, f1, vec2d(0.21 ,  0.5))) return 0;
    if (!compare_point(f, f1, vec2d(0.46 , 0.33))) return 0;
    if (!compare_point(f, f1, vec2d(0.2  , 0.2 ))) return 0;
    if (f.time() >= total_energy_check_time) {
      if (!compare(f.electric_energy_in_box(v.surroundings()),
                   f1.electric_energy_in_box(v.surroundings()),
                   "electric energy")) return 0;
      if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                   f1.magnetic_energy_in_box(v.surroundings()),
                   "magnetic energy")) return 0;
      if (!compare(f.total_energy(), f1.total_energy(),
                   "   total energy")) return 0;
      total_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_yperiodic_ymirror(double eps(const vec &)) {
  double a = 10.0;
  double ttot = 5.0;

  const volume v = voltwo(1.0, 1.0, a);
  the_center = v.center();
  const symmetry S = mirror(Y,v);
  structure s(v, eps, 0, S);
  structure s1(v, eps, 0, identity());
  s.set_output_directory(dirname);
  s1.set_output_directory(dirname);
  master_printf("Testing Y periodic with mirror symmetry...\n");

  fields f1(&s1);
  f1.use_bloch(vec2d(0.1*pi/2,0.0));
  //f1.add_point_source(Ex, 0.7, 2.5, 0.0, 4.0, vec2d(0.85 ,0.5));
  f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec2d(0.401,0.5));
  fields f(&s);
  f.use_bloch(vec2d(0.1*pi/2,0.0));
  //f.add_point_source(Ex, 0.7, 2.5, 0.0, 4.0, vec2d(0.85 ,0.5));
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec2d(0.401,0.5));
  double total_energy_check_time = 1.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec2d(0.951 ,  0.5))) return 0;
    if (!compare_point(f, f1, vec2d(0.01 ,  0.5))) return 0;
    if (!compare_point(f, f1, vec2d(0.21 ,  0.5))) return 0;
    if (!compare_point(f, f1, vec2d(0.46 , 0.33))) return 0;
    if (!compare_point(f, f1, vec2d(0.2  , 0.2 ))) return 0;
    if (f.time() >= total_energy_check_time) {
      if (!compare(f.electric_energy_in_box(v.surroundings()),
                   f1.electric_energy_in_box(v.surroundings()),
                   "electric energy")) {
        f.output_real_imaginary_slices("multi");
        f1.output_real_imaginary_slices("single");
        return 0;
      }
      if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                   f1.magnetic_energy_in_box(v.surroundings()),
                   "magnetic energy")) return 0;
      if (!compare(f.total_energy(), f1.total_energy(),
                   "   total energy")) return 0;
      total_energy_check_time += 1.0;
    }
  }
  return 1;
}

int test_metal_rot2y(double eps(const vec &)) {
  double a = 10.0;
  double ttot = 5.0;

  const volume v = voltwo(1.0, 1.0, a);
  the_center = v.center();
  const symmetry S = rotate2(Y,v);
  structure s(v, eps, 0, S);
  structure s1(v, eps, 0, identity());
  master_printf("Testing Y twofold rotational symmetry...\n");

  fields f1(&s1);
  f1.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec2d(0.25, 0.85), 1.0);
  f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec2d(0.25,0.4), 1.0);
  f1.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec2d(0.75, 0.85),-1.0);
  f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec2d(0.75,0.4),-1.0);
  fields f(&s);
  f.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec2d(0.25,0.85 ), 1.0);
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec2d(0.25,0.4), 1.0);
  f.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec2d(0.75,0.85 ),-1.0);
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec2d(0.75,0.4),-1.0);
  double total_energy_check_time = 1.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec2d(0.01 ,  0.5))) return 0;
    if (!compare_point(f, f1, vec2d(0.21 ,  0.5))) return 0;
    if (!compare_point(f, f1, vec2d(0.46 , 0.33))) return 0;
    if (!compare_point(f, f1, vec2d(0.2  , 0.2 ))) return 0;
    if (f.time() >= total_energy_check_time) {
      if (!compare(f.electric_energy_in_box(v.surroundings()),
                   f1.electric_energy_in_box(v.surroundings()),
                   "electric energy")) return 0;
      if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                   f1.magnetic_energy_in_box(v.surroundings()),
                   "magnetic energy")) return 0;
      if (!compare(f.total_energy(), f1.total_energy(),
                   "   total energy")) return 0;
      total_energy_check_time += 1.0;
    }
  }
  return 1;
}

int exact_metal_rot2y(double eps(const vec &)) {
  double a = 10.0;
  double ttot = 5.0;

  const volume v = voltwo(1.0, 1.5, a);
  the_center = v.center();
  const symmetry S = rotate2(Y,v);
  structure s(v, eps, 0, S);
  structure s1(v, eps, 0, identity());
  master_printf("Testing exact Y twofold rotational symmetry...\n");

  fields f1(&s1);
  f1.add_point_source(Ey, 0.7, 2.5, 0.0, 4.0, vec2d(0.5, 0.85));
  f1.add_point_source(Hy, 0.8, 0.6, 0.0, 4.0, vec2d(0.5,0.4));
  fields f(&s);
  f.add_point_source(Ey, 0.7, 2.5, 0.0, 4.0, vec2d(0.5, 0.85));
  f.add_point_source(Hy, 0.8, 0.6, 0.0, 4.0, vec2d(0.5,0.4));
  double total_energy_check_time = 1.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec2d(0.01 ,  0.5))) return 0;
    if (!compare_point(f, f1, vec2d(0.21 ,  0.5))) return 0;
    if (!compare_point(f, f1, vec2d(0.46 , 0.33))) return 0;
    if (!compare_point(f, f1, vec2d(0.2  , 0.2 ))) return 0;
    if (f.time() >= total_energy_check_time) {
      if (!compare(f.electric_energy_in_box(v.surroundings()),
                   f1.electric_energy_in_box(v.surroundings()),
                   "electric energy")) return 0;
      if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                   f1.magnetic_energy_in_box(v.surroundings()),
                   "magnetic energy")) return 0;
      if (!compare(f.total_energy(), f1.total_energy(),
                   "   total energy")) return 0;
      total_energy_check_time += 1.0;
    }
  }
  return 1;
}

int pml_twomirrors(double eps(const vec &)) {
  double a = 10.0;
  double ttot = 10.0;

  const volume v = voltwo(2.0, 2.0, a);
  the_center = v.center();
  const symmetry S = mirror(X,v) + mirror(Y,v);

  structure s_mm(v, eps, 0, S);
  structure s1(v, eps, 0, identity());
  structure ss[2] = { s1, s_mm };

  for (int i=0;i<2;i++)
    ss[i].use_pml_everywhere(0.5);
  master_printf("Testing two mirrors with PML...\n");

  fields fs[2] = { fields(&ss[0]), fields(&ss[1]) };
  for (int i=0;i<2;i++) {
    fs[i].add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec2d(1.0,1.0),-1.5);
    fs[i].add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec2d(0.75,0.75));
    fs[i].add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec2d(0.75,1.25));
    fs[i].add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec2d(1.25,0.75));
    fs[i].add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec2d(1.25,1.25));
  }
  double total_energy_check_time = 3.0;
  while (fs[0].time() < ttot) {
    for (int i=0;i<2;i++) fs[i].step();
    if (!compare_point(fs[1], fs[0], vec2d(0.01 ,  0.5))) return 0;
    if (!compare_point(fs[1], fs[0], vec2d(0.21 ,  0.5))) return 0;
    if (!compare_point(fs[1], fs[0], vec2d(0.46 , 0.33))) return 0;
    if (!compare_point(fs[1], fs[0], vec2d(0.2  , 0.2 ))) return 0;
    if (fs[0].time() >= total_energy_check_time) {
      if (!compare(fs[0].electric_energy_in_box(v.surroundings()),
                   fs[1].electric_energy_in_box(v.surroundings()),
                   "electric energy")) return 0;
      total_energy_check_time += 3.0;
    }
  }
  fs[0].eps_slices("mirror_single");
  fs[1].eps_slices("mirror_multi");
  return 1;
}

int exact_metal_rot4z(double eps(const vec &)) {
  double a = 10.0;
  double ttot = 5.0;

  const volume v = voltwo(1.0, 1.0, a);
  the_center = v.center();
  const symmetry S = rotate4(Z,v);

  structure s(v, eps, 0, S);
  structure s1(v, eps, 0, identity());
  master_printf("Testing Z fourfold rotational symmetry...\n");

  fields f1(&s1);
  f1.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec2d(0.5,0.5));
  f1.add_point_source(Hz, 0.8, 0.6, 0.0, 4.0, vec2d(0.5,0.5));
  fields f(&s);
  f.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec2d(0.5,0.5));
  f.add_point_source(Hz, 0.8, 0.6, 0.0, 4.0, vec2d(0.5,0.5));
  double total_energy_check_time = 1.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec2d(0.01 ,  0.5))) return 0;
    if (!compare_point(f, f1, vec2d(0.21 ,  0.5))) return 0;
    if (!compare_point(f, f1, vec2d(0.46 , 0.33))) return 0;
    if (!compare_point(f, f1, vec2d(0.2  , 0.2 ))) return 0;
    if (f.time() >= total_energy_check_time) {
      if (!compare(f.electric_energy_in_box(v.surroundings()),
                   f1.electric_energy_in_box(v.surroundings()),
                   "electric energy")) return 0;
      if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                   f1.magnetic_energy_in_box(v.surroundings()),
                   "magnetic energy")) return 0;
      if (!compare(f.total_energy(), f1.total_energy(),
                   "   total energy")) return 0;
      total_energy_check_time += 1.0;
    }
  }
  return 1;
}

int exact_metal_rot4z_nonlinear(double eps(const vec &)) {
  double a = 10.0;
  double ttot = 5.0;

  const volume v = voltwo(1.0, 1.0, a);
  the_center = v.center();
  const symmetry S = rotate4(Z,v);

  structure s(v, eps, 0, S);
  structure s1(v, eps, 0, identity());
  s.set_kerr(one);
  s1.set_kerr(one);
  master_printf("Testing nonlinear Z fourfold rotational symmetry...\n");

  fields f1(&s1);
  f1.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec2d(0.5,0.5));
  f1.add_point_source(Hz, 0.8, 0.6, 0.0, 4.0, vec2d(0.5,0.5));
  fields f(&s);
  f.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec2d(0.5,0.5));
  f.add_point_source(Hz, 0.8, 0.6, 0.0, 4.0, vec2d(0.5,0.5));
  double total_energy_check_time = 1.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec2d(0.01 ,  0.5))) return 0;
    if (!compare_point(f, f1, vec2d(0.21 ,  0.5))) return 0;
    if (!compare_point(f, f1, vec2d(0.46 , 0.33))) return 0;
    if (!compare_point(f, f1, vec2d(0.2  , 0.2 ))) return 0;
    if (f.time() >= total_energy_check_time) {
      if (!compare(f.electric_energy_in_box(v.surroundings()),
                   f1.electric_energy_in_box(v.surroundings()),
                   "electric energy")) return 0;
      if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                   f1.magnetic_energy_in_box(v.surroundings()),
                   "magnetic energy")) return 0;
      if (!compare(f.total_energy(), f1.total_energy(),
                   "   total energy")) return 0;
      total_energy_check_time += 1.0;
    }
  }
  return 1;
}

int exact_pml_rot2x_tm(double eps(const vec &)) {
  double a = 10.0;
  double ttot = 30.0;

  const volume v = voltwo(3.0, 3.0, a);
  the_center = v.center();
  const symmetry S = rotate2(X,v);

  structure s(v, eps, 0, S);
  structure s1(v, eps, 0, identity());
  s.set_output_directory(dirname);
  s1.set_output_directory(dirname);
  s.use_pml_everywhere(1.0);
  s1.use_pml_everywhere(1.0);
  master_printf("Testing X twofold rotational symmetry with PML...\n");

  fields f1(&s1);
  f1.add_point_source(Hx, 0.7, 2.5, 0.0, 4.0, vec2d(1.3,1.5));
  fields f(&s);
  f.add_point_source(Hx, 0.7, 2.5, 0.0, 4.0, vec2d(1.3,1.5));
  double total_energy_check_time = 1.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec2d(0.01 ,  1.5))) return 0;
    if (!compare_point(f, f1, vec2d(1.21 ,  1.5))) return 0;
    if (!compare_point(f, f1, vec2d(1.46 , 0.33))) return 0;
    if (!compare_point(f, f1, vec2d(1.2  , 1.2 ))) return 0;
    if (f.time() >= total_energy_check_time) {
      if (!compare(f.electric_energy_in_box(v.surroundings()),
                   f1.electric_energy_in_box(v.surroundings()),
                   "electric energy")) return 0;
      if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                   f1.magnetic_energy_in_box(v.surroundings()),
                   "magnetic energy")) return 0;
      if (!compare(f.total_energy(), f1.total_energy(),
                   "   total energy")) return 0;
      total_energy_check_time += 1.0;
    }
  }
  // Just to make sure slice output doesn't actually crash...
  f.output_real_imaginary_slices("hello");
  f1.output_real_imaginary_slices("world");
  f.eps_slices("hello");
  f1.eps_slices("world");
  return 1;
}

double polariton_ex(const volume &v, double eps(const vec &)) {
  const double ttot = 10.0;
  master_printf("Testing polariton in %s...\n", dimension_name(v.dim));
  the_center = v.center();
  const symmetry S = mirror(Z,v);
  structure s(v, eps);
  structure sS(v, eps, 0, S);
  s.add_polarizability(one, 0.3, 0.1, 7.63);
  sS.add_polarizability(one, 0.3, 0.1, 7.63);
  fields f(&s);
  f.add_point_source(Ex, 0.2, 3.0, 0.0, 2.0, v.center());
  fields fS(&s);
  fS.add_point_source(Ex, 0.2, 3.0, 0.0, 2.0, v.center());
  f.use_real_fields();
  fS.use_real_fields();
  f.use_bloch(zero_vec(v.dim));
  fS.use_bloch(zero_vec(v.dim));
  while (f.time() < ttot) {
    f.step();
    fS.step();
    if (!compare_point(fS, f, v.center())) return 0;
    if (!compare_point(fS, f, zero_vec(v.dim))) return 0;
    if (!compare_point(fS, f, v.center()*0.3)) return 0;
  }
  return 1;
}

double saturated_gain_ez(const volume &v, double eps(const vec &)) {
  const double ttot = 10.0;
  master_printf("Testing saturated gain in %s...\n", dimension_name(v.dim));
  the_center = v.center();
  const symmetry S = mirror(Z,v)*(-1);
  structure s(v, eps);
  structure sS(v, eps, 0, S);
  s.add_polarizability(one, 0.3, -0.1, 7.63, 0.5);
  sS.add_polarizability(one, 0.3, -0.1, 7.63, 0.5);
  fields f(&s);
  f.add_point_source(Ez, 0.2, 3.0, 0.0, 2.0, v.center());
  fields fS(&s);
  fS.add_point_source(Ez, 0.2, 3.0, 0.0, 2.0, v.center());
  f.use_real_fields();
  fS.use_real_fields();
  f.use_bloch(zero_vec(v.dim));
  fS.use_bloch(zero_vec(v.dim));
  while (f.time() < ttot) {
    f.step();
    fS.step();
    if (!compare_point(fS, f, v.center())) return 0;
    if (!compare_point(fS, f, zero_vec(v.dim))) return 0;
    if (!compare_point(fS, f, v.center()*0.3)) return 0;
  }
  return 1;
}

double saturated_gain_te(const volume &v, double eps(const vec &)) {
  const double ttot = 10.0;
  master_printf("Testing saturated gain in %s...\n", dimension_name(v.dim));
  the_center = v.center();
  const symmetry S = mirror(X,v)*(-1);
  structure s(v, eps);
  structure sS(v, eps, 0, S);
  s.add_polarizability(one, 0.3, -0.1, 7.63, 0.5);
  sS.add_polarizability(one, 0.3, -0.1, 7.63, 0.5);
  fields f(&s);
  f.add_point_source(Ex, 0.2, 3.0, 0.0, 2.0, v.center());
  fields fS(&s);
  fS.add_point_source(Ex, 0.2, 3.0, 0.0, 2.0, v.center());
  f.use_real_fields();
  fS.use_real_fields();
  f.use_bloch(zero_vec(v.dim));
  fS.use_bloch(zero_vec(v.dim));
  while (f.time() < ttot) {
    f.step();
    fS.step();
    if (!compare_point(fS, f, v.center())) return 0;
    if (!compare_point(fS, f, zero_vec(v.dim))) return 0;
    if (!compare_point(fS, f, v.center()*0.3)) return 0;
  }
  return 1;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  trash_output_directory(dirname);
  master_printf("Testing with various kinds of symmetry...\n");

  if (!polariton_ex(vol1d(1.0, 10.0), one))
    abort("error in 1D polariton vacuum\n");

  if (!polariton_ex(vol3d(1.0, 1.2, 0.8, 10.0), one))
    abort("error in 3D polariton vacuum\n");

  if (0) // FIXME: disable until divergence bug is fixed
  if (!saturated_gain_ez(vol3d(0.5, 1.2, 0.8, 10.0), one))
    abort("error in 3D saturated gain\n");

  if (0) // FIXME: disable until divergence bug is fixed
  if (!saturated_gain_ez(volcyl(0.5, 1.2, 10.0), one))
    abort("error in cylindrical saturated gain\n");

  if (0) // FIXME: disable until divergence bug is fixed
  if (!saturated_gain_te(vol2d(0.6, 1.2, 10.0), one))
    abort("error in 2D TE saturated gain\n");

  if (!test_1d_periodic_mirror(one))
    abort("error in test_1d_periodic_mirror vacuum\n");

  if (!test_cyl_metal_mirror(one))
    abort("error in test_cyl_metal_mirror vacuum\n");

  if (!test_cyl_metal_mirror(one))
    abort("error in test_cyl_metal_mirror nonlinear vacuum\n");

  if (!test_yperiodic_ymirror(one))
    abort("error in test_yperiodic_ymirror vacuum\n");

  if (!pml_twomirrors(one))
    abort("error in pml_twomirrors vacuum\n");

  if (!test_origin_shift(dirname))
    abort("error in test_origin_shift\n");

  if (!exact_pml_rot2x_tm(one))
    abort("error in exact_pml_rot2x_tm vacuum\n");

  if (!test_metal_xmirror(one))
    abort("error in test_metal_xmirror vacuum\n");
  if (!test_metal_xmirror(rods_2d))
    abort("error in test_metal_xmirror rods_2d\n");

  if (!test_metal_ymirror(one))
    abort("error in test_metal_ymirror vacuum\n");
  if (!test_metal_ymirror(rods_2d))
    abort("error in test_metal_ymirror rods_2d\n");

  if (!test_metal_rot2y(one))
    abort("error in test_metal_rot2y vacuum\n");
  if (!test_metal_rot2y(rods_2d))
    abort("error in test_metal_rot2y rods_2d\n");

  if (!exact_metal_rot2y(one))
    abort("error in exact_metal_rot2y vacuum\n");
  if (!exact_metal_rot2y(rods_2d))
    abort("error in exact_metal_rot2y rods_2d\n");

  if (!exact_metal_rot4z(one))
    abort("error in exact_metal_rot4z vacuum\n");
  if (!exact_metal_rot4z(rods_2d))
    abort("error in exact_metal_rot4z rods_2d\n");

  if (!exact_metal_rot4z_nonlinear(one))
    abort("error in exact_metal_rot4z nonlinear vacuum\n");
  if (!exact_metal_rot4z_nonlinear(rods_2d))
    abort("error in exact_metal_rot4z nonlinear rods_2d\n");

  if (!test_3D_metal_xmirror(one))
    abort("error in test_3D_metal_xmirror vacuum\n");

  if (!test_3D_metal_zmirror(one))
    abort("error in test_3D_metal_zmirror vacuum\n");

  if (!test_3D_metal_odd_zmirror(one))
    abort("error in test_3D_metal_odd_zmirror vacuum\n");

  if (!test_3D_metal_rot4z(one)) {
    all_wait();
    abort("error in test_3D_metal_rot4z vacuum\n");
  }

  if (!test_3D_metal_rot4z_mirror(one))
    abort("error in test_3D_metal_rot4z_mirror vacuum\n");

  return 0;
}
