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

#include "meep.h"

double eps(const vec &) { return 2.0; }

int compare(double a, double b, const char *n) {
  if (fabs(a-b) > fabs(b)*1e-15) {
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
      if (abs(v1 - v2) > 2e-13*abs(v2)) {
        master_printf("%s differs:  %lg %lg out of %lg %lg\n",
               component_name(c), real(v2-v1), imag(v2-v1), real(v2), imag(v2));
        master_printf("This comes out to a fractional error of %lg\n",
               abs(v1 - v2)/abs(v2));
        switch (p.dim) {
        case D2: master_printf("Right now I'm looking at %lg %lg, time %lg\n",
                               p.x(), p.y(), f1.time());
          break;
        case D1: master_printf("Right now I'm looking at %lg, time %lg\n",
                               p.z(), f1.time());
          break;
        }
        f1.output_real_imaginary_slices("new");
        f2.output_real_imaginary_slices("old");
        f1.eps_slices("new");
        f2.eps_slices("old");
        return 0;
      }
    }
  }
  return 1;
}

bool step_metal_1d(const char *dirname) {
  double a = 10.0;
  double ttot = 2.0;
  const volume v = volone(1.0, a);
  mat ma(v, eps);
  mat ma_old(v, eps);
  ma.set_output_directory(dirname);
  ma_old.set_output_directory(dirname);
  master_printf("Testing step algorithm in 1D...\n");

  fields f(&ma);
  f.use_metal_everywhere();
  f.add_point_source(Ex, 0.7, 2.5, 0.0, 4.0, vec(.3));
  fields f_old(&ma);
  f_old.use_metal_everywhere();
  f_old.add_point_source(Ex, 0.7, 2.5, 0.0, 4.0, vec(.3));
  while (f.time() < ttot) {
    f.step();
    f_old.step_old();
    if (!compare_point(f, f_old, vec(0.01))) return 0;
    if (!compare_point(f, f_old, vec(.301))) return 0;
    if (!compare_point(f, f_old, vec(.46 ))) return 0;
    if (!compare(f.electric_energy_in_box(v.surroundings()),
                 f_old.electric_energy_in_box(v.surroundings()),
                 "electric energy")) return 0;
  }
  return 1;
}

bool step_metal_1d_pml(const char *dirname) {
  double a = 10.0;
  double ttot = 10.0;
  const volume v = volone(2.0, a);
  mat ma(v, eps);
  mat ma_old(v, eps);
  ma.set_output_directory(dirname);
  ma_old.set_output_directory(dirname);
  ma.use_pml_everywhere(0.4);
  ma_old.use_pml_everywhere(0.4);
  master_printf("Testing step algorithm in 1D with PML...\n");

  fields f(&ma);
  f.add_point_source(Ex, 0.7, 2.5, 0.0, 4.0, vec(.6));
  fields f_old(&ma);
  f_old.add_point_source(Ex, 0.7, 2.5, 0.0, 4.0, vec(.6));
  while (f.time() < ttot) {
    f.step();
    f_old.step_old();
    if (!compare_point(f, f_old, vec(0.01))) return 0;
    if (!compare_point(f, f_old, vec(.301))) return 0;
    if (!compare_point(f, f_old, vec(.46 ))) return 0;
    if (!compare(f.electric_energy_in_box(v.surroundings()),
                 f_old.electric_energy_in_box(v.surroundings()),
                 "electric energy")) return 0;
  }
  return 1;
}

bool step_metal_2d_tm(const char *dirname) {
  double a = 10.0;
  double ttot = 2.0;

  const volume v = voltwo(1.0, 1.0, a);

  mat ma(v, eps);
  mat ma_old(v, eps);
  ma.set_output_directory(dirname);
  ma_old.set_output_directory(dirname);
  master_printf("Testing step algorithm in 2D TM...\n");

  fields f(&ma);
  f.use_metal_everywhere();
  f.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec2d(.3,.5));
  fields f_old(&ma);
  f_old.use_metal_everywhere();
  f_old.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec2d(.3,.5));
  while (f.time() < ttot) {
    f.step();
    f_old.step_old();
    if (!compare_point(f, f_old, vec2d(0.01 ,  0.5))) return 0;
    if (!compare_point(f, f_old, vec2d(.301 ,  .5))) return 0;
    if (!compare_point(f, f_old, vec2d(.46 , 0.33))) return 0;
    if (!compare_point(f, f_old, vec2d(.2  , .2 ))) return 0;
    if (!compare(f.electric_energy_in_box(v.surroundings()),
                 f_old.electric_energy_in_box(v.surroundings()),
                 "electric energy")) return 0;
    if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                 f_old.magnetic_energy_in_box(v.surroundings()),
                 "magnetic energy")) return 0;
    if (!compare(f.total_energy(), f_old.total_energy(),
                 "   total energy")) return 0;
  }
  return 1;
}

bool step_metal_2d_te(const char *dirname) {
  double a = 10.0;
  double ttot = 2.0;

  const volume v = voltwo(1.0, 1.0, a);

  mat ma(v, eps);
  mat ma_old(v, eps);
  ma.set_output_directory(dirname);
  ma_old.set_output_directory(dirname);
  master_printf("Testing step algorithm in 2D TE...\n");

  fields f(&ma);
  f.use_metal_everywhere();
  f.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec2d(.3,.5));
  fields f_old(&ma);
  f_old.use_metal_everywhere();
  f_old.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec2d(.3,.5));
  while (f.time() < ttot) {
    f.step();
    f_old.step_old();
    if (!compare_point(f, f_old, vec2d(0.01 ,  0.5))) return 0;
    if (!compare_point(f, f_old, vec2d(.301 ,  .5))) return 0;
    if (!compare_point(f, f_old, vec2d(.46 , 0.33))) return 0;
    if (!compare_point(f, f_old, vec2d(.2  , .2 ))) return 0;
    if (!compare(f.electric_energy_in_box(v.surroundings()),
                 f_old.electric_energy_in_box(v.surroundings()),
                 "electric energy")) return 0;
    if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                 f_old.magnetic_energy_in_box(v.surroundings()),
                 "magnetic energy")) return 0;
    if (!compare(f.total_energy(), f_old.total_energy(),
                 "   total energy")) return 0;
  }
  return 1;
}

bool step_pml_2d_tm(const char *dirname) {
  double a = 10.0;
  double ttot = 2.0;

  const volume v = voltwo(1.0, 1.0, a);

  mat ma(v, eps);
  mat ma_old(v, eps);
  ma.use_pml_everywhere(0.3);
  ma_old.use_pml_everywhere(0.3);
  ma.set_output_directory(dirname);
  ma_old.set_output_directory(dirname);
  master_printf("Testing step algorithm in 2D with TM PML...\n");

  fields f(&ma);
  f.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec2d(.4,.503));
  fields f_old(&ma);
  f_old.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec2d(.4,.503));
  while (f.time() < ttot) {
    f.step();
    f_old.step_old();
    if (!compare_point(f, f_old, vec2d(0.01 ,  0.5))) return 0;
    if (!compare_point(f, f_old, vec2d(.401 ,  .5))) return 0;
    if (!compare_point(f, f_old, vec2d(.46 , 0.33))) return 0;
    if (!compare_point(f, f_old, vec2d(.2  , .5 ))) return 0;
    if (!compare_point(f, f_old, vec2d(.2  , .2 ))) return 0;
    if (!compare(f.electric_energy_in_box(v.surroundings()),
                 f_old.electric_energy_in_box(v.surroundings()),
                 "electric energy")) return 0;
    if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                 f_old.magnetic_energy_in_box(v.surroundings()),
                 "magnetic energy")) return 0;
  }
  return true;
}

bool step_pml_2d_te(const char *dirname) {
  double a = 10.0;
  double ttot = 2.0;

  const volume v = voltwo(1.0, 1.0, a);

  mat ma(v, eps);
  mat ma_old(v, eps);
  ma.use_pml_everywhere(0.3);
  ma_old.use_pml_everywhere(0.3);
  ma.set_output_directory(dirname);
  ma_old.set_output_directory(dirname);
  master_printf("Testing step algorithm in 2D with TE PML...\n");

  fields f(&ma);
  f.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec2d(.4,.503));
  fields f_old(&ma);
  f_old.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec2d(.4,.503));
  while (f.time() < ttot) {
    f.step();
    f_old.step_old();
    if (!compare_point(f, f_old, vec2d(0.01 ,  0.5))) return 0;
    if (!compare_point(f, f_old, vec2d(.401 ,  .5))) return 0;
    if (!compare_point(f, f_old, vec2d(.29 , 0.53))) return 0;
    if (!compare_point(f, f_old, vec2d(.2  , .2 ))) return 0;
    if (!compare(f.electric_energy_in_box(v.surroundings()),
                 f_old.electric_energy_in_box(v.surroundings()),
                 "electric energy")) return 0;
    if (!compare(f.magnetic_energy_in_box(v.surroundings()),
                 f_old.magnetic_energy_in_box(v.surroundings()),
                 "magnetic energy")) return 0;
  }
  return 1;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  const char *dirname = "step_algorithm-out";
  trash_output_directory(dirname);

  if (!step_metal_1d(dirname))
    abort("error in step_metal_1d\n");

  if (!step_metal_1d_pml(dirname))
    abort("error in step_metal_1d_pml\n");

  if (!step_metal_2d_tm(dirname))
    abort("error in step_metal_2d_tm\n");

  if (!step_metal_2d_te(dirname))
    abort("error in step_metal_2d_te\n");

  if (!step_pml_2d_tm(dirname))
    abort("error in step_pml_2d_tm\n");

  if (!step_pml_2d_te(dirname))
    abort("error in step_pml_2d_te\n");

  return 0;
}

