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
double targets(const vec &v) {
  const double r = sqrt(v.x()*v.x() + v.y()*v.y());
  double dr = r;
  while (dr > 1) dr -= 1;
  if (dr > 0.7001) return 12.0;
  return 1.0;
}

int compare(double a, double b, const char *n) {
  if (fabs(a-b) > fabs(b)*2.0e-12) {
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
        master_printf("Right now I'm looking at %lg %lg %lg, time %lg\n",
                      p.x(), p.y(), p.z(), f1.time());
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

int approx_point(fields &f1, fields &f2, const vec &p) {
  monitor_point m1, m_test;
  f1.get_point(&m_test, p);
  f2.get_point(&m1, p);
  for (int i=0;i<10;i++) {
    component c = (component) i;
    if (f1.v.has_field(c)) {
      complex<double> v1 = m_test.get_component(c), v2 = m1.get_component(c);
      if (abs(v1 - v2) > 4e-15*(5.0 + abs(v2))) {
        master_printf("%s differs:  %lg %lg out of %lg %lg\n",
               component_name(c), real(v2-v1), imag(v2-v1), real(v2), imag(v2));
        master_printf("This comes out to a fractional error of %lg\n",
               abs(v1 - v2)/abs(v2));
        master_printf("Right now I'm looking at %lg %lg %lg, time %lg\n",
                      p.x(), p.y(), p.z(), f1.time());
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

int test_metal(double eps(const vec &), int splitting, const char *dirname) {
  double a = 10.0;
  double ttot = 17.0;

  volume v = vol3d(3.0, 2.0, 1.0, a);
  mat ma1(v, eps, 1);
  mat ma(v, eps, splitting);
  ma.set_output_directory(dirname);
  ma1.set_output_directory(dirname);

  master_printf("Metal test using %d chunks...\n", splitting);
  fields f(&ma);
  f.use_metal_everywhere();
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(1.299,1.299,0.401), 1.0);
  fields f1(&ma1);
  f1.use_metal_everywhere();
  f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(1.299,1.299,0.401), 1.0);
  double total_energy_check_time = 8.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.5  , 0.5  , 0.01))) return 0;
    if (!compare_point(f, f1, vec(0.46 , 0.33 , 0.33))) return 0;
    if (!compare_point(f, f1, vec(1.301  , 1.301  , 0.399 ))) return 0;
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

int test_periodic(double eps(const vec &), int splitting, const char *dirname) {
  double a = 10.0;
  double ttot = 17.0;

  volume v = vol3d(3.0, 2.0, 1.0, a);
  mat ma1(v, eps, 1);
  mat ma(v, eps, splitting);
  ma.set_output_directory(dirname);
  ma1.set_output_directory(dirname);

  master_printf("Periodic test using %d chunks...\n", splitting);
  fields f(&ma);
  f.use_bloch(vec(0.1,0.7,0.3));
  f.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.3,0.5,0.5), 1.0);
  fields f1(&ma1);
  f1.use_bloch(vec(0.1,0.7,0.3));
  f1.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.3,0.5,0.5), 1.0);
  double total_energy_check_time = 8.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.5  , 0.01, 0.5  ))) return 0;
    if (!compare_point(f, f1, vec(0.46 , 0.33, 0.2  ))) return 0;
    if (!compare_point(f, f1, vec(1.0  , 1.0 , 0.301))) return 0;
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

int test_pml(double eps(const vec &), const char *dirname) {
  double a = 10.0;

  volume v = vol3d(3.0, 2.0, 2.0, a);
  mat ma(v, eps, 0);
  ma.set_output_directory(dirname);
  ma.use_pml_everywhere(0.5);

  master_printf("Testing pml quality...\n");
  fields f(&ma);
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(1.299,1.299,0.801), 1.0);
  const double deltaT = 20.0;
  const double ttot = 3.1*deltaT;
  double total_energy_check_time = deltaT;

  while (f.time() < f.find_last_source()) f.step();

  double last_energy = f.total_energy();
  while (f.time() < ttot) {
    f.step();
    if (f.time() >= total_energy_check_time) {
      const double new_energy = f.total_energy();
      if (new_energy > last_energy*1e-5) { // FIXME: problem here? this is pretty slow...
        master_printf("Energy decaying too slowly: from %lg to %lg (%lg)\n",
                      last_energy, new_energy, new_energy/last_energy);
        return 0;
      } else {
        master_printf("Got newE/oldE of %lg\n", new_energy/last_energy);
      }
      total_energy_check_time += deltaT;
    }
  }
  return 1;
}

int test_pml_splitting(double eps(const vec &), int splitting, const char *dirname) {
  double a = 10.0;

  volume v = vol3d(3.0, 2.0, 2.0, a);
  mat ma1(v, eps, 1);
  mat ma(v, eps, splitting);
  ma.set_output_directory(dirname);
  ma1.set_output_directory(dirname);
  ma.use_pml_everywhere(0.3);
  ma1.use_pml_everywhere(0.3);

  master_printf("Testing pml while splitting into %d chunks...\n", splitting);
  fields f(&ma);
  f.add_point_source(Ez, 0.8, 1.6, 0.0, 4.0, vec(1.299,1.299,0.401), 1.0);
  fields f1(&ma1);
  f1.add_point_source(Ez, 0.8, 1.6, 0.0, 4.0, vec(1.299,1.299,0.401), 1.0);
  const double ttot = 31.0;

  double next_energy_time = 10.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    //f.output_real_imaginary_slices("multi");
    //f1.output_real_imaginary_slices("single");
    if (!approx_point(f, f1, vec(0.5  , 0.01 , 1.0 ))) return 0;
    if (!approx_point(f, f1, vec(0.46 , 0.33 , 0.33))) return 0;
    if (!approx_point(f, f1, vec(1.0  , 1.0  , 0.33))) return 0;
    if (!approx_point(f, f1, vec(1.3  , 1.3  , 0.15))) return 0;
    if (f.time() > next_energy_time) {
      if (!compare(f.total_energy(), f1.total_energy(),
                   "   total energy")) return 0;
      next_energy_time += 10.0;
    }
  }
  return 1;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  const char *dirname = "three_d-out";
  trash_output_directory(dirname);
  master_printf("Testing 3D...\n");

  if (!test_pml(one, dirname)) abort("error in test_pml vacuum\n");

  for (int s=2;s<7;s++)
    if (!test_periodic(targets, s, dirname))
      abort("error in test_periodic targets\n");

  for (int s=2;s<8;s++)
    if (!test_metal(one, s, dirname)) abort("error in test_metal vacuum\n");

  for (int s=2;s<5;s++)
    if (!test_pml_splitting(one, s, dirname))
      abort("error in test_pml_splitting vacuum\n");

  exit(0);
}

