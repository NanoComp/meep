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
double targets(const vec &pt) {
  const double r = sqrt(pt.x() * pt.x() + pt.y() * pt.y());
  double dr = r;
  while (dr > 1)
    dr -= 1;
  if (dr > 0.7001) return 12.0;
  return 1.0;
}

#if MEEP_SINGLE
static const double tol = 1e-3, thresh = 1e-5;
#else
static const double tol = 1e-11, thresh = 1e-5;
#endif

int compare(double a, double b, const char *n) {
  if (fabs(a - b) > fabs(b) * tol && fabs(b) > thresh) {
    master_printf("%s differs by\t%g out of\t%g\n", n, a - b, b);
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
      if (abs(v1 - v2) > tol * abs(v2) && abs(v2) > thresh) {
        master_printf("%s differs:  %g %g out of %g %g\n", component_name(c), real(v2 - v1),
                      imag(v2 - v1), real(v2), imag(v2));
        master_printf("This comes out to a fractional error of %g\n", abs(v1 - v2) / abs(v2));
        master_printf("Right now I'm looking at %g %g, time %g\n", p.x(), p.y(), f1.time());
        return 0;
      }
    }
  }
  return 1;
}

int test_metal(double eps(const vec &), int splitting) {
  double a = 10.0;
  double ttot = 17.0;

  grid_volume gv = voltwo(3.0, 2.0, a);
  structure s1(gv, eps);
  structure s(gv, eps, no_pml(), identity(), splitting);

  s.add_susceptibility(one, E_stuff, lorentzian_susceptibility(0.3, 0.1));
  s1.add_susceptibility(one, E_stuff, lorentzian_susceptibility(0.3, 0.1));

  master_printf("Metal+dispersion test using %d chunks...\n", splitting);
  fields f(&s);
  f.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec(0.3, 0.5), 1.0);
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(1.299, 0.401), 1.0);
  fields f1(&s1);
  f1.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec(0.3, 0.5), 1.0);
  f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(1.299, 0.401), 1.0);
  double field_energy_check_time = 8.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.5, 0.01))) return 0;
    if (!compare_point(f, f1, vec(0.46, 0.33))) return 0;
    if (!compare_point(f, f1, vec(1.0, 1.0))) return 0;
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
  return 1;
}

int test_periodic(double eps(const vec &), int splitting) {
  double a = 10.0;
  double ttot = 17.0;

  grid_volume gv = voltwo(3.0, 2.0, a);
  structure s1(gv, eps);
  structure s(gv, eps, no_pml(), identity(), splitting);

  master_printf("Periodic test using %d chunks...\n", splitting);
  fields f(&s);
  f.use_bloch(vec(0.1, 0.7));
  f.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec(0.3, 0.5), 1.0);
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(1.299, 0.401), 1.0);
  fields f1(&s1);
  f1.use_bloch(vec(0.1, 0.7));
  f1.add_point_source(Hz, 0.7, 2.5, 0.0, 4.0, vec(0.3, 0.5), 1.0);
  f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(1.299, 0.401), 1.0);
  double field_energy_check_time = 8.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.5, 0.01))) return 0;
    if (!compare_point(f, f1, vec(0.46, 0.33))) return 0;
    if (!compare_point(f, f1, vec(1.0, 1.0))) return 0;
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
  return 1;
}

int test_periodic_tm(double eps(const vec &), int splitting) {
  double a = 10.0;
  double ttot = 17.0;

  grid_volume gv = voltwo(3.0, 2.0, a);
  structure s1(gv, eps);
  structure s(gv, eps, no_pml(), identity(), splitting);

  master_printf("Periodic 2D TM test using %d chunks...\n", splitting);
  fields f(&s);
  f.use_bloch(vec(0.1, 0.7));
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(1.299, 0.401), 1.0);
  fields f1(&s1);
  f1.use_bloch(vec(0.1, 0.7));
  f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(1.299, 0.401), 1.0);
  double field_energy_check_time = 8.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.5, 0.01))) return 0;
    if (!compare_point(f, f1, vec(0.46, 0.33))) return 0;
    if (!compare_point(f, f1, vec(1.0, 1.0))) return 0;
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
  return 1;
}

int test_pml(double eps(const vec &), int splitting) {
  double a = 10.0;

  grid_volume gv = voltwo(3.0, 2.0, a);
  structure s1(gv, eps, pml(1.0, X) + pml(1.0, Y, High));
  structure s(gv, eps, pml(1.0, X) + pml(1.0, Y, High), identity(), splitting);

  master_printf("Testing pml while splitting into %d chunks...\n", splitting);
  fields f(&s);
  f.add_point_source(Hz, 0.7, 1.5, 0.0, 4.0, vec(1.5, 0.5), 1.0);
  f.add_point_source(Ez, 0.8, 1.6, 0.0, 4.0, vec(1.299, 0.401), 1.0);
  fields f1(&s1);
  f1.add_point_source(Hz, 0.7, 1.5, 0.0, 4.0, vec(1.5, 0.5), 1.0);
  f1.add_point_source(Ez, 0.8, 1.6, 0.0, 4.0, vec(1.299, 0.401), 1.0);
  const double deltaT = 100.0;
  const double ttot = 3.1 * deltaT;
  double field_energy_check_time = deltaT;

  while (f.time() < f.last_source_time())
    f.step();
  while (f1.time() < f1.last_source_time())
    f1.step();

  double last_energy = f.field_energy();
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (f.time() >= field_energy_check_time) {
      if (!compare_point(f, f1, vec(0.5, 0.01))) return 0;
      if (!compare_point(f, f1, vec(0.46, 0.33))) return 0;
      if (!compare_point(f, f1, vec(1.0, 1.0))) return 0;
      const double new_energy = f.field_energy();
      if (!compare(new_energy, f1.field_energy(), "   total energy")) return 0;
      if (new_energy > last_energy * 1e-6) {
        master_printf("Energy decaying too slowly: from %g to %g (%g)\n", last_energy, new_energy,
                      new_energy / last_energy);
        return 0;
      }
      else { master_printf("Got newE/oldE of %g\n", new_energy / last_energy); }
      field_energy_check_time += deltaT;
    }
  }
  return 1;
}

int test_pml_tm(double eps(const vec &), int splitting) {
  double a = 10.0;

  grid_volume gv = voltwo(3.0, 3.0, a);
  structure s1(gv, eps, pml(1.0));
  structure s(gv, eps, pml(1.0), identity(), splitting);

  master_printf("Testing TM pml while splitting into %d chunks...\n", splitting);
  fields f(&s);
  f.add_point_source(Ez, 0.8, 1.6, 0.0, 4.0, vec(1.299, 1.401), 1.0);
  fields f1(&s1);
  f1.add_point_source(Ez, 0.8, 1.6, 0.0, 4.0, vec(1.299, 1.401), 1.0);
  const double deltaT = 100.0;
  const double ttot = 3.1 * deltaT;
  double field_energy_check_time = deltaT;

  while (f.time() < f.last_source_time())
    f.step();
  while (f1.time() < f1.last_source_time())
    f1.step();

  double last_energy = f.field_energy();
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (f.time() >= field_energy_check_time) {
      if (!compare_point(f, f1, vec(0.5, 0.01))) return 0;
      if (!compare_point(f, f1, vec(0.46, 0.33))) return 0;
      if (!compare_point(f, f1, vec(1.0, 1.0))) return 0;
      const double new_energy = f.field_energy();
      if (!compare(new_energy, f1.field_energy(), "   total energy")) return 0;
      if (new_energy > last_energy * 4e-6) {
        master_printf("Energy decaying too slowly: from %g to %g (%g)\n", last_energy, new_energy,
                      new_energy / last_energy);
        return 0;
      }
      else { master_printf("Got newE/oldE of %g\n", new_energy / last_energy); }
      field_energy_check_time += deltaT;
    }
  }
  return 1;
}

int test_pml_te(double eps(const vec &), int splitting) {
  double a = 10.0;

  grid_volume gv = voltwo(3.0, 3.0, a);
  structure s1(gv, eps, pml(1.0));
  structure s(gv, eps, pml(1.0), identity(), splitting);

  master_printf("Testing TE pml while splitting into %d chunks...\n", splitting);
  fields f(&s);
  f.add_point_source(Hz, 0.7, 1.5, 0.0, 4.0, vec(1.5, 1.5), 1.0);
  f.add_point_source(Hz, 0.7, 1.5, 0.0, 4.0, vec(1.37, 1.27), 1.0);
  fields f1(&s1);
  f1.add_point_source(Hz, 0.7, 1.5, 0.0, 4.0, vec(1.5, 1.5), 1.0);
  f1.add_point_source(Hz, 0.7, 1.5, 0.0, 4.0, vec(1.37, 1.27), 1.0);
  const double deltaT = 100.0;
  const double ttot = 3.1 * deltaT;
  double field_energy_check_time = deltaT;

  while (f.time() < f.last_source_time())
    f.step();
  while (f1.time() < f1.last_source_time())
    f1.step();

  double last_energy = f.field_energy();
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (f.time() >= field_energy_check_time) {
      if (!compare_point(f, f1, vec(0.5, 0.01))) return 0;
      if (!compare_point(f, f1, vec(0.46, 0.33))) return 0;
      if (!compare_point(f, f1, vec(1.0, 1.0))) return 0;
      const double new_energy = f.field_energy();
      if (!compare(new_energy, f1.field_energy(), "   total energy")) return 0;
      if (new_energy > last_energy * 1.1e-6) {
        master_printf("Energy decaying too slowly: from %g to %g (%g)\n", last_energy, new_energy,
                      new_energy / last_energy);
        return 0;
      }
      else { master_printf("Got newE/oldE of %g\n", new_energy / last_energy); }
      field_energy_check_time += deltaT;
    }
  }
  return 1;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 0;
  master_printf("Testing 2D...\n");

  for (int s = 2; s < 4; s++)
    if (!test_pml(one, s)) meep::abort("error in test_pml vacuum\n");

  for (int s = 2; s < 4; s++)
    if (!test_pml_tm(one, s)) meep::abort("error in test_pml_tm vacuum\n");

  for (int s = 2; s < 4; s++)
    if (!test_pml_te(one, s)) meep::abort("error in test_pml_te vacuum\n");

  for (int s = 2; s < 4; s++)
    if (!test_metal(one, s)) meep::abort("error in test_metal vacuum\n");
  // if (!test_metal(one, 200)) meep::abort("error in test_metal vacuum\n");

  for (int s = 2; s < 5; s++)
    if (!test_metal(targets, s)) meep::abort("error in test_metal targets\n");
  // if (!test_metal(targets, 60)) meep::abort("error in test_metal targets\n");

  for (int s = 2; s < 5; s++)
    if (!test_periodic(targets, s)) meep::abort("error in test_periodic targets\n");
  // if (!test_periodic(one, 200))
  //  meep::abort("error in test_periodic targets\n");

  for (int s = 2; s < 4; s++)
    if (!test_periodic_tm(one, s)) meep::abort("error in test_periodic_tm vacuum\n");

  return 0;
}
