/* Copyright (C) 2005-2021 Massachusetts Institute of Technology
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

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>

#include <meep.hpp>
using namespace meep;
using std::complex;

double one(const vec &) { return 1.0; }
double targets(const vec &pt) {
  const double r = sqrt(pt.x() * pt.x() + pt.y() * pt.y());
  double dr = r;
  while (dr > 1) dr -= 1;
  if (dr > 0.7001) return 12.0;
  return 1.0;
}

#if MEEP_SINGLE
static const double tol = 1e-3, thresh = 1e-3;
#else
static const double tol = 1e-9, thresh = 1e-15;
#endif

int compare(double a, double b, const char *n) {
  if (fabs(a - b) > fabs(b) * tol && fabs(b) > thresh) {
    master_printf("%s differs by\t%g out of\t%g\n", n, a - b, b);
    master_printf("This gives a fractional error of %g\n",
                  fabs(a - b) / fabs(b));
    return 0;
  } else {
    return 1;
  }
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
        master_printf("%s differs:  %g %g out of %g %g\n", component_name(c),
                      real(v2 - v1), imag(v2 - v1), real(v2), imag(v2));
        master_printf("This comes out to a fractional error of %g\n",
                      abs(v1 - v2) / abs(v2));
        master_printf("Right now I'm looking at %g %g %g, time %g\n", p.x(),
                      p.y(), p.z(), f1.time());
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
  for (int i = 0; i < 10; i++) {
    component c = (component)i;
    if (f1.gv.has_field(c)) {
      complex<double> v1 = m_test.get_component(c), v2 = m1.get_component(c);
      if (abs(v1 - v2) > tol * abs(v2) && abs(v2) > thresh) {
        master_printf("%s differs:  %g %g out of %g %g\n", component_name(c),
                      real(v2 - v1), imag(v2 - v1), real(v2), imag(v2));
        master_printf("This comes out to a fractional error of %g\n",
                      abs(v1 - v2) / abs(v2));
        master_printf("Right now I'm looking at %g %g %g, time %g\n", p.x(),
                      p.y(), p.z(), f1.time());
        return 0;
      }
    }
  }
  return 1;
}

std::string structure_dump(structure *s, const std::string &filename_prefix,
                           const std::string &output_name) {
  std::string filename = filename_prefix + "-structure-" + output_name;
  s->dump(filename.c_str());
  master_printf("Dumping structure: %s\n", filename.c_str());
  return filename;
}

void structure_load(structure *s, const std::string &filename) {
  master_printf("Loading structure: %s\n", filename.c_str());
  s->load(filename.c_str());
}

std::string fields_dump(fields *f, const std::string &filename_prefix,
                        const std::string &output_name) {
  std::string filename = filename_prefix + "-fields-" + output_name;
  f->dump(filename.c_str());
  master_printf("Dumping fields: %s\n", filename.c_str());
  return filename;
}

void fields_load(fields *f, const std::string &filename) {
  master_printf("Loading fields: %s\n", filename.c_str());
  f->load(filename.c_str());
}

int test_metal(double eps(const vec &), int splitting, const char *tmpdir) {
  double a = 10.0;
  double ttot = 17.0;

  grid_volume gv = vol3d(1.5, 0.5, 1.0, a);
  structure s1(gv, eps);
  structure s(gv, eps, no_pml(), identity(), splitting);

  std::string filename_prefix =
      std::string(tmpdir) + "/test_metal_" + std::to_string(splitting);
  std::string structure_filename_original =
      structure_dump(&s, filename_prefix, "original");

  master_printf("Metal test using %d chunks...\n", splitting);
  fields f(&s);
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(1.299, 0.299, 0.401), 1.0);

  std::string fields_filename_original =
      fields_dump(&f, filename_prefix, "original");

  fields f1(&s1);
  f1.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(1.299, 0.299, 0.401), 1.0);

  double field_energy_check_time = 8.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.5, 0.5, 0.01))) return 0;
    if (!compare_point(f, f1, vec(0.46, 0.33, 0.33))) return 0;
    if (!compare_point(f, f1, vec(1.301, 0.301, 0.399))) return 0;
    if (f.time() >= field_energy_check_time) {
      if (!compare(f.field_energy(), f1.field_energy(), "   total energy"))
        return 0;
      if (!compare(f.electric_energy_in_box(gv.surroundings()),
                   f1.electric_energy_in_box(gv.surroundings()),
                   "electric energy"))
        return 0;
      if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                   f1.magnetic_energy_in_box(gv.surroundings()),
                   "magnetic energy"))
        return 0;
      field_energy_check_time += 5.0;
    }
  }

  std::string structure_filename_after_sim =
      structure_dump(&s, filename_prefix, "after-sim");

  std::string fields_filename_after_sim =
      fields_dump(&f, filename_prefix, "after_sim");

  structure s_load(gv, eps, no_pml(), identity(), splitting);
  structure_load(&s_load, structure_filename_after_sim);

  std::string structure_filename_dump_loaded =
      structure_dump(&s, filename_prefix, "dump-loaded");

  fields f_load(&s_load);
  f_load.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(1.299, 0.299, 0.401),
                          1.0);
  fields_load(&f_load, fields_filename_after_sim);

  std::string fields_filename_dump_loaded =
      fields_dump(&f_load, filename_prefix, "dump-loaded");
  return 1;
}

int test_periodic(double eps(const vec &), int splitting, const char *tmpdir) {
  double a = 10.0;
  double ttot = 17.0;

  grid_volume gv = vol3d(1.5, 0.5, 1.0, a);
  structure s1(gv, eps);
  structure s(gv, eps, no_pml(), identity(), splitting);

  std::string filename_prefix =
      std::string(tmpdir) + "/test_periodic_" + std::to_string(splitting);
  std::string structure_filename_original =
      structure_dump(&s, filename_prefix, "original");

  master_printf("Periodic test using %d chunks...\n", splitting);
  fields f(&s);
  f.use_bloch(vec(0.1, 0.7, 0.3));
  f.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.3, 0.25, 0.5), 1.0);
  std::string fields_filename_original =
      fields_dump(&f, filename_prefix, "original");

  fields f1(&s1);
  f1.use_bloch(vec(0.1, 0.7, 0.3));
  f1.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.3, 0.25, 0.5), 1.0);

  double field_energy_check_time = 8.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!compare_point(f, f1, vec(0.5, 0.01, 0.5))) return 0;
    if (!compare_point(f, f1, vec(0.46, 0.33, 0.2))) return 0;
    if (!compare_point(f, f1, vec(1.0, 0.25, 0.301))) return 0;
    if (f.time() >= field_energy_check_time) {
      if (!compare(f.field_energy(), f1.field_energy(), "   total energy"))
        return 0;
      if (!compare(f.electric_energy_in_box(gv.surroundings()),
                   f1.electric_energy_in_box(gv.surroundings()),
                   "electric energy"))
        return 0;
      if (!compare(f.magnetic_energy_in_box(gv.surroundings()),
                   f1.magnetic_energy_in_box(gv.surroundings()),
                   "magnetic energy"))
        return 0;
      field_energy_check_time += 5.0;
    }
  }

  std::string structure_filename_after_sim =
      structure_dump(&s, filename_prefix, "after-sim");

  std::string fields_filename_after_sim =
      fields_dump(&f, filename_prefix, "after_sim");

  structure s_load(gv, eps, no_pml(), identity(), splitting);
  structure_load(&s_load, structure_filename_after_sim);

  std::string structure_filename_dump_loaded =
      structure_dump(&s, filename_prefix, "dump-loaded");

  fields f_load(&s_load);
  f_load.use_bloch(vec(0.1, 0.7, 0.3));
  f_load.add_point_source(Ez, 0.7, 2.5, 0.0, 4.0, vec(0.3, 0.25, 0.5), 1.0);
  fields_load(&f_load, fields_filename_after_sim);

  std::string fields_filename_dump_loaded =
      fields_dump(&f_load, filename_prefix, "dump-loaded");
  return 1;
}

int test_pml(double eps(const vec &), const char *tmpdir) {
  master_printf("Testing pml quality...\n");

  double a = 10.0;
  grid_volume gv = vol3d(1.5, 1.0, 1.2, a);
  structure s(gv, eps, pml(0.401));

  std::string filename_prefix = std::string(tmpdir) + "/test_pml";
  std::string structure_filename_original =
      structure_dump(&s, filename_prefix, "original");

  fields f(&s);
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.751, 0.5, 0.601), 1.0);
  const double deltaT = 10.0;
  const double ttot = 3.1 * deltaT;
  double field_energy_check_time = deltaT;

  std::string fields_filename_original =
      fields_dump(&f, filename_prefix, "original");

  while (f.time() < f.last_source_time()) f.step();

  double last_energy = f.field_energy();
  while (f.time() < ttot) {
    f.step();
    if (f.time() >= field_energy_check_time) {
      const double new_energy = f.field_energy();
      master_printf("Got newE/oldE of %g\n", new_energy / last_energy);
      if (new_energy > last_energy * 4e-3) {
        master_printf("Energy decaying too slowly: from %g to %g (%g)\n",
                      last_energy, new_energy, new_energy / last_energy);
        return 0;
      }
      field_energy_check_time += deltaT;
    }
  }

  std::string structure_filename_after_sim =
      structure_dump(&s, filename_prefix, "after-sim");

  std::string fields_filename_after_sim =
      fields_dump(&f, filename_prefix, "after_sim");

  structure s_load(gv, eps, pml(0.401));
  structure_load(&s_load, structure_filename_after_sim);

  std::string structure_filename_dump_loaded =
      structure_dump(&s, filename_prefix, "dump-loaded");

  fields f_load(&s_load);
  fields_load(&f_load, fields_filename_after_sim);

  std::string fields_filename_dump_loaded =
      fields_dump(&f_load, filename_prefix, "dump-loaded");
  return 1;
}

int test_pml_splitting(double eps(const vec &), int splitting,
                       const char *tmpdir) {
  double a = 10.0;

  grid_volume gv = vol3d(1.5, 1.0, 1.2, a);
  structure s1(gv, eps, pml(0.3));
  structure s(gv, eps, pml(0.3), identity(), splitting);

  master_printf("Testing pml while splitting into %d chunks...\n", splitting);
  fields f(&s);
  f.add_point_source(Ez, 0.8, 1.6, 0.0, 4.0, vec(1.099, 0.499, 0.501), 1.0);
  fields f1(&s1);
  f1.add_point_source(Ez, 0.8, 1.6, 0.0, 4.0, vec(1.099, 0.499, 0.501), 1.0);
  const double ttot = 31.0;

  double next_energy_time = 10.0;
  while (f.time() < ttot) {
    f.step();
    f1.step();
    if (!approx_point(f, f1, vec(0.5, 0.01, 1.0))) return 0;
    if (!approx_point(f, f1, vec(0.46, 0.33, 0.33))) return 0;
    if (!approx_point(f, f1, vec(1.0, 1.0, 0.33))) return 0;
    if (!approx_point(f, f1, vec(1.3, 0.3, 0.15))) return 0;
    if (f.time() > next_energy_time) {
      if (!compare(f.field_energy(), f1.field_energy(), "   total energy"))
        return 0;
      next_energy_time += 10.0;
    }
  }
  return 1;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 0;

  std::unique_ptr<const char[]> temp_dir(make_output_directory());
  master_printf("Testing 3D dump/load: temp_dir = %s...\n", temp_dir.get());

  if (!test_pml(one, temp_dir.get())) abort("error in test_pml vacuum\n");

  for (int s = 2; s < 7; s++)
    if (!test_periodic(targets, s, temp_dir.get()))
      abort("error in test_periodic targets\n");

  for (int s = 2; s < 8; s++)
    if (!test_metal(one, s, temp_dir.get()))
      abort("error in test_metal vacuum\n");

  for (int s = 2; s < 4; s++)
    if (!test_pml_splitting(one, s, temp_dir.get()))
      abort("error in test_pml_splitting vacuum\n");

  // delete_directory(temp_dir.get());
  return 0;
}
