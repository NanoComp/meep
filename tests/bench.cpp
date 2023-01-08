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

double one(const vec &) { return 1.0; }
static double width = 20.0;
double bump(const vec &pt) { return (fabs(pt.z() - 50.0) > width) ? 1.0 : 12.0; }

struct bench {
  double time; // In seconds.
  double gridsteps;
};

bench bench_periodic(const double rmax, const double zmax, double eps(const vec &)) {
  const double a = 10.0;
  const double gridpts = (zmax == 0.0) ? a * rmax : a * a * rmax * zmax;
  const double ttot = 5.0 + 1e5 / gridpts;
  const int m = 0;

  grid_volume gv = volcyl(rmax, zmax, a);
  structure s(gv, eps);
  fields f(&s, m);
  f.use_bloch(0.0);
  f.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, veccyl(0.5, 0.4), 1.0);
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, veccyl(0.401, 0.301), 1.0);

  double start = wall_time();
  while (f.time() < ttot)
    f.step();
  bench b;
  b.time = (wall_time() - start);
  b.gridsteps = ttot * a * 2 * gridpts;
  // f.print_times();
  return b;
}

bench bench_flux_1d(const double zmax, double eps(const vec &)) {
  const double a = 10.0;
  const double gridpts = a * zmax;
  const double ttot = 10.0 + 1e5 / zmax;

  grid_volume gv = volone(zmax, a);
  structure s(gv, eps, pml(zmax / 6));

  fields f(&s);
  f.use_real_fields();
  f.add_point_source(Ex, 0.7, 2.5, 0.0, 3.0, vec(zmax / 2 + 0.3), 1.0);
  flux_vol *left = f.add_flux_plane(vec(zmax / 3.0), vec(zmax / 3.0));
  flux_vol *right = f.add_flux_plane(vec(zmax * 2.0 / 3.0), vec(zmax * 2.0 / 3.0));

  while (f.time() <= f.last_source_time())
    f.step();

  grid_volume mid = volone(zmax / 3, a);
  mid.set_origin(vec(zmax / 3));
  double flux_energy = 0.0;
  double start = wall_time();
  while (f.time() < ttot) {
    f.step();
    flux_energy += f.dt * (right->flux() - left->flux());
  }
  bench b;
  b.time = (wall_time() - start);
  b.gridsteps = ttot * a * 2 * gridpts;
  // f.print_times();
  return b;
}

bench bench_2d(const double xmax, const double ymax, double eps(const vec &)) {
  const double a = 10.0;
  const double gridpts = a * a * xmax * ymax;
  const double ttot = 5.0 + 1e5 / gridpts;

  grid_volume gv = voltwo(xmax, ymax, a);
  structure s(gv, eps);
  fields f(&s);
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.401, 0.301));
  f.add_point_source(Hz, 0.8, 0.7, 0.0, 4.0, vec(0.431, 0.2));

  while (f.time() < f.last_source_time())
    f.step();
  const double tend = f.time() + ttot;
  double start = wall_time();
  while (f.time() < tend)
    f.step();
  bench b;
  b.time = (wall_time() - start);
  b.gridsteps = ttot * a * 2 * gridpts;
  // f.print_times();
  return b;
}

const double te_tm_2d_time = 2e5;

bench bench_2d_tm_nonlinear(const double xmax, const double ymax, double eps(const vec &)) {
  const double a = 10.0;
  const double gridpts = a * a * xmax * ymax;
  const double ttot = 5.0 + te_tm_2d_time / gridpts;

  grid_volume gv = voltwo(xmax, ymax, a);
  structure s(gv, eps);
  s.set_chi3(eps);
  fields f(&s);
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.401, 0.301));

  while (f.time() < f.last_source_time())
    f.step();
  const double tend = f.time() + ttot;
  double start = wall_time();
  while (f.time() < tend)
    f.step();
  bench b;
  b.time = (wall_time() - start);
  b.gridsteps = ttot * a * 2 * gridpts;
  // f.print_times();
  return b;
}

bench bench_2d_tm(const double xmax, const double ymax, double eps(const vec &)) {
  const double a = 10.0;
  const double gridpts = a * a * xmax * ymax;
  const double ttot = 5.0 + te_tm_2d_time / gridpts;

  grid_volume gv = voltwo(xmax, ymax, a);
  structure s(gv, eps);
  fields f(&s);
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.401, 0.301));

  while (f.time() < f.last_source_time())
    f.step();
  const double tend = f.time() + ttot;
  double start = wall_time();
  while (f.time() < tend)
    f.step();
  bench b;
  b.time = (wall_time() - start);
  b.gridsteps = ttot * a * 2 * gridpts;
  // f.print_times();
  return b;
}

bench bench_2d_te(const double xmax, const double ymax, double eps(const vec &)) {
  const double a = 10.0;
  const double gridpts = a * a * xmax * ymax;
  const double ttot = 5.0 + te_tm_2d_time / gridpts;

  grid_volume gv = voltwo(xmax, ymax, a);
  structure s(gv, eps);
  fields f(&s);
  f.add_point_source(Ex, 0.8, 0.6, 0.0, 4.0, vec(0.401, 0.301));
  f.add_point_source(Hz, 0.6, 0.6, 0.0, 4.0, vec(0.7, 0.5));

  while (f.time() < f.last_source_time())
    f.step();
  const double tend = f.time() + ttot;
  double start = wall_time();
  while (f.time() < tend)
    f.step();
  bench b;
  b.time = (wall_time() - start);
  b.gridsteps = ttot * a * 2 * gridpts;
  // f.print_times();
  return b;
}

bench bench_2d_te_nonlinear(const double xmax, const double ymax, double eps(const vec &)) {
  const double a = 10.0;
  const double gridpts = a * a * xmax * ymax;
  const double ttot = 5.0 + te_tm_2d_time / gridpts;

  grid_volume gv = voltwo(xmax, ymax, a);
  structure s(gv, eps);
  s.set_chi3(eps);
  fields f(&s);
  f.add_point_source(Ex, 0.8, 0.6, 0.0, 4.0, vec(0.401, 0.301));
  f.add_point_source(Hz, 0.6, 0.6, 0.0, 4.0, vec(0.7, 0.5));

  while (f.time() < f.last_source_time())
    f.step();
  const double tend = f.time() + ttot;
  double start = wall_time();
  while (f.time() < tend)
    f.step();
  bench b;
  b.time = (wall_time() - start);
  b.gridsteps = ttot * a * 2 * gridpts;
  // f.print_times();
  return b;
}

#define showbench(name, bb)                                                                        \
  {                                                                                                \
    bench b = bb;                                                                                  \
    master_printf("bench:, %s, %g, %g\n", name, b.time, b.time * 1e6 / b.gridsteps);               \
  }

// 3D benchmarks:

inline double max(double a, double b) { return (a > b) ? a : b; }

bench bench_3d_periodic(const double xmax, const double ymax, const double zmax,
                        double eps(const vec &)) {
  const double a = 10.0;
  const double gridpts = a * a * a * max(xmax, 1 / a) * max(ymax, 1 / a) * max(zmax, 1 / a);
  const double ttot = 5.0 + 1e5 / gridpts;

  grid_volume gv = vol3d(xmax, ymax, zmax, a);
  structure s(gv, eps);
  fields f(&s);
  if (xmax == 0) f.use_bloch(X, 0.0);
  if (ymax == 0) f.use_bloch(Y, 0.0);
  if (ymax == 0) f.use_bloch(Z, 0.0);
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(xmax * .5, ymax * .5, zmax * .5));

  while (f.time() < f.last_source_time())
    f.step();
  const double tend = f.time() + ttot;
  double start = wall_time();
  while (f.time() < tend)
    f.step();
  bench b;
  b.time = (wall_time() - start);
  b.gridsteps = ttot * a * 2 * gridpts;
  // f.print_times();
  return b;
}

bench bench_3d(const double xmax, const double ymax, const double zmax, double eps(const vec &)) {
  const double a = 10.0;
  const double gridpts = a * a * a * xmax * ymax * zmax;
  const double ttot = 5.0 + 1e5 / gridpts;

  grid_volume gv = vol3d(xmax, ymax, zmax, a);
  structure s(gv, eps);
  fields f(&s);
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(xmax * .5, ymax * .5, zmax * .5));

  while (f.time() < f.last_source_time())
    f.step();
  const double tend = f.time() + ttot;
  double start = wall_time();
  while (f.time() < tend)
    f.step();
  bench b;
  b.time = (wall_time() - start);
  b.gridsteps = ttot * a * 2 * gridpts;
  // f.print_times();
  return b;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 0;
  master_printf("Benchmarking with %d processor%s...\n", count_processors(),
                count_processors() > 1 ? "s" : "");

  master_printf("bench:, test, total time (s), normalized time (s/Mgs)\n");

  showbench("Periodic 6x4 ", bench_periodic(6.0, 4.0, one));
  showbench("Periodic 12x1", bench_periodic(12.0, 1.0, one));
  showbench("Periodic 1x12", bench_periodic(1.0, 12.0, one));
  showbench("Periodic 12x0", bench_periodic(12.0, 0.0, one));
  showbench("Periodic 12x12", bench_periodic(12.0, 12.0, one));

  width = 20.0;
  showbench("Flux 1D 100", bench_flux_1d(100.0, bump));
  width = 10.0;
  showbench("Flux 1D 100", bench_flux_1d(100.0, bump));
  width = 300.0;
  showbench("Flux 1D 100", bench_flux_1d(100.0, bump));

  showbench("3D 1x1x10", bench_3d(1.0, 1.0, 10.0, one));
  showbench("3D 10x1x1", bench_3d(10.0, 1.0, 1.0, one));
  showbench("3D 1x1x1 ", bench_3d(1.0, 1.0, 1.0, one));
  showbench("3D 3x3x3 ", bench_3d(3.0, 3.0, 3.0, one));
  showbench("3D 10x3x0", bench_3d_periodic(10.0, 3.0, 0.0, one));
  showbench("3D 0x3x10", bench_3d_periodic(0.0, 3.0, 10.0, one));

  showbench("2D 6x4 ", bench_2d(6.0, 4.0, one));
  showbench("2D 12x12 ", bench_2d(12.0, 12.0, one));
  showbench("2D 12x12 ", bench_2d(12.0, 12.0, one));

  showbench("2D TM 6x4 nonlinear ", bench_2d_tm_nonlinear(6.0, 4.0, one));
  showbench("2D TM 6x4 ", bench_2d_tm(6.0, 4.0, one));
  showbench("2D TM 12x12 ", bench_2d_tm(12.0, 12.0, one));

  showbench("2D TE 2x2 nonlinear ", bench_2d_te_nonlinear(2.0, 2.0, one));
  showbench("2D TE 2x2 ", bench_2d_te(2.0, 2.0, one));
  showbench("2D TE 10x11 nonlinear ", bench_2d_te_nonlinear(10.0, 11.0, one));
  showbench("2D TE 10x11 ", bench_2d_te(10.0, 11.0, one));

  master_printf("\nnote: 1 Mgs = 1 million grid point time steps\n");

  return 0;
}
