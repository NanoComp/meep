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
#include <time.h>

#include "dactyl.h"

double one(const vec &) { return 1.0; }

struct bench {
  double time; // In seconds.
  double gridsteps;
};

bench bench_periodic(const double rmax, const double zmax,
                     double eps(const vec &)) {
  const double a = 10.0;
  const double gridpts = (zmax==0.0)?a*rmax:a*a*rmax*zmax;
  const double ttot = 5.0 + 1e5/gridpts;
  const int m = 0;

  volume v = volcyl(rmax,zmax,a);
  mat ma(v, eps);
  fields f(&ma, m);
  f.use_bloch(0.0);
  f.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, vec(0.5, 0.4), 1.0);
  f.add_point_source(Ez, 0.8, 0.6, 0.0, 4.0, vec(0.401, 0.301), 1.0);

  clock_t start = clock();
  while (f.time() < ttot) f.step();
  bench b;
  b.time = (clock()-start)*(1.0/CLOCKS_PER_SEC);
  b.gridsteps = ttot*a*2*gridpts;
  return b;
}

void showbench(const char *name, const bench &b) {
  master_printf("%s\n  total time:    \t%lg s\n  normalized time:\t%lg s/Mgs\n",
                name, b.time, b.time*1e6/b.gridsteps);
}

int main(int argc, char **argv) {
  initialize(argc, argv);
  master_printf("Benchmarking...\n");

  showbench("Periodic 6x4 ", bench_periodic(6.0, 4.0, one));
  showbench("Periodic 12x1", bench_periodic(12.0, 1.0, one));
  showbench("Periodic 1x12", bench_periodic(1.0, 12.0, one));
  showbench("Periodic 12x0", bench_periodic(12.0, 0.0, one));
  showbench("Periodic 12x12", bench_periodic(12.0, 12.0, one));

  master_printf("\nnote: 1 Mgs = 1 million grid point time steps\n");

  finished();
  exit(0);
}

