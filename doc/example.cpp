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

double guided_eps(const vec &x) {
  return 1.0;
  if (x.z() < 1.0) return 1.0;
  double rcore = 1.14;
  double rr = x.r() - rcore;
  if (rr > 6.001) return 1.0; // outside the entire waveguide
  if (rr < 0.0) return 1.0;   // in the core
  while (rr > 1.0) rr -= 1.0; // calculate (r - rcore) % 1
  if (rr < 0.3) return 21.16; // in the high dielectric
  return 1.6*1.6;             // in the low dielectric
}

complex<double> checkers(const vec &v) {
  const double ther = v.r() + 0.001; // Just to avoid roundoff issues.
  if (v.z() > 0.5) return 0.0;
  int z = (int) (v.z()*5.0);
  int r = (int) (ther*5.0);
  int zz = (int) (v.z()*10.0);
  int rr = (int) (ther*10.0);
  if ((r & 1) ^ (z & 1)) return -1.0;
  if ((rr & 1) ^ (zz & 1)) return 1.0;
  return 0.0;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  master_printf("I've got %d processors!\n", count_processors());
  deal_with_ctrl_c();
  master_printf("Running example program!\n");

  double a = 10;
  int m=1;
  double ttot = 100;
  
  mat ma(volcyl(4.0,3.6,a), guided_eps);
  const char *dirname = make_output_directory(__FILE__);
  ma.set_output_directory(dirname);
  //ma.use_pml_right(1.0);
  //ma.use_pml_left(1.0);
  //ma.use_pml_radial(1.0);
  for (m=0;m<1 && !interrupt;m++) {
    char m_str[10];
    snprintf(m_str, 10, "%d", m);
    master_printf("Working on m = %d with a=%lg...\n", m, a);
    fields f(&ma, m);
    f.use_bloch(0.0);
    f.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, vec(0.6, 2.2), 1.0);
    f.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, vec(0.6, 3.2), 1.0);
    f.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, vec(0.6, 0.2), 1.0);
    //f.initialize_field(Ep, checkers);

    double next_print = 0.0;
    while (f.time() < ttot && !interrupt) {
      if (f.time() >= next_print) {
        master_printf("%d is Working on time %lg...  ", my_rank(), f.time());
        f.eps_slices(m_str);
        master_printf("energy is %lg\n", f.field_energy());
        next_print += 10.0;
      }
      f.step();
      //f.step_right();
    }
  }
  delete[] dirname;
}
