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

int rad;

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
  if (v.z() > 0.5) return 0.0;
  int z = (int) (v.z()*5.0);
  int r = (int) (v.r()*5.0);
  int zz = (int) (v.z()*10.0);
  int rr = (int) (v.r()*10.0);
  if ((r & 1) ^ (z & 1)) return -1.0;
  if ((rr & 1) ^ (zz & 1)) return 1.0;
  return 0.0;
}

int main(int argc, char **argv) {
  deal_with_ctrl_c();
  printf("Running example program!\n");

  rad = 10;
  int m=1;
  double ttot = 100;
  
  mat ma(volcyl(3.0,3.5,rad), guided_eps);
  const char *dirname = make_output_directory(argv[0]);
  ma.set_output_directory(dirname);
  //ma.use_pml_right(1.0);
  //ma.use_pml_left(1.0);
  for (m=0;m<4 && !interrupt;m++) {
    char m_str[10];
    snprintf(m_str, 10, "%d", m);
    printf("Working on m = %d with a=%d...\n", m, rad);
    fields f(&ma, m);
    f.use_bloch(0.0);
    f.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, vec(0.6, 1.2), 1.0);
    //f.initialize_field(Ep, checkers);

    double next_print = 0.0;
    while (f.time() < ttot) {
      if (f.time() >= next_print) {
        printf("Working on time %lg...  ", f.time());
        f.eps_slices(m_str);
        printf("energy is %lg\n", f.total_energy());
        next_print += 10.0;
      }
      f.step();
      //f.step_right();
    }
  }
  delete[] dirname;
}

