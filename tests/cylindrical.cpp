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

int die(const char *msg) {
  puts(msg);
  exit(1);
  return 1;
}

double one(const vec &) {
  return 1.0;
}

int test_simple_periodic(double eps(const vec &), int splitting) {
  double a = 10.0;
  double ttot = 30.0;
  
  volume v = volcyl(6.0,4.2,a);
  mat ma1(v, eps, 1);
  mat ma(v, eps, splitting);
  for (int m=0;m<3;m++) {
    char m_str[10];
    snprintf(m_str, 10, "%d", m);
    printf("Trying with m = %d and a splitting into %d chunks...\n",
           m, splitting);
    fields f(&ma, m);
    f.use_bloch(0.0);
    f.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, vec(0.6, 1.2), 1.0);
    fields f1(&ma, m);
    f1.use_bloch(0.0);
    f1.add_point_source(Ep, 0.7, 2.5, 0.0, 4.0, vec(0.6, 1.2), 1.0);
    while (f.time() < ttot) {
      f.step();
      f1.step();
    }
  }
  return 1;
}

int main(int argc, char **argv) {
  printf("Testing cylindrical coords under different splittings...\n");

  for (int s=2;s<11;s++) {
    test_simple_periodic(one, s) || die("error in test_simple_periodic\n");
  }
}

