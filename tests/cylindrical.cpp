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

double one(const vec &) {
  return 1.0;
}

int test_simple_periodic(int splitting) {
  double a = 10;
  double ttot = 30;
  
  mat ma1(volcyl(2.0,1.5,a), one, 1);
  mat ma(volcyl(2.0,1.5,a), one, splitting);
  for (int m=0;m<4;m++) {
    char m_str[10];
    snprintf(m_str, 10, "%d", m);
    printf("Trying m = %d with a=%lg and a splitting into %d chunks...\n",
           m, a, splitting);
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
  return 0;
}

int main(int argc, char **argv) {
  printf("Testing cylindrical coords under different splittings...\n");

  return test_simple_periodic(3);
}

