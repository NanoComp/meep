/* Copyright (C) 2003 David Roundy  
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

#include "dactyl.h"

double eps(double r, double z) {
  return 1.0;
}

int main() {
  const char *name = "d_pml_m_1";
  printf("Running %s!\n", name);

  const int a=1;
  const int m=1;
  mat ma(eps, 10.0, 10.0, a);
  ma.use_pml(2,2);
  ma.output_slices(name);
  fields f(&ma, m);

  f.ep[0][3+(f.nz+1)*3] = 1;
  f.ez[0][3+(f.nz+1)*3] = 0.7;

  f.output_real_imaginary_slices(name);
  for (int t=0; t < 40; t++) {
    f.step();
    f.output_real_imaginary_slices(name);
  }
  printf("t is %d\n", f.t);
}

