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

int rad;

double guided_eps(double r, double z) {
  return 1.0;
  double rcore = 2.14;
  double rr = r - rcore;
  if (rr > 3.001) return 1.0; // outside the entire waveguide
  if (rr < 0.0) return 1.0;   // in the core
  while (rr > 1.0) rr -= 1.0; // calculate (r - rcore) % 1
  if (rr < 0.3) return 21.16; // in the high dielectric
  return 1.6*1.6;             // in the low dielectric
}

double source_sharp2(double r) {
  if (r == 0) return 0;
  double dr = r - 1.2; //.5;
  double sig = 0.15;
  return exp(-dr*dr/(sig*sig));
}

double source_sharp(double r) {
  if (r == 0) return 0;
  double dr = r - 0.3;
  double sig = 0.15;
  return exp(-dr*dr/(sig*sig));
}

int main() {
  printf("Running example program!\n");
  FILE *ban = fopen("bands", "w");

  rad = 10;
  int m=1;
  double k = 0.0;
  mat ma(guided_eps, 3.5, 3.0, rad);
  ma.use_pml(8,8);
  ma.output_slices("example");
  for (m=0;m<1;m++) {
    for (k=0.3;k<0.31;k+=0.1) {
      printf("Working on k of %g and m = %d with a=%d...\n", k, m, rad);
      fields f(&ma, m);
      //f.use_bloch(k);
      //f.ep[0][14] = 1;
      //f.add_er_source(0.35 , 0.3, 0.0, 3.0, 0.0*rad, source_sharp);
      //f.add_er_source(0.38 , 0.3, 0.0, 3.0, 0.0*rad, source_sharp2);
      //f.add_ep_source(0.13 , 0.5, 0.0, 5.0, 2.0*rad, source_sharp);
      //f.add_ep_source(0.12 , 0.5, 0.0, 5.0, 2.0*rad, source_sharp2);
      //f.add_ez_source(0.125, 0.5, 0.0, 5.0, 2.0*rad, source_sharp);
      f.add_ez_source(0.135, 0.51, 0.0, 5.0, 2.0*rad, source_sharp2);

      //f.set_frequency_range(0.0,0.4, ((double) f.a) / (c * ((double)ttot)) );
      //f.add_zfluxplane(0,(int)(5*rad),80);

      int ttot = 3000*rad;
      ttot = 140*rad;
      
      //f.prepare_for_bands(0, ttot, 3.0, 20);
      for (int t=0; t < ttot+1; t++) {
        if (t % /*(100*rad)*/1 == 0 && t > 137*rad) {
          printf("Working on time step %d...  ", t);
          f.output_slices("example");
          printf("energy is %lg\n", f.total_energy()/3.66e-6);
          //printf("flux is %lg\n", f.zflux(0,22,83));
        }
        //f.record_bands();
        f.step();
        f.dft_flux();
      }
      //f.fluxw_output(sumwflux, "sumwflux");
      //f.output_bands(ban, "band", 10);
    }
  }
  fclose(ban);
}

