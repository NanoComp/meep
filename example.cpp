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

static int stopnow = 0;
void handle_control_c(int i) {
  printf("Be patient, I'll stop as soon as it's convenient...\n");
  stopnow = 1;
}

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

static double rsource = 0.15;

complex<double> source_sharp(double r) {
  if (r == 0) return 0;
  double dr = r - 0.6;
  double sig = 0.15;
  return exp(-dr*dr/(sig*sig));
}

void add_clever_sources(fields &f, double fmin, double fmax, double r) {
  rsource = r;
  const double A = 0.5;
  f.add_ep_source(0.33*fmin, 0.13*fmin/fmax, 0.0, 8.0*(1+0*fmin/fmax*A), 0.0, source_sharp);
  f.add_ez_source(0.33*fmin, 0.13*fmin/fmax, 0.0, 8.0*(1+1*fmin/fmax*A), 0.0, source_sharp);
  f.add_hp_source(0.33*fmin, 0.13*fmin/fmax, 0.0, 8.0*(1+2*fmin/fmax*A), 0.0, source_sharp);
  f.add_hz_source(0.33*fmin, 0.13*fmin/fmax, 0.0, 8.0*(1+3*fmin/fmax*A), 0.0, source_sharp);
  f.use_real_sources();
}

int main(int argc, char **argv) {
  signal(SIGINT, handle_control_c);
  printf("Running example program!\n");

  rad = 10;
  int m=1;
  double k = 0.0;
  int ttot = 4000*rad;
  
  mat ma(guided_eps, 1.0, 0.0, rad);
  const char *dirname = make_output_directory(argv[0]);
  printf("Storing output in directory %s/\n", dirname);
  FILE *ban = create_output_file(dirname, "bands");
  //FILE *fluxf = create_output_file(dirname, "flux");
  ma.set_output_directory(dirname);
  //ma.use_pml(8,8);
  ma.output_slices("");
  for (m=0;m<4 && !stopnow;m++) {
    for (k=0.0;k<5.1 && !stopnow;k+=1.0) {
      char k_and_m[10];
      snprintf(k_and_m, 10, "%g-%d", k, m);
      printf("Working on k of %g and m = %d with a=%d...\n", k, m, rad);
      fields f(&ma, m);
      f.use_bloch(k);
      //f.ep[0][14] = 1;
      add_clever_sources(f, 0.30, 2.5, 0.15);
      add_clever_sources(f, 0.31, 2.5, 0.6);
      add_clever_sources(f, 0.29, 2.5, 0.45);
      f.prepare_for_bands(0, ttot, 2.5, 200);
//      f.add_ez_source(0.125, 0.02,  0.0, 8.0, 0.0*rad, source_sharp);
      //f.set_frequency_range(0.0,0.4, ((double) f.a) / (c * ((double)ttot)) );
      //f.add_zfluxplane(0,(int)(5*rad),80);
      
      for (int t=0; t < ttot+1 && !stopnow; t++) {
        if (t % (1000*rad) == 0 && t > 137*rad) {
          printf("Working on time step %d...  ", t);
          f.output_slices(k_and_m);
          printf("energy is %lg\n", f.total_energy()/3.66e-6);
          //printf("flux is %lg\n", f.zflux(0,22,83));
        }
        f.record_bands();
        //f.output_point(ban, 0, 0, "point");
        f.step();
        f.dft_flux();
      }
      //f.fluxw_output(fluxf, "sumwflux");
      f.output_bands(ban, "band", 20);
    }
  }
  fclose(ban);
}

