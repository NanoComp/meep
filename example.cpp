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
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
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

double source_sharp2(double r) {
  if (r == 0) return 0;
  double dr = r - 0.6;
  double sig = 0.15;
  return exp(-dr*dr/(sig*sig));
}

double source_sharp(double r) {
  if (r == 0) return 0;
  double dr = r - rsource;
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
  char directory[100], target[100], sumwflux_name[100], temp[100];
  char executable_name[100], source_name[100], target_name[100];
  char bands_name[100];
  char *lastptr, *currptr;
  DIR *dir;
  FILE *tempf, *targetf;
  int myc;

  signal(SIGINT, handle_control_c);
  printf("Running example program!\n");

  rad = 10;
  int m=1;
  double k = 0.0;
  int ttot = 8000*rad;

  sprintf(directory,"example-dir/");
  if ((dir = opendir(directory)) != NULL)
    closedir(dir);
  else
    mkdir(directory, 00777);
  sprintf(bands_name, "%s/bands", directory);
  FILE *ban = fopen(bands_name, "w");
  strcpy(executable_name, argv[0]);
  lastptr = strtok(executable_name, "/");
  do {
    currptr = strtok(NULL, "/"); 
    if (currptr != NULL)
      lastptr = currptr;
  } while (currptr != NULL);
  sprintf(source_name, "%s.cpp", lastptr);
  sprintf(target_name, "%s/%s", directory, source_name);
  if ((tempf = fopen(source_name,"r")) != NULL 
    && (targetf = fopen(target_name,"w")) != NULL) {
    while ( (myc=getc(tempf)) != EOF) fprintf(targetf, "%c", myc);
    fclose(tempf);
    fclose(targetf);
  } else
    printf("Error, couldn't open source file %s and/or target file %s.\n",
	   source_name, target_name);
  mat ma(guided_eps, 1.0, 0.0, rad);
  //ma.use_pml(8,8);
  ma.output_slices("example");
  for (m=1;m<2 && !stopnow;m++) {
    for (k=0.3;k<0.31 && !stopnow;k+=0.1) {
      strcpy(target,directory);
      sprintf(temp,"m%d",m);
      strcat(target,temp);
      strcpy(sumwflux_name,directory);
      sprintf(temp, "sumwflux-m%d.dat", m);
      strcat(sumwflux_name,temp);
      printf("Working on k of %g and m = %d with a=%d...\n", k, m, rad);
      fields f(&ma, m);
      f.use_bloch(k);
      //f.ep[0][14] = 1;
      add_clever_sources(f, 0.30, 2.5, 0.15);
      add_clever_sources(f, 0.31, 2.5, 0.6);
      add_clever_sources(f, 0.29, 2.5, 0.45);
      f.prepare_for_bands(0, ttot, 2.5, 20);
//f.add_ep_source(0.11 , 0.02,  0.0, 8.0, 0.0*rad, source_sharp);
//      f.add_ep_source(0.12 , 0.02,  0.0, 8.0, 0.0*rad, source_sharp2);
//      f.add_ez_source(0.125, 0.02,  0.0, 8.0, 0.0*rad, source_sharp);
//      f.add_ez_source(0.115, 0.021, 0.0, 8.0, 0.0*rad, source_sharp2);

      //f.set_frequency_range(0.0,0.4, ((double) f.a) / (c * ((double)ttot)) );
      //f.add_zfluxplane(0,(int)(5*rad),80);
      
      for (int t=0; t < ttot+1 && !stopnow; t++) {
        if (t % (1000*rad) == 0 && t > 137*rad) {
          printf("Working on time step %d...  ", t);
          f.output_slices("example");
          printf("energy is %lg\n", f.total_energy()/3.66e-6);
          //printf("flux is %lg\n", f.zflux(0,22,83));
        }
        f.record_bands();
        f.output_point(ban, 0, 0, "point");
        f.step();
        f.dft_flux();
      }
      //FILE *sumwflux = fopen(sumwflux_name, "w");
      //f.fluxw_output(sumwflux, "sumwflux");
      f.output_bands(ban, "band", 15);
    }
  }
  fclose(ban);
}

