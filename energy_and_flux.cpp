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
#include <math.h>

#include "dactyl.h"
#include "dactyl_internals.h"

// Flux stuff

double fields::zflux(int ri, int ro, int z) {
  double flux = 0;
  double rph, rtemp;
 
  if (ri > ro) SWAP(ri,ro)
  DOCMP
    for (int r=ri;r<=ro;r++) {
      rph   = ((double)r) + 0.5;
      flux += rph*CM(er,r,z)*CM(hp,r,z)-((double)r)*CM(ep,r,z)*CM(hr,r,z);
    }
  return flux;
}

double fields::rflux(int zl, int zu, int r) {
  double flux = 0;
  double rph, rtemp;

  if (zl > zu) SWAP(zl,zu)
  DOCMP
    for (int z=zl;z<=zu;z++)
      flux += CM(ep,r,z)*CM(hz,r,z) - CM(hp,r,z) * CM(ez,r,z);
  return sqrt(((double)r)*(((double)r)+0.5))*flux;
}

flux_plane::flux_plane(double the_ymin, double the_ymax, double the_xconst, int the_is_rflux, double a)  {
  verbosity = 1;
  if (the_ymin <= the_ymax) {
    ymin = the_ymin*a;
    ymax = the_ymax*a;
  }
  else {
    printf("Switching ymin and ymax. Should we change sign of flux?\n");
    ymin = the_ymax*a;
    ymax = the_ymin*a;
  }
  xconst = the_xconst*a;
  is_rflux = the_is_rflux;
  if (verbosity)  
    printf("Creating flux plane with ymin=%f\tymax=%f\txconst=%f\n", ymin, ymax, xconst);

  int i_xconst = (int)floor(xconst); 
  xpos[0] = i_xconst;
  xpos[1] = i_xconst+1;
  weights[0] = i_xconst + 1 - xconst;
  weights[1] = xconst - i_xconst;
  if (weights[1] == 0.0)
    num_wf = 1;
  else 
    num_wf = 2;
  
  if (verbosity) {
    printf("Going to use %d flux_plane(s)\n", num_wf);
    for (int j=0; j<num_wf; j++)
      printf("%d at %d\t with weight %f\n", j, xpos[j], weights[j]);
  }
  
  int i_ymin = (int)floor(ymin);
  int i_ymax = (int)ceil(ymax);
  double dy_min = ymin - i_ymin;
  double dy_max = i_ymax - ymax;

  if (i_ymin == i_ymax)
    printf("Flux plane too small!!! \n\n\n *****************");
  
  for (int j=0; j<num_wf; j++) {
    wf[j] = new weighted_flux_plane(i_ymin, i_ymax, xpos[j], dy_min, dy_max, is_rflux);
    wf[j]->verbosity = verbosity;
  }
}

flux_plane::flux_plane(const flux_plane &fp) {
  verbosity = fp.verbosity;
  ymin = fp.ymin;
  ymax = fp.ymax;
  xconst = fp.xconst;
  is_rflux = fp.is_rflux;
  num_wf = fp.num_wf;
  for (int j=0; j<num_wf; j++) {
    weights[j] = fp.weights[j];
    xpos[j] = fp.xpos[j];
    wf[j] = new weighted_flux_plane( * fp.wf[j] );
  }
}

flux_plane::~flux_plane() {
  for (int j=0; j<num_wf; j++)
    delete wf[j];
}

weighted_flux_plane::weighted_flux_plane(int the_ymin, int the_ymax, int the_xconst, 
	    double the_dy_min, double the_dy_max, int the_is_rflux) {
  ymin = the_ymin;
  ymax = the_ymax;
  xconst = the_xconst;
  dy_min = the_dy_min;
  dy_max = the_dy_max;
  is_rflux = the_is_rflux;
  if (verbosity)
    printf("Creating weighted flux plane with ymin=%d  ymax=%d  xc=%d\t%.3f %.3f\n", ymin, ymax, xconst, dy_min, dy_max);
}

complex<double> flux_plane::flux(fields *f) {
  complex<double> fl = 0.0;
  for (int j=0;j<num_wf;j++)
    fl += weights[j] * wf[j]->flux(f);
  return fl;
}

complex<double> weighted_flux_plane::flux(fields *f) {
  complex<double> fl = 0.0;
  int nz = f->nz;
  int y;
  if (is_rflux) {
    DOCMP { // + ep x hz 
      y = ymin;
      fl += max(0.0, 0.5-dy_min) *
	CM(f->ep, xconst, y)*0.5*(CM(f->hz, xconst-1, y)+CM(f->hz, xconst, y));
      y = ymax;
      fl += max(0.0, 0.5-dy_max) *
	CM(f->ep, xconst, y)*0.5*(CM(f->hz, xconst-1, y)+CM(f->hz, xconst, y));
      for (y=ymin+1; y<ymax; y++) {
	fl += CM(f->ep, xconst, y)*0.5*(CM(f->hz, xconst-1, y)+CM(f->hz, xconst, y));
      }
      for (y=ymin; y<ymax; y++) { // - ez x hp 
	fl += -CM(f->ez, xconst, y)*0.5*(CM(f->hp, xconst-1, y)+CM(f->hp, xconst, y));
      }
    
      y = ymin+1;
      fl -= max(0.0,dy_min-0.5)*CM(f->ep, xconst, y)*0.5*(CM(f->hz, xconst-1, y)+CM(f->hz, xconst, y));
      y = ymax-1;
      fl -= max(0.0,dy_max-0.5)*CM(f->ep, xconst, y)*0.5*(CM(f->hz, xconst-1, y)+CM(f->hz, xconst, y));
    
      y = ymin;
      fl -= -dy_min * CM(f->ez, xconst, y)*0.5*(CM(f->hp, xconst-1, y)+CM(f->hp, xconst, y));
      y = ymax-1;
      fl -= -dy_max * CM(f->ez, xconst, y)*0.5*(CM(f->hp, xconst-1, y)+CM(f->hp, xconst, y));
    }
    fl = fl * xconst;
  }
  else {
    DOCMP { // - ep x hr 
      y = ymin;
      fl += max(0.0, 0.5-dy_min) * y *
	(-1)*CM(f->ep, y, xconst)*0.5*(CM(f->hr, y, xconst-1)+CM(f->hr, y, xconst));
      y = ymax;
      fl += max(0.0, 0.5-dy_max) * y *
	(-1)*CM(f->ep, y, xconst)*0.5*(CM(f->hr, y, xconst-1)+CM(f->hr, y, xconst));
      for (y=ymin+1; y<ymax; y++) {
	fl += y*(-1)*CM(f->ep, y, xconst)*0.5*(CM(f->hr, y, xconst-1)+CM(f->hr, y, xconst));
      }
      for (y=ymin; y<ymax; y++) { // + er x hp 
	fl += (y+0.5)*CM(f->er, y, xconst)*0.5*(CM(f->hp, y, xconst-1)+CM(f->hp, y, xconst));
      }
          
      y = ymin+1;
      fl -= max(0.0,dy_min-0.5)*y*(-1)*CM(f->ep, y, xconst)*0.5*(CM(f->hr, y, xconst-1)+CM(f->hr, y, xconst));
      y = ymax-1;
      fl -= max(0.0,dy_max-0.5)*y*(-1)*CM(f->ep, y, xconst)*0.5*(CM(f->hr, y, xconst-1)+CM(f->hr, y, xconst));
    
      y = ymin;
      fl -= dy_min * (y+0.5)*CM(f->er, y, xconst)*0.5*(CM(f->hp, y, xconst-1)+CM(f->hp, y, xconst));
      y = ymax-1;
      fl -= dy_max * (y+0.5)*CM(f->er, y, xconst)*0.5*(CM(f->hp, y, xconst-1)+CM(f->hp, y, xconst));
    }
  }
  return fl;
}

flux_plane fields::create_rflux_plane(double zmin, double zmax, double rconst) {
  return flux_plane(zmin, zmax, rconst, 1, a);
}

flux_plane fields::create_zflux_plane(double rmin, double rmax, double zconst) {
  return flux_plane(rmin, rmax, zconst, 0, a);
}

complex<double> fields::get_flux(flux_plane *fp) {
  return fp->flux(this);
}


/* Energy calculation */

double fields::total_energy() {
  
  if (backup_hr[0] == NULL) {
    DOCMP {
      backup_hr[cmp] = new double[(nr+1)*(nz+1)];
      backup_hp[cmp] = new double[(nr+1)*(nz+1)];
      backup_hz[cmp] = new double[(nr+1)*(nz+1)];
    }
  }

  DOCMP {
    for (int j=0; j<(nr+1)*(nz+1); j++)
      backup_hr[cmp][j] = hr[cmp][j];
    for (int j=0; j<(nr+1)*(nz+1); j++)
      backup_hp[cmp][j] = hp[cmp][j];
    for (int j=0; j<(nr+1)*(nz+1); j++)
      backup_hz[cmp][j] = hz[cmp][j];
  }

  step_h_bulk();
  step_h_boundaries();
  double next_step_magnetic_energy = magnetic_energy();

  DOCMP {
    for (int j=0; j<(nr+1)*(nz+1); j++)
      hr[cmp][j] = backup_hr[cmp][j];
    for (int j=0; j<(nr+1)*(nz+1); j++)
      hp[cmp][j] = backup_hp[cmp][j];
    for (int j=0; j<(nr+1)*(nz+1); j++)
      hz[cmp][j] = backup_hz[cmp][j];
  }

 return electric_energy() + 0.5*next_step_magnetic_energy + 0.5*magnetic_energy();
}

double fields::electric_energy() {
  double energy = 0, norm = 0;
  DOCMP {
    for (int r=0;r<nr;r++) {
      double rph = r+0.5;
      for (int z=0;z<nz;z++) {
        energy += rph*(1./MA(ma->invepser,r,z))*CM(er,r,z)*CM(er,r,z);
        energy += r*(1./MA(ma->invepsep,r,z))*CM(ep,r,z)*CM(ep,r,z);
        energy += r*(1./MA(ma->invepsez,r,z))*CM(ez,r,z)*CM(ez,r,z);
        norm += (r + rph)/2;
      }
    }
  }
  return energy/norm/(8*pi);
}

double fields::magnetic_energy() {
  double energy = 0, norm = 0;
  DOCMP {
    for (int r=0;r<nr;r++) {
      double rph = r+0.5;
      for (int z=0;z<nz;z++) {
        energy += r*CM(hr,r,z)*CM(hr,r,z);
        energy += rph*CM(hp,r,z)*CM(hp,r,z);
        energy += rph*CM(hz,r,z)*CM(hz,r,z);
        norm += (r + rph)/2;
      }
    }
  }
  return energy/norm/(8*pi);
}

