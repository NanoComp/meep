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

/* Below are the monitor point routines. */

a_field::a_field() {
  next = NULL;
}

a_field::~a_field() {
  if (next) delete next;
}

a_field::a_field(complex<double> tr, complex<double> tp, complex<double> tz) {
  r = tr;
  p = tp;
  z = tz;
  next = NULL;
}

void a_field::another_field(complex<double> tr, complex<double> tp, complex<double> tz) {
  if (next) next->another_field(tr,tp,tz);
  else next = new a_field(tr,tp,tz);
}

monitor_point::monitor_point() {
  p = NULL;
  next = NULL;
}

monitor_point::~monitor_point() {
  if (p) delete p;
  if (next) delete next;
}

void fields::output_point(FILE *o, double r, double z, const char *name) {
  monitor_point tmp;
  get_point(&tmp, r, z);
  fprintf(o, "%s\t%8lg", name, t*inva);
  fprintf(o, "\t%8lg\t%8lg", real(tmp.e.r), imag(tmp.e.r));
  fprintf(o, "\t%8lg\t%8lg", real(tmp.e.p), imag(tmp.e.p));
  fprintf(o, "\t%8lg\t%8lg", real(tmp.e.z), imag(tmp.e.z));
  fprintf(o, "\t%8lg\t%8lg", real(tmp.h.r), imag(tmp.h.r));
  fprintf(o, "\t%8lg\t%8lg", real(tmp.h.p), imag(tmp.h.p));
  fprintf(o, "\t%8lg\t%8lg", real(tmp.h.z), imag(tmp.h.z));
  fprintf(o, "\n");
}

static inline complex<double> interpolate(double *f[2], double r, double z, 
                                          int nr, int nz,
                                          int odd_with_respect_to_r,
                                          complex<double> eiknz) {
  int rlo = (int) r, zlo = (int) z, rhi = rlo+1, zhi = zlo+1;
  double dr = r - (rlo + 0.5), dz = z - (zlo + 0.5);
  complex<double> phrlo = 1.0, phzlo = 1.0, phzhi = 1.0;
  if (rlo < 0) {
    rlo = 0;
    if (odd_with_respect_to_r) phrlo = -1;
  }
  if (rhi > nr-1) rhi = nr-1;
  if (zlo < 0) {
    zlo = nz-1;
    phzlo = 1.0/eiknz;
  }
  if (zhi > nz-1) {
    zhi = 0;
    phzhi = eiknz;
  }
  if (rlo > nr-1 || zlo > nz-1 || zhi < 0 || rhi < 0) {
    printf("interpolated point is out of range! (%lg,%lg)\n", r, z);
    exit(1);
  }
  complex<double> fsw = phrlo*phzlo*complex<double>(RE(f,rlo,zlo),IM(f,rlo,zlo));
  complex<double> fnw =       phzlo*complex<double>(RE(f,rhi,zlo),IM(f,rhi,zlo));
  complex<double> fse = phrlo*phzhi*complex<double>(RE(f,rlo,zhi),IM(f,rlo,zhi));
  complex<double> fne =       phzhi*complex<double>(RE(f,rhi,zhi),IM(f,rhi,zhi));
  complex<double> fmid = 0.25*(fsw+fse+fnw+fne);
  if (dr+dz > 0.0) {
    if (dr-dz > 0.0) { // North     
      return 2*dr*( (0.5+dz)*fne + (0.5-dz)*fnw ) + (1-2*dr)*fmid;
    } else { // East
      return 2*dz*( (0.5+dr)*fne + (0.5-dr)*fse ) + (1-2*dz)*fmid;
    }
  } else {
    if (dr-dz > 0.0) { // West
      return -2*dz*( (0.5+dr)*fnw + (0.5-dr)*fsw ) + (1+2*dz)*fmid;
    } else { // South
      return -2*dr*( (0.5+dz)*fse + (0.5-dz)*fsw ) + (1+2*dr)*fmid;
    }
  }
  return fmid;
}

void fields::get_point(monitor_point *pt, double r, double z) {
  if (pt == NULL) {
    printf("Error:  get_point passed a null pointer!\n");
    exit(1);
  }
  pt->r = r;
  pt->z = z;
  pt->t = t*inva;
  pt->e.r = interpolate(er,r*a-0.5,z*a    ,nr,nz,m==1,eiknz);
  pt->e.p = interpolate(ep,r*a    ,z*a    ,nr,nz,m==1,eiknz);
  pt->e.z = interpolate(ez,r*a    ,z*a-0.5,nr,nz,0   ,eiknz);
  pt->h.r = interpolate(hr,r*a    ,z*a-0.5,nr,nz,m==1,eiknz);
  pt->h.p = interpolate(hp,r*a-0.5,z*a-0.5,nr,nz,m==1,eiknz);
  pt->h.z = interpolate(hz,r*a-0.5,z*a    ,nr,nz,0   ,eiknz);
  polarization *p = pol;
  if (p) {
    pt->p = new a_field(interpolate(p->Pr,r*a-0.5,z*a    ,nr,nz,m==1,eiknz),
                        interpolate(p->Pp,r*a    ,z*a    ,nr,nz,m==1,eiknz),
                        interpolate(p->Pz,r*a    ,z*a-0.5,nr,nz,0   ,eiknz));
    p = p->next;
  }
  while (p) {
    pt->p->another_field(interpolate(p->Pr,r*a-0.5,z*a    ,nr,nz,m==1,eiknz),
                         interpolate(p->Pp,r*a    ,z*a    ,nr,nz,m==1,eiknz),
                         interpolate(p->Pz,r*a    ,z*a-0.5,nr,nz,0   ,eiknz));
    p = p->next;
  }
}

monitor_point *fields::get_new_point(double r, double z, monitor_point *the_list) {
  monitor_point *p = new monitor_point();
  get_point(p, r, z);
  p->next = the_list;
  return p;
}
