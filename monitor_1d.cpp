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

#include "tidod.h"
#include "tidod_internals.h"

/* Below are the monitor point routines. */

monitor_point_1d::monitor_point_1d() {
  next = NULL;
}

monitor_point_1d::~monitor_point_1d() {
  if (next) delete next;
}

void fields_1d::output_point(FILE *o, double z, const char *name) {
  monitor_point_1d tmp;
  get_point(&tmp, z);
  fprintf(o, "%s\t%8lg", name, t*inva);
  fprintf(o, "\t%8lg\t%8lg", real(tmp.ex), imag(tmp.ex));
  fprintf(o, "\t%8lg\t%8lg", real(tmp.hy), imag(tmp.hy));
  fprintf(o, "\n");
}

static inline complex<double> interpolate(double *f[2], double z, int nz,
                                          complex<double> eiknz) {
  int zlo = (int) z, zhi = zlo+1;
  double dzlo = z - zlo, dzhi = zhi - zlo;
  complex<double> phzlo = 1.0, phzhi = 1.0;
  if (zlo < 0) {
    zlo = nz-1;
    phzlo = 1.0/eiknz;
  }
  if (zhi > nz-1) {
    zhi = 0;
    phzhi = eiknz;
  }
  if (zlo > nz-1 || zhi < 0) {
    printf("interpolated point is out of range! %lg\n", z);
    exit(1);
  }
  complex<double> flo = phzlo*complex<double>(RE(f,zlo),IM(f,zlo));
  complex<double> fhi = phzhi*complex<double>(RE(f,zhi),IM(f,zhi));
  return dzlo*flo + dzhi*fhi;
}

void fields_1d::get_point(monitor_point_1d *pt, double z) {
  if (pt == NULL) {
    printf("Error:  get_point passed a null pointer!\n");
    exit(1);
  }
  pt->z = z;
  pt->t = t*inva*c;
  pt->ex = interpolate(ex, z*a    ,nz,eiknz);
  pt->hy = interpolate(ex, z*a-0.5,nz,eiknz);
}

monitor_point_1d *fields_1d::get_new_point(double z, monitor_point_1d *the_list) {
  monitor_point_1d *p = new monitor_point_1d();
  get_point(p, z);
  p->next = the_list;
  return p;
}

#include <fftw.h>

void monitor_point_1d::fourier_transform(component_1d w,
                                         complex<double> **a, complex<double> **f,
                                         int *numout, double fmin, double fmax,
                                         int maxbands) {
  int n = 1;
  monitor_point_1d *p = next;
  double tmax = t, tmin = t;
  while (p) {
    n++;
    if (p->t > tmax) tmax = p->t;
    if (p->t < tmin) tmin = p->t;
    p = p->next;
  }
  p = this;
  complex<double> *d = new complex<double>[n];
  for (int i=0;i<n;i++,p=p->next) {
    d[i] = p->get_component(w);
  }
  if (fmin > 0.0) { // Get rid of any static fields_chunk!
    complex<double> mean = 0.0;
    for (int i=0;i<n;i++) mean += d[i];
    mean /= n;
    for (int i=0;i<n;i++) d[i] -= mean;
  }
  if ((fmin > 0.0 || fmax > 0.0) && maxbands > 0) {
    *a = new complex<double>[maxbands];
    *f = new complex<double>[maxbands];
    *numout = maxbands;
    delete[] d;
    for (int i = 0;i<maxbands;i++) {
      double df = (maxbands == 1) ? 0.0 : (fmax-fmin)/(maxbands-1);
      (*f)[i] = fmin + i*df;
      (*a)[i] = 0.0;
      p = this;
      while (p) {
        double inside = 2*pi*real((*f)[i])*p->t;
        (*a)[i] += p->get_component(w)*complex<double>(cos(inside),sin(inside));
        p = p->next;
      }
      (*a)[i] /= (tmax-tmin);
    }
  } else {
    *numout = n;
    *a = new complex<double>[n];
    *f = d;
    fftw_complex *in = (fftw_complex *) d, *out = (fftw_complex *) *a;
    fftw_plan p;
    p = fftw_create_plan(n, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_one(p, in, out);
    fftw_destroy_plan(p);
    for (int i=0;i<n;i++) {
      (*f)[i] = i*(1.0/(tmax-tmin));
      if (real((*f)[i]) > 0.5*n/(tmax-tmin)) (*f)[i] -= n/(tmax-tmin);
      (*a)[i] *= (tmax-tmin)/n;
    }
  }
}

void monitor_point_1d::harminv(component_1d w,
                               complex<double> **a, double **f_re, double **f_im,
                               int *numout, double fmin, double fmax,
                               int maxbands) {
  int n = 1;
  monitor_point_1d *p = next;
  double tmax = t, tmin = t;
  while (p) {
    n++;
    if (p->t > tmax) tmax = p->t;
    if (p->t < tmin) tmin = p->t;
    p = p->next;
  }
  p = this;
  complex<double> *d = new complex<double>[n];
  for (int i=0;i<n;i++,p=p->next) {
    d[i] = p->get_component(w);
  }
  *a = new complex<double>[n];
  *f_re = new double[n];
  *f_im = new double[n];
  printf("using an a of %lg\n", (n-1)/(tmax-tmin)*c);
  *numout = do_harminv(d, n, 1, (n-1)/(tmax-tmin)*c, fmin, fmax, maxbands,
                       *a, *f_re, *f_im, NULL);
  delete[] d;
}
