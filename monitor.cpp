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

#include "meep.h"
#include "meep_internals.h"

/* Below are the monitor point routines. */

monitor_point::monitor_point() {
  next = NULL;
}

monitor_point::~monitor_point() {
  if (next) delete next;
}

inline complex<double> getcm(const double * const f[2], int i) {
  return complex<double>(f[0][i],f[1][i]);
}

static void dumbsort(complex<double> val[8]) {
  for (int i=0;i<7;i++) {
    int lowest = i;
    for (int j=i+1;j<8;j++) if (abs(val[j]) < abs(val[lowest])) lowest = j;
    complex<double> tmp = val[i];
    val[i] = val[lowest];
    val[lowest] = tmp;
  }
}

static void dumbsort(double val[8]) {
  for (int i=0;i<7;i++) {
    int lowest = i;
    for (int j=i+1;j<8;j++) if (abs(val[j]) < abs(val[lowest])) lowest = j;
    double tmp = val[i];
    val[i] = val[lowest];
    val[lowest] = tmp;
  }
}

void fields::get_point(monitor_point *pt, const vec &loc) const {
  if (pt == NULL) abort("Error:  get_point passed a null pointer!\n");
  for (int i=0;i<10;i++) pt->f[i] = 0.0;
  pt->loc = loc;
  pt->t = time();
  FOR_COMPONENTS(c)
    if (v.has_field(c))
      pt->f[c] = get_field(c,loc);
}

complex<double> fields::get_field(component c, const vec &loc) const {
  ivec ilocs[8];
  double w[8];
  complex<double> val[8];
  for (int i=0;i<8;i++) val[i] = 0.0;
  v.interpolate(c, loc, ilocs, w);
  for (int argh=0;argh<8&&w[argh];argh++)
    val[argh] = w[argh]*get_field(c,ilocs[argh]);
  dumbsort(val);
  complex<double> res = 0.0;
  for (int i=0;i<8;i++) res += val[i];
  return res;
}

complex<double> fields::get_field(component c, const ivec &origloc) const {
  ivec iloc = origloc;
  complex<double> kphase = 1.0;
  locate_point_in_user_volume(&iloc, &kphase);
  for (int sn=0;sn<S.multiplicity();sn++)
    for (int i=0;i<num_chunks;i++)
      if (chunks[i]->v.contains(S.transform(iloc,sn)))
        return S.phase_shift(c,sn)*kphase*
          chunks[i]->get_field(S.transform(c,sn),S.transform(iloc,sn));
}

complex<double> fields_chunk::get_field(component c, const ivec &iloc) const {
  complex<double> res = 0.0;
  if (is_mine()) {
    if (f[c][0] && f[c][1]) res = getcm(f[c], v.index(c, iloc));
    else if (f[c][0]) res = f[c][0][v.index(c,iloc)];
  }
  return broadcast(n_proc(), res);
}

double fields::get_eps(const ivec &origloc) const {
  ivec iloc = origloc;
  complex<double> aaack = 1.0;
  locate_point_in_user_volume(&iloc, &aaack);
  for (int sn=0;sn<S.multiplicity();sn++)
    for (int i=0;i<num_chunks;i++)
      if (chunks[i]->v.contains(S.transform(iloc,sn)))
        return chunks[i]->get_eps(S.transform(iloc,sn));
  return 0.0;
}

double fields_chunk::get_eps(const ivec &iloc) const {
  double res = 0.0;
  if (is_mine()) res = ma->eps[v.index(v.eps_component(), iloc)];
  return broadcast(n_proc(), res);
}

double fields::get_eps(const vec &loc) const {
  double theeps = 0.0;
  double val[8];
  for (int i=0;i<8;i++) val[i] = 0.0;
  for (int i=0;i<num_chunks;i++) {
    if (chunks[i]->v.contains(loc))
      chunks[i]->ma->interpolate_eps(loc, val);
      dumbsort(val);
      for (int i=0;i<8;i++) theeps += val[i];
    }
  return theeps;
}

double mat::get_eps(const vec &loc) const {
  double theeps = 0.0;
  double val[8];
  for (int i=0;i<8;i++) val[i] = 0.0;
  for (int i=0;i<num_chunks;i++) {
    if (chunks[i]->v.contains(loc))
      chunks[i]->interpolate_eps(loc, val);
      dumbsort(val);
      for (int i=0;i<8;i++) theeps += val[i];
    }
  return theeps;
}

void mat_chunk::interpolate_eps(const vec &loc, double val[8]) const {
  if (is_mine()) {
    int ind[8];
    double w[8];
    v.interpolate(v.eps_component(),loc,ind,w);
    int startingat = 0;
    for (int i=0;i<8 && val[i]!=0.0;i++) startingat = i+1;
    for (int i=0;i<8 && w[i] && (i+startingat<8);i++) {
      val[i+startingat] = w[i]*eps[ind[i]];
      if (val[i+startingat] == 0.0) startingat--;
    }
  }
  broadcast(n_proc(), val, 8);
}

monitor_point *fields::get_new_point(const vec &loc, monitor_point *the_list) const {
  monitor_point *p = new monitor_point();
  get_point(p, loc);
  p->next = the_list;
  return p;
}

complex<double> monitor_point::get_component(component w) {
  return f[w];
}

#include <fftw.h>

void monitor_point::fourier_transform(component w,
                                      complex<double> **a, complex<double> **f,
                                      int *numout, double fmin, double fmax,
                                      int maxbands) {
  int n = 1;
  monitor_point *p = next;
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

void monitor_point::harminv(component w,
                            complex<double> **a, double **f_re, double **f_im,
                            int *numout, double fmin, double fmax,
                            int maxbands) {
  int n = 1;
  monitor_point *p = next;
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
