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
#include <complex>

#include "dactyl.h"
#include "dactyl_internals.h"

// Cylindrical coordinates:

#include <gsl/gsl_sf_bessel.h>

static double J(int m, double kr) { return gsl_sf_bessel_Jn(m, kr); }
static double Jprime(int m, double kr) { 
  if (m) return 0.5*(J(m-1,kr)-J(m+1,kr));
  else return -J(1,kr);
}
static double Jroot(int m, int n) { return gsl_sf_bessel_zero_Jnu(m, n+1); }
static double Jmax(int m, int n) {
  double rlow, rhigh = Jroot(m,n), rtry;
  if (n == 0) rlow = 0;
  else rlow = Jroot(m, n-1);
  double jplow = Jprime(m,rlow), jptry;
  do {
    rtry = rlow + (rhigh - rlow)*0.5;
    jptry = Jprime(m,rtry);
    if (jplow*jptry < 0) rhigh = rtry;
    else rlow = rtry;
  } while (rhigh - rlow > rhigh*1e-15);
  return rtry;
}

static double ktrans, kax;
static int m_for_J;
static complex<double> JJ(const vec &v) {
  return polar(J(m_for_J, ktrans*v.r()),kax*v.r());
}
static complex<double> JP(const vec &v) {
  return polar(Jprime(m_for_J, ktrans*v.r()),kax*v.r());
}

void fields::initialize_with_nth_te(int np0) {
  if (v.dim == dcyl) {
    const int n = (m==0) ? np0 - 0 : np0 - 1;
    const double rmax = Jmax(m,n);
    ktrans = rmax/(v.nr()*inva);
    kax = k*2*pi*inva;
    m_for_J = m;
    initialize_field(Hz, JJ);
  } else {
    printf("Can't initialize with TE in this dimension.\n");
  }
}

void fields::initialize_with_nth_tm(int np1) {
  if (v.dim == dcyl) {
    const int n = np1 - 1;
    const double rroot = Jroot(m,n);
    ktrans = rroot/(v.nr()*inva);
    kax = k*2*pi*inva;
    m_for_J = m;
    initialize_field(Ez, JJ);
    initialize_field(Hp, JP);
  } else {
    printf("Can't initialize with TM in this dimension.\n");
  }
}

void fields::initialize_with_n_te(int ntot) {
  for (int n=0;n<ntot;n++) initialize_with_nth_te(n+1);
}

void fields::initialize_with_n_tm(int ntot) {
  for (int n=0;n<ntot;n++) initialize_with_nth_tm(n+1);
}

void fields::initialize_field(component c, complex<double> func(const vec &)) {
  for (int i=0;i<v.ntot();i++) {
    complex<double> val = func(v.loc(c,i));
    f[c][0][i] += real(val);
    f[c][1][i] += imag(val);
  }
  step_h_boundaries();
  step_e_boundaries();
}
