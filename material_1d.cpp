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

mat_1d::~mat_1d() {
  delete[] inveps;
  delete[] eps;

  delete[] Czex;
  delete[] Czhy;
  if (pb) delete pb;
}

static double sig(double r, double power);

double mat_1d::use_pml(double pmlz, double fmin) {
  return use_integer_pml((int)(0.5+pmlz*a), fmin);
}

static double minimize_badness(double sig[], int thickness, double eps, double fmin, int i);
inline void reverse(double sig[], int l) {
  for (int i=0;i<l/2;i++) {
    double temp = sig[i];
    sig[i] = sig[l-1-i];
    sig[l-1-i] = temp;
  }
}

double mat_1d::use_integer_pml(int numpmlz, double fmin) {
  pml_fmin = fmin;
  npmlz = numpmlz;
  if (pb) pb->use_integer_pml(npmlz);
  if (npmlz * 2 >= nz) {
    printf("Not enough room for the z PML. nz = %d\n", nz);
    exit(1);
  }
    
  // Delete any previously allocated conductivity arrays...
  delete[] Czex;
  delete[] Czhy;
  // Allocate the conductivity arrays:
  Czex = new double[npmlz];
  Czhy = new double[npmlz];
  // Initialize them... (for now I'm setting them to zero...)
  const double Cmax = 0.5;
  double meaneps = 0;
  for (int i=0;i<nz;i++) {
    meaneps += eps[i]; // really I should use the min eps within the pml area...
  }
  meaneps /= nz;
  double oldrefl = 1e6;
  double reflection = 0;
  for (int z=0;z<npmlz;z++) {
    double rr = (z)/(double)npmlz;
    double rp = (1+z)/(double)npmlz;
    Czhy[z] = Czex[z] = Cmax*sig(rp, 2.0);
  }
  reverse(Czhy, npmlz);
  reflection = 0;
  for (int i=0;i<npmlz*20;i++) {
    for (int z=0;z<npmlz;z++) {
      reflection = minimize_badness(Czhy, npmlz, meaneps, fmin/a*c, z);
      reflection = minimize_badness(Czhy, npmlz, meaneps, fmin/a*c, z);
      reflection = minimize_badness(Czhy, npmlz, meaneps, fmin/a*c, z);
      if (oldrefl == reflection) break;
      oldrefl = reflection;
    }
  }
  //reverse(Czhy, npmlz);
  for (int z=0;z<npmlz;z++) Czex[z] = Czhy[z];
  return reflection;
}

static double badness(double sig[], int thickness, double epsilon, double fmin) {
  if (thickness < 1) return 1;
  const double A = .0001/fmin*.1/fmin, K = 6.0/epsilon*2.25/epsilon;
  double sofar = 1.0;
  for (int i=0;i<thickness-1;i++) {
    double first_trans = exp(-K*sig[i+1]);
    double refl = A*fabs(sig[i]-sig[i+1])*fabs(sig[i]-sig[i+1]);
    double total_trans = exp(-K*sig[i])*first_trans;
    sofar = refl + (1-refl)*total_trans*sofar;
    if (sofar > 1.0) sofar = 1.0;
  }
  double last_refl = A*fabs(sig[thickness-1]);
  sofar = last_refl + (1-last_refl)*sofar;
  return sofar;
}

static double minimize_badness(double sig[], int thickness,
                               double epsilon, double fmin, int i) {
  double behind_reflection = badness(sig, i-1, epsilon, fmin);
  

  double now = badness(sig, thickness, epsilon, fmin);
  double tried = now;
  do {
    now = tried;
    sig[i] *= 1.001;
    tried = badness(sig, thickness, epsilon, fmin);
  } while (tried < now);
  sig[i] /= 1.001;
  tried = now = badness(sig, thickness, epsilon, fmin);
  do {
    now = tried;
    sig[i] /= 1.001;
    tried = badness(sig, thickness, epsilon, fmin);
  } while (tried < now);
  sig[i] *= 1.001;
  return badness(sig, thickness, epsilon, fmin);
}

static double sig(double r, double power) {
  return pow(r, power);
}

void mat_1d::mix_with(const mat_1d *n, double f) {
  for (int i=0;i<nz+1;i++) {
    eps[i] = 1.0/(1.0/eps[i] + f*(1.0/n->eps[i]-1.0/eps[i]));
    inveps[i] += f*(n->inveps[i] - inveps[i]);
  }
  // Mix in the polarizability...
  polarizability_1d *po = pb, *pn = n->pb;
  while (po && pn) {
    for (int i=0;i<nz+1;i++) {
      po->sigma[i] += f*(pn->sigma[i] - po->sigma[i]);
    }
    po = po->next;
    pn = pn->next;
  }
}

void mat_1d::make_average_eps() {
  double meaneps = 0;
  for (int i=0;i<nz;i++) {
    meaneps += eps[i];
  }
  meaneps /= nz;
  for (int i=0;i<nz+1;i++) {
    eps[i] = meaneps;
    inveps[i] = 1/meaneps;
  }
}

mat_1d::mat_1d(const mat_1d *o) {
  outdir = o->outdir;
  if (o->pb) pb = new polarizability_1d(o->pb);
  else pb = NULL;
  a = o->a;
  nz = o->nz;
  eps = new double[nz+1];
  if (eps == NULL) {
    printf("Out of memory!\n");
    exit(1);
  }
  for (int z=0;z<nz+1;z++) {
    MA(eps,z) = MA(o->eps,z);
  }
  inveps = new double[nz+1];
  for (int i=0;i<nz+1;i++) inveps[i] = o->inveps[i];
  // Allocate the conductivity arrays:
  Czex = Czhy = NULL;
  npmlz = 0;
  use_integer_pml(o->npmlz, o->pml_fmin);
}

mat_1d::mat_1d(double feps(double z),double zmax, double ta) {
  pml_fmin = 0.2;
  outdir = ".";
  pb = NULL;
  a = ta;
  nz = (int) (zmax*a);
  npmlz = 0;
  if (nz == 0) nz = 1;
  eps = new double[nz+1];
  if (eps == NULL) {
    printf("Out of memory!\n");
    exit(1);
  }
  for (int z=0;z<nz+1;z++)
    MA(eps,z) = (feps) ? feps(z/a) : 1.0; // Null feps means vacuum.
  MA(eps,nz) = MA(eps,0);
  inveps = new double[nz+1];
  for (int z=0;z<nz+1;z++) inveps[z] = 1/eps[z];
  // Allocate the conductivity arrays:
  Czex = Czhy = NULL;
}
