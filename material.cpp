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

mat::~mat() {
  delete[] invepser;
  delete[] invepsep;
  delete[] invepsez;
  delete[] eps;

  delete[] Crez;
  delete[] Crep;
  delete[] Crhz;
  delete[] Crhp;
  delete[] Cper;
  delete[] Cpez;
  delete[] Cphr;
  delete[] Cphz;

  delete[] Czer;
  delete[] Czep;
  delete[] Czhr;
  delete[] Czhp;
  if (pb) delete pb;
}

static double sig(double r, double power);

double mat::use_pml(double pmlr, double pmlz, double fmin) {
  return use_integer_pml((int)(0.5+pmlr*a), (int)(0.5+pmlz*a), fmin);
}

static double minimize_badness(double sig[], int thickness, double eps, double fmin, int i);
inline void reverse(double sig[], int l) {
  for (int i=0;i<l/2;i++) {
    double temp = sig[i];
    sig[i] = sig[l-1-i];
    sig[l-1-i] = temp;
  }
}

double mat::use_integer_pml(int numpmlr, int numpmlz, double fmin) {
  pml_fmin = fmin;
  npmlz = numpmlz;
  npmlr = numpmlr;
  if (pb) pb->use_integer_pml(npmlr,npmlz);
  if (npmlz * 2 >= nz) {
    printf("Not enough room for the z PML. nz = %d\n", nz);
    exit(1);
  }
  if (npmlr >= nr) {
    printf("Not enough room for the r PML. nr = %d\n", nr);
    exit(1);
  }
    
  // Delete any previously allocated conductivity arrays...
  delete[] Crez;
  delete[] Crep;
  delete[] Crhz;
  delete[] Crhp;
  delete[] Cper;
  delete[] Cpez;
  delete[] Cphr;
  delete[] Cphz;

  delete[] Czer;
  delete[] Czep;
  delete[] Czhr;
  delete[] Czhp;
  // Allocate the conductivity arrays:
  Crez = new double[npmlr];
  Crep = new double[npmlr];
  Crhz = new double[npmlr];
  Crhp = new double[npmlr];
  Cper = new double[npmlr];
  Cpez = new double[npmlr];
  Cphr = new double[npmlr];
  Cphz = new double[npmlr];

  Czer = new double[npmlz];
  Czep = new double[npmlz];
  Czhr = new double[npmlz];
  Czhp = new double[npmlz];
  // Initialize them... (for now I'm setting them to zero...)
  const double Cmax = 0.5;
  for (int r=0;r<npmlr;r++) {
    double rr = (r)/(double)npmlr;
    double rp = (1+r)/(double)npmlr;
    Crhz[r] = Crhp[r] = Cmax*sig(rp, 2.0);
  }
  double sigintegrated = 0.0;
  for (int r=0;r<npmlr;r++) {
    double rr = (r)/(double)npmlr;
    double rp = (1+r)/(double)npmlr;
    sigintegrated += Cmax*0.5*(sig(rp, 2.0)+sig(rr, 2.0));
    Cper[r] = Cphz[r] = 0.0; //sigintegrated/(nr-npmlr+r+0.5);
  }
  double meaneps = 0;
  for (int i=0;i<nr*(nz+1);i++) {
    meaneps += eps[i]; // really I should use the min eps within the pml area...
  }
  meaneps /= nr*(nz+1);
  reverse(Crhz, npmlr);
  double oldrefl = 1e6;
  double reflection = 0;
  for (int i=0;i<npmlr*20;i++) {
    for (int r=0;r<npmlr;r++) {
      reflection = minimize_badness(Crhz, npmlr, meaneps, fmin/a*c, r);
      reflection = minimize_badness(Crhz, npmlr, meaneps, fmin/a*c, r);
      reflection = minimize_badness(Crhz, npmlr, meaneps, fmin/a*c, r);
      if (oldrefl == reflection) break;
      oldrefl = reflection;
    }
  }
  reverse(Crhz, npmlr);
  for (int r=0;r<npmlr;r++) {
    Crhp[r] = Crhz[r];
    if (r==0) Crep[r] = Crez[r] = 0.5*Crhz[r];
    else Crep[r] = Crez[r] = 0.5*(Crhz[r]+Crhz[r-1]);
    //if (r==0) Cphr[r] = Cpez[r] = 0.5*Cphz[r];
    //else Cphr[r] = Cpez[r] = 0.5*(Cphz[r]+Crhz[r-1]);
    Cphr[r] = Cpez[r] = 0.0; // Avoid instability...
  }
  for (int z=0;z<npmlz;z++) {
    double rr = (z)/(double)npmlz;
    double rp = (1+z)/(double)npmlz;
    Czhr[z] = Czhp[z] = Cmax*sig(rp, 2.0);
  }
  reverse(Czhr, npmlz);
  reflection = 0;
  for (int i=0;i<npmlz*20;i++) {
    for (int z=0;z<npmlz;z++) {
      reflection = minimize_badness(Czhr, npmlz, meaneps, fmin/a*c, z);
      reflection = minimize_badness(Czhr, npmlz, meaneps, fmin/a*c, z);
      reflection = minimize_badness(Czhr, npmlz, meaneps, fmin/a*c, z);
      if (oldrefl == reflection) break;
      oldrefl = reflection;
    }
  }
  reverse(Czhr, npmlz);
  for (int z=0;z<npmlz;z++) {
    Czhp[z] = Czhr[z];
    if (z==0) Czer[z] = Czep[z] = Cmax*0.5*Czhr[z];
    else Czer[z] = Czep[z] = Cmax*0.5*(Czhr[z]+Czhr[z-1]);//(sig(rp, 2.0)+sig(rr, 2.0));
  }
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

void mat::mix_with(const mat *n, double f) {
  for (int i=0;i<nr*(nz+1);i++) {
    eps[i] = 1.0/(1.0/eps[i] + f*(1.0/n->eps[i]-1.0/eps[i]));
    invepser[i] += f*(n->invepser[i] - invepser[i]);
    invepsep[i] += f*(n->invepsep[i] - invepsep[i]);
    invepsez[i] += f*(n->invepsez[i] - invepsez[i]);
  }
  // Mix in the polarizability...
  polarizability *po = pb, *pn = n->pb;
  while (po && pn) {
    for (int i=0;i<nr*(nz+1);i++) {
      po->sr[i] += f*(pn->sr[i] - po->sr[i]);
      po->sp[i] += f*(pn->sp[i] - po->sp[i]);
      po->sz[i] += f*(pn->sz[i] - po->sz[i]);
    }
    po = po->next;
    pn = pn->next;
  }
}

void mat::make_average_eps() {
  double meaneps = 0;
  for (int i=0;i<nr*(nz+1);i++) {
    meaneps += eps[i];
  }
  meaneps /= nr*(nz+1);
  for (int i=0;i<nr*(nz+1);i++) {
    eps[i] = meaneps;
    invepser[i] = 1/meaneps;
    invepsep[i] = 1/meaneps;
    invepsez[i] = 1/meaneps;
  }
}

mat::mat(const mat *o) {
  outdir = o->outdir;
  if (o->pb) pb = new polarizability(o->pb);
  else pb = NULL;
  a = o->a;
  nr = o->nr;
  nz = o->nz;
  eps = new double[nr*(nz+1)];
  if (eps == NULL) {
    printf("Out of memory!\n");
    exit(1);
  }
  for (int r=0;r<nr;r++) {
    for (int z=0;z<nz+1;z++) {
      MA(eps,r,z) = MA(o->eps,r,z);
    }
  }
  invepser = new double[(nr+1)*(nz+1)];
  invepsep = new double[(nr+1)*(nz+1)];
  invepsez = new double[(nr+1)*(nz+1)];

  for (int i=0;i<(nr+1)*(nz+1);i++) {
    invepser[i] = o->invepser[i];
    invepsep[i] = o->invepsep[i];
    invepsez[i] = o->invepsez[i];
  }
  // Allocate the conductivity arrays:
  Crez = Crep = Crhz = Crhp = Cper = Cpez = Cphr = Cphz = NULL;
  Czer = Czep = Czhr = Czhp = NULL;
  npmlz = 0;
  npmlr = 0;
  use_integer_pml(o->npmlr,o->npmlz, o->pml_fmin);
}

mat::mat(double feps(double r, double z),
         double rmax, double zmax, double ta) {
  pml_fmin = 0.2;
  int r,z;
  outdir = ".";
  pb = NULL;
  a = ta;
  nr = (int) (rmax*a);
  nz = (int) (zmax*a);
  npmlz = 0;
  npmlr = 0;
  if (nz == 0) nz = 1;
  eps = new double[nr*(nz+1)];
  if (eps == NULL) {
    printf("Out of memory!\n");
    exit(1);
  }
  for (r=0;r<nr;r++) {
    for (z=0;z<nz+1;z++) {
      MA(eps,r,z) = (feps) ? feps((r+0.5)/a,(z+0.5)/a) : 1.0; // Null feps means vacuum.
    }
  }
  invepser = new double[(nr+1)*(nz+1)];
  invepsep = new double[(nr+1)*(nz+1)];
  invepsez = new double[(nr+1)*(nz+1)];

  // Initialize eps to 1;
  for (int i=0;i<(nr+1)*(nz+1);i++)
    invepser[i] = invepsep[i] = invepsez[i] = 1;
  for (r=1;r<nr;r++) {
    for (z=1;z<=nz;z++) {
      MA(invepser,r,z) = 2./(MA(eps,r,z)+MA(eps,r,z-1));
      MA(invepsep,r,z) = 4./(MA(eps,r,z)+MA(eps,r-1,z)+
                            MA(eps,r,z-1)+MA(eps,r-1,z-1));
      MA(invepsez,r,z) = 2./(MA(eps,r-1,z)+MA(eps,r,z));
    }
  }
  for (r=1;r<nr;r++) {
    {
      const int z=0;
      MA(invepser,r,z) = 2./(MA(eps,r,z)+MA(eps,r,nz-1));
      MA(invepsep,r,z) = 4./(MA(eps,r,z)+MA(eps,r-1,z)+
                          MA(eps,r,nz-1)+MA(eps,r-1,nz-1));
      MA(invepsez,r,z) = 2./(MA(eps,r-1,z)+MA(eps,r,z));
    }
    {
      const int z=nz;
      MA(invepser,r,z) = 2./(MA(eps,r,0)+MA(eps,r,z-1));
      MA(invepsep,r,z) = 4./(MA(eps,r,0)+MA(eps,r-1,0)+
                          MA(eps,r,z-1)+MA(eps,r-1,z-1));
    }
  }
  for (z=1;z<=nz;z++) {
    {
      const int r=0;
      MA(invepser,r,z) = 2./(MA(eps,r,z)+MA(eps,r,z-1));
      MA(invepsep,r,z) = 2./(MA(eps,r,z)+ MA(eps,r,z-1));
      MA(invepsez,r,z) = 1./MA(eps,r,z);
    }
  }
  {
    const int r=0,z=0;
    MA(invepser,r,z) = 2./(MA(eps,r,z)+MA(eps,r,nz-1));
    MA(invepsep,r,z) = 2./(MA(eps,r,z)+MA(eps,r,nz-1));
    MA(invepsez,r,z) = 1./MA(eps,r,z);
  }
  {
    const int r=0,z=nz;
    MA(invepser,r,z) = 2./(MA(eps,r,0)+MA(eps,r,nz-1));
    MA(invepsep,r,z) = 2./(MA(eps,r,0)+MA(eps,r,z-1));
    MA(invepsez,r,z) = 1./MA(eps,r,0);
  }
  // Allocate the conductivity arrays:
  Crez = Crep = Crhz = Crhp = Cper = Cpez = Cphr = Cphz = NULL;
  Czer = Czep = Czhr = Czhp = NULL;
}
