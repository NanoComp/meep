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
}

static double sig(double r);

void mat::use_pml(int numpmlr, int numpmlz) {
  npmlz = numpmlz;
  npmlr = numpmlr;
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
    Crez[r] = Crep[r] = Cmax*0.5*(sig(rp)+sig(rr));
    Crhz[r] = Crhp[r] = Cmax*sig(rp);
  }
  double sigintegrated = 0.0;
  for (int r=0;r<npmlr;r++) {
    double rr = (r)/(double)npmlr;
    double rp = (1+r)/(double)npmlr;
    sigintegrated += Cmax*0.5*(sig(rp)+sig(rr));
    Cper[r] = Cphz[r] = sigintegrated/(nr-npmlr+r+0.5);
    if (r==0) Cpez[r] = Cphr[r] = 0.5*Cper[r];
    else Cpez[r] = Cphr[r] = 0.5*(Cper[r]+Cper[r-1]);

    Cphz[r] = 0.0;
    Cper[r] = 0.0;
    Cphr[r] = 0.0;
    Cpez[r] = 0.0;
    printf("Big sig(Ez)[%d  ] is %10lg, little sig is %10lg\n", r, Crez[r], Cpez[r]);
    printf("Big sig(Hz)[%d  ] is %10lg, little sig is %10lg\n", r, Crhz[r], Cphz[r]);
  }
  for (int z=0;z<npmlz;z++) {
    double rr = (z)/(double)npmlz;
    double rp = (1+z)/(double)npmlz;
    Czer[z] = Czep[z] = Cmax*0.5*(sig(rp)+sig(rr));
    Czhr[z] = Czhp[z] = Cmax*sig(rp);
  }
}

static double sig(double r) {
  return pow(r, 2);
}

void mat::mix_with(const mat *n, double f) {
  for (int i=0;i<nr*(nz+1);i++) {
    eps[i] = 1.0/(1.0/eps[i] + f*(1.0/n->eps[i]-1.0/eps[i]));
    invepser[i] += f*(n->invepser[i] - invepser[i]);
    invepsep[i] += f*(n->invepsep[i] - invepsep[i]);
    invepsez[i] += f*(n->invepsez[i] - invepsez[i]);
  }
}

void mat::make_vacuum() {
  for (int i=0;i<nr*(nz+1);i++) {
    eps[i] = 1;
    invepser[i] = 1;
    invepsep[i] = 1;
    invepsez[i] = 1;
  }
}

mat::mat(const mat *o) {
  outdir = o->outdir;
  numpols = o->numpols;
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
  invepser = new double[nr*(nz+1)];
  invepsep = new double[nr*(nz+1)];
  invepsez = new double[nr*(nz+1)];

  for (int i=0;i<nr*(nz+1);i++) {
    invepser[i] = o->invepser[i];
    invepsep[i] = o->invepsep[i];
    invepsez[i] = o->invepsez[i];
  }
  // Allocate the conductivity arrays:
  Crez = Crep = Crhz = Crhp = Cper = Cpez = Cphr = Cphz = NULL;
  Czer = Czep = Czhr = Czhp = NULL;
  npmlz = 0;
  npmlr = 0;
  use_pml(npmlr,npmlz);
}

mat::mat(double feps(double r, double z),
         double rmax, double zmax, double ta) {
  int r,z;
  outdir = ".";
  numpols = 0;
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
      MA(eps,r,z) = (feps) ? feps((r+0.5)/a,(z+0.5)) : 1.0; // Null feps means vacuum.
    }
  }
  invepser = new double[nr*(nz+1)];
  invepsep = new double[nr*(nz+1)];
  invepsez = new double[nr*(nz+1)];

  // Initialize eps to NaNs so we'll notice if we don't set it properly.
  for (int i=0;i<nr*(nz+1);i++)
    invepser[i] = invepsep[i] = invepsez[i] = sqrt(-1);
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
