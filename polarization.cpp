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

polarization *polarization::set_up_polarizations(const mat *ma) {
  if (ma->pb == NULL) return NULL;
  return new polarization(ma->nr, ma->nz, ma->pb);
}

polarization::polarization(const int nr, const int nz, const polarizability *the_pb) {
  DOCMP {
    Pr[cmp] = new double[nr*(nz+1)];
    Pp[cmp] = new double[(nr+1)*(nz+1)];
    Pz[cmp] = new double[(nr+1)*(nz+1)];
  }
  if (Pz[1] == NULL) {
    printf("Allocation error in polarization!\n");
    exit(1);
  }
  pb = the_pb;
  if (pb->next == NULL) {
    next = NULL;
  } else {
    next = new polarization(nr, nz, pb->next);
  }
}

polarization::~polarization() {
  DOCMP {
    delete[] Pr[cmp];
    delete[] Pp[cmp];
    delete[] Pz[cmp];
  }
  if (next) delete next;
}

polarizability::polarizability(const polarizability *pb) {
  omeganot = pb->omeganot;
  gamma = pb->gamma;
  nr = pb->nr;
  nz = pb->nz;
  sr = new double[nr*(nz+1)];
  sp = new double[nr*(nz+1)];
  sz = new double[nr*(nz+1)];
  sigma = new double[nr*(nz+1)];
  for (int i=0;i<nr*(nz+1);i++) {
    sr[i] = pb->sr[i];
    sp[i] = pb->sp[i];
    sz[i] = pb->sz[i];
    sigma[i] = pb->sigma[i];
  }
}

polarizability::polarizability(const mat *ma, double sig(double,double), double om, double ga) {
  nr = ma->nr;
  nz = ma->nz;
  const double a = ma->a;
  omeganot = om;
  gamma = ga;

  sr = new double[nr*(nz+1)];
  sp = new double[nr*(nz+1)];
  sz = new double[nr*(nz+1)];
  sigma = new double[nr*(nz+1)];
  if (sigma == NULL) {
    printf("Out of memory!\n");
    exit(1);
  }
  for (int r=0;r<nr;r++) {
    for (int z=0;z<nz+1;z++) {
      MA(sigma,r,z) = (sig) ? sig((r+0.5)/a,(z+0.5)/a) : 0.0; // Null sig means vacuum.
    }
  }
  // Average out sigma over the grid...
  for (int i=0;i<nr*(nz+1);i++)
    sr[i] = sp[i] = sz[i] = sqrt(-1);
  for (int r=1;r<nr;r++) {
    for (int z=1;z<=nz;z++) {
      MA(sr,r,z) = 2./(MA(sigma,r,z)+MA(sigma,r,z-1));
      MA(sp,r,z) = 4./(MA(sigma,r,z)+MA(sigma,r-1,z)+
                            MA(sigma,r,z-1)+MA(sigma,r-1,z-1));
      MA(sz,r,z) = 2./(MA(sigma,r-1,z)+MA(sigma,r,z));
    }
  }
  for (int r=1;r<nr;r++) {
    {
      const int z=0;
      MA(sr,r,z) = 0.5*(MA(sigma,r,z)+MA(sigma,r,nz-1));
      MA(sp,r,z) = 0.25*(MA(sigma,r,z)+MA(sigma,r-1,z)+
                          MA(sigma,r,nz-1)+MA(sigma,r-1,nz-1));
      MA(sz,r,z) = 0.5*(MA(sigma,r-1,z)+MA(sigma,r,z));
    }
    {
      const int z=nz;
      MA(sr,r,z) = 0.5*(MA(sigma,r,0)+MA(sigma,r,z-1));
      MA(sp,r,z) = 0.25*(MA(sigma,r,0)+MA(sigma,r-1,0)+
                          MA(sigma,r,z-1)+MA(sigma,r-1,z-1));
    }
  }
  for (int z=1;z<=nz;z++) {
    {
      const int r=0;
      MA(sr,r,z) = 0.5*(MA(sigma,r,z)+MA(sigma,r,z-1));
      MA(sp,r,z) = 0.5*(MA(sigma,r,z)+ MA(sigma,r,z-1));
      MA(sz,r,z) = MA(sigma,r,z);
    }
  }
  {
    const int r=0,z=0;
    MA(sr,r,z) = 0.5*(MA(sigma,r,z)+MA(sigma,r,nz-1));
    MA(sp,r,z) = 0.5*(MA(sigma,r,z)+MA(sigma,r,nz-1));
    MA(sz,r,z) = MA(sigma,r,z);
  }
  {
    const int r=0,z=nz;
    MA(sr,r,z) = 0.5*(MA(sigma,r,0)+MA(sigma,r,nz-1));
    MA(sp,r,z) = 0.5*(MA(sigma,r,0)+MA(sigma,r,z-1));
    MA(sz,r,z) = MA(sigma,r,0);
  }
}

polarizability::~polarizability() {
  delete[] sr;
  delete[] sp;
  delete[] sz;
  delete[] sigma;
}
