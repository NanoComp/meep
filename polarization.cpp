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
    for (int i=0;i<nr*(nz+1);i++) Pr[cmp][i] = 0.0;
    for (int i=0;i<(nr+1)*(nz+1);i++) Pp[cmp][i] = 0.0;
    for (int i=0;i<(nr+1)*(nz+1);i++) Pz[cmp][i] = 0.0;
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
  next = NULL;
}

polarizability::polarizability(const mat *ma, double sig(double,double),
                               double om, double ga, double sigscale) {
  nr = ma->nr;
  nz = ma->nz;
  const double a = ma->a;
  omeganot = om;
  gamma = ga;
  next = NULL;

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
      MA(sigma,r,z) = (sig) ? sigscale*sig((r+0.5)/a,(z+0.5)/a) : 0.0; // Null sig means vacuum.
    }
  }
  // Average out sigma over the grid...
  //for (int i=0;i<nr*(nz+1);i++)
  //  sr[i] = sp[i] = sz[i] = sqrt(-1);
  for (int r=1;r<nr;r++) {
    for (int z=1;z<=nz;z++) {
      MA(sr,r,z) = 0.5*(MA(sigma,r,z)+MA(sigma,r,z-1));
      MA(sp,r,z) = 0.25*(MA(sigma,r,z)+MA(sigma,r-1,z)+
                            MA(sigma,r,z-1)+MA(sigma,r-1,z-1));
      MA(sz,r,z) = 0.5*(MA(sigma,r-1,z)+MA(sigma,r,z));
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

void mat::add_polarizability(double sigma(double,double),
                             double omega, double gamma, double delta_epsilon) {
  const double freq_conversion = 2*pi*c/a;
  double sigma_scale  = freq_conversion*freq_conversion*omega*omega*delta_epsilon;
  polarizability *npb = new polarizability(this, sigma,
                                           freq_conversion*omega,
                                           freq_conversion*gamma,
                                           sigma_scale);
  npb->next = pb;
  pb = npb;
}

void fields::step_polarization_itself(polarization *op, polarization *np) {
  if (op == NULL && np == NULL && olpol != NULL && pol != NULL) {
    // This is the initial call... so I should start running from olpol and pol.
    step_polarization_itself(olpol, pol);
    polarization *temp = olpol;
    olpol = pol;
    pol = temp; // They got switched....
  } else if (olpol != NULL && pol != NULL) {
    const double g = op->pb->gamma;
    const double om = op->pb->omeganot;
    const double funinv = 1.0/(1-0.5*g);
    const double *sr = np->pb->sr;
    const double *sp = np->pb->sp;
    const double *sz = np->pb->sz;
    DOCMP {
      for (int r=0;r<nr;r++) for (int z=0;z<=nz;z++)
        CM(op->Pr,r,z) = funinv*((2-om*om)*CM(np->Pr,r,z)+
                                 (0.5*g-1)*CM(op->Pr,r,z)+
                                 MA(sr,r,z)*CM(er,r,z));
      for (int r=0;r<nr;r++) for (int z=0;z<=nz;z++)
        CM(op->Pp,r,z) = funinv*((2-om*om)*CM(np->Pp,r,z)+
                                 (0.5*g-1)*CM(op->Pp,r,z)+
                                 MA(sp,r,z)*CM(ep,r,z));
      for (int r=0;r<nr;r++) for (int z=0;z<nz;z++)
        CM(op->Pz,r,z) = funinv*((2-om*om)*CM(np->Pz,r,z)+
                                 (0.5*g-1)*CM(op->Pz,r,z)+
                                 MA(sz,r,z)*CM(ez,r,z));
    }
    if (op->next && np->next) step_polarization_itself(op->next, np->next);
  }
}

inline double expi(int cmp, double x) {
  return (cmp) ? cos(x) : sin(x);
}

void fields::initialize_polarizations(polarization *op, polarization *np) {
  // Set up polarizations so we'll have them nicely excited, which should
  // give us a handy way of getting all the modes out of a polaritonic
  // material.
  if (op == NULL && np == NULL && olpol != NULL && pol != NULL) {
    initialize_polarizations(olpol, pol);
  } else {
    double omt = op->pb->omeganot;
    double amp_shift = exp(op->pb->gamma);
    double sinkz = sin(-omt);
    double coskz = cos(-omt);
    DOCMP {
      for (int r=0;r<nr;r++) for (int z=0;z<=nz;z++) {
        CM(op->Pr,r,z) = CM(er,r,z);
        CM(np->Pr,r,z) = EIKZ(er,r,z);
      }
      for (int r=0;r<nr;r++) for (int z=0;z<=nz;z++) {
        CM(op->Pp,r,z) = CM(ep,r,z);
        CM(np->Pp,r,z) = EIKZ(ep,r,z);
      }
      for (int r=0;r<nr;r++) for (int z=0;z<nz;z++) {
        CM(op->Pz,r,z) = CM(ez,r,z);
        CM(np->Pz,r,z) = EIKZ(ez,r,z);
      }
    }
    if (op->next && np->next) initialize_polarizations(op->next, np->next);
  }
}

void fields::step_e_polarization(polarization *op, polarization *np) {
  if (op == NULL && np == NULL && olpol != NULL && pol != NULL) {
    // This is the initial call... so I should start running from olpol and pol.
    step_e_polarization(olpol, pol);
  } else if (olpol != NULL && pol != NULL) {
    DOCMP {
      for (int r=0;r<nr;r++) for (int z=0;z<=nz;z++) {
        if (CM(op->Pr,r,z) != CM(op->Pr,r,z)) printf("nan in op pr\n");
        if (CM(np->Pr,r,z) != CM(np->Pr,r,z)) printf("nan in np pr\n");
        CM(er,r,z) -= MA(ma->invepser,r,z)*(CM(np->Pr,r,z)-CM(op->Pr,r,z));
      }
      for (int r=0;r<nr;r++) for (int z=0;z<=nz;z++)
        CM(ep,r,z) -= MA(ma->invepsep,r,z)*(CM(np->Pp,r,z)-CM(op->Pp,r,z));
      for (int r=0;r<nr;r++) for (int z=0;z<nz;z++)
        CM(ez,r,z) -= MA(ma->invepsez,r,z)*(CM(np->Pz,r,z)-CM(op->Pz,r,z));
    }
    if (op->next && np->next) step_e_polarization(op->next, np->next);
  }
}
