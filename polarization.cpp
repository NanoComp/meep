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

polarization *polarization::set_up_polarizations(const mat *ma) {
  if (ma->pb == NULL) return NULL;
  return new polarization(ma->pb);
}

polarization::polarization(const polarizability *the_pb) {
  const int nr = the_pb->nr, nz = the_pb->nz,
    npmlr = the_pb->npmlr, npmlz = the_pb->npmlz;
  DOCMP {
    Pr[cmp] = new double[nr*(nz+1)];
    Pp[cmp] = new double[(nr+1)*(nz+1)];
    Pz[cmp] = new double[(nr+1)*(nz+1)];
    for (int i=0;i<nr*(nz+1);i++) Pr[cmp][i] = 0.0;
    for (int i=0;i<(nr+1)*(nz+1);i++) Pp[cmp][i] = 0.0;
    for (int i=0;i<(nr+1)*(nz+1);i++) Pz[cmp][i] = 0.0;
    Prp[cmp] = new double[npmlr*(nz+1)];
    Ppz[cmp] = new double[npmlr*(nz+1)];
    Pzr[cmp] = new double[npmlr*(nz+1)];
    if (npmlr) for (int i=0;i<npmlr*(nz+1);i++) Prp[cmp][i] = 0.0;
    if (npmlr) for (int i=0;i<npmlr*(nz+1);i++) Ppz[cmp][i] = 0.0;
    if (npmlr) for (int i=0;i<npmlr*(nz+1);i++) Pzr[cmp][i] = 0.0;
    z_Prp[cmp][0] = new double[npmlz*(nr+1)];
    z_Prp[cmp][1] = new double[npmlz*(nr+1)];
    z_Ppz[cmp][0] = new double[npmlz*(nr+1)];
    z_Ppz[cmp][1] = new double[npmlz*(nr+1)];
    if (npmlz) for (int i=0;i<npmlz*(nr+1);i++) z_Prp[cmp][0][i] = 0.0;
    if (npmlz) for (int i=0;i<npmlz*(nr+1);i++) z_Prp[cmp][1][i] = 0.0;
    if (npmlz) for (int i=0;i<npmlz*(nr+1);i++) z_Ppz[cmp][0][i] = 0.0;
    if (npmlz) for (int i=0;i<npmlz*(nr+1);i++) z_Ppz[cmp][1][i] = 0.0;
  }
  if (Pz[1] == NULL) {
    printf("Allocation error in polarization!\n");
    exit(1);
  }
  pb = the_pb;
  if (pb->next == NULL) {
    next = NULL;
  } else {
    next = new polarization(pb->next);
  }
}

polarization::~polarization() {
  DOCMP {
    delete[] Pr[cmp];
    delete[] Pp[cmp];
    delete[] Pz[cmp];
    delete[] Prp[cmp];
    delete[] Ppz[cmp];
    delete[] Pzr[cmp];
    delete[] z_Prp[cmp][0];
    delete[] z_Prp[cmp][1];
    delete[] z_Ppz[cmp][0];
    delete[] z_Ppz[cmp][1];
  }
  if (next) delete next;
}

polarizability::polarizability(const polarizability *pb) {
  omeganot = pb->omeganot;
  gamma = pb->gamma;
  nr = pb->nr;
  nz = pb->nz;
  npmlr = pb->npmlr;
  npmlz = pb->npmlz;
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
  if (pb->next) next = new polarizability(pb->next);
  else next = NULL;
}

void polarizability::use_integer_pml(int new_npmlr, int new_npmlz) {
  npmlz = new_npmlz;
  npmlr = new_npmlr;
  if (next) next->use_integer_pml(npmlr, npmlz);
}

polarizability::polarizability(const mat *ma, double sig(double,double),
                               double om, double ga, double sigscale) {
  nr = ma->nr;
  nz = ma->nz;
  npmlr = ma->npmlr;
  npmlz = ma->npmlz;
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

inline double expi(int cmp, double x) {
  return (cmp) ? cos(x) : sin(x);
}

void fields::initialize_polarizations(polarization *op, polarization *np) {
  // Set up polarizations so we'll have them nicely excited, which should
  // give us a handy way of getting all the modes out of a polaritonic
  // material.
  if (op == NULL && np == NULL && olpol != NULL && pol != NULL) {
    initialize_polarizations(olpol, pol);
  } else if (olpol != NULL && pol != NULL) {
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

void fields::step_polarization_itself(polarization *op, polarization *np) {
  if (op == NULL && np == NULL && olpol != NULL && pol != NULL) {
    // This is the initial call... so I should start running from olpol and pol.
    step_polarization_itself(olpol, pol);
    // The two polarizations get switched during the pml step...
  } else if (olpol != NULL && pol != NULL) {
    const double g = op->pb->gamma;
    const double om = op->pb->omeganot;
    const double funinv = 1.0/(1+0.5*g);
    const double *sr = np->pb->sr;
    const double *sp = np->pb->sp;
    const double *sz = np->pb->sz;
    DOCMP {
      for (int r=rmin_bulk(m)-1;r<nr;r++) for (int z=0;z<=nz;z++)
        CM(op->Pr,r,z) = funinv*((2-om*om)*CM(np->Pr,r,z)+
                                 (0.5*g-1)*CM(op->Pr,r,z))+
                         MA(sr,r,z)*CM(er,r,z);
      for (int r=rmin_bulk(m)-1;r<nr;r++) for (int z=0;z<=nz;z++)
        CM(op->Pp,r,z) = funinv*((2-om*om)*CM(np->Pp,r,z)+
                                 (0.5*g-1)*CM(op->Pp,r,z))+
                         MA(sp,r,z)*CM(ep,r,z);
      for (int r=rmin_bulk(m)-1;r<nr;r++) for (int z=0;z<nz;z++)
        CM(op->Pz,r,z) = funinv*((2-om*om)*CM(np->Pz,r,z)+
                                 (0.5*g-1)*CM(op->Pz,r,z))+
                         MA(sz,r,z)*CM(ez,r,z);
    }
    if (op->next && np->next) step_polarization_itself(op->next, np->next);
  }
}

void fields::step_polarization_pml(polarization *op, polarization *np) {
  if (op == NULL && np == NULL && olpol != NULL && pol != NULL) {
    // This is the initial call... so I should start running from olpol and pol.
    step_polarization_pml(olpol, pol);
    polarization *temp = olpol;
    olpol = pol;
    pol = temp; // They got switched....
  } else if (olpol != NULL && pol != NULL) {
    const double g = op->pb->gamma;
    const double om = op->pb->omeganot;
    const double funinv = 1.0/(1+0.5*g);
    const double *sr = np->pb->sr;
    const double *sp = np->pb->sp;
    const double *sz = np->pb->sz;
    DOCMP {
      if (npmlz) {
        for (int r=rmin_bulk(m)-1;r<nr-npmlr;r++) {
          int z0 = npmlz;
          for (int lr=-1;lr<2;lr+=2,z0=nz-npmlz) {
            int z = z0;
            for (int iz=0;iz<npmlz;iz++,z+=lr) {
              PMLZ(op->z_Prp,r) = funinv*((2-om*om)*PMLZ(np->z_Prp,r)+
                                          (0.5*g-1)*PMLZ(op->z_Prp,r))+
                                  MA(sr,r,z)*PMLZ(z_erp,r);
              PMLZ(op->z_Ppz,r) = funinv*((2-om*om)*PMLZ(np->z_Ppz,r)+
                                          (0.5*g-1)*PMLZ(op->z_Ppz,r))+
                                  MA(sp,r,z)*PMLZ(z_epz,r);
            }
          }
        }
      }
      if (npmlr) {
        // update large r pml for all z (except actual boundary)...
        for (int r=nr-npmlr;r<nr;r++) {
          for (int z=1;z<nz;z++) {
            PMLR(op->Prp,r,z) = funinv*((2-om*om)*PMLR(np->Prp,r,z)+
                                        (0.5*g-1)*PMLR(op->Prp,r,z))+
                                MA(sr,r,z)*PMLR(erp,r,z);
            PMLR(op->Ppz,r,z) = funinv*((2-om*om)*PMLR(np->Ppz,r,z)+
                                        (0.5*g-1)*PMLR(op->Ppz,r,z))+
                                MA(sp,r,z)*PMLR(epz,r,z);
            PMLR(op->Pzr,r,z) = funinv*((2-om*om)*PMLR(np->Pzr,r,z)+
                                        (0.5*g-1)*PMLR(op->Pzr,r,z))+
                                MA(sz,r,z)*PMLR(ezr,r,z);
          }
        }
      }
    }
    if (op->next && np->next) step_polarization_pml(op->next, np->next);
  }
}

void fields::step_e_polarization(polarization *op, polarization *np) {
  if (op == NULL && np == NULL && olpol != NULL && pol != NULL) {
    // This is the initial call... so I should start running from olpol and pol.
    step_e_polarization(olpol, pol);
  } else if (olpol != NULL && pol != NULL) {
    DOCMP {
      for (int r=rmin_bulk(m)-1;r<nr-npmlr;r++) {
        for (int z=npmlz+1;z<nz-npmlz;z++) {
          CM(er,r,z) -= MA(ma->invepser,r,z)*(CM(np->Pr,r,z)-CM(op->Pr,r,z));
          CM(ep,r,z) -= MA(ma->invepsep,r,z)*(CM(np->Pp,r,z)-CM(op->Pp,r,z));
        }
        {
          const int z = 0;
          CM(er,r,z) -= MA(ma->invepser,r,z)*(CM(np->Pr,r,z)-CM(op->Pr,r,z));
          CM(ep,r,z) -= MA(ma->invepsep,r,z)*(CM(np->Pp,r,z)-CM(op->Pp,r,z));
        }
        {
          const int z = nz;
          CM(er,r,z) -= MA(ma->invepser,r,z)*(CM(np->Pr,r,z)-CM(op->Pr,r,z));
          CM(ep,r,z) -= MA(ma->invepsep,r,z)*(CM(np->Pp,r,z)-CM(op->Pp,r,z));
        }
      }
      for (int r=rmin_bulk(m)-1;r<nr-npmlr;r++) {
        for (int z=0;z<=nz;z++) {
          CM(ez,r,z) -= MA(ma->invepsez,r,z)*(CM(np->Pz,r,z)-CM(op->Pz,r,z));
        }
      }
    }
    if (op->next && np->next) step_e_polarization(op->next, np->next);
  }
}

void fields::step_e_pml_polarization(polarization *op, polarization *np) {
  if (op == NULL && np == NULL && olpol != NULL && pol != NULL) {
    // This is the initial call... so I should start running from olpol and pol.
    step_e_pml_polarization(olpol, pol);
  } else if (olpol != NULL && pol != NULL) {
    DOCMP {
      if (npmlz) {
        for (int r=rmin_bulk(m)-1;r<nr-npmlr;r++) {
          int z0 = npmlz;
          for (int lr=-1;lr<2;lr+=2,z0=nz-npmlz) {
            int z = z0;
            for (int iz=0;iz<npmlz;iz++,z+=lr) {
              double Czer = ma->Czer[iz];
              double Czep = ma->Czep[iz];

              double derp = -MA(ma->invepser,r,z)*(PMLZ(np->z_Prp,r)-PMLZ(op->z_Prp,r));
              double dPrz = (CM(np->Pr,r,z)-CM(op->Pr,r,z))-
                (PMLZ(np->z_Prp,r)-PMLZ(op->z_Prp,r));
              PMLZ(z_erp,r) += derp;
              CM(er,r,z) += derp - MA(ma->invepser,r,z)*dPrz/(1+.5*MA(ma->invepser,r,z)*Czer);
              
              double depz = -MA(ma->invepsep,r,z)*(PMLZ(np->z_Ppz,r)-PMLZ(op->z_Ppz,r))
                /(1+.5*MA(ma->invepsep,r,z)*Czep);
              double dPpr = (CM(np->Pp,r,z)-CM(op->Pp,r,z))-
                (PMLZ(np->z_Ppz,r)-PMLZ(op->z_Ppz,r));
              PMLZ(z_epz,r) += depz;
              CM(ep,r,z) += depz - MA(ma->invepsep,r,z)*dPpr;
            }
          }
        }
      }
      if (npmlr) {
        // update large r pml for all z (except actual boundary)...
        for (int r=nr-npmlr;r<nr;r++) {
          double Cper = ma->Cper[r-nr+npmlr];
          double Crep = ma->Crep[r-nr+npmlr];
          double Crez = ma->Crez[r-nr+npmlr];
          double Cpez = ma->Cpez[r-nr+npmlr];
          for (int z=1;z<nz;z++) {
            double Czep, Czer;
            if (z <= npmlz) {
              Czer = ma->Czer[npmlz - z];
              Czep = ma->Czep[npmlz - z];
            } else if (z >= nz - npmlz) {
              Czer = ma->Czer[z+npmlz-nz];
              Czep = ma->Czep[z+npmlz-nz];
            } else {
              Czer = 0;
              Czep = 0;
            }
            double derp = -MA(ma->invepser,r,z)*(PMLR(np->Prp,r,z)-PMLR(op->Prp,r,z))
              /(1+.5*MA(ma->invepser,r,z)*Cper);
            double dPrz = (CM(np->Pr,r,z)-CM(op->Pr,r,z))-(PMLR(np->Prp,r,z)-PMLR(op->Prp,r,z));
            PMLR(erp,r,z) += derp;
            CM(er,r,z) += derp - MA(ma->invepser,r,z)*dPrz/(1+.5*MA(ma->invepser,r,z)*Czer);

            double depz = -MA(ma->invepsep,r,z)*(PMLR(np->Ppz,r,z)-PMLR(op->Ppz,r,z))
              /(1+.5*MA(ma->invepsep,r,z)*Czep);
            double dPpr = (CM(np->Pp,r,z)-CM(op->Pp,r,z))-(PMLR(np->Ppz,r,z)-PMLR(op->Ppz,r,z));
            PMLR(epz,r,z) += depz;
            CM(ep,r,z) += depz - MA(ma->invepsep,r,z)*dPpr/(1+.5*MA(ma->invepsep,r,z)*Crep);

            double dezr = -MA(ma->invepsez,r,z)*(PMLR(np->Pzr,r,z)-PMLR(op->Pzr,r,z))
              /(1+.5*MA(ma->invepsez,r,z)*Crez);
            double dPzp = (CM(np->Pz,r,z)-CM(op->Pz,r,z))-(PMLR(np->Pzr,r,z)-PMLR(op->Pzr,r,z));
            PMLR(ezr,r,z) += dezr;
            CM(ez,r,z) += dezr - MA(ma->invepsez,r,z)*dPzp/(1+.5*MA(ma->invepsez,r,z)*Cpez);
          }
        }
      }
    }
    if (op->next && np->next) step_e_pml_polarization(op->next, np->next);
  }
}
