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

#include "tidod.h"
#include "tidod_internals.h"

polarization_1d *polarization_1d::set_up_polarizations(const mat_1d *ma) {
  if (ma->pb == NULL) return NULL;
  return new polarization_1d(ma->pb);
}

polarization_1d::polarization_1d(const polarizability_1d *the_pb) {
  const int nz = the_pb->nz, npmlz = the_pb->npmlz;
  DOCMP {
    Px[cmp] = new double[nz+1];
    for (int i=0;i<nz+1;i++) Px[cmp][i] = 0.0;
  }
  if (Px[1] == NULL) {
    printf("Allocation error in polarization!\n");
    exit(1);
  }
  pb = the_pb;
  if (pb->next == NULL) {
    next = NULL;
  } else {
    next = new polarization_1d(pb->next);
  }
}

polarization_1d::~polarization_1d() {
  DOCMP {
    delete[] Px[cmp];
  }
  if (next) delete next;
}

polarizability_1d::polarizability_1d(const polarizability_1d *pb) {
  omeganot = pb->omeganot;
  gamma = pb->gamma;
  nz = pb->nz;
  npmlz = pb->npmlz;
  sigma = new double[nz+1];
  for (int i=0;i<nz+1;i++) sigma[i] = pb->sigma[i];
  if (pb->next) next = new polarizability_1d(pb->next);
  else next = NULL;
}

void polarizability_1d::use_integer_pml(int new_npmlz) {
  npmlz = new_npmlz;
  if (next) next->use_integer_pml(npmlz);
}

polarizability_1d::polarizability_1d(const mat_1d *ma, double sig(double),
                                     double om, double ga, double sigscale) {
  nz = ma->nz;
  npmlz = ma->npmlz;
  const double a = ma->a;
  omeganot = om;
  gamma = ga;
  next = NULL;

  sigma = new double[nz+1];
  if (sigma == NULL) {
    printf("Out of memory!\n");
    exit(1);
  }
  for (int z=0;z<nz+1;z++) {
    MA(sigma,z) = (sig) ? sigscale*sig(z/a) : 0.0; // Null sig means vacuum.
  }
}

polarizability_1d::~polarizability_1d() {
  delete[] sigma;
}

void mat_1d::add_polarizability(double sigma(double),
                                double omega, double gamma, double delta_epsilon) {
  const double freq_conversion = 2*pi*c/a;
  double sigma_scale  = freq_conversion*freq_conversion*omega*omega*delta_epsilon;
  polarizability_1d *npb = new polarizability_1d(this, sigma,
                                                 freq_conversion*omega,
                                                 freq_conversion*gamma,
                                                 sigma_scale);
  npb->next = pb;
  pb = npb;
}

void fields_1d::initialize_polarizations(polarization_1d *op, polarization_1d *np) {
  // Set this up as a noop for now.
}

void fields_1d::step_polarization_itself(polarization_1d *op, polarization_1d *np) {
  if (op == NULL && np == NULL && olpol != NULL && pol != NULL) {
    // This is the initial call... so I should start running from olpol and pol.
    step_polarization_itself(olpol, pol);
    // The two polarizations get switched during the pml step...
  } else if (olpol != NULL && pol != NULL) {
    const double g = op->pb->gamma;
    const double om = op->pb->omeganot;
    const double funinv = 1.0/(1+0.5*g);
    const double *s = np->pb->sigma;
    DOCMP {
      for (int z=0;z<=nz;z++)
        CM(op->Px,z) = funinv*((2-om*om)*CM(np->Px,z)+ (0.5*g-1)*CM(op->Px,z))
          + MA(s,z)*CM(ex,z);
    }
    if (op->next && np->next) step_polarization_itself(op->next, np->next);
  }
}

void fields_1d::step_e_polarization(polarization_1d *op, polarization_1d *np) {
  if (op == NULL && np == NULL && olpol != NULL && pol != NULL) {
    // This is the initial call... so I should start running from olpol and pol.
    step_e_polarization(olpol, pol);
  } else if (olpol != NULL && pol != NULL) {
    DOCMP {
      for (int z=0;z<=nz;z++)
        CM(ex,z) -= MA(ma->inveps,z)*(CM(np->Px,z)-CM(op->Px,z));
    }
    if (op->next && np->next) step_e_polarization(op->next, np->next);
  }
}
