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

polarization *polarization::set_up_polarizations(const mat *ma, int is_r) {
  if (ma->pb == NULL) return NULL;
  return new polarization(ma->pb, is_r);
}

void polarization::use_real_fields() {
  is_real = 1;
  for (int c=0;c<10;c++) delete[] P[c][1];
  for (int c=0;c<10;c++) delete[] P_pml[c][1];
  for (int c=0;c<10;c++) P[c][1] = NULL;
  for (int c=0;c<10;c++) P_pml[c][1] = NULL;
  if (next) next->use_real_fields();
}

polarization::polarization(const polarizability *the_pb, int is_r) {
  const volume &v = the_pb->v;
  is_real = is_r;
  DOCMP {
    for (int c=0;c<10;c++)
      if (v.has_field((component)c) && is_electric((component)c)) {
        P[c][cmp] = new double[v.ntot()];
        for (int i=0;i<v.ntot();i++) P[c][cmp][i] = 0.0;
        // FIXME perhaps shouldn't allocate the PML split fields if we don't
        // have pml...
        P_pml[c][cmp] = new double[v.ntot()];
        if (P_pml[c][cmp] == NULL) {
          printf("Allocation error in polarization!\n");
          exit(1);
        }
        for (int i=0;i<v.ntot();i++) P_pml[c][cmp][i] = 0.0;
      } else {
        P[c][cmp] = NULL;
        P_pml[c][cmp] = NULL;
      }
  }
  for (int c=0;c<10;c++)
    if (v.has_field((component)c) && is_electric((component)c)) {
      energy[c] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) energy[c][i] = 0.0;
    } else {
      energy[c] = NULL;
    }
  pb = the_pb;
  // Initialize the s[] arrays that point to sigma.
  for (int c=0;c<10;c++)
    if (pb->energy_saturation != 0.0) {
      if (pb->s[c]) {
        s[c] = new double[v.ntot()];
        for (int i=0;i<v.ntot();i++) s[c][i] = pb->s[c][i];
      } else s[c] = NULL;
    } else s[c] = pb->s[c];
  // Deal with saturation stuff.
  if (pb->energy_saturation != 0.0) {
    saturation_factor = pb->saturated_sigma/pb->energy_saturation;
    const double isf = 1.0/fabs(saturation_factor);
    for (int c=0;c<10;c++)
      if (pb->s[c]) for (int i=0;i<v.ntot();i++) energy[c][i] = -isf*s[c][i];
  } else {
    saturation_factor = 0.0;
    for (int c=0;c<10;c++)
      if (energy[c]) for (int i=0;i<v.ntot();i++) energy[c][i] = 0.0;
  }
  if (pb->next == NULL) {
    next = NULL;
  } else {
    next = new polarization(pb->next, is_r);
  }
}

polarization::~polarization() {
  DOCMP {
    for (int c=0;c<10;c++) delete[] P[c][cmp];
    for (int c=0;c<10;c++) delete[] P_pml[c][cmp];
  }
  for (int c=0;c<10;c++) delete[] energy[c];
  if (saturation_factor != 0.0)
    for (int c=0;c<10;c++) delete[] s[c];
  if (next) delete next;
}

double polarization::total_energy(const volume &what) {
  const volume v = pb->v;
  double e = 0.0;
  for (int c=0;c<10;c++)
    if (energy[c])
      for (int i=0;i<v.ntot();i++)
        if (what.contains(v.loc((component)c,i)))
          e += v.dv((component)c,i)*energy[c][i];
  if (next) e += next->total_energy(what);
  return e;
}

polarizability::polarizability(const polarizability *pb) {
  omeganot = pb->omeganot;
  gamma = pb->gamma;
  v = pb->v;
  
  energy_saturation = pb->energy_saturation;
  saturated_sigma = pb->saturated_sigma;
  sigma = new double[v.ntot()];
  for (int i=0;i<v.ntot();i++) sigma[i] = pb->sigma[i];
  for (int c=0;c<10;c++)
    if (v.has_field((component)c) && is_electric((component)c)) {
      s[c] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) s[c][i] = pb->s[c][i];
    } else {
      s[c] = NULL;
    }
  if (pb->next) next = new polarizability(pb->next);
  else next = NULL;
}

void polarizability::use_pml() {
  // Dummy function for now...
  if (next) next->use_pml();
}

polarizability::polarizability(const mat *ma, double sig(const vec &),
                               double om, double ga, double sigscale,
                               double energy_sat) {
  v = ma->v;
  omeganot = om;
  gamma = ga;
  next = NULL;
  energy_saturation = energy_sat;
  saturated_sigma = sigscale;

  for (int c=0;c<10;c++)
    if (v.has_field((component)c) && is_electric((component)c)) {
      s[c] = new double[v.ntot()];
    } else {
      s[c] = NULL;
    }
  sigma = new double[v.ntot()];
  if (sigma == NULL) {
    printf("Out of memory in polarizability!\n");
    exit(1);
  }

  if (v.dim == dcyl) {
    for (int i=0;i<v.ntot();i++) sigma[i] = sigscale*sig(v.loc(Hp,i));
  } else if (v.dim == d1) {
    for (int i=0;i<v.ntot();i++) sigma[i] = sigscale*sig(v.loc(Ex,i));
  } else {
    printf("Unsupported dimensionality!\n");
    exit(1);
  }
  for (int c=0;c<10;c++) if (s[c])
    for (int i=0;i<v.ntot();i++) s[c][i] = 0.0;
  // Average out sigma over the grid...
  if (v.dim == dcyl) {
    const vec dr = v.dr()*0.5; // The distance between Yee field components
    const vec dz = v.dz()*0.5; // The distance between Yee field components
    for (int r=1;r<v.nr();r++) {
      const int ir = r*(v.nz()+1);
      const int irm1 = (r-1)*(v.nz()+1);
      for (int z=1;z<=v.nz();z++) {
        s[Er][z + ir] = 0.5*(sigma[z+ir] + sigma[z+ir-1]);
        s[Ep][z + ir] = 0.25*(sigma[z+ir] + sigma[z+ir-1] +
                           sigma[z+irm1] + sigma[z+irm1-1]);
        s[Ez][z + ir] = 0.5*(sigma[z+ir] + sigma[z+irm1]);
      }
    }
    for (int r=0;r<v.nr();r++) {
      const int ir = r*(v.nz()+1);
      const vec here = v.loc(Ep,ir);
      s[Er][ir] = 0.5*sigscale*(sig(here+dr+dz) + sig(here+dr-dz));
      s[Ep][ir] = 0.25*sigscale*(sig(here+dr+dz) + sig(here-dr+dz) +
                              sig(here+dr-dz) + sig(here-dr-dz));
      s[Ez][ir] = 0.5*sigscale*(sig(here+dr+dz) + sig(here-dr+dz));
    }
    for (int z=0;z<v.nz();z++) {
      const vec here = v.loc(Ep,z);
      s[Er][z] = 0.5*sigscale*(sig(here+dr+dz) + sig(here+dr-dz));
      s[Ep][z] = 0.25*sigscale*(sig(here+dr+dz) + sig(here-dr+dz) +
                             sig(here+dr-dz) + sig(here-dr-dz));
      s[Ez][z] = 0.5*sigscale*(sig(here+dr+dz) + sig(here-dr+dz));
    }
  } else if (v.dim == d1) {
    // There's just one field point...
    for (int i=0;i<v.ntot();i++) s[Ex][i] = sigma[i];
  } else {
    printf("Unsupported dimensionality!\n");
    exit(1);
  }
}

polarizability::~polarizability() {
  for (int c=0;c<10;c++) delete[] s[c];
  delete[] sigma;
}

void mat::add_polarizability(double sigma(const vec &),
                             double omega, double gamma, double delta_epsilon,
                             double energy_sat) {
  const double freq_conversion = 2*pi*c/a;
  double sigma_scale  = freq_conversion*freq_conversion*omega*omega*delta_epsilon;
  polarizability *npb = new polarizability(this, sigma,
                                           freq_conversion*omega,
                                           freq_conversion*gamma,
                                           sigma_scale, energy_sat);
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
      for (int c=0;c<10;c++)
        if (v.has_field((component)c) && is_electric((component)c))
          for (int i=0;i<v.ntot();i++) np->P[c][cmp][i] = op->P[c][cmp][i] = f[c][cmp][i];
    }
    if (op->next && np->next) initialize_polarizations(op->next, np->next);
  }
}

void fields::prepare_step_polarization_energy(polarization *op, polarization *np) {
  if (op == NULL && np == NULL && olpol != NULL && pol != NULL) {
    // This is the initial call... so I should start running from olpol and pol.
    prepare_step_polarization_energy(olpol, pol);
  } else if (op != NULL && np != NULL) {
    for (int c=0;c<10;c++)
      if (np->energy[c])
        for (int i=0;i<v.ntot();i++)
          np->energy[c][i] = op->energy[c][i];
    if (op->next && np->next) prepare_step_polarization_energy(op->next, np->next);
  }
}

void fields::half_step_polarization_energy(polarization *op, polarization *np) {
  if (op == NULL && np == NULL && olpol != NULL && pol != NULL) {
    // This is the initial call... so I should start running from olpol and pol.
    half_step_polarization_energy(olpol, pol);
  } else if (op != NULL && np != NULL) {
    DOCMP
      for (int c=0;c<10;c++)
        if (np->energy[c])
          for (int i=0;i<v.ntot();i++)
            np->energy[c][i] += 0.5*(np->P[c][cmp][i] - op->P[c][cmp][i])*f[c][cmp][i];
    if (op->next && np->next) half_step_polarization_energy(op->next, np->next);
  }
}

void fields::update_polarization_saturation(polarization *op, polarization *np) {
  if (op == NULL && np == NULL && olpol != NULL && pol != NULL) {
    // This is the initial call... so I should start running from olpol and pol.
    update_polarization_saturation(olpol, pol);
  } else if (olpol != NULL && pol != NULL) {
    if (np->saturation_factor != 0.0) {
      const volume v = np->pb->v;
      const double fac = np->saturation_factor;
      if (v.dim == d1) {
        if (fac > 0.0)
          for (int i=0;i<v.ntot();i++) {
            const double shere = -np->energy[Ex][i]*fac;
            if (shere < 0.0) np->s[Ex][i] = 0;
            else np->s[Ex][i] = shere;
          }
        else for (int i=0;i<v.ntot();i++)
          //np->s[Ex][i] = (np->energy[Ex][i] - f[Ex][0][i]*np->P[Ex][0][i]/(8*pi))*fac;
          np->s[Ex][i] = np->energy[Ex][i]*fac;
      } else {
        printf("I don't yet support saturation in this dimension.\n");
        exit(1);
      }
    }
    if (op->next && np->next) update_polarization_saturation(op->next, np->next);
  }
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
    const double funinv = 1.0/(1+0.5*g);
    DOCMP {
      for (int cc=0;cc<10;cc++)
        if (v.has_field((component)cc) && is_electric((component)cc)) {
          for (int i=0;i<v.ntot();i++)
            op->P[cc][cmp][i] = funinv*((2-om*om)*np->P[cc][cmp][i]+
                                        (0.5*g-1)*op->P[cc][cmp][i])+
              np->s[cc][i]*f[cc][cmp][i];
          if (f_pml[cc][cmp])
            for (int i=0;i<v.ntot();i++)
              op->P_pml[cc][cmp][i] = funinv*((2-om*om)*np->P_pml[cc][cmp][i]+
                                              (0.5*g-1)*op->P_pml[cc][cmp][i])+
                np->s[cc][i]*f_pml[cc][cmp][i];
        }
    }
    if (op->next && np->next) step_polarization_itself(op->next, np->next);
  }
}

void fields::step_e_polarization(polarization *op, polarization *np) {
  if (op == NULL && np == NULL && olpol != NULL && pol != NULL) {
    // This is the initial call... so I should start running from olpol and pol.
    step_e_polarization(olpol, pol);
  } else if (olpol != NULL && pol != NULL) {
    DOCMP {
      for (int cc=0;cc<10;cc++)
        if (op->P[cc][cmp]) {
          for (int i=0;i<v.ntot();i++)
            f[cc][cmp][i] -= ma->inveps[cc][i]*(np->P[cc][cmp][i]-op->P[cc][cmp][i]);
          if (f_pml[cc][cmp])
            for (int i=0;i<v.ntot();i++)
              f_pml[cc][cmp][i] -=
                ma->inveps[cc][i]*(np->P_pml[cc][cmp][i]-op->P_pml[cc][cmp][i]);
        }
    }
    if (op->next && np->next) step_e_polarization(op->next, np->next);
  }
}
