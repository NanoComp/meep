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

#include "meep.h"
#include "meep_internals.h"

polarization *polarization::set_up_polarizations(const mat_chunk *ma, int is_r) {
  if (ma->pb == NULL) return NULL;
  return new polarization(ma->pb, is_r);
}

void polarization::use_real_fields() {
  is_real = 1;
  FOR_COMPONENTS(c) delete[] P[c][1];
  FOR_COMPONENTS(c) delete[] P_p_pml[c][1];
  FOR_COMPONENTS(c) delete[] P_m_pml[c][1];
  FOR_COMPONENTS(c) P[c][1] = NULL;
  FOR_COMPONENTS(c) P_p_pml[c][1] = NULL;
  FOR_COMPONENTS(c) P_m_pml[c][1] = NULL;
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
        // FIXME perhaps shouldn't allocate the PML split fields_chunk if we don't
        // have pml...
        P_p_pml[c][cmp] = new double[v.ntot()];
        P_m_pml[c][cmp] = new double[v.ntot()];
        if (P_m_pml[c][cmp] == NULL)
          abort("Allocation error in polarization!\n");
        for (int i=0;i<v.ntot();i++) P_p_pml[c][cmp][i] = 0.0;
        for (int i=0;i<v.ntot();i++) P_m_pml[c][cmp][i] = 0.0;
      } else {
        P[c][cmp] = NULL;
        P_p_pml[c][cmp] = NULL;
        P_m_pml[c][cmp] = NULL;
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
    FOR_COMPONENTS(c) delete[] P[c][cmp];
    FOR_COMPONENTS(c) delete[] P_p_pml[c][cmp];
    FOR_COMPONENTS(c) delete[] P_m_pml[c][cmp];
  }
  for (int c=0;c<10;c++) delete[] energy[c];
  if (saturation_factor != 0.0)
    for (int c=0;c<10;c++) delete[] s[c];
  if (next) delete next;
}

double polarization::total_energy(const geometric_volume &what) {
  const volume v = pb->v;
  double e = 0.0;
  FOR_ELECTRIC_COMPONENTS(c)
    if (energy[c])
      for (int i=0;i<v.ntot();i++)
        e += what.intersect_with(v.dV(c,i)).full_volume()*energy[c][i];
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

polarizability::polarizability(const mat_chunk *ma, double sig(const vec &),
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
  if (sigma == NULL) abort("Out of memory in polarizability!\n");

  for (int i=0;i<v.ntot();i++)
    sigma[i] = sigscale*sig(v.loc(v.eps_component(),i));
  for (int c=0;c<10;c++) if (s[c])
    for (int i=0;i<v.ntot();i++) s[c][i] = 0.0;
  // Average out sigma over the grid...
  if (v.dim == Dcyl) {
    const vec dr = v.dr()*0.5; // The distance between Yee field components
    const vec dz = v.dz()*0.5; // The distance between Yee field components
    for (int i=0;i<v.ntot();i++) {
      const vec here = v.loc(Ep,i);
      s[Er][i] = 0.5*sigscale*(sig(here+dr+dz) + sig(here+dr-dz));
      s[Ep][i] = 0.25*sigscale*(sig(here+dr+dz) + sig(here-dr+dz) +
                                sig(here+dr-dz) + sig(here-dr-dz));
      s[Ez][i] = 0.5*sigscale*(sig(here+dr+dz) + sig(here-dr+dz));
    }
  } else if (v.dim == D1) {
    // There's just one field point...
    for (int i=0;i<v.ntot();i++) s[Ex][i] = sigma[i];
  } else {
    abort("Unsupported dimensionality!\n");
  }
}

polarizability::~polarizability() {
  for (int c=0;c<10;c++) delete[] s[c];
  delete[] sigma;
}

complex<double> polarization::analytic_epsilon(double freq, const vec &p) const {
  const complex<double> I = complex<double>(0,1);
  double w[8];
  int in[8];
  pb->v.interpolate(pb->v.eps_component(), p, in, w);
  complex<double> epsi = 0.0;
  for (int i=0;i<8 && w[i];i++)
    epsi += w[i]*pb->sigma[in[i]]/
      (pb->omeganot*pb->omeganot - freq*freq - freq*pb->gamma*I);
  if (next) epsi += next->analytic_epsilon(freq, p);
  return epsi;
}

complex<double> fields::analytic_epsilon(double f, const vec &p) const {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->v.contains(p))
      return chunks[i]->analytic_epsilon(f,p);
  return 0.0;
}

complex<double> fields_chunk::analytic_epsilon(double f, const vec &p) const {
  complex<double> epsi = 0.0;
  if (is_mine()) {
    const double freq_conversion = 2*pi*c/a;
    double freq = f*freq_conversion;
    const component c = v.eps_component();
    int in[8];
    double w[8];
    v.interpolate(c,p,in,w);
    for (int i=0;i<8 && w[i];i++)
      epsi += ma->eps[in[i]]*w[i];
    if (pol) epsi += pol->analytic_epsilon(freq, p);
  }
  return broadcast(n_proc(), epsi);
}

void mat::add_polarizability(double sigma(const vec &), double omega, double gamma,
                             double delta_epsilon, double energy_saturation) {
  for (int i=0;i<num_chunks;i++)
    chunks[i]->add_polarizability(sigma, omega, gamma, delta_epsilon, energy_saturation);
}

void mat_chunk::add_polarizability(double sigma(const vec &),
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

void fields::initialize_polarizations() {
  for (int i=0;i<num_chunks;i++)
    chunks[i]->initialize_polarizations();
}

void fields_chunk::initialize_polarizations(polarization *op, polarization *np) {
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

void fields::prepare_step_polarization_energy() {
  for (int i=0;i<num_chunks;i++)
    chunks[i]->prepare_step_polarization_energy();
}

void fields_chunk::prepare_step_polarization_energy(polarization *op, polarization *np) {
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

void fields::half_step_polarization_energy() {
  for (int i=0;i<num_chunks;i++)
    chunks[i]->half_step_polarization_energy();
}

void fields_chunk::half_step_polarization_energy(polarization *op, polarization *np) {
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

void fields::update_polarization_saturation() {
  for (int i=0;i<num_chunks;i++)
    chunks[i]->update_polarization_saturation();
}

void fields_chunk::update_polarization_saturation(polarization *op, polarization *np) {
  if (op == NULL && np == NULL && olpol != NULL && pol != NULL) {
    // This is the initial call... so I should start running from olpol and pol.
    update_polarization_saturation(olpol, pol);
  } else if (olpol != NULL && pol != NULL) {
    if (np->saturation_factor != 0.0) {
      const volume v = np->pb->v;
      const double fac = np->saturation_factor;
      if (v.dim == D1) {
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
        abort("I don't yet support saturation in this dimension.\n");
      }
    }
    if (op->next && np->next) update_polarization_saturation(op->next, np->next);
  }
}

void fields::step_polarization_itself() {
  for (int i=0;i<num_chunks;i++)
    chunks[i]->step_polarization_itself();
}

void fields_chunk::step_polarization_itself(polarization *op, polarization *np) {
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
          if (f_p_pml[cc][cmp])
            for (int i=0;i<v.ntot();i++)
              op->P_p_pml[cc][cmp][i] = funinv*((2-om*om)*np->P_p_pml[cc][cmp][i]+
                                              (0.5*g-1)*op->P_p_pml[cc][cmp][i])+
                np->s[cc][i]*f_p_pml[cc][cmp][i];
          if (f_m_pml[cc][cmp])
            for (int i=0;i<v.ntot();i++)
              op->P_m_pml[cc][cmp][i] = funinv*((2-om*om)*np->P_m_pml[cc][cmp][i]+
                                              (0.5*g-1)*op->P_m_pml[cc][cmp][i])+
                np->s[cc][i]*f_m_pml[cc][cmp][i];
        }
    }
    if (op->next && np->next) step_polarization_itself(op->next, np->next);
  }
}

void fields::step_e_polarization() {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->step_e_polarization();
}

void fields_chunk::step_e_polarization(polarization *op, polarization *np) {
  if (op == NULL && np == NULL && olpol != NULL && pol != NULL) {
    // This is the initial call... so I should start running from olpol and pol.
    step_e_polarization(olpol, pol);
  } else if (olpol != NULL && pol != NULL) {
    DOCMP {
      FOR_ELECTRIC_COMPONENTS(cc)
        if (op->P[cc][cmp] && f[cc][cmp]) {
          for (int i=0;i<v.ntot();i++)
            f[cc][cmp][i] -= ma->inveps[cc][component_direction(cc)][i]*
              (np->P[cc][cmp][i]-op->P[cc][cmp][i]);
          if (f_p_pml[cc][cmp])
            for (int i=0;i<v.ntot();i++)
              f_p_pml[cc][cmp][i] -=
                ma->inveps[cc][component_direction(cc)][i]*
                (np->P_p_pml[cc][cmp][i]-op->P_p_pml[cc][cmp][i]);
          if (f_m_pml[cc][cmp])
            for (int i=0;i<v.ntot();i++)
              f_m_pml[cc][cmp][i] -=
                ma->inveps[cc][component_direction(cc)][i]*
                (np->P_m_pml[cc][cmp][i]-op->P_m_pml[cc][cmp][i]);
        }
    }
    if (op->next && np->next) step_e_polarization(op->next, np->next);
  }
}
