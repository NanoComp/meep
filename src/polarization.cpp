/* Copyright (C) 2005-2009 Massachusetts Institute of Technology
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

#include "meep.hpp"
#include "meep_internals.hpp"
#include "config.h"

namespace meep {

void polarization::set_up_polarizations(polarization *pols[NUM_FIELD_TYPES], const structure_chunk *sc, int is_r, bool store_enrgy) {
  if (sc->is_mine()) FOR_FIELD_TYPES(ft) {
    const polarizability *pb = sc->pb;
    while (pb && pb->ft != ft) pb = pb->next;
    if (pb) {
      pols[ft] = new polarization(pb, is_r, store_enrgy);
      polarization *pol = pols[ft];
      for (pb = pb->next; pb; pb = pb->next) if (pb->ft == ft)
	pol = (pol->next = new polarization(pb, is_r, store_enrgy));
    }
  }
}

void polarization::use_real_fields() {
  is_real = 1;
  FOR_COMPONENTS(c) {
    delete[] P[c][1];
    P[c][1] = NULL;
  }
  if (next) next->use_real_fields();
}

void polarization::zero_fields() {
  const volume &v = pb->v;
  DOCMP FOR_COMPONENTS(c) if (P[c][cmp])
    for (int i=0;i<v.ntot();i++) P[c][cmp][i] = 0;
  if (pb->energy_saturation != 0.0) {
    FOR_COMPONENTS(c) if (s[c] && pb->s[c])
      for (int i=0;i<v.ntot();i++) s[c][i] = pb->s[c][i];
    double num_components = 0.0;
    FOR_FT_COMPONENTS(pb->ft,c) if (s[c]) num_components += 1.0;
    const double isf = 1.0/saturation_factor/num_components;
    FOR_COMPONENTS(c) if (energy[c] && s[c])
      for (int i=0;i<v.ntot();i++) energy[c][i] = isf*s[c][i];
  }
  else
    FOR_COMPONENTS(c) if (energy[c])
      for (int i=0;i<v.ntot();i++) energy[c][i] = 0;
  if (next) next->zero_fields();
}

polarization::polarization(const polarizability *the_pb, 
			   int is_r, bool store_enrgy) {
  const volume &v = the_pb->v;
  is_real = is_r;
  store_energy = store_enrgy;
  DOCMP2 FOR_COMPONENTS(c) P[c][cmp] = NULL;
  DOCMP FOR_FT_COMPONENTS(the_pb->ft, c) if (v.has_field(c))
    P[c][cmp] = new realnum[v.ntot()];
  FOR_COMPONENTS(c) energy[c] = NULL;
  FOR_FT_COMPONENTS(the_pb->ft, c) if (v.has_field(c) && store_energy)
    energy[c] = new realnum[v.ntot()];
  pb = the_pb;
#ifndef WITH_SATURABLE_ABSORBERS
  if (pb->energy_saturation != 0.0)
    abort("saturable absorber, but not configured --with-saturable-absorbers");
#endif
  // Initialize the s[] arrays that point to sigma.
  FOR_COMPONENTS(c) {
    if (pb->energy_saturation != 0.0) {
      if (pb->s[c])
        s[c] = new realnum[v.ntot()];
      else 
	s[c] = NULL;
    } 
    else 
      s[c] = pb->s[c];
  }
  // Deal with saturation stuff.
  if (pb->energy_saturation != 0.0)
    saturation_factor = pb->saturated_sigma/pb->energy_saturation;
  else
    saturation_factor = 0.0;

  next = NULL;
  zero_fields();
}

polarization::~polarization() {
  DOCMP FOR_COMPONENTS(c) delete[] P[c][cmp];
  FOR_COMPONENTS(c) delete[] energy[c];
  if (pb->energy_saturation != 0.0)
    FOR_COMPONENTS(c) delete[] s[c];
  if (next) delete next;
}

double polarization::local_energy(const ivec &iloc) {
  if (pb->v.dim != D1) abort("Can't do local_energy in these dims.\n");
  double res = 0.0;
  FOR_COMPONENTS(c) if (energy[c])
    res += energy[c][pb->v.index(c,iloc)];
  return res;
}

polarizability::polarizability(const polarizability *pb) {
  omeganot = pb->omeganot;
  gamma = pb->gamma;
  v = pb->v;
  ft = pb->ft;
  energy_saturation = pb->energy_saturation;
  saturated_sigma = pb->saturated_sigma;
  is_it_mine = pb->is_it_mine;
  FOR_COMPONENTS(c) s[c] = NULL;
  if (is_mine()) {
    FOR_COMPONENTS(c) if (v.has_field(c) && pb->s[c]) {
      s[c] = new realnum[v.ntot()];
      for (int i=0;i<v.ntot();i++) s[c][i] = pb->s[c][i];
    }
  }
  if (pb->next) next = new polarizability(pb->next);
  else next = NULL;
}

polarizability::polarizability(const structure_chunk *sc, material_function &sig,
                               field_type ft_, double om, double ga, 
			       vec sigscale,
                               double energy_sat, bool mine) {
  v = sc->v;
  is_it_mine = mine;
  ft = ft_;
  omeganot = om;
  gamma = ga;
  next = NULL;
  energy_saturation = energy_sat;
  {
    saturated_sigma = 0.0;
    int numdims = 0;
    LOOP_OVER_DIRECTIONS(v.dim, d) {
      saturated_sigma += sigscale.in_direction(d);
      ++numdims;
    }
    saturated_sigma /= numdims;
  }

  sig.set_volume(sc->v.pad().surroundings());
  FOR_COMPONENTS(c) s[c] = NULL;
  if (is_mine()) {
    FOR_FT_COMPONENTS(ft,c) if (v.has_field(c))
      s[c] = new realnum[v.ntot()];

    FOR_COMPONENTS(c) if (s[c]) for (int i=0;i<v.ntot();i++) s[c][i] = 0.0;

    // TODO:  should we be doing some kind of subpixel averaging here?
    FOR_FT_COMPONENTS(ft,c) if (s[c]) {
      double sigscale_c = sigscale.in_direction(component_direction(c));
      LOOP_OVER_VOL(v, c, i) {
	IVEC_LOOP_LOC(v, here);
	s[c][i] = sigscale_c*sig.sigma(here);
      }
    }
  } else { // Not mine, don't store arrays...
    FOR_COMPONENTS(c) s[c] = 0;
  }
}

polarizability::~polarizability() {
  FOR_COMPONENTS(c) delete[] s[c];
}

complex<double> polarization::analytic_chi1(component c, double freq, const vec &p) const {
  const complex<double> I = complex<double>(0,1);
  double w[8];
  int in[8];
  pb->v.interpolate(c, p, in, w);
  complex<double> epsi = 0.0;
  if (pb->s[c])
    for (int i=0;i<8 && w[i];i++)
      epsi += w[i]*pb->s[c][in[i]]/
	(pb->omeganot*pb->omeganot - freq*freq - freq*pb->gamma*I);
  if (next) epsi += next->analytic_chi1(c, freq, p);
  return epsi;
}

complex<double> fields::analytic_chi1(component c, double f, const vec &p) const {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->v.contains(p))
      return chunks[i]->analytic_chi1(c,f,p) + get_eps(p);
  return 0.0;
}

  complex<double> fields_chunk::analytic_chi1(component c, double f, const vec &p) const {
  complex<double> epsi = 0.0;
  if (is_mine() && pols[type(c)])
    epsi += pols[type(c)]->analytic_chi1(c, f * (2*pi*dt), p);
  return broadcast(n_proc(), epsi);
}

polarizability_identifier structure::add_polarizability(material_function &sigma,
						  field_type ft,
						  double omega, double gamma,
                                                  double delta_epsilon,
                                                  double energy_saturation) {
  return add_polarizability(sigma, ft, omega, gamma,
			    one_vec(v.dim) * delta_epsilon, energy_saturation);
}

polarizability_identifier structure::add_polarizability(material_function &sigma,
						  field_type ft,
						  double omega, double gamma,
                                                  vec delta_epsilon,
                                                  double energy_saturation) {
  changing_chunks();
  for (int i=0;i<num_chunks;i++)
    chunks[i]->add_polarizability(sigma, ft,omega, gamma, delta_epsilon, energy_saturation);
  return chunks[0]->pb->get_identifier();
}

polarizability_identifier structure::add_polarizability(double sigma(const vec &),
							field_type ft,
                                                  double omega, double gamma,
                                                  double delta_epsilon,
                                                  double energy_saturation) {
  return add_polarizability(sigma, ft, omega, gamma,
			    one_vec(v.dim) * delta_epsilon, energy_saturation);
}

polarizability_identifier structure::add_polarizability(double sigma(const vec &),
							field_type ft,
                                                  double omega, double gamma,
                                                  vec delta_epsilon,
                                                  double energy_saturation) {
  simple_material_function sig(sigma);
  return add_polarizability(sig,ft,omega,gamma,delta_epsilon,energy_saturation);
}

polarizability_identifier polarizability::get_identifier() const {
  polarizability_identifier pi;
  pi.ft = ft;
  pi.gamma = gamma;
  pi.omeganot = omeganot;
  pi.energy_saturation = energy_saturation;
  pi.saturated_sigma = saturated_sigma;
  return pi;
}

bool polarizability_identifier::operator==(const polarizability_identifier &a) {
  return ft == a.ft && gamma == a.gamma && omeganot == a.omeganot &&
    energy_saturation == a.energy_saturation && saturated_sigma == a.saturated_sigma;
}

void structure_chunk::add_polarizability(material_function &sigma,
					 field_type ft, double omega,
					 double gamma, vec delta_epsilon,
					 double energy_sat) {
  sigma.set_polarizability(ft, omega, gamma, delta_epsilon, energy_sat);
  const double freq_conversion = 2*pi*dt;
  vec sigma_scale  = delta_epsilon * (freq_conversion*freq_conversion*omega*omega); 
  polarizability *npb = new polarizability(this, sigma,
                                           ft, freq_conversion*omega,
                                           freq_conversion*gamma,
                                           sigma_scale, energy_sat,
                                           is_mine());
  npb->next = pb;
  pb = npb;
}

void structure_chunk::add_polarizability(double sigma(const vec &),
	      field_type ft, double omega, double gamma, vec delta_epsilon,
              double energy_sat) {
  simple_material_function sig(sigma);
  add_polarizability(sig, ft, omega,gamma,delta_epsilon,energy_sat);
}

} // namespace meep
