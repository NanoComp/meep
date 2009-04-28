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
  FOR_COMPONENTS(c) s[c] = pb->s[c];
  next = NULL;
  zero_fields();
}

polarization::~polarization() {
  DOCMP FOR_COMPONENTS(c) delete[] P[c][cmp];
  FOR_COMPONENTS(c) delete[] energy[c];
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
			       double sigscale, bool mine) {
  v = sc->v;
  is_it_mine = mine;
  ft = ft_;
  omeganot = om;
  gamma = ga;
  next = NULL;

  sig.set_volume(sc->v.pad().surroundings());
  FOR_COMPONENTS(c) s[c] = NULL;
  if (is_mine()) {
    FOR_FT_COMPONENTS(ft,c) if (v.has_field(c))
      s[c] = new realnum[v.ntot()];

    FOR_COMPONENTS(c) if (s[c]) for (int i=0;i<v.ntot();i++) s[c][i] = 0.0;

    // TODO:  should we be doing some kind of subpixel averaging here?
    FOR_FT_COMPONENTS(ft,c) if (s[c]) {
      double sigrow[3];
      int ic = component_index(c);
      LOOP_OVER_VOL(v, c, i) {
	IVEC_LOOP_LOC(v, here);
	sig.sigma_row(c, sigrow, here);
	s[c][i] = sigscale*sigrow[ic];
	if (sigrow[(ic+1)%3] != 0.0 || sigrow[(ic+2)%3] != 0.0)
	  abort("non-diagonal polarizabilities are not yet supported");
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
						  double omega, double gamma) {
  changing_chunks();
  for (int i=0;i<num_chunks;i++)
    chunks[i]->add_polarizability(sigma, ft,omega, gamma);
  return chunks[0]->pb->get_identifier();
}

polarizability_identifier structure::add_polarizability(
				  double sigma(const vec &),
				  field_type ft, double omega, double gamma) {
  simple_material_function sig(sigma);
  return add_polarizability(sig, ft, omega, gamma);
}

polarizability_identifier polarizability::get_identifier() const {
  polarizability_identifier pi;
  pi.ft = ft;
  pi.gamma = gamma;
  pi.omeganot = omeganot;
  return pi;
}

bool polarizability_identifier::operator==(const polarizability_identifier &a) {
  return ft == a.ft && gamma == a.gamma && omeganot == a.omeganot;
}

void structure_chunk::add_polarizability(material_function &sigma,
					 field_type ft, double omega,
					 double gamma) {
  sigma.set_polarizability(ft, omega, gamma);
  const double freq_conversion = 2*pi*dt;
  double sigma_scale  = freq_conversion*freq_conversion*omega*omega; 
  polarizability *npb = new polarizability(this, sigma,
                                           ft, freq_conversion*omega,
                                           freq_conversion*gamma,
                                           sigma_scale, is_mine());
  npb->next = pb;
  pb = npb;
}

} // namespace meep
