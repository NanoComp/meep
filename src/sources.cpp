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

namespace meep {

complex<double> src::get_dPdt_at_time(double time, double dt) const {
  return (get_dipole_at_time(time) - get_dipole_at_time(time - dt))/dt;
}

complex<double> src::get_dipole_at_time(double time) const {
  double envelope = get_envelope_at_time(time);
  if (envelope == 0.0)
    return 0.0;
  double tt = time - peaktime;
  return polar(1.0,-2*pi*freq*tt)*envelope;
}

void src::update_dipole(double time) {
  pol_now = get_dipole_at_time(time);
}

double src::get_envelope_at_time(double time) const {
  double tt = time - peaktime;
  if (is_continuous && tt > 0) {
    return 1.0;
  } else if (fabs(tt) > cutoff) {
    return 0.0;
  } else {
    return exp(-tt*tt/(2*width*width));
    return 2.0/(exp(-tt/width) + exp(tt/width));
  }
}

src::src() {
  next = NULL;
}

src::~src() {
  delete next;
}

void fields::add_point_source(component whichf, double freq,
                              double width, double peaktime,
                              double cutoff, const vec &p,
                              complex<double> amp, int is_c) {
  ivec ilocs[8];
  double w[8];
  v.interpolate(whichf, p, ilocs, w);
  for (int argh=0;argh<8&&w[argh];argh++)
    add_point_source(whichf, freq, width, peaktime,
                     cutoff, ilocs[argh], w[argh]*amp, is_c);
}

void fields::add_point_source(component whichf, double freq,
                              double width, double peaktime,
                              double cutoff, const ivec &p,
                              complex<double> amp, int is_c) {
  const double invmul = 1.0/S.multiplicity();
  int need_to_connect = 0;
  ivec iloc = p;
  complex<double> kphase = 1.0;
  locate_point_in_user_volume(&iloc, &kphase);
  for (int sn=0;sn<S.multiplicity();sn++) {
    component cc = S.transform(whichf,sn);
    complex<double> ph = S.phase_shift(whichf,sn);
    const ivec pp = S.transform(iloc,sn);
    for (int i=0;i<num_chunks;i++)
      if (chunks[i]->is_mine())
        need_to_connect += 
          chunks[i]->add_point_source(cc, freq, width, peaktime,
                                      cutoff, pp, invmul*ph*amp/kphase, is_c,
                                      time());
  }
  if (sum_to_all(need_to_connect)) connect_chunks();
}

int fields_chunk::add_point_source(component whichf, double freq,
                                   double width, double peaktime,
                                   double cutoff, const ivec &p,
                                   complex<double> amp, int is_c, double tim) {
  if (p.dim != v.dim)
    abort("Error:  source doesn't have right dimensions! %s %s\n",
	  dimension_name(p.dim), dimension_name(v.dim));
  if (!v.has_field(whichf))
    abort("Error:  source component %s is invalid.\n", component_name(whichf));
  // Allocate fields if they haven't already been allocated:
  int need_reconnection = 0;
  if (!f[whichf][0]) {
    alloc_f(whichf);
    need_reconnection = 1;
  }
  if (amp == 0.0) return need_reconnection; // No source here...
  double prefac = 1.0;
  switch (v.dim) {
  case Dcyl: prefac = a; break;
  case D3: prefac = a*sqrt(a); break;
  case D2: prefac = a; break; // FIXME: verify that this works right.
  case D1: prefac = 1; break;
  }
  if (v.owns(p))
    add_indexed_source(whichf, freq, width, peaktime, cutoff, v.index(whichf,p),
                       amp*prefac, is_c, tim);
  return need_reconnection;
}

void fields_chunk::add_indexed_source(component whichf, double freq, double width,
                                      double peaktime, double cutoff, int theindex, 
                                      complex<double> amp, int is_c, double time) {
  if (theindex >= v.ntot() || theindex < 0)
    abort("Error:  source is outside of cell! (%d)\n", theindex);
  src tmp;
  tmp.freq = freq;
  tmp.width = width/tmp.freq; // this is now time width
  for (int com=0;com<10;com++) tmp.A[com] = 0;
  tmp.A[whichf] = amp;
  tmp.i = theindex;
  tmp.is_continuous = is_c;
  tmp.cutoff = inva+ cutoff*tmp.width;
  while (exp(-tmp.cutoff*tmp.cutoff/(2*tmp.width*tmp.width)) == 0.0)
    tmp.cutoff *= 0.9;
  tmp.peaktime = peaktime;
  if (peaktime <= 0.0) tmp.peaktime = time+tmp.cutoff;
  // Apply a shift so that we won't end up with a static polarization when
  // the source is gone:  (FIXME: is there a bug here?)
  if (is_magnetic(whichf)) {
    h_sources = tmp.add_to(h_sources);
  } else {
    e_sources = tmp.add_to(e_sources);
  }
}

src *src::add_to(src *others) const {
  if (!others) {
    src *t = new src(*this);
    t->next = NULL;
    return t;
  }
  if (others->i == i &&
      others->is_continuous == is_continuous &&
      others->cutoff == cutoff &&
      others->peaktime == peaktime) {
    for (int com=0;com<10;com++)
      others->A[com] += A[com];
    return others;
  } else {
    others->next = add_to(others->next);
    return others;
  }
}

} // namespace meep
