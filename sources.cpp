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

complex<double> src::get_amplitude_at_time(double time) const {
  double envelope = get_envelope_at_time(time);
  if (envelope == 0.0)
    return 0.0;
  double tt = time - peaktime;
  return (polar(1.0,-2*pi*freq*tt) - amp_shift)*envelope;
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

static double integrate_envelope(const src *s, const double inva) {
  if (s == NULL) abort("Bad arg to integrate_envelope!\n");
  double sofar = 0.0;
  for (int t=(int)((s->peaktime-s->cutoff)/inva);t<(1<<30);t++) {
    double e = s->get_envelope_at_time(t*inva);
    sofar += e;
    if (e == 0) break; // Bug here if there is a source that starts late,
                       // or a source that never stops.
  }
  return sofar*inva;
}

static complex<double> integrate_source(const src *s, const double inva) {
  if (s == NULL) abort("Bad arg to integrate_source!\n");
  complex<double> sofar = 0.0;
  for (int t=(int)((s->peaktime-s->cutoff)/inva);t<(1<<30);t++) {
    complex<double> A = s->get_amplitude_at_time(t*inva);
    sofar += A;
    if (A == 0.0) break; // Bug here if there is a source that starts late,
                         // or a source that never stops.
  }
  return sofar*inva;
}

void fields::add_point_source(component whichf, double freq,
                              double width, double peaktime,
                              double cutoff, const vec &p,
                              complex<double> amp, int is_c) {
  const double invmul = 1.0/S.multiplicity();
  int need_to_connect = 0;
  for (int sn=0;sn<S.multiplicity();sn++) {
    component cc = S.transform(whichf,sn);
    complex<double> ph = S.phase_shift(whichf,sn);
    const vec pp = S.transform(p,sn);
    for (int i=0;i<num_chunks;i++)
      if (chunks[i]->is_mine())
        need_to_connect += 
          chunks[i]->add_point_source(cc, freq, width, peaktime,
                                      cutoff, pp, invmul*ph*amp, is_c,
                                      time());
  }
  if (sum_to_all(need_to_connect)) connect_chunks();
}

int fields_chunk::add_point_source(component whichf, double freq,
                                   double width, double peaktime,
                                   double cutoff, const vec &p,
                                   complex<double> amp, int is_c, double tim) {
  if (p.dim != v.dim)
    abort("Error:  source doesn't have right dimensions!\n");
  if (!v.has_field(whichf))
    abort("Error:  source component %s is invalid.\n", component_name(whichf));
  // Allocate fields if they haven't already been allocated:
  int need_reconnection = 0;
  if (v.dim == D2 && !f[whichf][0]) {
    switch (whichf) {
    case Ex: case Ey: case Hz:
      alloc_f(Ex); alloc_f(Ey); alloc_f(Hz); break;
    case Hx: case Hy: case Ez:
      alloc_f(Hx); alloc_f(Hy); alloc_f(Ez); break;
    }
    need_reconnection = 1;
  }
  int ind[8];
  double w[8];
  v.interpolate(whichf, p, ind, w);
  if (w[0] == 0.0) return need_reconnection; // No source here...
  double prefac = 1.0;
  switch (v.dim) {
  case Dcyl: prefac = a; break;
  case D2: prefac = a; break; // FIXME: verify that this works right.
  case D1: prefac = 1; break;
  }
  for (int i=0;i<8 && w[i];i++)
    add_indexed_source(whichf, freq, width, peaktime, cutoff, ind[i],
                       amp*prefac*w[i], is_c, tim);
  return need_reconnection;
}

void fields::add_plane_source(double freq, double width, double peaktime,
                              double cutoff, double envelope (const vec &),
                              const vec &p, const vec &norm,
                              int is_c) {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->add_plane_source(freq, width, peaktime,
                                  cutoff, envelope, p, norm, is_c, time());
}

void fields_chunk::add_plane_source(double freq, double width, double peaktime,
                                    double cutoff, double envelope (const vec &),
                                    const vec &p, const vec &norm,
                                    int is_c, double time) {
  if (v.dim == Dcyl) {
    // We ignore norm in this case...
    if (m != 1) abort("Can only use plane source with m == 1!\n");
    const complex<double> I = complex<double>(0,1);
    const double z = p.z();
    const double eps = sqrt(ma->eps[(int)(z+0.5)]);
    for (int ir=0;ir<v.nr();ir++) {
      {
        const double r = ir*inva;
        // E_phi
        add_point_source(Ep, freq, width, peaktime, cutoff, vec(r,z),
                         envelope(vec(r,z)), is_c, time);
        // iH_r = d(rH_phi)/dr
        const double slope = ((r+0.5)*envelope(vec(r+0.5*inva,z)) -
                              (r-0.5)*envelope(vec(r-0.5*inva,z)))*a;
        add_point_source(Hr, freq, width, peaktime, cutoff, vec(r,z),
                         -eps*slope, is_c, time);
      }
      {
        const double r = (ir+0.5)*inva;
        const double sc = (ir == 0)?0.5:1.0;
        // iE_r = d(rE_phi)/dr
        const double slope = ((r+0.5)*envelope(vec(r+0.5*inva,z)) -
                              (r-0.5)*envelope(vec(r-0.5*inva,z)))*a;
        add_point_source(Er, freq, width, peaktime, cutoff, vec(r,z),
                         -I*sc*slope, is_c, time);
        // H_phi
        add_point_source(Hp, freq, width, peaktime, cutoff, vec(r,z),
                         -I*eps*sc*envelope(vec(r,z)), is_c, time);
      }
    }
  } else if (v.dim == D1) {
    const double z = p.z();
    const double eps = sqrt(ma->eps[(int)(z+0.5)]);
    add_point_source(Ex, freq, width, peaktime, cutoff, vec(z),
                     envelope(vec(z)), is_c, time);
    add_point_source(Hy, freq, width, peaktime, cutoff, vec(z),
                     envelope(vec(z))*eps, is_c, time);
  } else {
    abort("Can't use plane source in this number of dimensions.\n");
  }
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
  tmp.amp_shift = 0.0;
  tmp.is_continuous = is_c;
  tmp.cutoff = inva+ cutoff*tmp.width;
  while (exp(-tmp.cutoff*tmp.cutoff/(2*tmp.width*tmp.width)) == 0.0)
    tmp.cutoff *= 0.9;
  tmp.peaktime = peaktime;
  if (peaktime <= 0.0) tmp.peaktime = time+tmp.cutoff;
  // Apply a shift so that we won't end up with a static polarization when
  // the source is gone:  (FIXME: is there a bug here?)
  if (is_c) tmp.amp_shift = 0.0;
  else tmp.amp_shift = integrate_source(&tmp, inva)/integrate_envelope(&tmp, inva);
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
