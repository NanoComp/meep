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

complex<double> src::get_amplitude_at_time(int t) const {
  double envelope = get_envelope_at_time(t);
  if (envelope == 0.0)
    return 0.0;
  double tt = t - peaktime;
  return (polar(1.0,-2*pi*freq*tt) - amp_shift)*envelope;
}

double src::get_envelope_at_time(int t) const {
  double tt = t - peaktime;
  if (is_continuous && tt > 0) {
    return 1.0;
  } else if (fabs(tt) > cutoff) {
    return 0.0;
  } else {
    return exp(-tt*tt/(2*width*width));
    return 2.0/(exp(-tt/width) + exp(tt/width));
  }
}

src::~src() {
  delete next;
}

static double integrate_envelope(const src *s) {
  if (s == NULL) {
    printf("Bad arg to integrate_envelope!\n");
    exit(1);
  }
  double sofar = 0.0;
  for (int t=(int)s->peaktime-s->cutoff;t<(1<<30);t++) {
    double e = s->get_envelope_at_time(t);
    sofar += e;
    if (e == 0) break; // Bug here if there is a source that starts late,
                       // or a source that never stops.
  }
  return sofar;
}

static complex<double> integrate_source(const src *s) {
  if (s == NULL) {
    printf("Bad arg to integrate_source!\n");
    exit(1);
  }
  complex<double> sofar = 0.0;
  for (int t=0;1<<30;t++) {
    complex<double> A = s->get_amplitude_at_time(t);
    sofar += A;
    if (A == 0) break; // Bug here if there is a source that starts late,
                       // or a source that never stops.
  }
  return sofar;
}

void fields::add_point_source(component whichf, double freq,
                              double width, double peaktime,
                              double cutoff, const vec &p,
                              complex<double> amp, int is_c) {
  for (int i=0;i<num_chunks;i++)
    chunks[i]->add_point_source(whichf, freq, width, peaktime,
                                cutoff, p, amp, is_c);
}

void fields_chunk::add_point_source(component whichf, double freq,
                                    double width, double peaktime,
                                    double cutoff, const vec &p,
                                    complex<double> amp, int is_c) {
  // FIXME this really should call an interpolation routine...
  if (p.dim != v.dim) {
    printf("Error:  source doesn't have right dimensions!\n");
    exit(1);
  } else if (!v.has_field(whichf)) {
    printf("Error:  source component %s is invalid.\n", component_name(whichf));
    exit(1);
  }
  int ind[8];
  double w[8];
  v.interpolate(whichf, p, ind, w);
  double prefac = 1.0;
  switch (v.dim) {
  case dcyl: prefac = a; break;
  case d1: prefac = 1; break;
  }
  for (int i=0;i<8 && w[i];i++)
    add_indexed_source(whichf, freq, width, peaktime, cutoff, ind[i], amp*prefac*w[i], is_c);
}

void fields::add_plane_source(double freq, double width, double peaktime,
                              double cutoff, double envelope (const vec &),
                              const vec &p, const vec &norm,
                              int is_c) {
  for (int i=0;i<num_chunks;i++)
    chunks[i]->add_plane_source(freq, width, peaktime,
                                cutoff, envelope, p, norm, is_c);
}

void fields_chunk::add_plane_source(double freq, double width, double peaktime,
                                    double cutoff, double envelope (const vec &),
                                    const vec &p, const vec &norm,
                                    int is_c) {
  if (v.dim == dcyl) {
    // We ignore norm in this case...
    if (m != 1) {
      printf("Can only use plane source with m == 1!\n");
      exit(1);
    }
    const complex<double> I = complex<double>(0,1);
    const double z = p.z();
    const double eps = sqrt(ma->eps[(int)(z+0.5)]);
    for (int ir=0;ir<v.nr();ir++) {
      {
        const double r = ir*inva;
        // E_phi
        add_point_source(Ep, freq, width, peaktime, cutoff, vec(r,z),
                         envelope(vec(r,z)), is_c);        
        // iH_r = d(rH_phi)/dr
        const double slope = ((r+0.5)*envelope(vec(r+0.5*inva,z)) -
                              (r-0.5)*envelope(vec(r-0.5*inva,z)))*a;
        add_point_source(Hr, freq, width, peaktime, cutoff, vec(r,z), -eps*slope, is_c);
      }
      {
        const double r = (ir+0.5)*inva;
        const double sc = (ir == 0)?0.5:1.0;
        // iE_r = d(rE_phi)/dr
        const double slope = ((r+0.5)*envelope(vec(r+0.5*inva,z)) -
                              (r-0.5)*envelope(vec(r-0.5*inva,z)))*a;
        add_point_source(Er, freq, width, peaktime, cutoff, vec(r,z), -I*sc*slope, is_c);
        // H_phi
        add_point_source(Hp, freq, width, peaktime, cutoff, vec(r,z),
                         -I*eps*sc*envelope(vec(r,z)), is_c);
      }
    }
  } else if (v.dim == d1) {
    const double z = p.z();
    const double eps = sqrt(ma->eps[(int)(z+0.5)]);
    add_point_source(Ex, freq, width, peaktime, cutoff, vec(z), envelope(vec(z)), is_c);
    add_point_source(Hy, freq, width, peaktime, cutoff, vec(z), envelope(vec(z))*eps, is_c);
  } else {
    printf("Can't use plane source in this number of dimensions.\n");
    exit(1);
  }
}

void fields_chunk::add_indexed_source(component whichf, double freq, double width,
                                double peaktime, int cutoff, int theindex, 
                                complex<double> amp, int is_c) {
  if (theindex >= v.ntot() || theindex < 0) {
    printf("Error:  source is outside of cell! (%d)\n", theindex);
    exit(1);
  }
  src *tmp = new src;
  tmp->freq = freq*c*inva;
  tmp->width = width/tmp->freq; // this is now time width
  for (int com=0;com<10;com++) tmp->A[com] = 0;
  tmp->A[whichf] = amp;
  tmp->i = theindex;
  tmp->amp_shift = 0.0;
  tmp->is_continuous = is_c;
  if (is_magnetic(whichf)) {
    tmp->next = h_sources;
    h_sources = tmp;
  } else {
    tmp->next = e_sources;
    e_sources = tmp;
  }
  tmp->cutoff = 1+ (int)(cutoff*tmp->width);
  tmp->peaktime = peaktime*a/c;
  if (peaktime <= 0.0) tmp->peaktime = t+tmp->cutoff;
  // Apply a shift so that we won't end up with a static polarization when
  // the source is gone:
  if (is_c) tmp->amp_shift = 0.0;
  else tmp->amp_shift = integrate_source(tmp)/integrate_envelope(tmp);
}
