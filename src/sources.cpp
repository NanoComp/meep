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

/*********************************************************************/

src_time *src_time::add_to(src_time *others) const
{
  src_time *t = clone();
  t->next = others;
  return t;
}

double src_time::last_time_max(double after)
{
  after = max(last_time(), after);
  if (next)
    return next->last_time_max(after);
  else
    return after;
}

gaussian_src_time::gaussian_src_time(double f, double w, double st, double et)
{
  freq = f;
  width = w;
  peak_time = 0.5 * (st + et);
  cutoff = (et - st) * 0.5;

  // TODO: why bother with this?
  while (exp(-cutoff*cutoff / (2*width*width)) == 0.0)
    cutoff *= 0.9;
}

complex<double> gaussian_src_time::current(double time) const
{
  double tt = time - peak_time;
  if (fabs(tt) > cutoff)
    return 0.0;
  return exp(-tt*tt / (2*width*width)) * polar(1.0, -2*pi*freq*tt);
}

continuous_src_time::continuous_src_time(double f, double w, double st, double et, double s)
{
  freq = f;
  width = w == 0.0 ? 1e-20 : w; // hack to prevent NaN in current(t), below
  start_time = st;
  end_time = et;
  slowness = s;
}

complex<double> continuous_src_time::current(double time) const
{
  if (time < start_time || time > end_time)
    return 0.0;

  double ts = (time - start_time) / width - slowness;
  double te = (end_time - time) / width - slowness;

  return polar(1.0, -2*pi*freq*time) 
    * (1.0 + tanh(ts))  // goes from 0 to 2
    * (1.0 + tanh(te))  // goes from 2 to 0
    * 0.25;
}

/*********************************************************************/

src_pt *src_pt::add_to(src_pt *others) const {
  src_pt *t = new src_pt(*this);
  t->next = others;
  return t;
}

/*********************************************************************/

void fields::add_point_source(component whichf, double freq,
                              double width, double peaktime,
                              double cutoff, const vec &p,
                              complex<double> amp, int is_c) {
  width /= freq;

  if (is_c) { // TODO: don't ignore peaktime?
    continuous_src_time src(freq, width, time(), infinity, cutoff);
    add_point_source(whichf, src, p, amp);
  }
  else {
    cutoff = inva + cutoff * width;
    if (peaktime <= 0.0)
      peaktime = time() + cutoff;
  
    gaussian_src_time src(freq, width,
			     peaktime - cutoff, peaktime + cutoff);
    add_point_source(whichf, src, p, amp);
  }
}

void fields::add_point_source(component whichf, const src_time &src,
			      const vec &p, complex<double> amp) {
  ivec ilocs[8];
  double w[8];
  v.interpolate(whichf, p, ilocs, w);
  for (int argh=0;argh<8&&w[argh];argh++)
    add_point_source(whichf, src, ilocs[argh], w[argh]*amp);
}

void fields::add_point_source(component whichf, const src_time &src,
                              const ivec &p, complex<double> amp) {
  int added_src = 0;
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
      if (chunks[i]->is_mine()) {
	if (!added_src) {
	  sources = src.add_to(sources);
	  added_src = 1;
	}
        need_to_connect += 
	    chunks[i]->add_point_source(cc, sources, pp, invmul*ph*amp/kphase);
      }
  }
  if (sum_to_all(need_to_connect)) connect_chunks();
}

int fields_chunk::add_point_source(component whichf, src_time *src,
                                   const ivec &p, complex<double> amp) {
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
    add_indexed_source(whichf, src, v.index(whichf,p), amp*prefac);
  return need_reconnection;
}

void fields_chunk::add_indexed_source(component whichf, src_time *src,
                                      int theindex, complex<double> amp) {
  if (theindex >= v.ntot() || theindex < 0)
    abort("Error:  source is outside of cell! (%d)\n", theindex);
  src_pt tmp(src);
  tmp.A = amp * (inva * c); // multiply by dt for time-stepping
  tmp.c = whichf;
  tmp.i = theindex;
  if (is_magnetic(whichf)) {
    h_sources = tmp.add_to(h_sources);
  } else {
    e_sources = tmp.add_to(e_sources);
  }
}

} // namespace meep
