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

// this function is necessary to make equality commutative ... ugh
bool src_times_equal(const src_time &t1, const src_time &t2)
{
     return t1.is_equal(t2) && t2.is_equal(t1);
}

src_time *src_time::add_to(src_time *others, src_time **added) const
{
     if (others) {
	  if (src_times_equal(*this, *others))
	       *added = others;
	  else
	       add_to(others->next, added);
	  return others;
     }
     else {
	  src_time *t = clone();
	  t->next = others;
	  *added = t;
	  return t;
     }
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

complex<double> gaussian_src_time::dipole(double time) const
{
  double tt = time - peak_time;
  if (fabs(tt) > cutoff)
    return 0.0;
  return exp(-tt*tt / (2*width*width)) * polar(1.0, -2*pi*freq*tt);
}

bool gaussian_src_time::is_equal(const src_time &t) const
{
     const gaussian_src_time *tp = dynamic_cast<const gaussian_src_time*>(&t);
     if (tp)
	  return(tp->freq == freq && tp->width == width &&
		 tp->peak_time == peak_time && tp->cutoff == cutoff);
     else
	  return 0;
}

continuous_src_time::continuous_src_time(double f, double w, double st, double et, double s)
{
  freq = f;
  width = w == 0.0 ? 1e-20 : w; // hack to prevent NaN in dipole(t), below
  start_time = st;
  end_time = et;
  slowness = s;
}

complex<double> continuous_src_time::dipole(double time) const
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

bool continuous_src_time::is_equal(const src_time &t) const
{
     const continuous_src_time *tp = 
	  dynamic_cast<const continuous_src_time*>(&t);
     if (tp)
	  return(tp->freq == freq && tp->width == width &&
		 tp->start_time == start_time && tp->end_time == end_time &&
		 tp->slowness == slowness);
     else
	  return 0;
}

/*********************************************************************/

src_vol::src_vol(component cc, src_time *st, int n, const int *ind, const complex<double> *amps) {
  c = cc;
  t = st; next = NULL;
  npts = n;
  index = new int[npts];
  A = new complex<double>[npts];
  for (int j=0; j<npts; j++) {
    index[j] = ind[j];
    A[j] = amps[j];
  }
}

src_vol::src_vol(const src_vol &sv) {
  c = sv.c;
  t = sv.t;
  npts = sv.npts;
  index = new int[npts];
  A = new complex<double>[npts];
  for (int j=0; j<npts; j++) {
    index[j] = sv.index[j];
    A[j] = sv.A[j];
  }
  if (sv.next)
    next = new src_vol(*sv.next);
  else
    next = NULL;
}

src_vol *src_vol::add_to(src_vol *others) const {
  if (others) {
    if (*this == *others) {
      if (npts != others->npts)
        abort("Cannot add volume sources with different number of points\n");
      /* Compare all of the indices...if this ever becomes too slow,
	 we can just compare the first and last indices. */
      for (int j=0; j<npts; j++) {
        if (others->index[j] != index[j])
          abort("Different indices\n");
        others->A[j] += A[j];
      }
    }
    else
      others->next = add_to(others->next);
    return others;
  }
  else {
    src_vol *sv = new src_vol(*this);
    sv->next = others;
    return sv;
  }
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

complex<double> one(const vec &v) {(void) v; return 1.0;}
void fields::add_point_source(component whichf, const src_time &src,
			      const vec &p, complex<double> amp) {
  add_volume_source(whichf, src, p, p, one, amp);
}

void fields::add_volume_source(component whichf, const src_time &src,
                               const vec &p1, const vec &p2,
                               complex<double> A(const vec &), complex<double> amp) {
  src_time *newsrc;
  sources = src.add_to(sources, &newsrc);

  amp = 8.0*(0.125*amp); // this puts three 0's at the end to prevent floating point errors below 
  const double invmul = 1.0/S.multiplicity();
  int need_to_connect = 0;
  
  LOOP_OVER_DIRECTIONS(v.dim, d) 
    if (abs((p1-p2).in_direction(d)) > 
        user_volume.boundary_location(High, d) - user_volume.boundary_location(Low, d)) 
      abort("Cannot accept source width larger than cell width in %s direction!\n", direction_name(d));
  vec newp1[8], newp2[8];
  complex<double> kphase[8];
  int ncopies;
  // possibly make copies of volume source if periodic boundary conditions are used
  locate_volume_source_in_user_volume(p1, p2, newp1, newp2, kphase, ncopies);
  
  for (int j=0; j<ncopies; j++)
    for (int sn=0;sn<S.multiplicity();sn++) {
    component cc = S.transform(whichf,sn);
    complex<double> ph = S.phase_shift(whichf,sn);
    const vec rotated1 = S.transform(newp1[j],sn);
    const vec rotated2 = S.transform(newp2[j],sn);
    for (int i=0;i<num_chunks;i++)
      if (chunks[i]->is_mine()) {
        need_to_connect += 
	    chunks[i]->add_volume_source(cc, newsrc, rotated1, rotated2, A, invmul*ph*amp/kphase[j], S, sn);
      }
    }
  if (sum_to_all(need_to_connect)) connect_chunks();
}

int fields_chunk::add_volume_source(component whichf, src_time *src, const vec &p1, const vec &p2,
                                    complex<double> A(const vec &), complex<double> amp, symmetry S, int sn) {
  // Allocate fields if they haven't already been allocated:
  int need_reconnection = 0;
  if (!f[whichf][0]) {
    alloc_f(whichf);
    need_reconnection = 1;
  }

  geometric_volume vbig = (v.pad()).surroundings();
  if (vbig.intersects(geometric_volume(p1, p2))) {
    geometric_volume Vinters = vbig.intersect_with(geometric_volume(p1, p2));
    int minindex = v.ntot(), maxindex = -1;
    int npts = 0;
    int nmax = 1;
    LOOP_OVER_DIRECTIONS(v.dim, d)
      nmax *= (2 + (int)ceil(v.a*abs((p1-p2).in_direction(d))));
    if (nmax > v.ntot()) nmax = v.ntot();
    int *index_array = new int[nmax];
    complex<double> *amps_array = new complex<double>[nmax];
    for (int index=0; index<v.ntot(); index++) {
      vec here = v.loc(whichf, index);
      double weight = 1.0;
      LOOP_OVER_DIRECTIONS(v.dim, d)
        if (p1.in_direction(d) == p2.in_direction(d)) { // delta function in this direction
          // using 1.0 below gives floating point errors in tests
          if (abs(here.in_direction(d) - p1.in_direction(d)) * v.a < 1.0 - 3e-15) 
            // multiply by v.a below to have radiated power indpendent of a
            weight *= (1 - v.a*abs(here.in_direction(d) - p1.in_direction(d))) * v.a;
          else
            weight = 0;
        }
        else { // finite thickness
          double dhere = abs((here - (p1+p2)*0.5).in_direction(d));
          double half_width = abs(((p1-p2)*0.5).in_direction(d));
          if (dhere < 0.5*inva + half_width) {
            if (dhere > half_width - 0.5*inva)
              weight *= 1 - v.a*(dhere - (half_width - 0.5*inva)); // border
            else
              weight *= 1.0; // completely inside
          }
          else
            weight = 0; // outside
        }
      if (weight > 0) {
        if (minindex > index) minindex = index;
        if (maxindex < index) maxindex = index;
        index_array[npts] = index;
        vec original_vec = S.transform(here, -sn); // inv(S,sn) * here
        amps_array[npts] = weight * amp * A(original_vec);  
        npts++;
      }
    }
    src_vol tmp(whichf, src, npts, index_array, amps_array);
    if (is_magnetic(whichf)) {
      h_sources = tmp.add_to(h_sources);
    } else {
      e_sources = tmp.add_to(e_sources);
    }
    delete[] index_array;
    delete[] amps_array;
  }  

  return need_reconnection;
}

} // namespace meep
