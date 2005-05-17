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

#include "meep.hpp"
#include "meep_internals.hpp"

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

src_vol::src_vol(component cc, src_time *st, int n, int *ind, complex<double> *amps) {
  c = cc;
  if (is_D(c)) c = direction_component(Ex, component_direction(c));
  t = st; next = NULL;
  npts = n;
  index = ind;
  A = amps;
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

src_vol *src_vol::add_to(src_vol *others) {
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
    next = others;
    return this;
  }
}

/*********************************************************************/

void fields::add_point_source(component c, double freq,
                              double width, double peaktime,
                              double cutoff, const vec &p,
                              complex<double> amp, int is_c) {
  width /= freq;

  if (is_c) { // TODO: don't ignore peaktime?
    continuous_src_time src(freq, width, time(), infinity, cutoff);
    add_point_source(c, src, p, amp);
  }
  else {
    cutoff = v.inva + cutoff * width;
    if (peaktime <= 0.0)
      peaktime = time() + cutoff;
  
    gaussian_src_time src(freq, width,
			     peaktime - cutoff, peaktime + cutoff);
    add_point_source(c, src, p, amp);
  }
}

void fields::add_point_source(component c, const src_time &src,
			      const vec &p, complex<double> amp) {
  add_volume_source(c, src, geometric_volume(p, p), amp);
}

complex<double> one(const vec &v) {(void) v; return 1.0;}
void fields::add_volume_source(component c, const src_time &src,
                               const geometric_volume &where,
			       complex<double> amp) {
  add_volume_source(c, src, where, one, amp);
}

struct src_vol_chunkloop_data {
  complex<double> (*A)(const vec &);
  complex<double> amp;
  src_time *src;
  vec center;
};

/* Adding source volumes can be treated as a kind of "integration"
   problem, since we need to loop over all the chunks that intersect
   the source volume, with appropriate interpolation weights at the
   boundaries so that the integral of the current is fixed regardless
   of resolution.  Unlike most uses of fields::loop_in_chunks, however, we
   set use_symmetry=false: we only find the intersection of the volume
   with the untransformed chunks (since the transformed versions are
   implicit). */
static void src_vol_chunkloop(fields_chunk *fc, component c,
			      ivec is, ivec ie,
			      vec s0, vec s1, vec e0, vec e1,
			      double dV0, double dV1,
			      ivec shift, complex<double> shift_phase, 
			      const symmetry &S, int sn,
			      void *data_)
{
  src_vol_chunkloop_data *data = (src_vol_chunkloop_data *) data_;
  
  (void) S; (void) sn; // these should be the identity
  (void) dV0; (void) dV1; // volume weighting is included in data->amp

  int npts = 1;
  LOOP_OVER_DIRECTIONS(is.dim, d)
    npts *= (ie.in_direction(d) - is.in_direction(d)) / 2 + 1;
  int *index_array = new int[npts];
  complex<double> *amps_array = new complex<double>[npts];

  complex<double> amp = data->amp * conj(shift_phase);

  direction cd = component_direction(c);

  double inva = fc->v.inva;
  int idx_vol = 0;
  LOOP_OVER_IVECS(fc->v, is, ie, idx) {
    IVEC_LOOP_LOC(fc->v, loc);
    loc += shift * (0.5*inva) - data->center;

    amps_array[idx_vol] = IVEC_LOOP_WEIGHT(s0,s1,e0,e1,1) * amp * data->A(loc);

    /* for "D" sources, multiply by epsilon.  FIXME: this is not quite
       right because it doesn't handle non-diagonal inveps! */
    if (is_D(c) && fc->s->inveps[c][cd]) 
      amps_array[idx_vol] /= fc->s->inveps[c][cd][idx];

    index_array[idx_vol++] = idx;
  }

  if (idx_vol != npts)
    abort("add_volume_source: computed wrong npts (%d vs. %d)", npts, idx_vol);

  src_vol *tmp = new src_vol(c, data->src, npts, index_array, amps_array);
  if (is_magnetic(c))
    fc->h_sources = tmp->add_to(fc->h_sources);
  else
    fc->e_sources = tmp->add_to(fc->e_sources);
}

void fields::require_component(component c) {
  if (!v.has_field(c))
    abort("cannot require a %s component in a %s grid",
	  component_name(c), dimension_name(v.dim));
  // allocate fields if they haven't been allocated yet for this component
  int need_to_reconnect = 0;
  for (int i = 0; i < num_chunks; ++i)
    if (chunks[i]->is_mine() && !chunks[i]->f[c][0]) {
      chunks[i]->alloc_f(c);
      need_to_reconnect++;
    }
  if (chunk_connections_valid && sum_to_all(need_to_reconnect)) 
    chunk_connections_valid = false;
}

void fields::add_volume_source(component c, const src_time &src,
                               const geometric_volume &where,
                               complex<double> A(const vec &), 
			       complex<double> amp) {
  if (v.dim != where.dim)
    abort("incorrect source volume dimensionality in add_volume_source");
  LOOP_OVER_DIRECTIONS(v.dim, d) 
    if (where.in_direction(d) > (user_volume.boundary_location(High, d) 
				 - user_volume.boundary_location(Low, d)))
      abort("Source width > cell width in %s direction!\n", direction_name(d));

  src_vol_chunkloop_data data;
  data.A = A;
  data.amp = amp;
  LOOP_OVER_DIRECTIONS(v.dim, d)
    if (where.in_direction(d) == 0.0) // delta-function direction
      data.amp *= v.a; // correct units for J delta-function amplitude
  sources = src.add_to(sources, &data.src);
  data.center = (where.get_min_corner() + where.get_max_corner()) * 0.5;
  loop_in_chunks(src_vol_chunkloop, (void *) &data, where, c, false);
  require_component(c);
}

} // namespace meep
