/* Copyright (C) 2005-2014 Massachusetts Institute of Technology
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

// For sine-integral function used in the band-source
#define HAVE_LIBGSL_EXPINT

#ifdef HAVE_LIBGSL_EXPINT
#  include <gsl/gsl_sf_expint.h>
#endif
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
	       others->next = add_to(others->next, added);
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

gaussian_src_time::gaussian_src_time(double f, double fwidth, double s)
{
  freq = f;
  width = 1.0 / fwidth;
  peak_time = width * s;
  cutoff = width * s * 2;

  // this is to make last_source_time as small as possible
  while (exp(-cutoff*cutoff / (2*width*width)) < 1e-100)
    cutoff *= 0.9;
  cutoff = float(cutoff); // don't make cutoff sensitive to roundoff error
}

gaussian_src_time::gaussian_src_time(double f, double w, double st, double et)
{
  freq = f;
  width = w;
  peak_time = 0.5 * (st + et);
  cutoff = (et - st) * 0.5;

  // this is to make last_source_time as small as possible
  while (exp(-cutoff*cutoff / (2*width*width)) < 1e-100)
    cutoff *= 0.9;
  cutoff = float(cutoff); // don't make cutoff sensitive to roundoff error
}

complex<double> gaussian_src_time::dipole(double time) const
{
  double tt = time - peak_time;
  if (float(fabs(tt)) > cutoff)
    return 0.0;

  // correction factor so that current amplitude (= d(dipole)/dt) is
  // ~ 1 near the peak of the Gaussian.
  complex<double> amp = 1.0 / complex<double>(0,-2*pi*freq);

  return exp(-tt*tt / (2*width*width)) * polar(1.0, -2*pi*freq*tt) * amp;
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
/// Band source
band_src_time::band_src_time(double f, double fwidth, double user_cutoff)
{
  freq = f;
  width = 1.0 / fwidth;
  peak_time = user_cutoff/2;
  cutoff = user_cutoff;
#ifdef HAVE_LIBGSL_EXPINT
  master_printf("Initializing band source for time_domain (peak_time = %g, cutoff = %g)\n", peak_time, cutoff);
  master_printf("\tExperimental: using GSL si() function for  flat-top spectrum\n");
#else
  master_printf("Initializing band source for time_domain (peak_time = %g, cutoff = %g)\n", peak_time, cutoff);
  master_printf("\tWarning: not compiled with GSL, the source spectrum will not be flat-top\n", peak_time, cutoff);
#endif
  cutoff = float(cutoff); // don't make cutoff sensitive to roundoff error
}

complex<double> band_src_time::dipole(double time) const
{
  double tt = time - peak_time;

  // The function that introduces a rectangular band in the spectrum (centered around zero frequency)
#ifdef HAVE_LIBGSL_EXPINT
  // The emitted field is the derivative of the dipole, so if possible, we use the sine integral: 
  // Si(x) = \int_0^x dt sin(t)/t
  complex<double> func        = gsl_sf_Si(tt*2*pi*(freq-.5/width)) - gsl_sf_Si(tt*2*pi*(freq+.5/width));
#else
  // If GSL not available, a reasonable approximation is the sinc(t) function (but has not flat top)
  complex<double> func        = sin(tt*2*pi / width/2)/(tt*2*pi * width) * polar(1.0, -2*pi*freq*tt); 
#endif

  // The envelope that suppresses ringing and side lobes of the rectangle
  double wnd_BlackmanN= (0.3635819 + 0.4891775*cos(tt/(cutoff)*pi*2) + 
		  0.1365995*cos(tt/(cutoff)*pi*4)+ 0.0106411*cos(tt/(cutoff)*pi*6));
  //double wnd_Hann     = (.5 + .5*cos(tt/(cutoff)*pi*2));
  //double wnd_Hamming  = (0.53836+0.46164*cos(tt/(cutoff)*pi*2))     ;
  //double wnd_Blackman = (0.42659 + 0.49656*cos(tt/(cutoff)*pi*2) + 0.076849*cos(tt/(cutoff)*pi*4))     ;
  //double wnd_Gauss    = exp(-((tt)/(cutoff)*3)**2)       // Gaussian window, char. width = 1/3;
  //double wnd_ContRect = // should crop sinc when passing zero: rect(tt, 0, (trunc((cutoff)/8/width))*8*width) ;
  
  // correction factor so that current amplitude (= d(dipole)/dt) is
  // ~ 1 near the peak of the band.
  complex<double> amp = 1.0 / complex<double>(0,-2*pi*freq);  // TODO

  // The needed time-domain source amplitude is the sinc function constrained by window and shifted
  // in frequency by complex exponential
  if (abs(tt) < cutoff/2) 
	  //return sinc;
	  return func * wnd_BlackmanN; 
	  //return func * wnd_BlackmanN  * amp; 
  else
	  return 0;

}

bool band_src_time::is_equal(const src_time &t) const
{
     const band_src_time *tp = dynamic_cast<const band_src_time*>(&t);
     if (tp)
	  return(tp->freq == freq && tp->width == width &&
		 tp->peak_time == peak_time && tp->cutoff == cutoff);
     else
	  return 0;
}

complex<double> continuous_src_time::dipole(double time) const
{
  float rtime = float(time);
  if (rtime < start_time || rtime > end_time)
    return 0.0;

  // correction factor so that current amplitude (= d(dipole)/dt) is 1.
  complex<double> amp = 1.0 / (complex<double>(0,-1.0) * (2*pi)*freq);

  if (width == 0.0)
    return exp(complex<double>(0,-1.0) * (2*pi)*freq*time) * amp;
  else {
    double ts = (time - start_time) / width - slowness;
    double te = (end_time - time) / width - slowness;
    
    return exp(complex<double>(0,-1.0) * (2*pi)*freq*time) * amp
      * (1.0 + tanh(ts))  // goes from 0 to 2
      * (1.0 + tanh(te))  // goes from 2 to 0
      * 0.25;
  }
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

bool custom_src_time::is_equal(const src_time &t) const
{
     const custom_src_time *tp = dynamic_cast<const custom_src_time*>(&t);
     if (tp)
	  return(tp->start_time == start_time && tp->end_time == end_time &&
		 tp->func == func && tp->data == data);
     else
	  return 0;
}

/*********************************************************************/

src_vol::src_vol(component cc, src_time *st, int n, int *ind, complex<double> *amps) {
  c = cc;
  if (is_D(c)) c = direction_component(Ex, component_direction(c));
  if (is_B(c)) c = direction_component(Hx, component_direction(c));
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
        abort("Cannot add grid_volume sources with different number of points\n");
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

// THIS VARIANT IS FOR BACKWARDS COMPATIBILITY, and is DEPRECATED:
void fields::add_point_source(component c, double freq,
                              double width, double peaktime,
                              double cutoff, const vec &p,
                              complex<double> amp, int is_c) {
  width /= freq;

  if (is_c) { // TODO: don't ignore peaktime?
    continuous_src_time src(freq, width, time(), infinity, cutoff);
    if (is_magnetic(c)) src.is_integrated = false;
    add_point_source(c, src, p, amp);
  }
  else {
    cutoff = gv.inva + cutoff * width;
    if (peaktime <= 0.0)
      peaktime = time() + cutoff;
    
    // backward compatibility (slight phase shift in old Meep version)
    peaktime += is_magnetic(c) ? -dt*0.5 : dt;
  
    gaussian_src_time src(freq, width,
			     peaktime - cutoff, peaktime + cutoff);
    if (is_magnetic(c)) src.is_integrated = false;
    add_point_source(c, src, p, amp);
  }
}

void fields::add_point_source(component c, const src_time &src,
			      const vec &p, complex<double> amp) {
  add_volume_source(c, src, volume(p, p), amp);
}

static complex<double> one(const vec &pt) {(void) pt; return 1.0;}
void fields::add_volume_source(component c, const src_time &src,
                               const volume &where,
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
   the source grid_volume, with appropriate interpolation weights at the
   boundaries so that the integral of the current is fixed regardless
   of resolution.  Unlike most uses of fields::loop_in_chunks, however, we
   set use_symmetry=false: we only find the intersection of the grid_volume
   with the untransformed chunks (since the transformed versions are
   implicit). */
static void src_vol_chunkloop(fields_chunk *fc, int ichunk, component c,
			      ivec is, ivec ie,
			      vec s0, vec s1, vec e0, vec e1,
			      double dV0, double dV1,
			      ivec shift, complex<double> shift_phase, 
			      const symmetry &S, int sn,
			      void *data_)
{
  src_vol_chunkloop_data *data = (src_vol_chunkloop_data *) data_;
  
  (void) S; (void) sn; // these should be the identity
  (void) dV0; (void) dV1; // grid_volume weighting is included in data->amp
  (void) ichunk;

  int npts = 1;
  LOOP_OVER_DIRECTIONS(is.dim, d)
    npts *= (ie.in_direction(d) - is.in_direction(d)) / 2 + 1;
  int *index_array = new int[npts];
  complex<double> *amps_array = new complex<double>[npts];

  complex<double> amp = data->amp * conj(shift_phase);

  direction cd = component_direction(c);

  double inva = fc->gv.inva;
  int idx_vol = 0;
  LOOP_OVER_IVECS(fc->gv, is, ie, idx) {
    IVEC_LOOP_LOC(fc->gv, loc);
    loc += shift * (0.5*inva) - data->center;

    amps_array[idx_vol] = IVEC_LOOP_WEIGHT(s0,s1,e0,e1,1) * amp * data->A(loc);

    /* for "D" sources, multiply by epsilon.  FIXME: this is not quite
       right because it doesn't handle non-diagonal chi1inv! 
       similarly, for "B" sources, multiply by mu. */
    if (is_D(c) && fc->s->chi1inv[c-Dx+Ex][cd]) 
      amps_array[idx_vol] /= fc->s->chi1inv[c-Dx+Ex][cd][idx];
    if (is_B(c) && fc->s->chi1inv[c-Bx+Hx][cd]) 
      amps_array[idx_vol] /= fc->s->chi1inv[c-Bx+Hx][cd][idx];

    index_array[idx_vol++] = idx;
  }

  if (idx_vol != npts)
    abort("add_volume_source: computed wrong npts (%d vs. %d)", npts, idx_vol);

  src_vol *tmp = new src_vol(c, data->src, npts, index_array, amps_array);
  field_type ft = is_magnetic(c) ? B_stuff : D_stuff;
  fc->sources[ft] = tmp->add_to(fc->sources[ft]);
}

void fields::add_volume_source(component c, const src_time &src,
                               const volume &where_,
                               complex<double> A(const vec &), 
			       complex<double> amp) {
  volume where(where_); // make a copy to adjust size if necessary
  if (gv.dim != where.dim)
    abort("incorrect source grid_volume dimensionality in add_volume_source");
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    double w = user_volume.boundary_location(High, d)
      - user_volume.boundary_location(Low, d);
    if (where.in_direction(d) > w + gv.inva)
      abort("Source width > cell width in %s direction!\n", direction_name(d));
    else if (where.in_direction(d) > w) { // difference is less than 1 pixel
      double dw = where.in_direction(d) - w;
      where.set_direction_min(d, where.in_direction_min(d) - dw * 0.5);
      where.set_direction_max(d, where.in_direction_min(d) + w);
    }
  }

  src_vol_chunkloop_data data;
  data.A = A ? A : one;
  data.amp = amp;
  LOOP_OVER_DIRECTIONS(gv.dim, d)
    if (where.in_direction(d) == 0.0 && !nosize_direction(d)) // delta-fun
      data.amp *= gv.a; // correct units for J delta-function amplitude
  sources = src.add_to(sources, &data.src);
  data.center = (where.get_min_corner() + where.get_max_corner()) * 0.5;
  loop_in_chunks(src_vol_chunkloop, (void *) &data, where, c, false);
  require_component(c);
}

} // namespace meep
