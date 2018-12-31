/* Copyright (C) 2005-2015 Massachusetts Institute of Technology
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

using namespace std;

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
  cutoff = width * s;

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

src_vol::src_vol(component cc, src_time *st, size_t n, ptrdiff_t *ind, complex<double> *amps) {
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
  index = new ptrdiff_t[npts];
  A = new complex<double>[npts];
  for (size_t j=0; j<npts; j++) {
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
      for (size_t j=0; j<npts; j++) {
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

// implementation of pure virtual destructor for source_indicator
source_indicator::~source_indicator() {}

static complex<double> one(const vec &pt) {(void) pt; return 1.0;}
void fields::add_volume_source(component c, const src_time &src,
                               const volume &where,
			       complex<double> amp,
			       source_indicator *indicator) {
  add_volume_source(c, src, where, one, amp, indicator);
}

struct src_vol_chunkloop_data {
  complex<double> (*A)(const vec &);
  complex<double> amp;
  src_time *src;
  vec center;
  source_indicator *indicator;
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
  
  size_t npts = 1;
  LOOP_OVER_DIRECTIONS(is.dim, d)
    npts *= (ie.in_direction(d) - is.in_direction(d)) / 2 + 1;
  ptrdiff_t *index_array = new ptrdiff_t[npts];
  complex<double> *amps_array = new complex<double>[npts];

  complex<double> amp = data->amp * conj(shift_phase);

  direction cd = component_direction(c);

  double inva = fc->gv.inva;
  size_t idx_vol = 0;
  LOOP_OVER_IVECS(fc->gv, is, ie, idx) {
    IVEC_LOOP_ILOC(fc->gv, iloc);
    if (!fc->gv.owns(iloc)) continue;

    IVEC_LOOP_LOC(fc->gv, loc);
    loc += shift * (0.5*inva);

    if (data->indicator && data->indicator->point_in_source(loc)==false) continue;

    vec rel_loc = loc - data->center;
    amps_array[idx_vol] = IVEC_LOOP_WEIGHT(s0,s1,e0,e1,1) * amp * data->A(rel_loc);

    /* for "D" sources, multiply by epsilon.  FIXME: this is not quite
       right because it doesn't handle non-diagonal chi1inv!
       similarly, for "B" sources, multiply by mu. */
    if (is_D(c) && fc->s->chi1inv[c-Dx+Ex][cd])
      amps_array[idx_vol] /= fc->s->chi1inv[c-Dx+Ex][cd][idx];
    if (is_B(c) && fc->s->chi1inv[c-Bx+Hx][cd])
      amps_array[idx_vol] /= fc->s->chi1inv[c-Bx+Hx][cd][idx];

    index_array[idx_vol++] = idx;
  }

  if (idx_vol > npts)
    abort("add_volume_source: computed wrong npts (%zd vs. %zd)", npts, idx_vol);

  src_vol *tmp = new src_vol(c, data->src, idx_vol, index_array, amps_array);
  field_type ft = is_magnetic(c) ? B_stuff : D_stuff;
  fc->sources[ft] = tmp->add_to(fc->sources[ft]);
}

static realnum *amp_func_data_re = NULL;
static realnum *amp_func_data_im = NULL;
static const volume *amp_func_vol = NULL;
static size_t amp_file_dims[3];

complex<double> amp_file_func(const vec &p) {
  double x_size = 0, y_size = 0, z_size = 0;

  switch (amp_func_vol->dim) {
  case D1:
    z_size = amp_func_vol->in_direction(Z);
    break;
  case D2:
    x_size = amp_func_vol->in_direction(X);
    y_size = amp_func_vol->in_direction(Y);
    break;
  case D3:
    x_size = amp_func_vol->in_direction(X);
    y_size = amp_func_vol->in_direction(Y);
    z_size = amp_func_vol->in_direction(Z);
    break;
  case Dcyl:
    x_size = amp_func_vol->in_direction(X);
    z_size = amp_func_vol->in_direction(Z);
    break;
  }

  double rx = x_size == 0 ? 0 : 0.5 + p.x() / x_size;
  double ry = y_size == 0 ? 0 : 0.5 + p.y() / y_size;
  double rz = z_size == 0 ? 0 : 0.5 + p.z() / z_size;

  complex<double> res;
  res.real(linear_interpolate(rx, ry, rz, amp_func_data_re,
                              amp_file_dims[0], amp_file_dims[1], amp_file_dims[2], 1));
  res.imag(linear_interpolate(rx, ry, rz, amp_func_data_im,
                              amp_file_dims[0], amp_file_dims[1], amp_file_dims[2], 1));
  return res;
}

void fields::add_volume_source(component c, const src_time &src, const volume &where_,
                               complex<double> *arr, size_t dim1, size_t dim2, size_t dim3,
                               complex<double> amp, source_indicator *indicator) {

  amp_func_vol = &where_;

  amp_file_dims[0] = dim1;
  amp_file_dims[1] = dim2;
  amp_file_dims[2] = dim3;

  size_t total_size = dim1 * dim2 * dim3;
  amp_func_data_re = new double[total_size];
  amp_func_data_im = new double[total_size];

  for (size_t i = 0; i < total_size; ++i) {
    amp_func_data_re[i] = real(arr[i]);
    amp_func_data_im[i] = imag(arr[i]);
  }

  add_volume_source(c, src, where_, amp_file_func, amp, indicator);

  delete[] amp_func_data_re;
  delete[] amp_func_data_im;
}

// Reads amplitude function data from h5file 'filename.' Assumes real and imaginary components
// of 'dataset' exist with '.re' and '.im' extensions.
void fields::add_volume_source(component c, const src_time &src, const volume &where_,
                               const char *filename, const char *dataset,
                               complex<double> amp, source_indicator *indicator) {

  meep::h5file eps_file(filename, meep::h5file::READONLY, false);
  int rank;
  std::string dataset_re = std::string(dataset) + ".re";
  std::string dataset_im = std::string(dataset) + ".im";

  size_t re_dims[] = {1, 1, 1};
  double *real_data = eps_file.read(dataset_re.c_str(), &rank, re_dims, 3);
  master_printf("read in %zdx%zdx%zd amplitude function file \"%s:%s\"\n",
                re_dims[0], re_dims[1], re_dims[2], filename, dataset_re.c_str());

  size_t im_dims[] = {1, 1, 1};
  double *imag_data = eps_file.read(dataset_im.c_str(), &rank, im_dims, 3);
  master_printf("read in %zdx%zdx%zd amplitude function file \"%s:%s\"\n",
                im_dims[0], im_dims[1], im_dims[2], filename, dataset_im.c_str());

  for (int i = 0; i < 3; ++i) {
    if (re_dims[i] != im_dims[i]) {
      abort("Imaginary and real datasets have different dimensions");
    }
  }

  size_t total_size = re_dims[0] * re_dims[1] * re_dims[2];
  complex<double> *amp_data = new complex<double>[total_size];
  for (size_t i = 0; i < total_size; ++i) {
    amp_data[i] = complex<double>(real_data[i], imag_data[i]);
  }

  add_volume_source(c, src, where_, amp_data, re_dims[0], re_dims[1], re_dims[2], amp, indicator);

  delete[] real_data;
  delete[] imag_data;
  delete[] amp_data;
}

void fields::add_volume_source(component c, const src_time &src,
                               const volume &where_,
                               complex<double> A(const vec &),
			       complex<double> amp,
			       source_indicator *indicator)
{
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
  data.indicator = indicator;
  loop_in_chunks(src_vol_chunkloop, (void *) &data, where, c, false);
  require_component(c);
}

} // namespace meep
