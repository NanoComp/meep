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
#include <math.h>

#include "meep.h"

namespace meep {

static inline int min(int a, int b) { return (a<b)?a:b; }
static inline int max(int a, int b) { return (a>b)?a:b; }
static inline int abs(int a) { return a < 0 ? -a : a; }

typedef struct {
  h5file *file;
  ivec min_corner, max_corner;
  int num_chunks;
  component c;
  int reim;
  double *buf;
  int bufsz;
  int rank;
  direction ds[3];
} h5_output_data;

static void update_datasize(h5_output_data *data, fields_chunk *fc,
			    const ivec is, const ivec ie,
			    const ivec shift, const symmetry &S, int sn)
{
    ivec isS = S.transform(is, sn) + shift;
    ivec ieS = S.transform(ie, sn) + shift;
    data->min_corner = min(data->min_corner, min(isS, ieS));
    data->max_corner = max(data->max_corner, max(isS, ieS));
    data->num_chunks++;
    int bufsz = 1;
    LOOP_OVER_DIRECTIONS(fc->v.dim, d)
      bufsz *= (ie.in_direction(d) - is.in_direction(d)) / 2 + 1;
    data->bufsz = max(data->bufsz, bufsz);
}

static void h5_findsize_integrand(fields_chunk *fc, component cgrid,
				  ivec is, ivec ie,
				  vec s0, vec s1, vec e0, vec e1,
				  double dV0, double dV1,
				  ivec shift, complex<double> shift_phase,
				  const symmetry &S, int sn,
				  void *data_)
{
  h5_output_data *data = (h5_output_data *) data_;
  component cS = S.transform(data->c, -sn);
  double *f = cS == Dielectric ? fc->s->eps : fc->f[cS][data->reim];
  if (f) update_datasize(data, fc, is, ie, shift, S, sn);
}

static void get_output_dimensions(int start[3], int count[3],
				  int offset[3], int stride[3],
				  h5_output_data *data, fields_chunk *fc,
				  const ivec &is, const ivec &ie,
				  const ivec &shift, const symmetry &S, int sn)
{
    ivec isS = S.transform(is, sn) + shift;
    ivec ieS = S.transform(ie, sn) + shift;

    for (int i = 0; i < 3; ++i) {
      start[i] = offset[i] = 0;
      count[i] = stride[i] = 1;
    }

    // figure out what yucky_directions (in LOOP_OVER_IVECS)
    // correspond to what directions in the transformed vectors (in output).
    ivec permute(zero_ivec(fc->v.dim));
    for (int i = 0; i < 3; ++i) 
      permute.set_direction(fc->v.yucky_direction(i), i);
    permute = S.transform_unshifted(permute, sn);
    LOOP_OVER_DIRECTIONS(permute.dim, d)
      permute.set_direction(d, abs(permute.in_direction(d)));

    // compute the size of the chunk to output, and its strides etc.
    for (int i = 0; i < data->rank; ++i) {
      direction d = data->ds[i];
      int isd = isS.in_direction(d), ied = ieS.in_direction(d);
      start[i] = (min(isd, ied) - data->min_corner.in_direction(d)) / 2;
      count[i] = abs(ied - isd) / 2 + 1;
      int j = permute.in_direction(d);
      if (ied < isd) offset[permute.in_direction(d)] = count[i] - 1;
    }
    for (int i = 0; i < data->rank; ++i) {
      direction d = data->ds[i];
      int j = permute.in_direction(d);
      for (int k = i + 1; k < data->rank; ++k) stride[j] *= count[k];
      offset[j] *= stride[j];
      if (offset[j]) stride[j] *= -1;
    }
}

static void h5_output_integrand(fields_chunk *fc, component cgrid,
				ivec is, ivec ie,
				vec s0, vec s1, vec e0, vec e1,
				double dV0, double dV1,
				ivec shift, complex<double> shift_phase,
				const symmetry &S, int sn,
				void *data_)
{
  h5_output_data *data = (h5_output_data *) data_;
  component cS = S.transform(data->c, -sn);
  double *f = cS == Dielectric ? fc->s->eps : fc->f[cS][data->reim];
  if (f) {
    int start[3], count[3], offset[3], stride[3];
    get_output_dimensions(start, count, offset, stride,
			  data, fc, is, ie, shift, S, sn);

    // cgrid is Dielectric grid, so we need to average onto it:
    int o1, o2;
    fc->v.yee2diel_offsets(cS, o1, o2);

    shift_phase *= S.phase_shift(cS, sn); // vector component may flip
    
    // Copy data to buffer, taking shift_phase into account:
    if (cS == Dielectric) { // no phase
      LOOP_OVER_IVECS(fc->v, is, ie, idx) {
	int idx2 = ((((offset[0] + offset[1] + offset[2])
		      + loop_i1 * stride[0]) 
		     + loop_i2 * stride[1]) + loop_i3 * stride[2]);
	data->buf[idx2] = f[idx];
      }
    }
    else if (imag(shift_phase) == 0.0) { // real phase (possibly real field)
      LOOP_OVER_IVECS(fc->v, is, ie, idx) {
	int idx2 = ((((offset[0] + offset[1] + offset[2])
		      + loop_i1 * stride[0]) 
		     + loop_i2 * stride[1]) + loop_i3 * stride[2]);
	data->buf[idx2] = (f[idx] + f[idx+o1] + f[idx+o2] + f[idx+o1+o2])
	  * (0.25 * real(shift_phase));
      }
    }
    else { // complex phase: do complex multiplication with complex field
      double *fr = fc->f[cS][0], *fi = fc->f[cS][1];
      if (!fi) abort("complex Bloch boundary condition with real field!");
      LOOP_OVER_IVECS(fc->v, is, ie, idx) {
	int idx2 = ((((offset[0] + offset[1] + offset[2])
		      + loop_i1 * stride[0])
		     + loop_i2 * stride[1]) + loop_i3 * stride[2]);
	double re = 0.25 * (fr[idx] + fr[idx+o1] + fr[idx+o2] + fr[idx+o1+o2]);
	double im = 0.25 * (fi[idx] + fi[idx+o1] + fi[idx+o2] + fi[idx+o1+o2]);
	data->buf[idx2] = data->reim 
	  ? re * imag(shift_phase) + im * real(shift_phase)
	  : re * real(shift_phase) - im * imag(shift_phase);
      }
    }
    data->file->write_chunk(data->rank, start, count, data->buf);
  }
}

void fields::output_hdf5(h5file *file, const char *dataname,
			 component c, int reim,
			 const geometric_volume &where,
			 bool append_data,
                         bool single_precision) {
  h5_output_data data;

  data.file = file;
  data.min_corner = v.round_vec(where.get_max_corner()) + one_ivec(v.dim);
  data.max_corner = v.round_vec(where.get_min_corner()) - one_ivec(v.dim);
  data.num_chunks = 0;
  data.bufsz = 0;
  data.c = c; data.reim = reim;

  integrate(h5_findsize_integrand, (void *) &data, 
	    where, Dielectric, true, true);

  data.max_corner = max_to_all(data.max_corner);
  data.min_corner = -max_to_all(-data.min_corner); // i.e., min_to_all
  data.num_chunks = sum_to_all(data.num_chunks);
  if (data.num_chunks == 0 || !(data.min_corner <= data.max_corner))
    return; // no data to write;

  int rank = 0, dims[3];
  LOOP_OVER_DIRECTIONS(v.dim, d) {
    if (rank >= 3) abort("too many dimensions in output_hdf5");
    int n = (data.max_corner.in_direction(d)
	     - data.min_corner.in_direction(d)) / 2 + 1;
    if (n > 1) {
      data.ds[rank] = d;
      dims[rank++] = n;
    }
  }
  data.rank = rank;

  file->create_or_extend_data(dataname, rank, dims,
                              append_data, single_precision);

  data.buf = new double[data.bufsz];

  integrate(h5_output_integrand, (void *) &data, 
	    where, Dielectric, true, true);

  delete[] data.buf;
  file->done_writing_chunks();
}

void fields::output_hdf5(h5file *file, component c,
			 const geometric_volume &where,
			 bool append_data,
                         bool single_precision) {
  char dataname[256];
  bool has_imag = !is_real && c != Dielectric;

  snprintf(dataname, 256, "%s%s", component_name(c), has_imag ? ".r" : "");
  output_hdf5(file, dataname, c, 0, where, append_data, single_precision);
  if (has_imag) {
    snprintf(dataname, 256, "%s.i", component_name(c));
    output_hdf5(file, dataname, c, 1, where, append_data,single_precision);
  }
}

void fields::output_hdf5(component c,
			 const geometric_volume &where,
			 bool single_precision, 
			 const char *prefix) {
  h5file *file = open_h5file(component_name(c), h5file::WRITE,
			     prefix, true);
  output_hdf5(file, c, where, false, single_precision);
  delete file;
}

h5file *fields::open_h5file(const char *name, h5file::access_mode mode, 
			    const char *prefix, bool timestamp)
{
  const int buflen = 1024;
  char filename[buflen];
  char time_step_string[32] = "";

  if (timestamp) snprintf(time_step_string, 32, "-%09.2f", time());

  snprintf(filename, buflen, "%s/" "%s%s" "%s" "%s" ".h5",
	   outdir,
	   prefix ? prefix : "", prefix && prefix[0] ? "-" : "",
	   name, time_step_string);
  return new h5file(filename, mode, true);
}

} // meep
