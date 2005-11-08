/* Copyright (C) 2005 Massachusetts Institute of Technology
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

/* HDF5 output of fields and arbitrary functions thereof.  Works
   very similarly to integrate.cpp (using fields::loop_in_chunks). */

#include <stdio.h>
#include <math.h>

#include "meep_internals.hpp"

namespace meep {

/***************************************************************************/

typedef struct {
  // information related to the HDF5 dataset (its size, etcetera)
  h5file *file;
  ivec min_corner, max_corner;
  int num_chunks;
  double *buf;
  int bufsz;
  int rank;
  direction ds[3];

  int reim; // whether to output the real or imaginary part

  // the function to output and related info (offsets for averaging, etc.)
  int num_fields;
  const component *components;
  component *cS;
  complex<double> *ph;
  complex<double> *fields;
  int *offsets;
  int ninveps;
  component inveps_cs[3];
  direction inveps_ds[3];
  int inveps_offsets[6];
  complex<long double> sum;
  double maxabs;
  field_function fun;
  void *fun_data_;
} h5_output_data;

static void h5_findsize_chunkloop(fields_chunk *fc, int ichnk, component cgrid,
				  ivec is, ivec ie,
				  vec s0, vec s1, vec e0, vec e1,
				  double dV0, double dV1,
				  ivec shift, complex<double> shift_phase,
				  const symmetry &S, int sn,
				  void *data_)
{
  h5_output_data *data = (h5_output_data *) data_;
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

static void h5_output_chunkloop(fields_chunk *fc, int ichnk, component cgrid,
				ivec is, ivec ie,
				vec s0, vec s1, vec e0, vec e1,
				double dV0, double dV1,
				ivec shift, complex<double> shift_phase,
				const symmetry &S, int sn,
				void *data_)
{
  h5_output_data *data = (h5_output_data *) data_;

  //-----------------------------------------------------------------------//
  // Find output chunk dimensions and strides, etc.

  int start[3]={0,0,0}, count[3]={1,1,1};
  int offset[3]={0,0,0}, stride[3]={1,1,1};

  ivec isS = S.transform(is, sn) + shift;
  ivec ieS = S.transform(ie, sn) + shift;
  
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
  
  //-----------------------------------------------------------------------//
  // Compute the function to output, exactly as in fields::integrate,
  // except that here we store its values in a buffer instead of integrating.

  int *off = data->offsets;
  component *cS = data->cS;
  complex<double> *fields = data->fields, *ph = data->ph;
  complex<long double> sum = 0.0;
  double maxabs = 0;
  const component *iecs = data->inveps_cs;
  const direction *ieds = data->inveps_ds;
  int *ieos = data->inveps_offsets;

  for (int i = 0; i < data->num_fields; ++i) {
    cS[i] = S.transform(data->components[i], -sn);
    if (cS[i] == Dielectric)
      ph[i] = 1.0;
    else {
      fc->v.yee2diel_offsets(cS[i], off[2*i], off[2*i+1]);
      ph[i] = shift_phase * S.phase_shift(cS[i], sn);
    }
  }
  for (int k = 0; k < data->ninveps; ++k)
    fc->v.yee2diel_offsets(iecs[k], ieos[2*k], ieos[2*k+1]);

  vec rshift(shift * (0.5*fc->v.inva));
  LOOP_OVER_IVECS(fc->v, is, ie, idx) {
    IVEC_LOOP_LOC(fc->v, loc);
    loc = S.transform(loc, sn) + rshift;

    for (int i = 0; i < data->num_fields; ++i) {
      if (cS[i] == Dielectric) {
	double tr = 0.0;
	for (int k = 0; k < data->ninveps; ++k)
	  tr += (fc->s->inveps[iecs[k]][ieds[k]][idx]
		 + fc->s->inveps[iecs[k]][ieds[k]][idx+ieos[2*k]]
		 + fc->s->inveps[iecs[k]][ieds[k]][idx+ieos[1+2*k]]
		 + fc->s->inveps[iecs[k]][ieds[k]][idx+ieos[2*k]+ieos[1+2*k]]);
	fields[i] = (4 * data->ninveps) / tr;
      }
      else {
	double f[2];
	for (int k = 0; k < 2; ++k)
	  if (fc->f[cS[i]][k])
	    f[k] = 0.25 * (fc->f[cS[i]][k][idx]
			   + fc->f[cS[i]][k][idx+off[2*i]]
			   + fc->f[cS[i]][k][idx+off[2*i+1]]
			   + fc->f[cS[i]][k][idx+off[2*i]+off[2*i+1]]);
	  else
	    f[k] = 0;
	fields[i] = complex<double>(f[0], f[1]) * ph[i];
      }
    }

    complex<double> fun = data->fun(fields, loc, data->fun_data_);
    int idx2 = ((((offset[0] + offset[1] + offset[2])
		  + loop_i1 * stride[0]) 
		 + loop_i2 * stride[1]) + loop_i3 * stride[2]);
    data->buf[idx2] = data->reim ? imag(fun) : real(fun);
  }

  //-----------------------------------------------------------------------//

  data->file->write_chunk(data->rank, start, count, data->buf);
}

void fields::output_hdf5(h5file *file, const char *dataname,
			 int num_fields, const component *components,
			 field_function fun, int reim,
			 const geometric_volume &where,
			 void *fun_data_,
			 bool append_data,
                         bool single_precision) {
  am_now_working_on(FieldOutput);
  h5_output_data data;

  data.file = file;
  data.min_corner = v.round_vec(where.get_max_corner()) + one_ivec(v.dim);
  data.max_corner = v.round_vec(where.get_min_corner()) - one_ivec(v.dim);
  data.num_chunks = 0;
  data.bufsz = 0;
  data.reim = reim;

  loop_in_chunks(h5_findsize_chunkloop, (void *) &data, 
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

  data.num_fields = num_fields;
  data.components = components;
  data.cS = new component[num_fields];
  data.ph = new complex<double>[num_fields];
  data.fields = new complex<double>[num_fields];
  data.sum = 0;
  data.maxabs = 0;
  data.fun = fun;
  data.fun_data_ = fun_data_;

  /* compute inverse-epsilon directions for computing Dielectric fields */
  data.ninveps = 0;
  bool needs_dielectric = false;
  for (int i = 0; i < num_fields; ++i)
    if (components[i] == Dielectric) { needs_dielectric = true; break; }
  if (needs_dielectric) 
    FOR_ELECTRIC_COMPONENTS(c) if (v.has_field(c)) {
      if (data.ninveps == 3) abort("more than 3 field components??");
      data.inveps_cs[data.ninveps] = c;
      data.inveps_ds[data.ninveps] = component_direction(c);
      ++data.ninveps;
    }
  
  data.offsets = new int[2 * num_fields];
  for (int i = 0; i < 2 * num_fields; ++i)
    data.offsets[i] = 0;
  
  loop_in_chunks(h5_output_chunkloop, (void *) &data, 
		 where, Dielectric, true, true);

  delete[] data.offsets;
  delete[] data.fields;
  delete[] data.ph;
  delete[] data.cS;
  delete[] data.buf;
  file->done_writing_chunks();
  finished_working();
}

/***************************************************************************/

void fields::output_hdf5(const char *dataname,
                         int num_fields, const component *components,
                         field_function fun,
                         const geometric_volume &where,
                         void *fun_data_,
			 h5file *file,
                         bool append_data,
                         bool single_precision,
			 const char *prefix)
{
  bool delete_file;
  if ((delete_file = !file))
    file = open_h5file(dataname, h5file::WRITE, prefix, true);

  int len = strlen(dataname) + 5;
  char *dataname2;

  dataname2 = new char[len];
  snprintf(dataname2, len, "%s%s", dataname, ".r");
  output_hdf5(file, dataname, num_fields, components, fun, 0, where,
              fun_data_, append_data, single_precision);
  snprintf(dataname2, len, "%s%s", dataname, ".i");
  output_hdf5(file, dataname, num_fields, components, fun, 1, where,
              fun_data_, append_data, single_precision);
  delete[] dataname2;

  if (delete_file) delete file;
}

/***************************************************************************/

typedef struct {
  field_rfunction fun;
  void *fun_data_;
} rintegrand_data;

static complex<double> rintegrand_fun(const complex<double> *fields,
                                     const vec &loc,
                                     void *data_)
{
  rintegrand_data *data = (rintegrand_data *) data_;
  return data->fun(fields, loc, data->fun_data_);
}

void fields::output_hdf5(const char *dataname,
                         int num_fields, const component *components,
                         field_rfunction fun,
                         const geometric_volume &where,
                         void *fun_data_,
			 h5file *file,
                         bool append_data,
                         bool single_precision,
			 const char *prefix)
{
  bool delete_file;
  if ((delete_file = !file))
    file = open_h5file(dataname, h5file::WRITE, prefix, true);

  rintegrand_data data; data.fun = fun; data.fun_data_ = fun_data_;
  output_hdf5(file, dataname, num_fields, components, rintegrand_fun, 0, where,
	      (void *) &data, append_data, single_precision);

  if (delete_file) delete file;
}

/***************************************************************************/

static complex<double> component_fun(const complex<double> *fields,
				     const vec &loc,
				     void *data_)
{
     (void) loc; // unused
     (void) data_; // unused
     return fields[0];
}

void fields::output_hdf5(component c,
			 const geometric_volume &where,
			 h5file *file,
			 bool append_data,
                         bool single_precision,
			 const char *prefix) {
  if (is_derived(int(c))) {
    output_hdf5(derived_component(c), 
		where, file, append_data, single_precision, prefix);
    return;
  }

  if (coordinate_mismatch(v.dim, c)) return;

  char dataname[256];
  bool has_imag = !is_real && c != Dielectric;

  bool delete_file;
  if ((delete_file = !file))
    file = open_h5file(component_name(c), h5file::WRITE, prefix, true);

  snprintf(dataname, 256, "%s%s", component_name(c), has_imag ? ".r" : "");
  output_hdf5(file, dataname, 1, &c, component_fun, 0, where, 0,
	      append_data, single_precision);
  if (has_imag) {
    snprintf(dataname, 256, "%s.i", component_name(c));
    output_hdf5(file, dataname, 1, &c, component_fun, 1, where, 0,
		append_data, single_precision);
  }

  if (delete_file) delete file;
}

/***************************************************************************/

static double poynting_fun(const complex<double> *fields,
		       const vec &loc, void *data_)
{
     (void) loc; // unused
     (void) data_; // unused
     return (real(conj(fields[0]) * fields[1])
	     - real(conj(fields[2])*fields[3]));
}

static double energy_fun(const complex<double> *fields,
		       const vec &loc, void *data_)
{
     (void) loc; // unused
     int nfields = *((int *) data_) / 2;
     double sum = 0;
     for (int k = 0; k < nfields; ++k)
       sum += real(conj(fields[2*k]) * fields[2*k+1]);
     return sum * 0.5;
}

void fields::output_hdf5(derived_component c,
			 const geometric_volume &where,
			 h5file *file,
			 bool append_data,
                         bool single_precision,
			 const char *prefix) {
  if (!is_derived(int(c))) {
    output_hdf5(component(c), 
		where, file, append_data, single_precision, prefix);
    return;
  }

  if (coordinate_mismatch(v.dim, c)) return;

  field_rfunction fun;
  int nfields;
  component cs[6];

  switch (c) {
  case Sx: case Sy: case Sz: case Sr: case Sp:
    switch (c) {
    case Sx: cs[0] = Ey; cs[1] = Hz; break;
    case Sy: cs[0] = Ez; cs[1] = Hx; break;
    case Sz: cs[0] = Ex; cs[1] = Hy; break;
    case Sr: cs[0] = Ep; cs[1] = Hz; break;
    case Sp: cs[0] = Ez; cs[1] = Hr; break;
    }
    nfields = 4;
    cs[2] = direction_component(Ex, component_direction(cs[1]));
    cs[3] = direction_component(Hx, component_direction(cs[0]));
    fun = poynting_fun;
    break;
  case EnergyDensity: case D_EnergyDensity: case H_EnergyDensity:
    nfields = 0;
    fun = energy_fun;
    if (c != H_EnergyDensity)
      FOR_ELECTRIC_COMPONENTS(c0) if (v.has_field(c0)) {
	cs[nfields++] = c0;
	cs[nfields++] = direction_component(Dx, component_direction(c0));
      }
    if (c != D_EnergyDensity)
      FOR_MAGNETIC_COMPONENTS(c0) if (v.has_field(c0)) {
	cs[nfields++] = c0;
	cs[nfields++] = c0;
      }
    if (nfields > 6) abort("too many field components");
    break;
  default: 
    abort("unknown derived_component in output_hdf5");
  }

  output_hdf5(component_name(c), nfields, cs, fun, where, &nfields,
	      file, append_data, single_precision, prefix);
}

/***************************************************************************/

const char *fields::h5file_name(const char *name, 
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
  return filename;
}

h5file *fields::open_h5file(const char *name, h5file::access_mode mode, 
			    const char *prefix, bool timestamp)
{
  const char *filename = h5file_name(name, prefix, timestamp);
  if (!quiet && mode == h5file::WRITE)
    master_printf("creating output file \"%s\"...\n", filename);
  return new h5file(filename, mode, true);
}

} // namespace meep
