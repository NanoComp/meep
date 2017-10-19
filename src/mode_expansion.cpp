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

/* given an arbitrary field configuration, compute the coefficients
   in an expansion of F as a linear combination of normal modes
   of the geometry as computed by mpb.
*/


#include <stdio.h>
#include <string.h>
#include <math.h>

#include "meep_internals.hpp"

using namespace std;

typedef complex<double> cdouble;

namespace meep {

/***************************************************************************/

typedef struct {

  // information related to the volume covered by the
  // array slice (its size, etcetera) 
  // these fields are filled in by get_array_slice_dimensions
  // if the data parameter is non-null
  ivec min_corner, max_corner;
  int num_chunks;
  int rank;
  direction ds[3];
  int slice_size;

  // the function to output and related info (offsets for averaging, etc.)
  // note: either fun *or* rfun should be non-NULL (not both)
  field_function fun;
  field_rfunction rfun;
  void *fun_data;
  std::vector<component> components;

  void *vslice;

  // temporary internal storage buffers
  component *cS;
  cdouble *ph;
  cdouble *fields;
  int *offsets;

  int ninveps;
  component inveps_cs[3];
  direction inveps_ds[3];

  int ninvmu;
  component invmu_cs[3];
  direction invmu_ds[3];

} mode_projection_data;

#define UNUSED(x) (void) x // silence compiler warnings

/***************************************************************/
/* callback function passed to loop_in_chunks to evaluate the  */
/* projection of the user-specified fields onto one or more    */
/* normal modes as computed by mpb                             */
/***************************************************************/
static void mode_projection_chunkloop(fields_chunk *fc, 
                                      int ichnk, component cgrid,
				      ivec is, ivec ie,
				      vec s0, vec s1, vec e0, vec e1,
				      double dV0, double dV1,
				      ivec shift, complex<double> shift_phase,
				      const symmetry &S, int sn,
				      void *data_)
{
  UNUSED(ichnk);UNUSED(cgrid);UNUSED(s0);UNUSED(s1);UNUSED(e0);UNUSED(e1);
  UNUSED(dV0);UNUSED(dV1);UNUSED(shift_phase); UNUSED(fc);
  mode_projection_data *data = (mode_projection_data *) data_;

  ivec isS = S.transform(is, sn) + shift;
  ivec ieS = S.transform(ie, sn) + shift;
  data->min_corner = min(data->min_corner, min(isS, ieS));
  data->max_corner = max(data->max_corner, max(isS, ieS));
  data->num_chunks++;

}

/***************************************************************/
/* callback function passed to loop_in_chunks to fill array slice */
/***************************************************************/
static void get_array_slice_chunkloop(fields_chunk *fc, int ichnk, component cgrid,
	      			  ivec is, ivec ie, vec s0, vec s1, vec e0, vec e1,
				  double dV0, double dV1,
	  			  ivec shift, complex<double> shift_phase,
				  const symmetry &S, int sn, void *data_)
{
  UNUSED(ichnk);UNUSED(cgrid);UNUSED(s0);UNUSED(s1);UNUSED(e0);UNUSED(e1);
  UNUSED(dV0);UNUSED(dV1);
  array_slice_data *data = (array_slice_data *) data_;

  //-----------------------------------------------------------------------//
  // Find output chunk dimensions and strides, etc.

  int count[3]={1,1,1}, offset[3]={0,0,0}, stride[3]={1,1,1};

  ivec isS = S.transform(is, sn) + shift;
  ivec ieS = S.transform(ie, sn) + shift;
  
  // figure out what yucky_directions (in LOOP_OVER_IVECS)
  // correspond to what directions in the transformed vectors (in output).
  ivec permute(zero_ivec(fc->gv.dim));
  for (int i = 0; i < 3; ++i) 
    permute.set_direction(fc->gv.yucky_direction(i), i);
  permute = S.transform_unshifted(permute, sn);
  LOOP_OVER_DIRECTIONS(permute.dim, d)
    permute.set_direction(d, abs(permute.in_direction(d)));
  
  // compute the size of the chunk to output, and its strides etc.
  for (int i = 0; i < data->rank; ++i) {
    direction d = data->ds[i];
    int isd = isS.in_direction(d), ied = ieS.in_direction(d);
    count[i] = abs(ied - isd) / 2 + 1;
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
  // Compute the function to output, exactly as in fields::integrate.
  int *off = data->offsets;
  component *cS = data->cS;
  complex<double> *fields = data->fields, *ph = data->ph;
  const component *iecs = data->inveps_cs;
  const direction *ieds = data->inveps_ds;
  int ieos[6];
  const component *imcs = data->invmu_cs;
  const direction *imds = data->invmu_ds;
  int imos[6];
  int num_components=data->components.size();
  
  double *slice=0;
  cdouble *zslice=0;
  bool complex_data = (data->rfun==0);
  if (complex_data)
   zslice = (cdouble *)data->vslice;
  else
   slice = (double *)data->vslice;

  for (int i = 0; i < num_components; ++i) {
    cS[i] = S.transform(data->components[i], -sn);
    if (cS[i] == Dielectric || cS[i] == Permeability)
      ph[i] = 1.0;
    else {
      fc->gv.yee2cent_offsets(cS[i], off[2*i], off[2*i+1]);
      ph[i] = shift_phase * S.phase_shift(cS[i], sn);
    }
  }
  for (int k = 0; k < data->ninveps; ++k)
    fc->gv.yee2cent_offsets(iecs[k], ieos[2*k], ieos[2*k+1]);
  for (int k = 0; k < data->ninvmu; ++k)
    fc->gv.yee2cent_offsets(imcs[k], imos[2*k], imos[2*k+1]);

  vec rshift(shift * (0.5*fc->gv.inva));
  LOOP_OVER_IVECS(fc->gv, is, ie, idx) {
    IVEC_LOOP_LOC(fc->gv, loc);
    loc = S.transform(loc, sn) + rshift;

    for (int i = 0; i < num_components; ++i) {
      if (cS[i] == Dielectric) {
	double tr = 0.0;
	for (int k = 0; k < data->ninveps; ++k) {
	  const realnum *ie = fc->s->chi1inv[iecs[k]][ieds[k]];
	  if (ie) tr += (ie[idx] + ie[idx+ieos[2*k]] + ie[idx+ieos[1+2*k]]
			 + ie[idx+ieos[2*k]+ieos[1+2*k]]);
	  else tr += 4; // default inveps == 1
	}
	fields[i] = (4 * data->ninveps) / tr;
      }
      else if (cS[i] == Permeability) {
	double tr = 0.0;
	for (int k = 0; k < data->ninvmu; ++k) {
	  const realnum *im = fc->s->chi1inv[imcs[k]][imds[k]];
	  if (im) tr += (im[idx] + im[idx+imos[2*k]] + im[idx+imos[1+2*k]]
			 + im[idx+imos[2*k]+imos[1+2*k]]);
	  else tr += 4; // default invmu == 1
	}
	fields[i] = (4 * data->ninvmu) / tr;
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
    int idx2 = ((((offset[0] + offset[1] + offset[2])
                   + loop_i1 * stride[0])
                   + loop_i2 * stride[1])
                   + loop_i3 * stride[2]);

    if (complex_data)
     zslice[idx2] = data->fun(fields, loc, data->fun_data);
    else
     slice[idx2]  = data->rfun(fields, loc, data->fun_data);

  };

}

/***************************************************************/
/* given a volume, fill in the dims[] and directions[] arrays  */
/* describing the array slice needed to store field data for   */
/* all grid points in the volume.                              */
/*                                                             */
/* return value is rank of array slice.                        */
/*                                                             */
/* if caller_data is non-NULL, it should point to a            */
/* caller-allocated array_slice_data structure which will be   */
/* initialized appopriately for subsequent use in              */
/* get_array_slice.                                            */
/***************************************************************/
int fields::get_mode_expansion_coefficients(dft_flux flux,
                                            int
{
  am_now_working_on(FieldOutput);

  // use a local data structure if the caller didn't provide one
  array_slice_data local_data;
  array_slice_data *data=(array_slice_data *)caller_data;
  if (data==0)
   data=&local_data;

  data->min_corner = gv.round_vec(where.get_max_corner()) + one_ivec(gv.dim);
  data->max_corner = gv.round_vec(where.get_min_corner()) - one_ivec(gv.dim);
  data->num_chunks = 0;

  loop_in_chunks(get_array_slice_dimensions_chunkloop,
                 (void *) data, where, Centered, true, true);

  data->max_corner = max_to_all(data->max_corner);
  data->min_corner = -max_to_all(-data->min_corner); // i.e., min_to_all
  data->num_chunks = sum_to_all(data->num_chunks);
  if (data->num_chunks == 0 || !(data->min_corner <= data->max_corner))
    return 0; // no data to write;

  int rank=0, slice_size=1;
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    if (rank >= 3) abort("too many dimensions in array_slice");
    int n = (data->max_corner.in_direction(d)
	     - data->min_corner.in_direction(d)) / 2 + 1;
    if (n > 1) {
      data->ds[rank] = d;
      dims[rank++] = n;
      slice_size *= n;
    }
  }
  data->rank=rank;
  data->slice_size=slice_size;
  finished_working();

  return rank;
}

/***************************************************************/
/* precisely one of fun, rfun should be non-NULL               */
/***************************************************************/
void *fields::do_get_array_slice(const volume &where,
                                 std::vector<component> components,
                                 field_function fun,
                                 field_rfunction rfun,
                                 void *fun_data,
                                 void *vslice) {

  am_now_working_on(FieldOutput);

  /***************************************************************/
  /* call get_array_slice_dimensions to get slice dimensions and */
  /* partially initialze an array_slice_data struct              */
  /***************************************************************/
  int dims[3];
  array_slice_data data;
  int rank=get_array_slice_dimensions(where, dims, &data);
  int slice_size=data.slice_size;
  if (rank==0 || slice_size==0) return 0; // no data to write

  bool complex_data = (rfun==0);
  cdouble *zslice;
  double *slice;
  if (vslice==0)
   { if (complex_data)
      { zslice = new cdouble[slice_size];
        vslice = (void *)zslice;
      } 
     else
      { slice  = new double[slice_size];
        vslice = (void *)slice;
      };
   };
   
  data.vslice     = vslice;
  data.fun        = fun;
  data.rfun       = rfun;
  data.fun_data   = fun_data;
  data.components = components;

  int num_components = components.size();

  data.cS      = new component[num_components];
  data.ph      = new cdouble[num_components];
  data.fields  = new cdouble[num_components];
  
  data.offsets = new int[2 * num_components];
  for (int i = 0; i < 2 * num_components; ++i)
    data.offsets[i] = 0;

  /* compute inverse-epsilon directions for computing Dielectric fields */
  data.ninveps = 0;
  bool needs_dielectric = false;
  for (int i = 0; i < num_components; ++i)
    if (components[i] == Dielectric) { needs_dielectric = true; break; }
  if (needs_dielectric) 
    FOR_ELECTRIC_COMPONENTS(c) if (gv.has_field(c)) {
      if (data.ninveps == 3) abort("more than 3 field components??");
      data.inveps_cs[data.ninveps] = c;
      data.inveps_ds[data.ninveps] = component_direction(c);
      ++data.ninveps;
    }
  
  /* compute inverse-mu directions for computing Permeability fields */
  data.ninvmu = 0;
  bool needs_permeability = false;
  for (int i = 0; i < num_components; ++i)
    if (components[i] == Permeability) { needs_permeability = true; break; }
  if (needs_permeability) 
    FOR_MAGNETIC_COMPONENTS(c) if (gv.has_field(c)) {
      if (data.ninvmu == 3) abort("more than 3 field components??");
      data.invmu_cs[data.ninvmu] = c;
      data.invmu_ds[data.ninvmu] = component_direction(c);
      ++data.ninvmu;
    }
 
  loop_in_chunks(get_array_slice_chunkloop, (void *) &data,
		 where, Centered, true, true);

  /***************************************************************/
  /* repeatedly call sum_to_all to consolidate full array slice  */
  /* on all cores                                                */
  /***************************************************************/
#define BUFSIZE 1<<16 // use 64k buffer
  if (complex_data)
   { cdouble *buffer = new cdouble[BUFSIZE];
     cdouble *slice = (cdouble *)vslice;
     int offset=0, remaining=slice_size;
     while(remaining!=0)
      { 
        int size = (remaining > BUFSIZE ? BUFSIZE : remaining);
        sum_to_all(slice + offset, buffer, size);
        memcpy(slice+offset, buffer, size*sizeof(cdouble));
        remaining-=size;
        offset+=size;
      };
     delete[] buffer;
   }
  else
   { double *buffer = new double[BUFSIZE];
     double *slice = (double *)vslice;
     int offset=0, remaining=slice_size;
     while(remaining!=0)
      { int size = (remaining > BUFSIZE ? BUFSIZE : remaining);
        sum_to_all(slice + offset, buffer, size);
        memcpy(slice+offset, buffer, size*sizeof(double));
        remaining-=size;
        offset+=size;
      };
     delete[] buffer;
   };

  delete[] data.offsets;
  delete[] data.fields;
  delete[] data.ph;
  delete[] data.cS;
  finished_working();

  return vslice;
}

/***************************************************************/
/* entry points to get_array_slice                             */
/***************************************************************/
double *fields::get_array_slice(const volume &where,
                                std::vector<component> components,
                                field_rfunction rfun, void *fun_data,
                                double *slice)
{
  return (double *)do_get_array_slice(where, components,
                                      0, rfun, fun_data,
                                      (void *)slice);
} 

cdouble *fields::get_complex_array_slice(const volume &where,
                                         std::vector<component> components,
                                         field_function fun, void *fun_data,
                                         cdouble *slice)
{
  return (cdouble *)do_get_array_slice(where, components,
                                       fun, 0, fun_data,
                                       (void *)slice);
}

double *fields::get_array_slice(const volume &where, component c,
                                double *slice, int slice_length)
{
  (void) slice_length;
  std::vector<component> components(1);
  components[0]=c;
  return (double *)do_get_array_slice(where, components,
                                      0, default_field_rfunc, 0,
                                      (void *)slice);
}

double *fields::get_array_slice(const volume &where,
                                derived_component c,
                                double *slice, int slice_length)
{
  (void) slice_length;
  int nfields;
  component carray[12];
  field_rfunction rfun = derived_component_func(c, gv, nfields, carray);
  std::vector<component> cs(carray, carray+nfields);
  return (double *)do_get_array_slice(where, cs,
                                      0, rfun, &nfields,
                                      (void *)slice);
}

cdouble *fields::get_complex_array_slice(const volume &where, component c,
                                         cdouble *slice, int slice_length)
{
  (void) slice_length;
  std::vector<component> components(1);
  components[0]=c;
  return (cdouble *)do_get_array_slice(where, components,
                                       default_field_func, 0, 0,
                                       (void *)slice);
}

} // namespace meep
