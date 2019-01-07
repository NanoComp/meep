/* Copyright (C) 2005-2019 Massachusetts Institute of Technology
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

/* create and return arrays of field components on user-specified
   spatial slices. Uses fields::loop_in_chunks analogous to
   h5fields.cpp
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
  size_t slice_size;

  // the function to output and related info (offsets for averaging, etc.)
  // note: either fun *or* rfun should be non-NULL (not both)
  field_function fun;
  field_rfunction rfun;
  void *fun_data;
  std::vector<component> components;
  component source_slice_component;
  bool get_source_slice;

  void *vslice;

  // temporary internal storage buffers
  component *cS;
  cdouble *ph;
  cdouble *fields;
  ptrdiff_t *offsets;

  int ninveps;
  component inveps_cs[3];
  direction inveps_ds[3];

  int ninvmu;
  component invmu_cs[3];
  direction invmu_ds[3];

} array_slice_data;

#define UNUSED(x) (void) x // silence compiler warnings

/* passthrough field function equivalent to component_fun in h5fields.cpp */
static cdouble default_field_func(const cdouble *fields,
				  const vec &loc, void *data_)
{
  (void) loc; // unused
  (void) data_; // unused
  return fields[0];
}

static double default_field_rfunc(const cdouble *fields,
				  const vec &loc, void *data_)
{
  (void) loc; // unused
  (void) data_; // unused
  return real(fields[0]);
}

/***************************************************************/
/* callback function passed to loop_in_chunks to compute       */
/* dimensions of array slice                                   */
/***************************************************************/
static void get_array_slice_dimensions_chunkloop(fields_chunk *fc, int ichnk, component cgrid,
				  ivec is, ivec ie,
				  vec s0, vec s1, vec e0, vec e1,
				  double dV0, double dV1,
				  ivec shift, complex<double> shift_phase,
				  const symmetry &S, int sn,
				  void *data_)
{
  UNUSED(ichnk);UNUSED(cgrid);UNUSED(s0);UNUSED(s1);UNUSED(e0);UNUSED(e1);
  UNUSED(dV0);UNUSED(dV1);UNUSED(shift_phase); UNUSED(fc);
  array_slice_data *data = (array_slice_data *) data_;
  ivec isS = S.transform(is, sn) + shift;
  ivec ieS = S.transform(ie, sn) + shift;
  data->min_corner = min(data->min_corner, min(isS, ieS));
  data->max_corner = max(data->max_corner, max(isS, ieS));
  data->num_chunks++;
}

/*****************************************************************/
/* populate the array slice with information about sources with  */
/* the component specified in source_slice_component.            */
/*****************************************************************/
void fill_chunk_source_slice(fields_chunk *fc, array_slice_data *data)
{
  ndim dim=fc->gv.dim;
  if (dim==Dcyl)
   { fprintf(stderr,"warning: source slices not implemented for cylindrical coordinates; array will be all zeros\n");
     return;
   }

  component slice_component = data->source_slice_component;
  cdouble *slice = (cdouble *)data->vslice;
  ivec slice_imin=data->min_corner, slice_imax=data->max_corner;

  ptrdiff_t NY=1, NZ=1;
  if (has_direction(fc->gv.dim,Z))
   NZ = ((slice_imax-slice_imin).in_direction(Z)/2) + 1;
  if (has_direction(fc->gv.dim,Y))
   NY = ((slice_imax-slice_imin).in_direction(Y)/2) + 1;

  for(int ft=0; ft<NUM_FIELD_TYPES; ft++)
   for(src_vol *s=fc->sources[ft]; s; s=s->next)
    {
      if (slice_component!=s->c)
       continue;

      // loop over point sources in this src_vol. for each point source,
      // the src_vol stores the amplitude and the global index of the grid point,
      // from which we need to compute the local index of the grid point within the
      // slice (that is, if it even lies within the slice)
      for(size_t npt=0; npt<s->npts; npt++)
       { cdouble amp = s->A[npt];
         ptrdiff_t chunk_index = s->index[npt];
         ivec iloc = fc->gv.iloc(Dielectric, chunk_index);
         if (iloc<slice_imin || iloc>slice_imax) continue; // source point outside slice
         ivec slice_offset = iloc-slice_imin;
         ptrdiff_t slice_index=0;
         if (has_direction(dim,Z))
          slice_index += slice_offset.in_direction(Z)/2;
         if (has_direction(dim,Y))
          slice_index += NZ*slice_offset.in_direction(Y)/2;
         if (has_direction(dim,X))
          slice_index += NY*NZ*slice_offset.in_direction(X)/2;

         slice[slice_index] = amp;
       }
    }
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
  //- If we're fetching a 'source slice,' branch off here to handle that.  //
  //-----------------------------------------------------------------------//
  if (data->get_source_slice)
   { fill_chunk_source_slice(fc,data);
     return;
   }

  //-----------------------------------------------------------------------//
  // Find output chunk dimensions and strides, etc.

  int start[3]={0,0,0}, count[3]={1,1,1};
  ptrdiff_t offset[3]={0,0,0};

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
  size_t slice_size=1;
  for (int i = 0; i < data->rank; ++i) {
    direction d = data->ds[i];
    int isd = isS.in_direction(d), ied = ieS.in_direction(d);
    start[i] = (min(isd, ied) - data->min_corner.in_direction(d)) / 2;
    count[i] = abs(ied - isd) / 2 + 1;
    slice_size *= count[i];
    if (ied < isd) offset[permute.in_direction(d)] = count[i] - 1;
  }

  // slightly confusing: for array_slice, in contrast to
  // h5fields, strides are computed using the dimensions of
  // the full array slice, not the dimensions of the chunk.
  size_t dims[3]={1,1,1};
  for (int i = 0; i<data->rank; i++) {
    direction d = data->ds[i];
    dims[i]= (data->max_corner.in_direction(d)
	     - data->min_corner.in_direction(d)) / 2 + 1;
   }

  ptrdiff_t stride[3]={1,1,1};
  for (int i = 0; i < data->rank; ++i) {
    direction d = data->ds[i];
    int j = permute.in_direction(d);
    for (int k = i + 1; k < data->rank; ++k) stride[j] *= dims[k];
    offset[j] *= stride[j];
    if (offset[j]) stride[j] *= -1;
  }

  // sco="slice chunk offset"
  ptrdiff_t sco=start[0]*dims[1]*dims[2] + start[1]*dims[2] + start[2];

  //-----------------------------------------------------------------------//
  // Otherwise proceed to compute the function of field components to be   //
  // tabulated on the slice, exactly as in fields::integrate.              //
  //-----------------------------------------------------------------------//
  double *slice=0;
  cdouble *zslice=0;
  bool complex_data = (data->rfun==0);
  if (complex_data)
   zslice = (cdouble *)data->vslice;
  else
   slice = (double *)data->vslice;

  ptrdiff_t *off = data->offsets;
  component *cS = data->cS;
  complex<double> *fields = data->fields, *ph = data->ph;
  const component *iecs = data->inveps_cs;
  const direction *ieds = data->inveps_ds;
  ptrdiff_t ieos[6];
  const component *imcs = data->invmu_cs;
  const direction *imds = data->invmu_ds;
  ptrdiff_t imos[6];
  int num_components=data->components.size();

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
  // main loop over all grid points owned by this field chunk.
  LOOP_OVER_IVECS(fc->gv, is, ie, idx) {

    // get real-space coordinates of grid point, taking into
    // account the complications of symmetries.
    IVEC_LOOP_LOC(fc->gv, loc);
    loc = S.transform(loc, sn) + rshift;

    // interpolate fields at the four nearest grid points
    // to get the value of the field component for this point
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

    // compute the index into the array for this grid point and store the result of the computation
    ptrdiff_t idx2 = sco + ((((offset[0] + offset[1] + offset[2])
                            + loop_i1 * stride[0])
                            + loop_i2 * stride[1])
                            + loop_i3 * stride[2]);

    if (complex_data)
     zslice[idx2] = data->fun(fields, loc, data->fun_data);
    else
     slice[idx2]  = data->rfun(fields, loc, data->fun_data);

  } // LOOP_OVER_IVECS

}

/***************************************************************/
/* given a volume, fill in the dims[] and dirs[] arrays        */
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
int fields::get_array_slice_dimensions(const volume &where, size_t dims[3],
                                       direction dirs[3], bool collapse_empty_dimensions,
                                       void *caller_data)
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

  int rank=0;
  size_t slice_size=1;
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    if (rank >= 3) abort("too many dimensions in array_slice");
    size_t n = (data->max_corner.in_direction(d)
	     - data->min_corner.in_direction(d)) / 2 + 1;
    if (where.in_direction(d)==0.0 && collapse_empty_dimensions)
     n=1;
    if (n > 1) {
      data->ds[rank] = d;
      dims[rank++] = n;
      slice_size *= n;
    }
  }
  for(int r=0; r<rank; r++)
      dirs[r]=(meep::direction)(data->ds[r] - X);
  data->rank=rank;
  data->slice_size=slice_size;
  finished_working();

  return rank;
}

/**********************************************************************/
/* precisely one of fun, rfun, source_slice should be non-NULL / true */
/**********************************************************************/
void *fields::do_get_array_slice(const volume &where,
                                 std::vector<component> components,
                                 field_function fun,
                                 field_rfunction rfun,
                                 void *fun_data,
                                 void *vslice,
                                 component source_slice_component,
                                 bool get_source_slice) {
  am_now_working_on(FieldOutput);

  /***************************************************************/
  /* call get_array_slice_dimensions to get slice dimensions and */
  /* partially initialze an array_slice_data struct              */
  /***************************************************************/
  // by tradition, empty dimensions in time-domain field arrays are *not* collapsed;
  // TODO make this a caller-specifiable parameter to get_array_slice()?
  bool collapse=false;
  size_t dims[3];
  direction dirs[3];
  array_slice_data data;
  int rank=get_array_slice_dimensions(where, dims, dirs, collapse, &data);
  size_t slice_size=data.slice_size;
  if (rank==0 || slice_size==0) return 0; // no data to write

  bool complex_data = (rfun==0);
  cdouble *zslice;
  double *slice;
  if (vslice==0)
   { if (complex_data)
      { zslice = new cdouble[slice_size];
        memset(zslice,0,slice_size*sizeof(cdouble));
        vslice = (void *)zslice;
      }
     else
      { slice  = new double[slice_size];
        memset(slice,0,slice_size*sizeof(double));
        vslice = (void *)slice;
      }
   }

  data.vslice       = vslice;
  data.fun          = fun;
  data.rfun         = rfun;
  data.fun_data     = fun_data;
  data.components   = components;

  int num_components = components.size();

  data.cS      = new component[num_components];
  data.ph      = new cdouble[num_components];
  data.fields  = new cdouble[num_components];

  data.offsets = new ptrdiff_t[2 * num_components];
  memset(data.offsets, 0, 2*num_components*sizeof(ptrdiff_t));

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
     ptrdiff_t offset=0;
     size_t remaining=slice_size;
     while(remaining!=0)
      {
        size_t size = (remaining > BUFSIZE ? BUFSIZE : remaining);
        sum_to_all(slice + offset, buffer, size);
        memcpy(slice+offset, buffer, size*sizeof(cdouble));
        remaining-=size;
        offset+=size;
      }
     delete[] buffer;
   }
  else
   { double *buffer = new double[BUFSIZE];
     double *slice = (double *)vslice;
     ptrdiff_t offset=0;
     size_t remaining=slice_size;
     while(remaining!=0)
      { size_t size = (remaining > BUFSIZE ? BUFSIZE : remaining);
        sum_to_all(slice + offset, buffer, size);
        memcpy(slice+offset, buffer, size*sizeof(double));
        remaining-=size;
        offset+=size;
      }
     delete[] buffer;
   }

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
                                double *slice)
{
  std::vector<component> components(1);
  components[0]=c;
  return (double *)do_get_array_slice(where, components,
                                      0, default_field_rfunc, 0,
                                      (void *)slice);
}

double *fields::get_array_slice(const volume &where,
                                derived_component c,
                                double *slice)
{
  int nfields;
  component carray[12];
  field_rfunction rfun = derived_component_func(c, gv, nfields, carray);
  std::vector<component> cs(carray, carray+nfields);
  return (double *)do_get_array_slice(where, cs,
                                      0, rfun, &nfields,
                                      (void *)slice);
}

cdouble *fields::get_complex_array_slice(const volume &where, component c,
                                         cdouble *slice)
{
  std::vector<component> components(1);
  components[0]=c;
  return (cdouble *)do_get_array_slice(where, components,
                                       default_field_func, 0, 0,
                                       (void *)slice);
}


cdouble *fields::get_source_slice(const volume &where, component source_slice_component, cdouble *slice) {
  vector<component> cs; // empty
  return (cdouble *)do_get_array_slice(where,cs,0,0,0,(void *)slice,source_slice_component, true);
}

} // namespace meep
