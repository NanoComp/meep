/* Copyright (C) 2005-2020 Massachusetts Institute of Technology
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

  // if non-null, min_max_loc[0,1] are filled in by get_array_slice_dimensions_chunkloop
  // with the (coordinate-wise) minimum and maximum grid points encountered
  // in looping over the slice region.
  vec *min_max_loc;

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
  ptrdiff_t *offsets;

  double omega;

  int ninveps;
  component inveps_cs[3];
  direction inveps_ds[3];

  int ninvmu;
  component invmu_cs[3];
  direction invmu_ds[3];

} array_slice_data;

#define UNUSED(x) (void)x // silence compiler warnings

/* passthrough field function equivalent to component_fun in h5fields.cpp */
static cdouble default_field_func(const cdouble *fields, const vec &loc, void *data_) {
  (void)loc;   // unused
  (void)data_; // unused
  return fields[0];
}

static double default_field_rfunc(const cdouble *fields, const vec &loc, void *data_) {
  (void)loc;   // unused
  (void)data_; // unused
  return real(fields[0]);
}

/***************************************************************/
/* callback function passed to loop_in_chunks to compute       */
/* dimensions of array slice                                   */
/***************************************************************/
static void get_array_slice_dimensions_chunkloop(fields_chunk *fc, int ichnk, component cgrid,
                                                 ivec is, ivec ie, vec s0, vec s1, vec e0, vec e1,
                                                 double dV0, double dV1, ivec shift,
                                                 complex<double> shift_phase, const symmetry &S,
                                                 int sn, void *data_) {
  UNUSED(ichnk);
  UNUSED(cgrid);
  UNUSED(s0);
  UNUSED(s1);
  UNUSED(e0);
  UNUSED(e1);
  UNUSED(dV0);
  UNUSED(dV1);
  UNUSED(shift_phase);

  array_slice_data *data = (array_slice_data *)data_;
  ivec isS = S.transform(is, sn) + shift;
  ivec ieS = S.transform(ie, sn) + shift;
  data->min_corner = min(data->min_corner, min(isS, ieS));
  data->max_corner = max(data->max_corner, max(isS, ieS));
  data->num_chunks++;

  if (data->min_max_loc) {
    vec rshift(shift * (0.5 * fc->gv.inva));
    LOOP_OVER_IVECS(fc->gv, is, ie, idx) {
      IVEC_LOOP_LOC(fc->gv, loc);
      loc = S.transform(loc, sn) + rshift;
      data->min_max_loc[0] = min(loc, data->min_max_loc[0]);
      data->min_max_loc[1] = max(loc, data->min_max_loc[1]);
    }
  }
}

/***************************************************************/
/* chunkloop for get_source_slice ******************************/
/***************************************************************/
typedef struct {
  component source_component;
  ivec slice_imin, slice_imax;
  cdouble *slice;
} source_slice_data;

bool in_range(int imin, int i, int imax) { return (imin <= i && i <= imax); }

bool in_subgrid(ivec ivmin, ivec iv, ivec ivmax) {
  LOOP_OVER_DIRECTIONS(iv.dim, d)
  if (!in_range(ivmin.in_direction(d), iv.in_direction(d), ivmax.in_direction(d))) return false;
  return true;
}

static void get_source_slice_chunkloop(fields_chunk *fc, int ichnk, component cgrid, ivec is,
                                       ivec ie, vec s0, vec s1, vec e0, vec e1, double dV0,
                                       double dV1, ivec shift, complex<double> shift_phase,
                                       const symmetry &S, int sn, void *data_) {

  UNUSED(ichnk);
  UNUSED(cgrid);
  UNUSED(is);
  UNUSED(ie);
  UNUSED(s0);
  UNUSED(s1);
  UNUSED(e0);
  UNUSED(e1);
  UNUSED(dV0);
  UNUSED(dV1);
  UNUSED(shift_phase);

  source_slice_data *data = (source_slice_data *)data_;
  ivec slice_imin = data->slice_imin, slice_imax = data->slice_imax;
  ndim dim = fc->gv.dim;

  // the following works in all cases except cylindrical coordinates
  ptrdiff_t NY = 1, NZ = 1;
  if (has_direction(dim, Z)) NZ = ((slice_imax - slice_imin).in_direction(Z) / 2) + 1;
  if (has_direction(dim, Y)) NY = ((slice_imax - slice_imin).in_direction(Y) / 2) + 1;

  for (int ft = 0; ft < NUM_FIELD_TYPES; ft++)
    for (src_vol *s = fc->sources[ft]; s; s = s->next) {
      component cS = S.transform(data->source_component, -sn);
      if (s->c != cS) continue;

      // loop over point sources in this src_vol. for each point source,
      // the src_vol stores the amplitude and the global index of the
      // symmetry-parent grid point, from which we need to compute the
      // local index of the symmetry-child grid point within this
      // slice (that is, if it even lies within the slice)
      for (size_t npt = 0; npt < s->npts; npt++) {
        cdouble amp = s->A[npt];
        ptrdiff_t chunk_index = s->index[npt];
        ivec iloc_parent = fc->gv.iloc(Dielectric, chunk_index);
        ivec iloc_child = S.transform(iloc_parent, sn) + shift;
        if (!in_subgrid(slice_imin, iloc_child, slice_imax)) continue; // source point outside slice
        ivec slice_offset = iloc_child - slice_imin;
        ptrdiff_t slice_index = 0;
        // the following works to set the slice_index in all cases except cylindrical coordinates
        if (has_direction(dim, Z)) slice_index += slice_offset.in_direction(Z) / 2;
        if (has_direction(dim, Y)) slice_index += NZ * slice_offset.in_direction(Y) / 2;
        if (has_direction(dim, X)) slice_index += NY * NZ * slice_offset.in_direction(X) / 2;
        data->slice[slice_index] = amp;
      }
    }
}

/***************************************************************/
/* callback function passed to loop_in_chunks to fill array slice */
/***************************************************************/
static void get_array_slice_chunkloop(fields_chunk *fc, int ichnk, component cgrid, ivec is,
                                      ivec ie, vec s0, vec s1, vec e0, vec e1, double dV0,
                                      double dV1, ivec shift, complex<double> shift_phase,
                                      const symmetry &S, int sn, void *data_) {
  UNUSED(ichnk);
  UNUSED(cgrid);
  UNUSED(s0);
  UNUSED(s1);
  UNUSED(e0);
  UNUSED(e1);
  UNUSED(dV0);
  UNUSED(dV1);
  array_slice_data *data = (array_slice_data *)data_;

  //-----------------------------------------------------------------------//
  // Find output chunk dimensions and strides, etc.
  //-----------------------------------------------------------------------//

  int start[3] = {0, 0, 0}, count[3] = {1, 1, 1};
  ptrdiff_t offset[3] = {0, 0, 0};

  ivec isS = S.transform(is, sn) + shift;
  ivec ieS = S.transform(ie, sn) + shift;

  // figure out what yucky_directions (in LOOP_OVER_IVECS)
  // correspond to what directions in the transformed vectors (in output).
  ivec permute(zero_ivec(fc->gv.dim));
  for (int i = 0; i < 3; ++i)
    permute.set_direction(fc->gv.yucky_direction(i), i);
  permute = S.transform_unshifted(permute, sn);
  LOOP_OVER_DIRECTIONS(permute.dim, d) { permute.set_direction(d, abs(permute.in_direction(d))); }

  // compute the size of the chunk to output, and its strides etc.
  size_t slice_size = 1;
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
  size_t dims[3] = {1, 1, 1};
  for (int i = 0; i < data->rank; i++) {
    direction d = data->ds[i];
    dims[i] = (data->max_corner.in_direction(d) - data->min_corner.in_direction(d)) / 2 + 1;
  }

  ptrdiff_t stride[3] = {1, 1, 1};
  for (int i = 0; i < data->rank; ++i) {
    direction d = data->ds[i];
    int j = permute.in_direction(d);
    for (int k = i + 1; k < data->rank; ++k)
      stride[j] *= dims[k];
    offset[j] *= stride[j];
    if (offset[j]) stride[j] *= -1;
  }

  // sco="slice chunk offset"
  ptrdiff_t sco = start[0] * dims[1] * dims[2] + start[1] * dims[2] + start[2];

  //-----------------------------------------------------------------------//
  // Otherwise proceed to compute the function of field components to be   //
  // tabulated on the slice, exactly as in fields::integrate.              //
  //-----------------------------------------------------------------------//
  double *slice = 0;
  cdouble *zslice = 0;
  bool complex_data = (data->rfun == 0);
  if (complex_data)
    zslice = (cdouble *)data->vslice;
  else
    slice = (double *)data->vslice;

  ptrdiff_t *off = data->offsets;
  component *cS = data->cS;
  double omega = data->omega;
  complex<double> *fields = data->fields, *ph = data->ph;
  const component *iecs = data->inveps_cs;
  const direction *ieds = data->inveps_ds;
  ptrdiff_t ieos[6];
  const component *imcs = data->invmu_cs;
  const direction *imds = data->invmu_ds;
  ptrdiff_t imos[6];
  int num_components = data->components.size();

  for (int i = 0; i < num_components; ++i) {
    cS[i] = S.transform(data->components[i], -sn);
    if (cS[i] == Dielectric || cS[i] == Permeability)
      ph[i] = 1.0;
    else {
      fc->gv.yee2cent_offsets(cS[i], off[2 * i], off[2 * i + 1]);
      ph[i] = shift_phase * S.phase_shift(cS[i], sn);
    }
  }
  for (int k = 0; k < data->ninveps; ++k)
    fc->gv.yee2cent_offsets(iecs[k], ieos[2 * k], ieos[2 * k + 1]);
  for (int k = 0; k < data->ninvmu; ++k)
    fc->gv.yee2cent_offsets(imcs[k], imos[2 * k], imos[2 * k + 1]);

  vec rshift(shift * (0.5 * fc->gv.inva));
  // main loop over all grid points owned by this field chunk.
  LOOP_OVER_IVECS(fc->gv, is, ie, idx) {

    // get real-space coordinates of grid point, taking into
    // account the complications of symmetries.
    IVEC_LOOP_LOC(fc->gv, loc);
    loc = S.transform(loc, sn) + rshift;

    // interpolate fields at the four nearest grid points
    // to get the value of the field component for this point
    for (int i = 0; i < num_components; ++i) {
      // special case for fetching grid point coordinates and weights
      if (cS[i] == NO_COMPONENT) {
        fields[i] = IVEC_LOOP_WEIGHT(s0, s1, e0, e1, dV0 + dV1 * loop_i2);
      }
      else if (cS[i] == Dielectric) {
        complex<double> tr(0.0,0.0);
        for (int k = 0; k < data->ninveps; ++k) {
          tr += (fc->s->get_chi1inv_at_pt(iecs[k], ieds[k], idx, omega) +
                 fc->s->get_chi1inv_at_pt(iecs[k], ieds[k], idx + ieos[2 * k], omega) +
                 fc->s->get_chi1inv_at_pt(iecs[k], ieds[k], idx + ieos[1 + 2 * k], omega) +
                 fc->s->get_chi1inv_at_pt(iecs[k], ieds[k], idx + ieos[2 * k] + ieos[1 + 2 * k],
                                          omega));
          if (abs(tr) == 0.0) tr += 4.0; // default inveps == 1
        }
        fields[i] = (4.0 * data->ninveps) / tr;
      }
      else if (cS[i] == Permeability) {
        complex<double> tr(0.0,0.0);
        for (int k = 0; k < data->ninvmu; ++k) {
          tr += (fc->s->get_chi1inv_at_pt(imcs[k], imds[k], idx, omega) +
                 fc->s->get_chi1inv_at_pt(imcs[k], imds[k], idx + imos[2 * k], omega) +
                 fc->s->get_chi1inv_at_pt(imcs[k], imds[k], idx + imos[1 + 2 * k], omega) +
                 fc->s->get_chi1inv_at_pt(imcs[k], imds[k], idx + imos[2 * k] + imos[1 + 2 * k],
                                          omega));
          if (abs(tr) == 0.0) tr += 4.0; // default invmu == 1
        }
        fields[i] = (4.0 * data->ninvmu) / tr;
      }
      else {
        double f[2];
        for (int k = 0; k < 2; ++k)
          if (fc->f[cS[i]][k])
            f[k] = 0.25 * (fc->f[cS[i]][k][idx] + fc->f[cS[i]][k][idx + off[2 * i]] +
                           fc->f[cS[i]][k][idx + off[2 * i + 1]] +
                           fc->f[cS[i]][k][idx + off[2 * i] + off[2 * i + 1]]);
          else
            f[k] = 0;
        fields[i] = complex<double>(f[0], f[1]) * ph[i];
      }
    }

    // compute the index into the array for this grid point and store the result of the computation
    ptrdiff_t idx2 =
        sco + ((((offset[0] + offset[1] + offset[2]) + loop_i1 * stride[0]) + loop_i2 * stride[1]) +
               loop_i3 * stride[2]);

    if (complex_data)
      zslice[idx2] = data->fun(fields, loc, data->fun_data);
    else
      slice[idx2] = data->rfun(fields, loc, data->fun_data);

  } // LOOP_OVER_IVECS
}

/***************************************************************/
/* repeatedly call sum_to_all to consolidate all entries of    */
/* an array on all cores.                                      */
/***************************************************************/
#define BUFSIZE 1 << 16 // use 64k buffer
double *array_to_all(double *array, size_t array_size) {
  double *buffer = new double[BUFSIZE];
  ptrdiff_t offset = 0;
  size_t remaining = array_size;
  while (remaining != 0) {
    size_t xfer_size = (remaining > BUFSIZE ? BUFSIZE : remaining);
    sum_to_all(array + offset, buffer, xfer_size);
    memcpy(array + offset, buffer, xfer_size * sizeof(double));
    remaining -= xfer_size;
    offset += xfer_size;
  }
  delete[] buffer;
  return array;
}

cdouble *array_to_all(cdouble *array, size_t array_size) {
  return (cdouble *)array_to_all((double *)array, 2 * array_size);
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
int fields::get_array_slice_dimensions(const volume &where, size_t dims[3], direction dirs[3],
                                       bool collapse_empty_dimensions, bool snap_empty_dimensions,
                                       vec *min_max_loc, void *caller_data) {
  am_now_working_on(FieldOutput);

  // use a local data structure if the caller didn't provide one
  array_slice_data local_data;
  array_slice_data *data = (array_slice_data *)caller_data;
  if (data == 0) data = &local_data;

  data->min_corner = gv.round_vec(where.get_max_corner()) + one_ivec(gv.dim);
  data->max_corner = gv.round_vec(where.get_min_corner()) - one_ivec(gv.dim);
  data->num_chunks = 0;

  data->min_max_loc = min_max_loc;
  vec *min_loc = 0, *max_loc = 0;
  if (min_max_loc) {
    min_loc = min_max_loc + 0;
    max_loc = min_max_loc + 1;
    LOOP_OVER_DIRECTIONS(gv.dim, d) {
      min_loc->set_direction(d, +infinity);
      max_loc->set_direction(d, -infinity);
    }
  }

  bool use_symmetry = true;
  loop_in_chunks(get_array_slice_dimensions_chunkloop, (void *)data, where, Centered, use_symmetry,
                 snap_empty_dimensions);

  data->min_corner = -max_to_all(-data->min_corner); // i.e., min_to_all
  data->max_corner = max_to_all(data->max_corner);
  if (min_max_loc) LOOP_OVER_DIRECTIONS(gv.dim, d) {
      min_loc->set_direction(d, -1.0 * max_to_all(-1.0 * min_loc->in_direction(d)));
      max_loc->set_direction(d, max_to_all(max_loc->in_direction(d)));
    }
  data->num_chunks = sum_to_all(data->num_chunks);
  if (data->num_chunks == 0 || !(data->min_corner <= data->max_corner))
    return 0; // no data to write;

  int rank = 0;
  size_t slice_size = 1;
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    if (rank >= 3) abort("too many dimensions in array_slice");
    size_t n = (data->max_corner.in_direction(d) - data->min_corner.in_direction(d)) / 2 + 1;
    if (where.in_direction(d) == 0.0 && collapse_empty_dimensions) n = 1;
    if (n > 1) {
      data->ds[rank] = d;
      dims[rank++] = n;
      slice_size *= n;
    }
  }
  for (int r = 0; r < rank; r++)
    dirs[r] = data->ds[r];
  data->rank = rank;
  data->slice_size = slice_size;
  finished_working();

  return rank;
}

/**********************************************************************/
/* precisely one of fun, rfun, should be non-NULL                     */
/**********************************************************************/
void *fields::do_get_array_slice(const volume &where, std::vector<component> components,
                                 field_function fun, field_rfunction rfun, void *fun_data,
                                 void *vslice, double omega) {
  am_now_working_on(FieldOutput);

  /***************************************************************/
  /* call get_array_slice_dimensions to get slice dimensions and */
  /* partially initialze an array_slice_data struct              */
  /***************************************************************/
  // by tradition, empty dimensions in time-domain field arrays are *not* collapsed
  // TODO make this a caller-specifiable parameter to get_array_slice()?
  bool collapse = false, snap = true;
  size_t dims[3];
  direction dirs[3];
  array_slice_data data;
  get_array_slice_dimensions(where, dims, dirs, collapse, snap, 0, &data);
  size_t slice_size = data.slice_size;

  bool complex_data = (rfun == 0);
  cdouble *zslice;
  double *slice;
  if (vslice == 0) {
    if (complex_data) {
      zslice = new cdouble[slice_size];
      memset(zslice, 0, slice_size * sizeof(cdouble));
      vslice = (void *)zslice;
    }
    else {
      slice = new double[slice_size];
      memset(slice, 0, slice_size * sizeof(double));
      vslice = (void *)slice;
    }
  }

  data.vslice = vslice;
  data.fun = fun;
  data.rfun = rfun;
  data.fun_data = fun_data;
  data.components = components;
  data.omega = omega;
  int num_components = components.size();
  data.cS = new component[num_components];
  data.ph = new cdouble[num_components];
  data.fields = new cdouble[num_components];
  data.offsets = new ptrdiff_t[2 * num_components];
  memset(data.offsets, 0, 2 * num_components * sizeof(ptrdiff_t));

  /* compute inverse-epsilon directions for computing Dielectric fields */
  data.ninveps = 0;
  bool needs_dielectric = false;
  for (int i = 0; i < num_components; ++i)
    if (components[i] == Dielectric) {
      needs_dielectric = true;
      break;
    }
  if (needs_dielectric) FOR_ELECTRIC_COMPONENTS(c) if (gv.has_field(c)) {
      if (data.ninveps == 3) abort("more than 3 field components??");
      data.inveps_cs[data.ninveps] = c;
      data.inveps_ds[data.ninveps] = component_direction(c);
      ++data.ninveps;
    }

  /* compute inverse-mu directions for computing Permeability fields */
  data.ninvmu = 0;
  bool needs_permeability = false;
  for (int i = 0; i < num_components; ++i)
    if (components[i] == Permeability) {
      needs_permeability = true;
      break;
    }
  if (needs_permeability) FOR_MAGNETIC_COMPONENTS(c) if (gv.has_field(c)) {
      if (data.ninvmu == 3) abort("more than 3 field components??");
      data.invmu_cs[data.ninvmu] = c;
      data.invmu_ds[data.ninvmu] = component_direction(c);
      ++data.ninvmu;
    }

  loop_in_chunks(get_array_slice_chunkloop, (void *)&data, where, Centered, true, true);

  /***************************************************************/
  /*consolidate full array on all cores                          */
  /***************************************************************/
  vslice = (void *)array_to_all((double *)vslice, (complex_data ? 2 : 1) * slice_size);

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
double *fields::get_array_slice(const volume &where, std::vector<component> components,
                                field_rfunction rfun, void *fun_data, double *slice, double omega) {
  return (double *)do_get_array_slice(where, components, 0, rfun, fun_data, (void *)slice, omega);
}

cdouble *fields::get_complex_array_slice(const volume &where, std::vector<component> components,
                                         field_function fun, void *fun_data, cdouble *slice,
                                         double omega) {
  return (cdouble *)do_get_array_slice(where, components, fun, 0, fun_data, (void *)slice, omega);
}

double *fields::get_array_slice(const volume &where, component c, double *slice, double omega) {
  std::vector<component> components(1);
  components[0] = c;
  return (double *)do_get_array_slice(where, components, 0, default_field_rfunc, 0, (void *)slice,
                                      omega);
}

double *fields::get_array_slice(const volume &where, derived_component c, double *slice,
                                double omega) {
  int nfields;
  component carray[12];
  field_rfunction rfun = derived_component_func(c, gv, nfields, carray);
  std::vector<component> cs(carray, carray + nfields);
  return (double *)do_get_array_slice(where, cs, 0, rfun, &nfields, (void *)slice, omega);
}

cdouble *fields::get_complex_array_slice(const volume &where, component c, cdouble *slice,
                                         double omega) {
  std::vector<component> components(1);
  components[0] = c;
  return (cdouble *)do_get_array_slice(where, components, default_field_func, 0, 0, (void *)slice,
                                       omega);
}

cdouble *fields::get_source_slice(const volume &where, component source_slice_component,
                                  cdouble *slice) {

  size_t dims[3];
  direction dirs[3];
  vec min_max_loc[2];
  bool collapse = false, snap = false;
  int rank = get_array_slice_dimensions(where, dims, dirs, collapse, snap, min_max_loc);
  size_t slice_size = dims[0] * (rank >= 2 ? dims[1] : 1) * (rank == 3 ? dims[2] : 1);

  source_slice_data data;
  data.source_component = source_slice_component;
  data.slice_imin = gv.round_vec(min_max_loc[0]);
  data.slice_imax = gv.round_vec(min_max_loc[1]);
  data.slice = slice ? slice : new cdouble[slice_size];
  if (!data.slice) abort("%s:%i: out of memory (%zu)", __FILE__, __LINE__, slice_size);

  loop_in_chunks(get_source_slice_chunkloop, (void *)&data, where, Centered, true, true);
  return array_to_all(data.slice, slice_size);
}

/**********************************************************************/
/* increment a multi-index <n_1, n_2, ..., n_r> in which n_i runs over*/
/* the range 0 <= n_i < nMax_i. return true if this is the increment  */
/* that brings the multiindex to all zeros.                           */
/**********************************************************************/
bool increment(size_t *n, size_t *nMax, int rank) {
  for (rank--; rank >= 0; rank--)
    if (++n[rank] < nMax[rank])
      return false;
    else
      n[rank] = 0;
  return true;
}

// data_size = 1,2 for real,complex-valued array
double *collapse_array(double *array, int *rank, size_t dims[3], direction dirs[3], volume where,
                       int data_size = 1) {

  /*--------------------------------------------------------------*/
  /*- detect empty dimensions and compute rank and strides for    */
  /*- collapsed array                                             */
  /*--------------------------------------------------------------*/
  int full_rank = *rank;
  if (full_rank == 0) return array;

  int reduced_rank = 0;
  size_t reduced_dims[3], reduced_stride[3] = {1, 1, 1};
  direction reduced_dirs[3];
  for (int r = 0; r < full_rank; r++) {
    if (where.in_direction(dirs[r]) == 0.0)
      reduced_stride[r] = 0; // degenerate dimension, to be collapsed
    else {
      reduced_dirs[reduced_rank] = dirs[r];
      reduced_dims[reduced_rank++] = dims[r];
    }
  }
  if (reduced_rank == full_rank) return array; // nothing to collapse

  /*--------------------------------------------------------------*/
  /*- set up strides into full and reduced arrays                 */
  /*--------------------------------------------------------------*/
  size_t stride[3] = {1, 1, 1}; // non-reduced array strides
  if (full_rank == 2)
    stride[0] = dims[1]; // rstride is already all set in this case
  else if (full_rank == 3) {
    stride[0] = dims[1] * dims[2];
    stride[1] = dims[2];
    if (reduced_stride[0] != 0)
      reduced_stride[0] = reduced_dims[1];
    else if (reduced_stride[1] != 0)
      reduced_stride[1] = reduced_dims[1];
    // else: two degenerate dimensions->reduced array is 1-diml, no strides needed
  }

  /*--------------------------------------------------------------*/
  /*- allocate reduced array and compress full array into it     -*/
  /*--------------------------------------------------------------*/
  size_t reduced_grid_size = reduced_dims[0] * (reduced_rank == 2 ? reduced_dims[1] : 1);
  size_t reduced_array_size = data_size * reduced_grid_size;
  double *reduced_array = new double[reduced_array_size];
  if (!reduced_array) abort("%s:%i: out of memory (%zu)", __FILE__, __LINE__, reduced_array_size);
  memset(reduced_array, 0, reduced_array_size * sizeof(double));

  size_t n[3] = {0, 0, 0};
  do {
    size_t index = n[0] * stride[0] + n[1] * stride[1] + n[2] * stride[2];
    size_t rindex = n[0] * reduced_stride[0] + n[1] * reduced_stride[1] + n[2] * reduced_stride[2];
    for (int i = 0; i < data_size; i++)
      reduced_array[data_size * rindex + i] += array[data_size * index + i];
  } while (!increment(n, dims, full_rank));

  *rank = reduced_rank;
  for (int r = 0; r < reduced_rank; r++) {
    dims[r] = reduced_dims[r];
    dirs[r] = reduced_dirs[r];
  }
  delete[] array;
  return reduced_array;
}

cdouble *collapse_array(cdouble *array, int *rank, size_t dims[3], direction dirs[3],
                        volume where) {
  return (cdouble *)collapse_array((double *)array, rank, dims, dirs, where, 2);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
std::vector<double> fields::get_array_metadata(const volume &where, bool collapse_empty_dimensions,
                                               bool snap_empty_dimensions) {

  /* get extremal corners of subgrid and array of weights, collapsed if necessary */
  size_t dims[3];
  direction dirs[3];
  vec min_max_loc[2]; // extremal points in subgrid
  int rank = get_array_slice_dimensions(where, dims, dirs, false /*collapse_empty_dimensions*/,
                                        snap_empty_dimensions, min_max_loc);
  int full_rank = rank;
  direction full_dirs[3];
  for (int fr = 0; fr < rank; fr++)
    full_dirs[fr] = dirs[fr];

  double *weights = get_array_slice(where, NO_COMPONENT);
  if (collapse_empty_dimensions) weights = collapse_array(weights, &rank, dims, dirs, where);

  /* get length and endpoints of x,y,z tics arrays */
  size_t nxyz[3] = {1, 1, 1};
  double xyzmin[3] = {0.0, 0.0, 0.0}, xyzmax[3] = {0.0, 0.0, 0.0};
  for (int fr = 0, rr = 0; fr < full_rank; fr++) {
    direction d = full_dirs[fr];
    int nd = d - X;
    if (where.in_direction(d) == 0.0 && collapse_empty_dimensions) {
      xyzmin[nd] = xyzmax[nd] = where.in_direction_min(d);
      nxyz[nd] = 1;
    }
    else {
      nxyz[nd] = dims[rr++];
      xyzmin[nd] = min_max_loc[0].in_direction(d);
      xyzmax[nd] = min_max_loc[1].in_direction(d);
    }
  }

  /* pack all data into a single vector with each tics array preceded by its */
  /* length: [ NX, xtics[:], NY, ytics[:], NZ, ztics[:], weights[:] ]        */
  std::vector<double> xyzw;
  for (int nd = 0; nd < 3; nd++) {
    xyzw.push_back((double)nxyz[nd]);
    for (size_t n = 0; n < nxyz[nd]; n++)
      xyzw.push_back(xyzmin[nd] + n * gv.inva);
  }
  for (unsigned nw = 0; nw < (nxyz[0] * nxyz[1] * nxyz[2]); nw++)
    xyzw.push_back(weights[nw]);

  return xyzw;

} // get_array_metadata

} // namespace meep
