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

static double get_reim(complex<double> x, int reim)
{
  return reim ? imag(x) : real(x);
}

void fields::output_hdf5(h5file *file, const char *dataname,
			 component c, int reim,
			 const geometric_volume &where, double res,
			 bool append_data,
                         bool single_precision) {
  geometric_volume vout(where); // FIXME: intersect with computational cell?
  vec loc0(vout.dim);

  // First, get total dimensions of output HDF5 file:
  int rank = 0, dims[5], start0[5];
  direction ds[5];
  LOOP_OVER_DIRECTIONS(vout.dim, d) {
    int minpt = int(ceil(vout.in_direction_min(d) * res));
    int maxpt = int(floor(vout.in_direction_max(d) * res));
    if (minpt < maxpt) {
      ds[rank] = d;
      start0[rank] = minpt;
      dims[rank++] = maxpt - minpt + 1;
    }
    else
      loc0.set_direction(d, 0.5 * (vout.in_direction_min(d) + 
				   vout.in_direction_max(d)));
  }

  // Next, determine number of parallel chunks to write, and required storage
  int parallel_chunks = 0;
  int nalloc = 1;
  for (int sn = 0; sn < S.multiplicity(); ++sn) {
    component cs = S.transform(c, -sn);
    for (int i = 0; i < num_chunks; ++i) {
      geometric_volume gvS = S.transform(chunks[i]->get_field_gv(cs), sn);
      if (!chunks[i]->is_mine() || !chunks[i]->have_component(cs, reim == 1))
	continue;

      // figure out range of lattice shifts for which gvS intersects vout:
      ivec min_ishift(gvS.dim), max_ishift(gvS.dim);
      LOOP_OVER_DIRECTIONS(gvS.dim, d) {
	if (boundaries[High][d] == Periodic) {
	  min_ishift.set_direction(d, int(floor((vout.in_direction_min(d) 
						 - gvS.in_direction_max(d))
                                       / lattice_vector(d).in_direction(d))));
	  max_ishift.set_direction(d, int(ceil((vout.in_direction_max(d) 
						- gvS.in_direction_min(d))
                                       / lattice_vector(d).in_direction(d))));
	}
	else {
	  min_ishift.set_direction(d, 0);
	  max_ishift.set_direction(d, 0);
	}
      }

      ivec ishift(min_ishift);
      do {
	vec shift(gvS.dim);
	LOOP_OVER_DIRECTIONS(gvS.dim, d)
	  shift.set_direction(d,
                lattice_vector(d).in_direction(d) * ishift.in_direction(d));
	geometric_volume gvSs(gvS + shift);
	if (gvSs && vout) {
	  geometric_volume cgv = gvSs & vout;
	  int nvol = 1;
	  for (int j = 0; j < rank; ++j) {
	    int minpt = int(ceil(cgv.in_direction_min(ds[j]) * res));
	    int maxpt = int(floor(cgv.in_direction_max(ds[j]) * res));
	    if (minpt > maxpt) nvol = 0;
	    else nvol *= maxpt - minpt + 1;
	  }
	  if (nvol > nalloc) nalloc = nvol;
	  if (nvol > 0)
	    ++parallel_chunks;
	}
	LOOP_OVER_DIRECTIONS(gvS.dim, d) {
	  if (ishift.in_direction(d) + 1 <= max_ishift.in_direction(d)) {
	    ishift.set_direction(d, ishift.in_direction(d) + 1);
	    break;
	  }
	  ishift.set_direction(d, min_ishift.in_direction(d));
	}
      } while (ishift != min_ishift);
    }
  }
  parallel_chunks = max_to_all(parallel_chunks);
  if (parallel_chunks == 0) return; // no data to write

  file->create_or_extend_data(dataname, rank, dims,
			      append_data, single_precision);

  double *data = new double[nalloc];

  double resinv = 1.0 / res;
  int start[5], count[5] = {1,1,1,1,1};

  // Finally, fetch the data from each chunk and write it to the file:
  for (int sn = 0; sn < S.multiplicity(); ++sn) {
    component cs = S.transform(c, -sn);
    complex<double> phS = S.phase_shift(cs, sn);
    for (int i = 0; i < num_chunks; ++i) {
      geometric_volume gvS = S.transform(chunks[i]->get_field_gv(cs), sn);
      if (!chunks[i]->is_mine() || !chunks[i]->have_component(cs, reim == 1))
	continue;

      // figure out range of lattice shifts for which gvS intersects vout:
      ivec min_ishift(gvS.dim), max_ishift(gvS.dim);
      LOOP_OVER_DIRECTIONS(gvS.dim, d) {
	if (boundaries[High][d] == Periodic) {
	  min_ishift.set_direction(d, int(floor((vout.in_direction_min(d) 
						 - gvS.in_direction_max(d))
                                       / lattice_vector(d).in_direction(d))));
	  max_ishift.set_direction(d, int(ceil((vout.in_direction_max(d) 
						- gvS.in_direction_min(d))
                                       / lattice_vector(d).in_direction(d))));
	}
	else {
	  min_ishift.set_direction(d, 0);
	  max_ishift.set_direction(d, 0);
	}
      }

      ivec ishift(min_ishift);
      do {
	complex<double> ph = phS;
	vec shift(gvS.dim);
	LOOP_OVER_DIRECTIONS(gvS.dim, d) {
	  shift.set_direction(d,
                lattice_vector(d).in_direction(d) * ishift.in_direction(d));
	  ph *= pow(eikna[d], ishift.in_direction(d));
	}
	if (c == Dielectric) ph = 1.0;
	geometric_volume gvSs(gvS + shift);
	if (gvSs && vout) {
	  geometric_volume cgv = gvSs & vout;
	
	  int j;
	  for (j = 0; j < rank; ++j) {
	    start[j] = int(ceil(cgv.in_direction_min(ds[j]) * res));
	    count[j] = 
	      int(floor(cgv.in_direction_max(ds[j]) * res)) - start[j]+1;
	    loc0.set_direction(ds[j], start[j] * resinv);
	    start[j] -= start0[j];
	    if (count[j] <= 0)
	      break;
	  }
	  if (j < rank)
	    goto next_shift;

	  loc0 -= shift;
	  
	  switch (rank) {
	  case 0:
	    data[0] =
	      get_reim(chunks[i]->get_field(cs, S.transform(loc0,-sn))*ph,
		       reim);
	    break;
	  case 1: {
	    vec loc = loc0;
	    for (int i0 = 0; i0 < count[0]; ++i0) {
	      loc.set_direction(ds[0], loc0.in_direction(ds[0]) + i0 * resinv);
	      data[i0] =
		get_reim(chunks[i]->get_field(cs, S.transform(loc,-sn))*ph,
			 reim);
	    }
	    break;
	  }
	  case 2: {
	    vec loc = loc0;
	    for (int i0 = 0; i0 < count[0]; ++i0) {
	      loc.set_direction(ds[0], loc0.in_direction(ds[0]) + i0 * resinv);
	      for (int i1 = 0; i1 < count[1]; ++i1) {
		loc.set_direction(ds[1], loc0.in_direction(ds[1])
				  + i1 * resinv);
		data[i0 * count[1] + i1] =
		  get_reim(chunks[i]->get_field(cs, S.transform(loc,-sn))*ph,
			   reim);
	      }
	    }
	    break;
	  }
	  case 3: {
	    vec loc = loc0;
	    for (int i0 = 0; i0 < count[0]; ++i0) {
	      loc.set_direction(ds[0], loc0.in_direction(ds[0]) + i0 * resinv);
	      for (int i1 = 0; i1 < count[1]; ++i1) {
		loc.set_direction(ds[1], loc0.in_direction(ds[1])
				  + i1 * resinv);
		for (int i2 = 0; i2 < count[2]; ++i2) {
		  loc.set_direction(ds[2], loc0.in_direction(ds[2])
				    + i2 * resinv);
		  data[(i0 * count[1] + i1) * count[2] + i2] =
		    get_reim(chunks[i]->get_field(cs, S.transform(loc,-sn))*ph,
			     reim);
		}
	      }
	    }
	    break;
	  }
	  default:
	    abort("unexpected dimensionality > 3 of HDF5 output data");
	  }
	  
	  file->write_chunk(rank, start, count, data);
	}
      next_shift:
	LOOP_OVER_DIRECTIONS(gvS.dim, d) {
	  if (ishift.in_direction(d) + 1 <= max_ishift.in_direction(d)) {
	    ishift.set_direction(d, ishift.in_direction(d) + 1);
	    break;
	  }
	  ishift.set_direction(d, min_ishift.in_direction(d));
	}
      } while (ishift != min_ishift);
    }
  }

  delete[] data;
}

void fields::output_hdf5(h5file *file, component c,
			 const geometric_volume &where, double res,
			 bool append_data,
                         bool single_precision) {
  char dataname[256];
  bool has_imag = !is_real && c != Dielectric;

  snprintf(dataname, 256, "%s%s", component_name(c), has_imag ? ".r" : "");
  output_hdf5(file, dataname, c, 0, where,res, append_data, single_precision);
  if (has_imag) {
    snprintf(dataname, 256, "%s.i", component_name(c));
    output_hdf5(file, dataname, c, 1, where,res, append_data,single_precision);
  }
}

void fields::output_hdf5(component c,
			 const geometric_volume &where, double res,
			 bool append_data, 
			 bool single_precision, 
			 const char *prefix) {
  h5file *file = open_h5file(component_name(c), h5file::WRITE,
			     prefix, !append_data);
  output_hdf5(file, c, where, res, append_data, single_precision);
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
