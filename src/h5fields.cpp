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
#include "h5io.h"

namespace meep {

static double get_reim(complex<double> x, int reim)
{
  return reim ? imag(x) : real(x);
}

void fields::output_hdf5(const char *filename,
			 const geometric_volume &where, double res,
			 component c, int reim,
			 bool append_data, int dindex,
                         bool single_precision, bool append_file) {
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
      geometric_volume gvs = S.transform(chunks[i]->get_field_gv(cs), sn);
      if (chunks[i]->is_mine() && chunks[i]->f[cs][reim] && (gvs && vout)) {
	geometric_volume cgv = gvs & vout;
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
    }
  }
  parallel_chunks = max_to_all(parallel_chunks);
  if (parallel_chunks == 0) return; // no data to write

  double *data = new double[nalloc];

  double resinv = 1.0 / res;
  int start[5], count[5] = {1,1,1,1,1};
  int chunks_written = 0;

  char dataname[256];
  snprintf(dataname, 256, "%s.%s", component_name(c), reim ? "i" : "r");

  // Finally, fetch the data from each chunk and write it to the file:
  for (int sn = 0; sn < S.multiplicity(); ++sn) {
    component cs = S.transform(c, -sn);
    complex<double> ph = S.phase_shift(cs, sn);
    for (int i = 0; i < num_chunks; ++i) {
      geometric_volume gvs = S.transform(chunks[i]->get_field_gv(cs), sn);
      if (chunks[i]->is_mine() && chunks[i]->f[cs][reim] && (gvs && vout)) {
	geometric_volume cgv = gvs & vout;
	
	int j;
	for (j = 0; j < rank; ++j) {
	  start[j] = int(ceil(cgv.in_direction_min(ds[j]) * res));
	  count[j] = int(floor(cgv.in_direction_max(ds[j]) * res))- start[j]+1;
	  loc0.set_direction(ds[j], start[j] * resinv);
	  start[j] -= start0[j];
	  if (count[j] <= 0)
	    break;
	}
	if (j < rank)
	  continue;

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
	
	h5io::write_chunk(filename, dataname,
			  rank, dims, data, start, count,
			  append_data, dindex,
			  true, !chunks_written && (!append_data || !dindex),
			  single_precision, append_file);
	++chunks_written;
      }
    }
  }

  delete[] data;

  /* All processes need to call write_chunk in parallel, even if
     some processes have nothing to write. */
  for (int j = 0; j < rank; ++j) start[j] = count[j] = 0; count[0] = 0;
  while (chunks_written < parallel_chunks) {
      h5io::write_chunk(filename, dataname,
			rank, dims, data, start, count,
			append_data, dindex,
			true, !chunks_written && (!append_data || !dindex),
			single_precision, append_file);
      ++chunks_written;
  }
}

void fields::output_hdf5(const char *filename,
			 const geometric_volume &where, double res,
			 component c,
			 bool append_data, int dindex,
                         bool single_precision, bool append_file) {
  output_hdf5(filename, where, res, c, 0,  // real part
	      append_data, dindex, single_precision, append_file);
  if (!is_real)
    output_hdf5(filename, where, res, c, 1,  // imaginary part
		append_data, dindex, single_precision, true);
}

} // meep
