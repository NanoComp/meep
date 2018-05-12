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

// Dump/load raw structure data to/from an HDF5 file.  Only
// works if the number of processors/chunks is the same.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "meep.hpp"
#include "meep_internals.hpp"

using namespace std;

namespace meep {

void structure::dump(const char *filename) {

  // make/save a num_chunks x NUM_FIELD_COMPONENTS x 5 array counting
  // the number of entries in the chi1inv array for each chunk.
  size_t *num_chi1inv_ = new size_t[num_chunks*NUM_FIELD_COMPONENTS*5];
  memset(num_chi1inv_, 0, sizeof(size_t)*size_t(num_chunks*NUM_FIELD_COMPONENTS*5));
  size_t my_ntot = 0;
  for (int i=0; i<num_chunks; i++)
    if (chunks[i]->is_mine()) {
      size_t ntot = chunks[i]->gv.ntot();
      for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c)
        for (int d = 0; d < 5; ++d)
          if (chunks[i]->chi1inv[c][d])
            my_ntot += (num_chi1inv_[(i*NUM_FIELD_COMPONENTS + c)*5 + d] = ntot);
    }
  size_t *num_chi1inv = new size_t[num_chunks*NUM_FIELD_COMPONENTS*5];
  sum_to_master(num_chi1inv_, num_chi1inv, num_chunks*NUM_FIELD_COMPONENTS*5);
  delete[] num_chi1inv_;

  // determine total dataset size and offset of this process's data
  size_t my_start = partial_sum_to_all(my_ntot) - my_ntot;
  size_t ntotal = sum_to_all(my_ntot);

  h5file file(filename, h5file::WRITE, true);
  size_t dims[3] = {(size_t)num_chunks, NUM_FIELD_COMPONENTS, 5};
  size_t start[3] = {0, 0, 0};
  file.create_data("num_chi1inv", 3, dims);
  if (am_master())
    file.write_chunk(3, start, dims, num_chi1inv);
  delete[] num_chi1inv;

  // write the data
  file.create_data("chi1inv", 1, &ntotal);
  for (int i=0; i<num_chunks; i++)
    if (chunks[i]->is_mine()) {
        size_t ntot = chunks[i]->gv.ntot();
        for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c)
          for (int d = 0; d < 5; ++d)
            if (chunks[i]->chi1inv[c][d]) {
              file.write_chunk(1, &my_start, &ntot, chunks[i]->chi1inv[c][d]);
              my_start += ntot;
            }
    }
}

void structure::load(const char *filename) {
  h5file file(filename, h5file::READONLY, true);

  // make/save a num_chunks x NUM_FIELD_COMPONENTS x 5 array counting
  // the number of entries in the chi1inv array for each chunk.
  size_t *num_chi1inv = new size_t[num_chunks*NUM_FIELD_COMPONENTS*5];
  int rank;
  size_t dims[3], _dims[3] = {(size_t)num_chunks, NUM_FIELD_COMPONENTS, 5};
  size_t start[3] = {0, 0, 0};
  file.read_size("num_chi1inv", &rank, dims, 3);
  if (rank != 3 || _dims[0] != dims[0] || _dims[1] != dims[1] || _dims[2] != dims[2])
    abort("chunk mismatch in structure::load");
  if (am_master())
    file.read_chunk(3, start, dims, num_chi1inv);

  file.prevent_deadlock();
  broadcast(0, num_chi1inv, dims[0]*dims[1]*dims[2]);

  changing_chunks();

  // allocate data as needed and check sizes
  size_t my_ntot = 0;
  for (int i=0; i<num_chunks; i++)
    if (chunks[i]->is_mine()) {
      size_t ntot = chunks[i]->gv.ntot();
      for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c)
        for (int d = 0; d < 5; ++d) {
          size_t n = num_chi1inv[(i*NUM_FIELD_COMPONENTS + c)*5 + d];
          if (n == 0) {
            delete[] chunks[i]->chi1inv[c][d];
            chunks[i]->chi1inv[c][d] = NULL;
          }
          else {
            if (n != ntot)
              abort("grid size mismatch %zd vs %zd in structure::load", n, ntot);
            chunks[i]->chi1inv[c][d] = new realnum[ntot];
            my_ntot += ntot;
          }
        }
    }

  // determine total dataset size and offset of this process's data
  size_t my_start = partial_sum_to_all(my_ntot) - my_ntot;
  size_t ntotal = sum_to_all(my_ntot);

  // read the data
  file.read_size("chi1inv", &rank, dims, 1);
  if (rank != 1 || dims[0] != ntotal)
    abort("inconsistent data size in structure::load");
  for (int i=0; i<num_chunks; i++)
    if (chunks[i]->is_mine()) {
        size_t ntot = chunks[i]->gv.ntot();
        for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c)
          for (int d = 0; d < 5; ++d)
            if (chunks[i]->chi1inv[c][d]) {
              file.read_chunk(1, &my_start, &ntot, chunks[i]->chi1inv[c][d]);
              my_start += ntot;
            }
    }
}

}
