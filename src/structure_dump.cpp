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

  // dump the susceptibilities
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine()) {
      for (int ft = 0; ft < NUM_FIELD_TYPES; ++ft) {
        if (chunks[i]->chiP[ft]) {
          susceptibility *susc = chunks[i]->chiP[ft];
          int n = 0;
          while (susc) {
            char prefix[12];
            snprintf(prefix, 12, "%d_%d_%d_", i, ft, n);
            susc->dump(&file, prefix);
            susc = susc->next;
            n++;
          }
        }
      }
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

  // load the susceptibilities
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine()) {
      for (int ft = 0; ft < NUM_FIELD_TYPES; ++ft) {
        int n = 0;
        char params_dset[25];
        snprintf(params_dset, 25, "%d_%d_%d_", i, ft, 0);
        strcat(params_dset, "params");
        susceptibility *susc_ptr;

        while (file.dataset_exists(params_dset)) {
          int params_rank;
          size_t params_start[1] = {0};
          size_t params_dims[1] = {0};
          file.read_size(params_dset, &params_rank, params_dims, 1);

          if (params_rank != 1 || params_dims[0] > 4 || params_dims[0] < 3)
            abort("inconsistent data size in structure::load");

          if (params_dims[0] == 3) {
            realnum params[3] = {0};

            if (am_master())
              file.read_chunk(1, params_start, params_dims, params);
            file.prevent_deadlock();
            broadcast(0, params, 3);

            if (n == 0) {
              chunks[i]->chiP[ft] = new lorentzian_susceptibility(params[0], params[1], (bool)params[2]);
              susc_ptr = chunks[i]->chiP[ft];
            }
            else
              susc_ptr->next = new lorentzian_susceptibility(params[0], params[1], (bool)params[2]);
          }
          else {
            realnum params[4] = {0};

            if (am_master())
              file.read_chunk(1, params_start, params_dims, params);
            file.prevent_deadlock();
            broadcast(0, params, 4);

            if (n == 0) {
              chunks[i]->chiP[ft] = new noisy_lorentzian_susceptibility(params[0], params[1], params[2],
                                                                   (bool)params[3]);
              susc_ptr = chunks[i]->chiP[ft];
            }
            else
              susc_ptr->next = new noisy_lorentzian_susceptibility(params[0], params[1], params[2],
                                                                   (bool)params[3]);
          }

          size_t ntot = chunks[i]->gv.ntot();
          if (n == 0)
            susc_ptr->ntot = ntot;
          else
            susc_ptr->next->ntot = ntot;
          char prefix[12];
          snprintf(prefix, 12, "%d_%d_%d_", i, ft, n);
          if (n == 0)
            susc_ptr->load(&file, prefix);
          else
            susc_ptr->next->load(&file, prefix);

          if (n != 0)
            susc_ptr = susc_ptr->next;
          n++;
          snprintf(params_dset, 25, "%d_%d_%d_", i, ft, n);
          strcat(params_dset, "params");
        }
      }
    }
  }
}
}
