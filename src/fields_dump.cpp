/* Copyright (C) 2005-2021 Massachusetts Institute of Technology
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

// Dump/load raw fields data to/from an HDF5 file.  Only
// works if the number of processors/chunks is the same.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cassert>

#include "meep.hpp"
#include "meep_internals.hpp"

namespace meep {

#define DUMP_FIELD(field_name)                                              \
  {                                                                         \
    /*                                                                      \
     * make/save a num_chunks x NUM_FIELD_COMPONENTS x 2 array counting     \
     * the number of entries in the 'field_name' array for each chunk.      \
     */                                                                     \
    size_t *num_f_ = new size_t[num_chunks * NUM_FIELD_COMPONENTS * 2];     \
    memset(num_f_, 0,                                                       \
           sizeof(size_t) * size_t(num_chunks * NUM_FIELD_COMPONENTS * 2)); \
    size_t my_ntot = 0;                                                     \
    for (int i = 0; i < num_chunks; i++) {                                  \
      if (chunks[i]->is_mine()) {                                           \
        size_t ntot = chunks[i]->gv.ntot();                                 \
        for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c) {                    \
          for (int d = 0; d < 2; ++d) {                                     \
            if (chunks[i]->field_name[c][d]) {                              \
              my_ntot +=                                                    \
                  (num_f_[(i * NUM_FIELD_COMPONENTS + c) * 2 + d] = ntot);  \
            }                                                               \
          }                                                                 \
        }                                                                   \
      }                                                                     \
    }                                                                       \
    size_t *num_f = new size_t[num_chunks * NUM_FIELD_COMPONENTS * 2];      \
    sum_to_master(num_f_, num_f, num_chunks *NUM_FIELD_COMPONENTS * 2);     \
    delete[] num_f_;                                                        \
                                                                            \
    /* determine total dataset size and offset of this process's data */    \
    size_t my_start = partial_sum_to_all(my_ntot) - my_ntot;                \
    size_t ntotal = sum_to_all(my_ntot);                                    \
                                                                            \
    size_t dims[3] = {(size_t)num_chunks, NUM_FIELD_COMPONENTS, 2};         \
    size_t start[3] = {0, 0, 0};                                            \
    file.create_data("num_" #field_name, 3, dims);                          \
    if (am_master()) file.write_chunk(3, start, dims, num_f);               \
    delete[] num_f;                                                         \
                                                                            \
    /* write the data */                                                    \
    file.create_data(#field_name, 1, &ntotal, false /* append_data */,      \
                     false /* single_precision */);                         \
    for (int i = 0; i < num_chunks; i++) {                                  \
      if (chunks[i]->is_mine()) {                                           \
        size_t ntot = chunks[i]->gv.ntot();                                 \
        for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c) {                    \
          for (int d = 0; d < 2; ++d)                                       \
            if (chunks[i]->field_name[c][d]) {                              \
              file.write_chunk(1, &my_start, &ntot,                         \
                               chunks[i]->field_name[c][d]);                \
              my_start += ntot;                                             \
            }                                                               \
        }                                                                   \
      }                                                                     \
    }                                                                       \
  }

void fields::dump(const char *filename) {
  if (verbosity > 0) {
    master_printf("creating fields output file \"%s\"...\n", filename);
  }

  h5file file(filename, h5file::WRITE, true);
  DUMP_FIELD(f);
  DUMP_FIELD(f_u);
  DUMP_FIELD(f_w);
  DUMP_FIELD(f_cond);
}
#undef DUMP_FIELD

#undef LOAD_FIELD
#define LOAD_FIELD(field_name)                                                \
  {                                                                           \
    /*                                                                        \
     * make/save a num_chunks x NUM_FIELD_COMPONENTS x 2 array counting       \
     * the number of entries in the 'field_name' array for each chunk.        \
     */                                                                       \
    size_t *num_f = new size_t[num_chunks * NUM_FIELD_COMPONENTS * 2];        \
    int rank;                                                                 \
    size_t dims[3], _dims[3] = {(size_t)num_chunks, NUM_FIELD_COMPONENTS, 2}; \
    size_t start[3] = {0, 0, 0};                                              \
    file.read_size("num_" #field_name, &rank, dims, 3);                       \
    if (rank != 3 || _dims[0] != dims[0] || _dims[1] != dims[1] ||            \
        _dims[2] != dims[2]) {                                                \
      meep::abort("chunk mismatch in fields::load");                          \
    }                                                                         \
    if (am_master()) file.read_chunk(3, start, dims, num_f);                  \
                                                                              \
    file.prevent_deadlock();                                                  \
    broadcast(0, num_f, dims[0] * dims[1] * dims[2]);                         \
                                                                              \
    /* allocate data as needed and check sizes */                             \
    size_t my_ntot = 0;                                                       \
    for (int i = 0; i < num_chunks; i++) {                                    \
      if (chunks[i]->is_mine()) {                                             \
        size_t ntot = chunks[i]->gv.ntot();                                   \
        for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c) {                      \
          for (int d = 0; d < 2; ++d) {                                       \
            size_t n = num_f[(i * NUM_FIELD_COMPONENTS + c) * 2 + d];         \
            if (n == 0) {                                                     \
              delete[] chunks[i]->field_name[c][d];                           \
              chunks[i]->field_name[c][d] = NULL;                             \
            } else {                                                          \
              if (n != ntot)                                                  \
                meep::abort("grid size mismatch %zd vs %zd in fields::load",  \
                            n, ntot);                                         \
              chunks[i]->field_name[c][d] = new realnum[ntot];                \
              my_ntot += ntot;                                                \
            }                                                                 \
          }                                                                   \
        }                                                                     \
      }                                                                       \
    }                                                                         \
                                                                              \
    /* determine total dataset size and offset of this process's data */      \
    size_t my_start = partial_sum_to_all(my_ntot) - my_ntot;                  \
    size_t ntotal = sum_to_all(my_ntot);                                      \
                                                                              \
    /* read the data */                                                       \
    file.read_size(#field_name, &rank, dims, 1);                              \
    if (rank != 1 || dims[0] != ntotal) {                                     \
      meep::abort("inconsistent data size in fields::load");                  \
    }                                                                         \
    for (int i = 0; i < num_chunks; i++) {                                    \
      if (chunks[i]->is_mine()) {                                             \
        size_t ntot = chunks[i]->gv.ntot();                                   \
        for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c) {                      \
          for (int d = 0; d < 2; ++d) {                                       \
            if (chunks[i]->field_name[c][d]) {                                \
              file.read_chunk(1, &my_start, &ntot,                            \
                              chunks[i]->field_name[c][d]);                   \
              my_start += ntot;                                               \
            }                                                                 \
          }                                                                   \
        }                                                                     \
      }                                                                       \
    }                                                                         \
  }

void fields::load(const char *filename) {
  if (verbosity > 0)
    master_printf("reading fields from file \"%s\"...\n", filename);

  h5file file(filename, h5file::READONLY, true);
  LOAD_FIELD(f);
  LOAD_FIELD(f_u);
  LOAD_FIELD(f_w);
  LOAD_FIELD(f_cond);
}

#undef LOAD_FIELD

}  // namespace meep
