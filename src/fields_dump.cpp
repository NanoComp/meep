/* Copyright (C) 2005-2023 Massachusetts Institute of Technology
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

void fields::dump_fields_chunk_field(h5file *h5f, bool single_parallel_file,
                                     const std::string &field_name,
                                     FieldPtrGetter field_ptr_getter) {
  /*
   * make/save a num_chunks x NUM_FIELD_COMPONENTS x 2 array counting
   * the number of entries in the 'field_name' array for each chunk.
   *
   * When 'single_parallel_file' is true, we are creating a single block of data
   * for ALL chunks (that are merged using MPI). Otherwise, we are just
   * making a copy of just the chunks that are ours.
   */
  int my_num_chunks = 0;
  for (int i = 0; i < num_chunks; i++) {
    my_num_chunks += (single_parallel_file || chunks[i]->is_mine());
  }
  size_t num_f_size = my_num_chunks * NUM_FIELD_COMPONENTS * 2;
  std::vector<size_t> num_f_(num_f_size);
  size_t my_ntot = 0;
  for (int i = 0, chunk_i = 0; i < num_chunks; i++) {
    if (chunks[i]->is_mine()) {
      FOR_FIELD_TYPES(ft) {
        if (chunks[i]->pol[ft] != NULL)
          meep::abort("non-null polarization_state in fields::dump (unsupported)");
      }
      size_t ntot = chunks[i]->gv.ntot();
      for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c) {
        for (int d = 0; d < 2; ++d) {
          realnum **f = field_ptr_getter(chunks[i], c, d);
          if (*f) { my_ntot += (num_f_[(chunk_i * NUM_FIELD_COMPONENTS + c) * 2 + d] = ntot); }
        }
      }
    }
    chunk_i += (chunks[i]->is_mine() || single_parallel_file);
  }

  std::vector<size_t> num_f;
  if (single_parallel_file) {
    num_f.resize(num_f_size);
    sum_to_master(num_f_.data(), num_f.data(), num_f_size);
  }
  else { num_f = std::move(num_f_); }

  /* determine total dataset size and offset of this process's data */
  size_t my_start = 0;
  size_t ntotal = my_ntot;
  if (single_parallel_file) {
    my_start = partial_sum_to_all(my_ntot) - my_ntot;
    ntotal = sum_to_all(my_ntot);
  }

  size_t dims[3] = {(size_t)my_num_chunks, NUM_FIELD_COMPONENTS, 2};
  size_t start[3] = {0, 0, 0};
  std::string num_f_name = std::string("num_") + field_name;
  h5f->create_data(num_f_name.c_str(), 3, dims);
  if (am_master() || !single_parallel_file) { h5f->write_chunk(3, start, dims, num_f.data()); }

  /* write the data */
  h5f->create_data(field_name.c_str(), 1, &ntotal, false /* append_data */,
                   false /* single_precision */);
  for (int i = 0; i < num_chunks; i++) {
    if (chunks[i]->is_mine()) {
      size_t ntot = chunks[i]->gv.ntot();
      for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c) {
        for (int d = 0; d < 2; ++d) {
          realnum **f = field_ptr_getter(chunks[i], c, d);
          if (*f) {
            h5f->write_chunk(1, &my_start, &ntot, *f);
            my_start += ntot;
          }
        }
      }
    }
  }
}

void fields::dump(const char *filename, bool single_parallel_file) {
  if (verbosity > 0) {
    printf("creating fields output file \"%s\" (%d)...\n", filename, single_parallel_file);
  }

  h5file file(filename, h5file::WRITE, single_parallel_file, !single_parallel_file);

  // Write out the current time 't'
  size_t dims[1] = {1};
  size_t start[1] = {0};
  size_t _t[1] = {(size_t)t};
  file.create_data("t", 1, dims);
  if (am_master() || !single_parallel_file) file.write_chunk(1, start, dims, _t);

  dump_fields_chunk_field(&file, single_parallel_file, "f",
                          [](fields_chunk *chunk, int c, int d) { return &(chunk->f[c][d]); });
  dump_fields_chunk_field(&file, single_parallel_file, "f_u",
                          [](fields_chunk *chunk, int c, int d) { return &(chunk->f_u[c][d]); });
  dump_fields_chunk_field(&file, single_parallel_file, "f_w",
                          [](fields_chunk *chunk, int c, int d) { return &(chunk->f_w[c][d]); });
  dump_fields_chunk_field(&file, single_parallel_file, "f_cond",
                          [](fields_chunk *chunk, int c, int d) { return &(chunk->f_cond[c][d]); });
  dump_fields_chunk_field(
      &file, single_parallel_file, "f_w_prev",
      [](fields_chunk *chunk, int c, int d) { return &(chunk->f_w_prev[c][d]); });

  // Dump DFT chunks.
  for (int i = 0; i < num_chunks; i++) {
    if (single_parallel_file || chunks[i]->is_mine()) {
      char dataname[1024];
      snprintf(dataname, 1024, "chunk%02d", i);
      save_dft_hdf5(chunks[i]->dft_chunks, dataname, &file, 0, single_parallel_file);
    }
  }
}

void fields::load_fields_chunk_field(h5file *h5f, bool single_parallel_file,
                                     const std::string &field_name,
                                     FieldPtrGetter field_ptr_getter) {
  int my_num_chunks = 0;
  for (int i = 0; i < num_chunks; i++) {
    my_num_chunks += (single_parallel_file || chunks[i]->is_mine());
  }
  size_t num_f_size = my_num_chunks * NUM_FIELD_COMPONENTS * 2;
  std::vector<size_t> num_f(num_f_size);

  int rank;
  size_t dims[3], _dims[3] = {(size_t)my_num_chunks, NUM_FIELD_COMPONENTS, 2};
  size_t start[3] = {0, 0, 0};

  std::string num_f_name = std::string("num_") + field_name;
  h5f->read_size(num_f_name.c_str(), &rank, dims, 3);
  if (rank != 3 || _dims[0] != dims[0] || _dims[1] != dims[1] || _dims[2] != dims[2])
    meep::abort("chunk mismatch in fields::load");
  if (am_master() || !single_parallel_file) h5f->read_chunk(3, start, dims, num_f.data());

  if (single_parallel_file) {
    h5f->prevent_deadlock();
    broadcast(0, num_f.data(), dims[0] * dims[1] * dims[2]);
  }

  /* allocate data as needed and check sizes */
  size_t my_ntot = 0;
  for (int i = 0, chunk_i = 0; i < num_chunks; i++) {
    if (chunks[i]->is_mine()) {
      size_t ntot = chunks[i]->gv.ntot();
      for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c) {
        for (int d = 0; d < 2; ++d) {
          size_t n = num_f[(chunk_i * NUM_FIELD_COMPONENTS + c) * 2 + d];
          realnum **f = field_ptr_getter(chunks[i], c, d);
          if (n == 0) {
            delete[] * f;
            *f = NULL;
          }
          else {
            if (n != ntot) meep::abort("grid size mismatch %zd vs %zd in fields::load", n, ntot);
            // here we need to allocate the fields array for H in the PML region
            // because of H = B in fields_chunk::alloc_f whereby H is lazily
            // allocated in fields_chunk::update_eh during the first timestep
            const direction d_c = component_direction(c);
            if (!(*f) || (*f && is_magnetic(component(c)) && chunks[i]->s->sigsize[d_c] > 1))
              *f = new realnum[ntot];
            my_ntot += ntot;
          }
        }
      }
    }
    chunk_i += (chunks[i]->is_mine() || single_parallel_file);
  }

  /* determine total dataset size and offset of this process's data */
  size_t my_start = 0;
  size_t ntotal = my_ntot;
  if (single_parallel_file) {
    my_start = partial_sum_to_all(my_ntot) - my_ntot;
    ntotal = sum_to_all(my_ntot);
  }

  /* read the data */
  h5f->read_size(field_name.c_str(), &rank, dims, 1);
  if (rank != 1 || dims[0] != ntotal) {
    meep::abort("inconsistent data size for '%s' in fields::load (rank, dims[0]): "
                "(%d, %zu) != (1, %zu)",
                field_name.c_str(), rank, dims[0], ntotal);
  }
  for (int i = 0; i < num_chunks; i++) {
    if (chunks[i]->is_mine()) {
      size_t ntot = chunks[i]->gv.ntot();
      for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c) {
        for (int d = 0; d < 2; ++d) {
          realnum **f = field_ptr_getter(chunks[i], c, d);
          if (*f) {
            h5f->read_chunk(1, &my_start, &ntot, *f);
            my_start += ntot;
          }
        }
      }
    }
  }
}

void fields::load(const char *filename, bool single_parallel_file) {
  if (verbosity > 0)
    printf("reading fields from file \"%s\" (%d)...\n", filename, single_parallel_file);

  h5file file(filename, h5file::READONLY, single_parallel_file, !single_parallel_file);

  // Read in the current time 't'
  int rank;
  size_t dims[1] = {1};
  size_t start[1] = {0};
  size_t _t[1];
  file.read_size("t", &rank, dims, 1);
  if (rank != 1 || dims[0] != 1) meep::abort("time size mismatch in fields::load");
  if (am_master() || !single_parallel_file) file.read_chunk(1, start, dims, _t);

  if (single_parallel_file) {
    file.prevent_deadlock();
    broadcast(0, _t, dims[0]);
  }

  t = static_cast<int>(_t[0]);
  calc_sources(time());

  load_fields_chunk_field(&file, single_parallel_file, "f",
                          [](fields_chunk *chunk, int c, int d) { return &(chunk->f[c][d]); });
  load_fields_chunk_field(&file, single_parallel_file, "f_u",
                          [](fields_chunk *chunk, int c, int d) { return &(chunk->f_u[c][d]); });
  load_fields_chunk_field(&file, single_parallel_file, "f_w",
                          [](fields_chunk *chunk, int c, int d) { return &(chunk->f_w[c][d]); });
  load_fields_chunk_field(&file, single_parallel_file, "f_cond",
                          [](fields_chunk *chunk, int c, int d) { return &(chunk->f_cond[c][d]); });
  load_fields_chunk_field(
      &file, single_parallel_file, "f_w_prev",
      [](fields_chunk *chunk, int c, int d) { return &(chunk->f_w_prev[c][d]); });

  // Load DFT chunks.
  for (int i = 0; i < num_chunks; i++) {
    if (single_parallel_file || chunks[i]->is_mine()) {
      char dataname[1024];
      snprintf(dataname, 1024, "chunk%02d", i);
      load_dft_hdf5(chunks[i]->dft_chunks, dataname, &file, 0, single_parallel_file);
    }
  }
}

} // namespace meep
