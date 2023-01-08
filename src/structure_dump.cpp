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

// Dump/load raw structure data to/from an HDF5 file.  Only
// works if the number of processors/chunks is the same.

#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "meep.hpp"
#include "meep_internals.hpp"

using namespace std;

namespace meep {

// Write the parameters required to reconstruct the susceptibility (id, noise_amp (for noisy),
// omega_0, gamma, no_omega_0_denominator)
void structure::write_susceptibility_params(h5file *file, bool single_parallel_file,
                                            const char *dname, int EorH) {
  // Get number of susceptibility params from first chunk, since all chunks will have
  // the same susceptibility list.
  size_t params_ntotal = 0;
  susceptibility *sus = chunks[0]->chiP[EorH];
  while (sus) {
    params_ntotal += sus->get_num_params() + 1;
    sus = sus->next;
  }

  // Write params
  size_t params_start = 0;
  file->create_data(dname, 1, &params_ntotal, false /* append_data */,
                    sizeof(realnum) == sizeof(float) /* single_precision */);
  if (am_master() || !single_parallel_file) {
    susceptibility *sus = chunks[0]->chiP[EorH];
    while (sus) {
      sus->dump_params(file, &params_start);
      sus = sus->next;
    }
  }
}

void structure::dump_chunk_layout(const char *filename) {
  // Write grid_volume info for each chunk so we can reconstruct chunk division from split_by_cost
  size_t sz = num_chunks * 3;
  double *origins = new double[sz];
  size_t *nums = new size_t[sz];
  memset(nums, 0, sizeof(size_t) * sz);
  memset(origins, 0, sizeof(double) * sz);
  for (int i = 0; i < num_chunks; ++i) {
    int idx = i * 3;
    LOOP_OVER_DIRECTIONS(gv.dim, d) {
      origins[idx++] = chunks[i]->gv.origin_in_direction(d);
      nums[i * 3 + ((int)d % 3)] = chunks[i]->gv.num_direction(d);
    }
  }
  h5file file(filename, h5file::WRITE, true);
  file.create_data("gv_origins", 1, &sz, false /* append_data */, false /* single_precision */);
  if (am_master()) {
    size_t gv_origins_start = 0;
    file.write_chunk(1, &gv_origins_start, &sz, origins);
  }
  file.create_data("gv_nums", 1, &sz);
  if (am_master()) {
    size_t nums_start = 0;
    file.write_chunk(1, &nums_start, &sz, nums);
  }
  delete[] origins;
  delete[] nums;
}

void structure::dump(const char *filename, bool single_parallel_file) {
  if (verbosity > 0)
    printf("creating epsilon from file \"%s\" (%d)...\n", filename, single_parallel_file);

  /*
   * make/save a num_chunks x NUM_FIELD_COMPONENTS x 5 array counting
   * the number of entries in the chi1inv array for each chunk.
   *
   * When 'single_parallel_file' is true, we are creating a single block of data
   * for ALL chunks (that are merged using MPI). Otherwise, we are just
   * making a copy of just the chunks that are ours.
   */
  int my_num_chunks = 0;
  for (int i = 0; i < num_chunks; i++) {
    my_num_chunks += (single_parallel_file || chunks[i]->is_mine());
  }
  size_t num_chi1inv_size = my_num_chunks * NUM_FIELD_COMPONENTS * 5;
  std::vector<size_t> num_chi1inv_(num_chi1inv_size);
  size_t my_ntot = 0;
  for (int i = 0, chunk_i = 0; i < num_chunks; i++) {
    if (chunks[i]->is_mine()) {
      size_t ntot = chunks[i]->gv.ntot();
      for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c) {
        for (int d = 0; d < 5; ++d) {
          if (chunks[i]->chi1inv[c][d])
            my_ntot += (num_chi1inv_[(chunk_i * NUM_FIELD_COMPONENTS + c) * 5 + d] = ntot);
        }
      }
    }
    chunk_i += (chunks[i]->is_mine() || single_parallel_file);
  }

  std::vector<size_t> num_chi1inv;
  if (single_parallel_file) {
    num_chi1inv.resize(num_chi1inv_size);
    sum_to_master(num_chi1inv_.data(), num_chi1inv.data(), num_chi1inv_size);
  }
  else { num_chi1inv = std::move(num_chi1inv_); }

  // determine total dataset size and offset of this process's data
  size_t my_start = 0;
  size_t ntotal = my_ntot;
  if (single_parallel_file) {
    my_start = partial_sum_to_all(my_ntot) - my_ntot;
    ntotal = sum_to_all(my_ntot);
  }

  h5file file(filename, h5file::WRITE, single_parallel_file, !single_parallel_file);
  size_t dims[3] = {(size_t)my_num_chunks, NUM_FIELD_COMPONENTS, 5};
  size_t start[3] = {0, 0, 0};
  file.create_data("num_chi1inv", 3, dims);
  if (am_master() || !single_parallel_file) {
    file.write_chunk(3, start, dims, num_chi1inv.data());
  }

  // write the data
  file.create_data("chi1inv", 1, &ntotal, false /* append_data */, false /* single_precision */);
  for (int i = 0; i < num_chunks; i++)
    if (chunks[i]->is_mine()) {
      size_t ntot = chunks[i]->gv.ntot();
      for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c)
        for (int d = 0; d < 5; ++d)
          if (chunks[i]->chi1inv[c][d]) {
            file.write_chunk(1, &my_start, &ntot, chunks[i]->chi1inv[c][d]);
            my_start += ntot;
          }
    }

  // Get the sizes of susceptibility lists for chiP[E_stuff] and chiP[H_stuff]
  // Since this info is copied to each chunk, we can just get it from the first chunk
  size_t num_sus[2] = {0, 0};
  susceptibility *Esus = chunks[0]->chiP[E_stuff];
  while (Esus) {
    num_sus[E_stuff] += 1;
    Esus = Esus->next;
  }
  susceptibility *Hsus = chunks[0]->chiP[H_stuff];
  while (Hsus) {
    num_sus[H_stuff] += 1;
    Hsus = Hsus->next;
  }

  {
    // Write the number of susceptibilites
    size_t len = 2;
    file.create_data("num_sus", 1, &len);
    if (am_master() || !single_parallel_file) {
      size_t start = 0;
      size_t ntot = 2;
      file.write_chunk(1, &start, &ntot, num_sus);
    }
  }

  // Get number of non-null sigma entries for each chiP in each chunk.
  // Assumes each susceptibility in the chiP[E_stuff] list has the
  // same number of non-null sigma elements. Likewise for chiP[H_stuff]
  std::vector<size_t> my_num_sigmas[2];
  std::vector<size_t> num_sigmas[2];
  for (int ft = 0; ft < 2; ++ft) {
    my_num_sigmas[ft].resize(num_chunks);
    num_sigmas[ft].resize(num_chunks);
  }

  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine()) {
      for (int ft = 0; ft < 2; ++ft) {
        susceptibility *sus = chunks[i]->chiP[ft];
        if (sus) {
          for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c) {
            for (int d = 0; d < 5; ++d) {
              if (sus->sigma[c][d]) { my_num_sigmas[ft][i] += 1; }
            }
          }
        }
      }
    }
  }

  // Write num_sigmas data.
  {
    size_t ntot = num_chunks * 2;
    file.create_data("num_sigmas", 1, &ntot);

    for (int i = 0; i < num_chunks; ++i) {
      if (chunks[i]->is_mine()) {
        for (int ft = 0; ft < 2; ++ft) {
          size_t start = ft * num_chunks + i;
          size_t count = 1;
          file.write_chunk(1, &start, &count, &my_num_sigmas[ft][i]);
        }
      }
    }
  }

  file.prevent_deadlock();
  for (int ft = 0; ft < 2; ++ft) {
    sum_to_all(my_num_sigmas[ft].data(), num_sigmas[ft].data(), num_chunks);
  }

  size_t num_E_sigmas = 0;
  size_t num_H_sigmas = 0;
  for (int i = 0; i < num_chunks; ++i) {
    if (num_sigmas[E_stuff][i] != 0) { num_E_sigmas = num_sigmas[E_stuff][i]; }
    if (num_sigmas[H_stuff][i] != 0) { num_H_sigmas = num_sigmas[H_stuff][i]; }
  }

  // Allocate space for component and direction of non-null sigmas
  std::vector<size_t> my_sigma_cd[2];
  my_sigma_cd[E_stuff].resize(num_E_sigmas * 2);
  my_sigma_cd[H_stuff].resize(num_H_sigmas * 2);

  std::vector<size_t> sigma_cd[2];
  sigma_cd[E_stuff].resize(num_E_sigmas * 2);
  sigma_cd[H_stuff].resize(num_H_sigmas * 2);

  // Find component and direction of non-null sigmas
  {
    for (int ft = 0; ft < 2; ++ft) {
      int j = 0;
      bool done = false;
      for (int i = 0; i < num_chunks; ++i) {
        if (done) { break; }
        if (chunks[i]->is_mine()) {
          susceptibility *sus = chunks[i]->chiP[ft];
          if (sus) {
            for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c) {
              for (int d = 0; d < 5; ++d) {
                if (sus->sigma[c][d]) {
                  my_sigma_cd[ft][j] = c;
                  my_sigma_cd[ft][j + 1] = d;
                  j += 2;
                  done = true;
                }
              }
            }
          }
        }
      }
    }
    bw_or_to_all(my_sigma_cd[E_stuff].data(), sigma_cd[E_stuff].data(), num_E_sigmas * 2);
    bw_or_to_all(my_sigma_cd[H_stuff].data(), sigma_cd[H_stuff].data(), num_H_sigmas * 2);
  }

  // Write location (component and direction) data of non-null sigmas (sigma[c][d])
  {
    size_t len = (num_E_sigmas + num_H_sigmas) * 2;
    file.create_data("sigma_cd", 1, &len);
    size_t start = 0;
    for (int ft = 0; ft < 2; ++ft) {
      if (am_master() || !single_parallel_file) {
        size_t count = sigma_cd[ft].size();
        file.write_chunk(1, &start, &count, sigma_cd[ft].data());
        start += count;
      }
    }
  }

  // Write the actual data in a particular non-null sigma[c][d] for each susceptibility in this
  // chunk's chiP lists.
  size_t nsig[2] = {num_E_sigmas, num_H_sigmas};
  for (int i = 0; i < num_chunks; ++i) {
    for (int ft = 0; ft < 2; ++ft) {
      if (nsig[ft] != 0 && num_sigmas[ft][i]) {
        for (size_t j = 0; j < nsig[ft] * 2; j += 2) {
          char dname[20];
          int c = sigma_cd[ft][j];
          int d = sigma_cd[ft][j + 1];
          snprintf(dname, 20, "%c_%d_sigma_%d_%d", ft == 0 ? 'E' : 'H', i, c, d);
          size_t ntot = chunks[i]->gv.ntot() * num_sus[ft];
          file.create_data(dname, 1, &ntot);
          if (chunks[i]->is_mine()) {
            susceptibility *sus = chunks[i]->chiP[ft];
            size_t start = 0;
            while (sus) {
              size_t count = chunks[i]->gv.ntot();
              file.write_chunk(1, &start, &count, sus->sigma[c][d]);
              sus = sus->next;
              start += count;
            }
          }
        }
      }
    }
  }

  write_susceptibility_params(&file, single_parallel_file, "E_params", E_stuff);
  write_susceptibility_params(&file, single_parallel_file, "H_params", H_stuff);
}

// Reconstruct the chiP lists of susceptibilities from the params hdf5 data
susceptibility *make_sus_list_from_params(h5file *file, int rank, size_t dims[3], size_t ntot) {
  susceptibility *sus = NULL;
  susceptibility *res = NULL;
  size_t start = 0;

  while (start < dims[0]) {
    size_t num_params_dims[3] = {1, 0, 0};
    realnum num_params;
    file->read_chunk(rank, &start, num_params_dims, &num_params);
    start += num_params_dims[0];
    if (num_params == 4) {
      // This is a lorentzian_susceptibility and the next 4 values in the dataset
      // are id, omega_0, gamma, and no_omega_0_denominator.
      size_t lorentzian_dims[3] = {4, 0, 0};
      realnum lorentzian_params[4];
      file->read_chunk(rank, &start, lorentzian_dims, lorentzian_params);
      start += lorentzian_dims[0];

      int id = (int)lorentzian_params[0];
      double omega_0 = lorentzian_params[1];
      double gamma = lorentzian_params[2];
      bool no_omega_0_denominator = (bool)lorentzian_params[3];
      if (sus) {
        sus->next = new lorentzian_susceptibility(omega_0, gamma, no_omega_0_denominator);
        sus->next->ntot = ntot;
        sus->next->set_id(id);
      }
      else {
        sus = new lorentzian_susceptibility(omega_0, gamma, no_omega_0_denominator);
        sus->ntot = ntot;
        sus->set_id(id);
        res = sus;
      }
      if (sus->next) sus = sus->next;
    }
    else if (num_params == 5) {
      // This is a noisy_lorentzian_susceptibility and the next 5 values in the dataset
      // are id, noise_amp, omega_0, gamma, and no_omega_0_denominator.
      size_t noisy_lorentzian_dims[3] = {5, 0, 0};
      realnum noisy_lorentzian_params[5];
      file->read_chunk(rank, &start, noisy_lorentzian_dims, noisy_lorentzian_params);
      start += noisy_lorentzian_dims[0];

      int id = (int)noisy_lorentzian_params[0];
      double noise_amp = noisy_lorentzian_params[1];
      double omega_0 = noisy_lorentzian_params[2];
      double gamma = noisy_lorentzian_params[3];
      bool no_omega_0_denominator = (bool)noisy_lorentzian_params[4];
      if (sus) {
        sus->next =
            new noisy_lorentzian_susceptibility(noise_amp, omega_0, gamma, no_omega_0_denominator);
        sus->next->ntot = ntot;
        sus->next->set_id(id);
      }
      else {
        sus =
            new noisy_lorentzian_susceptibility(noise_amp, omega_0, gamma, no_omega_0_denominator);
        sus->ntot = ntot;
        sus->set_id(id);
        res = sus;
      }
      if (sus->next) sus = sus->next;
    }
    else if (num_params == 8) {
      // This is a gyrotropic_susceptibility and the next 8 values in the dataset
      // are id, bias.x, bias.y, bias.z, omega_0, gamma, alpha, and model.
      size_t gyro_susc_dims[3] = {8, 0, 0};
      realnum gyro_susc_params[8];
      file->read_chunk(rank, &start, gyro_susc_dims, gyro_susc_params);
      start += gyro_susc_dims[0];

      int id = (int)gyro_susc_params[0];
      const vec bias(gyro_susc_params[1], gyro_susc_params[2], gyro_susc_params[3]);
      const double omega_0 = gyro_susc_params[4];
      const double gamma = gyro_susc_params[5];
      const double alpha = gyro_susc_params[6];
      const gyrotropy_model model = (gyrotropy_model)gyro_susc_params[7];
      if (sus) {
        sus->next = new gyrotropic_susceptibility(bias, omega_0, gamma, alpha, model);
        sus->next->ntot = ntot;
        sus->next->set_id(id);
      }
      else {
        sus = new gyrotropic_susceptibility(bias, omega_0, gamma, alpha, model);
        sus->ntot = ntot;
        sus->set_id(id);
        res = sus;
      }
      if (sus->next) sus = sus->next;
    }
    else { meep::abort("Invalid number of susceptibility parameters in structure::load"); }
  }
  return res;
}

void structure::set_chiP_from_file(h5file *file, const char *dataset, field_type ft) {
  int rank = 0;
  size_t dims[3] = {0, 0, 0};

  file->read_size(dataset, &rank, dims, 1);
  if (rank != 1) meep::abort("inconsistent data size in structure::set_chiP_from_file");

  if (dims[0] != 0) {
    for (int i = 0; i < num_chunks; ++i) {
      chunks[i]->chiP[ft] = make_sus_list_from_params(file, rank, dims, chunks[i]->gv.ntot());
    }
  }
}

binary_partition::binary_partition(int _proc_id) : proc_id(_proc_id), plane{NO_DIRECTION, 0.0} {}

binary_partition::binary_partition(const split_plane &_split_plane,
                                   std::unique_ptr<binary_partition> &&left_tree,
                                   std::unique_ptr<binary_partition> &&right_tree)
    : proc_id(-1), plane(_split_plane), left(std::move(left_tree)), right(std::move(right_tree)) {
  if (!left || !right) { meep::abort("Binary partition tree is required to be full"); }
}

binary_partition::binary_partition(const binary_partition &other)
    : proc_id(other.proc_id), plane(other.plane) {
  if (!other.is_leaf()) {
    left.reset(new binary_partition(*other.left));
    right.reset(new binary_partition(*other.right));
  }
}

bool binary_partition::is_leaf() const { return !left && !right; }

int binary_partition::get_proc_id() const {
  assert(is_leaf());
  return proc_id;
}

const split_plane &binary_partition::get_plane() const {
  assert(!is_leaf());
  return plane;
}

const binary_partition *binary_partition::left_tree() const {
  assert(!is_leaf());
  return left.get();
}

const binary_partition *binary_partition::right_tree() const {
  assert(!is_leaf());
  return right.get();
}

void split_by_binarytree(grid_volume gvol, std::vector<grid_volume> &result_gvs,
                         std::vector<int> &result_ids, const binary_partition *bp) {
  // reached a leaf
  if (bp->is_leaf()) {
    result_gvs.push_back(gvol);
    result_ids.push_back(bp->get_proc_id());
    return;
  }

  const auto &plane = bp->get_plane();
  int split_point = static_cast<int>((plane.pos - gvol.surroundings().in_direction_min(plane.dir)) /
                                         gvol.surroundings().in_direction(plane.dir) *
                                         gvol.num_direction(plane.dir) +
                                     0.5);
  // traverse left branch
  grid_volume left_gvol = gvol.split_at_fraction(false, split_point, plane.dir);
  split_by_binarytree(left_gvol, result_gvs, result_ids, bp->left_tree());

  // traverse right branch
  grid_volume right_gvol = gvol.split_at_fraction(true, split_point, plane.dir);
  split_by_binarytree(right_gvol, result_gvs, result_ids, bp->right_tree());
}

void structure::load_chunk_layout(const char *filename, boundary_region &br) {
  // Load chunk grid_volumes from a file
  h5file file(filename, h5file::READONLY, true);
  size_t sz = num_chunks * 3;
  double *origins = new double[sz];
  memset(origins, 0, sizeof(double) * sz);
  size_t *nums = new size_t[sz];
  memset(nums, 0, sizeof(size_t) * sz);

  int origins_rank;
  size_t origins_dims;
  file.read_size("gv_origins", &origins_rank, &origins_dims, 1);
  if (origins_rank != 1 || origins_dims != sz) { meep::abort("chunk mismatch in structure::load"); }
  if (am_master()) {
    size_t gv_origins_start = 0;
    file.read_chunk(1, &gv_origins_start, &origins_dims, origins);
  }
  file.prevent_deadlock();
  broadcast(0, origins, sz);

  int nums_rank;
  size_t nums_dims;
  file.read_size("gv_nums", &nums_rank, &nums_dims, 1);
  if (nums_rank != 1 || nums_dims != sz) { meep::abort("chunk mismatch in structure::load"); }
  if (am_master()) {
    size_t gv_nums_start = 0;
    file.read_chunk(1, &gv_nums_start, &nums_dims, nums);
  }
  file.prevent_deadlock();
  broadcast(0, nums, sz);

  std::vector<grid_volume> gvs;
  std::vector<int> ids;
  // Populate a vector with the new grid_volumes
  for (int i = 0; i < num_chunks; ++i) {
    int idx = i * 3;
    grid_volume new_gv = gv;
    vec new_origin(new_gv.dim);
    LOOP_OVER_DIRECTIONS(gv.dim, d) {
      new_origin.set_direction(d, origins[idx++]);
      new_gv.set_num_direction(d, nums[i * 3 + ((int)d % 3)]);
    }
    new_gv.set_origin(new_origin);
    gvs.push_back(new_gv);
    ids.push_back(i * count_processors() / num_chunks);
  }

  load_chunk_layout(gvs, ids, br);

  delete[] origins;
  delete[] nums;
}

void structure::load_chunk_layout(const std::vector<grid_volume> &gvs, const std::vector<int> &ids,
                                  boundary_region &br) {
  if (gvs.size() != size_t(num_chunks)) meep::abort("load_chunk_layout: wrong number of chunks.");
  // Recreate the chunks with the new grid_volumes
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->refcount-- <= 1) delete chunks[i];
    chunks[i] = new structure_chunk(gvs[i], v, Courant, ids[i] % count_processors());
    br.apply(this, chunks[i]);
  }
  check_chunks();
}

void structure::load(const char *filename, bool single_parallel_file) {
  h5file file(filename, h5file::READONLY, single_parallel_file, !single_parallel_file);

  if (verbosity > 0)
    printf("reading epsilon from file \"%s\" (%d)...\n", filename, single_parallel_file);

  /*
   * make/save a num_chunks x NUM_FIELD_COMPONENTS x 5 array counting
   * the number of entries in the chi1inv array for each chunk.
   *
   * When 'single_parallel_file' is true, we are creating a single block of data
   * for ALL chunks (that are merged using MPI). Otherwise, we are just
   * making a copy of just the chunks that are ours.
   */
  int my_num_chunks = 0;
  for (int i = 0; i < num_chunks; i++) {
    my_num_chunks += (single_parallel_file || chunks[i]->is_mine());
  }
  size_t num_chi1inv_size = my_num_chunks * NUM_FIELD_COMPONENTS * 5;
  std::vector<size_t> num_chi1inv(num_chi1inv_size);

  int rank;
  size_t dims[3], _dims[3] = {(size_t)my_num_chunks, NUM_FIELD_COMPONENTS, 5};
  size_t start[3] = {0, 0, 0};
  file.read_size("num_chi1inv", &rank, dims, 3);
  if (rank != 3 || _dims[0] != dims[0] || _dims[1] != dims[1] || _dims[2] != dims[2])
    meep::abort("chunk mismatch in structure::load");
  if (am_master() || !single_parallel_file) file.read_chunk(3, start, dims, num_chi1inv.data());

  if (single_parallel_file) {
    file.prevent_deadlock();
    broadcast(0, num_chi1inv.data(), dims[0] * dims[1] * dims[2]);
  }

  changing_chunks();

  // allocate data as needed and check sizes
  size_t my_ntot = 0;
  for (int i = 0, chunk_i = 0; i < num_chunks; i++) {
    if (chunks[i]->is_mine()) {
      size_t ntot = chunks[i]->gv.ntot();
      for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c)
        for (int d = 0; d < 5; ++d) {
          size_t n = num_chi1inv[(chunk_i * NUM_FIELD_COMPONENTS + c) * 5 + d];
          if (n == 0) {
            delete[] chunks[i]->chi1inv[c][d];
            chunks[i]->chi1inv[c][d] = NULL;
          }
          else {
            if (n != ntot) meep::abort("grid size mismatch %zd vs %zd in structure::load", n, ntot);
            if (!chunks[i]->chi1inv[c][d]) chunks[i]->chi1inv[c][d] = new realnum[ntot];
            my_ntot += ntot;
          }
        }
    }
    chunk_i += (chunks[i]->is_mine() || single_parallel_file);
  }

  // determine total dataset size and offset of this process's data
  size_t my_start = 0;
  size_t ntotal = my_ntot;
  if (single_parallel_file) {
    my_start = partial_sum_to_all(my_ntot) - my_ntot;
    ntotal = sum_to_all(my_ntot);
  }

  // read the data
  file.read_size("chi1inv", &rank, dims, 1);
  if (rank != 1 || dims[0] != ntotal) {
    meep::abort("inconsistent data size for chi1inv in structure::load (rank, dims[0]): "
                "(%d, %zu) != (1, %zu)",
                rank, dims[0], ntotal);
  }
  for (int i = 0; i < num_chunks; i++)
    if (chunks[i]->is_mine()) {
      size_t ntot = chunks[i]->gv.ntot();
      for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c)
        for (int d = 0; d < 5; ++d)
          if (chunks[i]->chi1inv[c][d]) {
            file.read_chunk(1, &my_start, &ntot, chunks[i]->chi1inv[c][d]);
            my_start += ntot;
          }
    }
  // Create susceptibilites from params datasets
  set_chiP_from_file(&file, "E_params", E_stuff);
  set_chiP_from_file(&file, "H_params", H_stuff);

  // Read the number of susceptibilities in the chiP lists.
  size_t num_sus[] = {0, 0};
  {
    int rank = 0;
    size_t dims[] = {0, 0, 0};
    file.read_size("num_sus", &rank, dims, 1);
    if (dims[0] > 0) {
      if (am_master() || !single_parallel_file) {
        size_t start = 0;
        size_t count = 2;
        file.read_chunk(rank, &start, &count, num_sus);
      }
    }
    if (single_parallel_file) {
      file.prevent_deadlock();
      broadcast(0, num_sus, 2);
    }
  }

  // Allocate and read non-null sigma entry data
  std::vector<size_t> num_sigmas[2];
  num_sigmas[E_stuff].resize(num_chunks);
  num_sigmas[H_stuff].resize(num_chunks);

  // Read num_sigmas data
  {
    int rank = 0;
    size_t dims[] = {0, 0, 0};
    file.read_size("num_sigmas", &rank, dims, 1);
    if (dims[0] != (size_t)num_chunks * 2) {
      meep::abort("inconsistent data size for num_sigmas in structure::load");
    }
    if (am_master() || !single_parallel_file) {
      size_t start = 0;
      size_t count = num_chunks;
      file.read_chunk(rank, &start, &count, num_sigmas[E_stuff].data());
      start += count;
      file.read_chunk(rank, &start, &count, num_sigmas[H_stuff].data());
    }
    if (single_parallel_file) {
      file.prevent_deadlock();
      broadcast(0, num_sigmas[E_stuff].data(), num_chunks);
      broadcast(0, num_sigmas[H_stuff].data(), num_chunks);
    }
  }

  // Allocate space for component and direction data of the non-null susceptibilities
  size_t nsig[2] = {0, 0};
  for (int ft = 0; ft < 2; ++ft) {
    for (int i = 0; i < num_chunks; ++i) {
      if (num_sigmas[ft][i] != 0) { nsig[ft] = num_sigmas[ft][i]; }
    }
  }

  std::vector<size_t> sigma_cd[2];
  sigma_cd[E_stuff].resize(nsig[E_stuff] * 2);
  sigma_cd[H_stuff].resize(nsig[H_stuff] * 2);

  // Read the component/direction data
  {
    int rank;
    size_t dims[] = {0, 0, 0};
    file.read_size("sigma_cd", &rank, dims, 1);
    if (dims[0] != 2 * (nsig[E_stuff] + nsig[H_stuff])) {
      meep::abort("inconsistent data size for sigma_cd in structure::load");
    }

    if (am_master() || !single_parallel_file) {
      size_t start = 0;
      for (int ft = 0; ft < 2; ++ft) {
        size_t count = sigma_cd[ft].size();
        file.read_chunk(rank, &start, &count, sigma_cd[ft].data());
        start += count;
      }
    }
    if (single_parallel_file) {
      file.prevent_deadlock();
      broadcast(0, sigma_cd[E_stuff].data(), nsig[E_stuff] * 2);
      broadcast(0, sigma_cd[H_stuff].data(), nsig[H_stuff] * 2);
    }
  }

  for (int ft = 0; ft < 2; ++ft) {
    for (int i = 0; i < num_chunks; ++i) {
      for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c) {
        for (int d = 0; d < 5; ++d) {
          susceptibility *sus = chunks[i]->chiP[ft];
          while (sus) {
            for (size_t j = 0; j < nsig[ft] * 2; j += 2) {
              int _c = sigma_cd[ft][j];
              int _d = sigma_cd[ft][j + 1];
              sus->trivial_sigma[_c][_d] = false;
            }
            sus = sus->next;
          }
        }
      }
    }
  }

  // Load sigma data
  for (int ft = 0; ft < 2; ++ft) {
    for (int i = 0; i < num_chunks; ++i) {
      for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c) {
        for (int d = 0; d < 5; ++d) {
          char dname[20];
          snprintf(dname, 20, "%c_%d_sigma_%d_%d", ft == 0 ? 'E' : 'H', i, c, d);
          if (file.dataset_exists(dname)) {
            int rank;
            size_t dims[3] = {0, 0, 0};
            file.read_size(dname, &rank, dims, 1);
            if (num_sigmas[ft][i] && chunks[i]->is_mine()) {
              susceptibility *sus = chunks[i]->chiP[ft];
              size_t start = 0;
              while (sus) {
                size_t count = chunks[i]->gv.ntot();
                if (sus->sigma[c][d]) { delete[] sus->sigma[c][d]; }
                sus->sigma[c][d] = new realnum[count];
                sus->trivial_sigma[c][d] = false;
                file.read_chunk(rank, &start, &count, sus->sigma[c][d]);
                sus = sus->next;
                start += count;
              }
            }
          }
        }
      }
    }
  }
}
} // namespace meep
