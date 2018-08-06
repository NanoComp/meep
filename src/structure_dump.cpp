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
  file.prevent_deadlock();

  // Get sizes of susceptibility lists for chiP[E_stuff] and chiP[H_stuff]
  size_t my_num_sus[2] = {0, 0};
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine()) {
      susceptibility *Esus = chunks[i]->chiP[E_stuff];
      while (Esus) {
        my_num_sus[0] += 1;
        Esus = Esus->next;
      }
      susceptibility *Hsus = chunks[i]->chiP[H_stuff];
      while (Hsus) {
        my_num_sus[1] += 1;
        Hsus = Hsus->next;
      }
    }
  }

  // Write susceptibility list sizes
  size_t E_sus_len = my_num_sus[0] == 0 ? 0 : num_chunks;
  file.create_data("num_E_sus", 1, &E_sus_len);
  for (size_t i = 0; i < E_sus_len; i++)
    if (chunks[i]->is_mine()) {
      size_t my_start = i;
      size_t ntot = 1;
      file.write_chunk(1, &my_start, &ntot, &my_num_sus[0]);
    }

  size_t H_sus_len = my_num_sus[1] == 0 ? 0 : num_chunks;
  file.create_data("num_H_sus", 1, &H_sus_len);
  for (size_t i = 0; i < H_sus_len; i++)
    if (chunks[i]->is_mine()) {
      size_t my_start = i;
      size_t ntot = 1;
      file.write_chunk(1, &my_start, &ntot, &my_num_sus[1]);
    }
  file.prevent_deadlock();

  // Get number of non-null sigma entries for each susceptibility of each chiP
  size_t *num_sigmas[2];
  num_sigmas[0] = new size_t[my_num_sus[0]];
  num_sigmas[1] = new size_t[my_num_sus[1]];
  memset(num_sigmas[0], 0, sizeof(size_t) * my_num_sus[0]);
  memset(num_sigmas[1], 0, sizeof(size_t) * my_num_sus[1]);

  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine()) {
      for (int ft = 0; ft < 2; ++ft) {
        susceptibility *sus = chunks[i]->chiP[ft];
        size_t n = 0;
        while (sus) {
          for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c) {
            for (int d = 0; d < 5; ++d) {
              if (sus->sigma[c][d]) {
                num_sigmas[ft][n] += 1;
              }
            }
          }
          n++;
          sus = sus->next;
        }
      }
    }
  }

  // Write num_sigmas data
  size_t sus_ntot[2] = {0, 0};
  sum_to_all(my_num_sus, sus_ntot, 2);
  size_t my_E_sigma_start = partial_sum_to_all(my_num_sus[0]) - my_num_sus[0];
  size_t my_H_sigma_start = partial_sum_to_all(my_num_sus[1]) - my_num_sus[1];

  file.create_data("num_E_sigmas", 1, &sus_ntot[0]);
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine()) {
      file.write_chunk(1, &my_E_sigma_start, &my_num_sus[0], num_sigmas[0]);
    }
  }

  file.create_data("num_H_sigmas", 1, &sus_ntot[1]);
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine()) {
      file.write_chunk(1, &my_H_sigma_start, &my_num_sus[1], num_sigmas[1]);
    }
  }
  file.prevent_deadlock();

  // Get component and direction of non-null sigmas
  size_t **sigma_cd[2] = {NULL, NULL};
  sigma_cd[0] = new size_t*[my_num_sus[0]];
  sigma_cd[1] = new size_t*[my_num_sus[1]];

  for (int ft = 0; ft < 2; ++ft) {
    for (size_t i = 0; i < my_num_sus[ft]; ++i) {
      sigma_cd[ft][i] = new size_t[num_sigmas[ft][i] * 2];
      memset(sigma_cd[ft][i], 0, sizeof(size_t) * num_sigmas[ft][i] * 2);
    }
  }

  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine()) {
      for (int ft = 0; ft < 2; ++ft) {
        susceptibility *sus = chunks[i]->chiP[ft];
        size_t n = 0;
        while (sus) {
          int j = 0;
          for (int c = 0; c < NUM_FIELD_COMPONENTS; ++c) {
            for (int d = 0; d < 5; ++d) {
              if (sus->sigma[c][d]) {
                sigma_cd[ft][n][j] = c;
                sigma_cd[ft][n][j + 1] = d;
                j += 2;
              }
            }
          }
          sus = sus->next;
          n++;
        }
      }
    }
  }

  size_t my_sigma_cd_E_ntot = 0;
  for (size_t i = 0; i < my_num_sus[E_stuff]; ++i) {
    my_sigma_cd_E_ntot += num_sigmas[E_stuff][i] * 2;
  }

  size_t my_sigma_cd_H_ntot = 0;
  for (size_t i = 0; i < my_num_sus[H_stuff]; ++i) {
    my_sigma_cd_H_ntot += num_sigmas[H_stuff][i] * 2;
  }

  size_t sigma_cd_E_ntot = sum_to_all(my_sigma_cd_E_ntot);
  size_t sigma_cd_H_ntot = sum_to_all(my_sigma_cd_H_ntot);

  size_t my_sigma_cd_E_start = partial_sum_to_all(my_sigma_cd_E_ntot) - my_sigma_cd_E_ntot;
  size_t my_sigma_cd_H_start = partial_sum_to_all(my_sigma_cd_H_ntot) - my_sigma_cd_H_ntot;

  // Write location (component and direction) of each non-null sigma
  file.create_data("sigma_cd_E", 1, &sigma_cd_E_ntot);
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine() && chunks[i]->chiP[E_stuff]) {
      for (size_t j = 0; j < my_num_sus[E_stuff]; ++j) {
        size_t cd_count = num_sigmas[E_stuff][j] * 2;
        file.write_chunk(1, &my_sigma_cd_E_start, &cd_count, sigma_cd[E_stuff][j]);
        my_sigma_cd_E_start += cd_count;
      }
    }
  }

  file.create_data("sigma_cd_H", 1, &sigma_cd_H_ntot);
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine() && chunks[i]->chiP[H_stuff]) {
      for (size_t j = 0; j < my_num_sus[H_stuff]; ++j) {
        size_t cd_count = num_sigmas[H_stuff][j];
        file.write_chunk(1, &my_sigma_cd_H_start, &cd_count, sigma_cd[H_stuff][j]);
        my_sigma_cd_H_start += cd_count;
      }
    }
  }
  file.prevent_deadlock();

  size_t my_tot_E_sus_points = 0;
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine() && chunks[i]->chiP[E_stuff]) {
      for (size_t j = 0; j < my_num_sus[E_stuff]; ++j) {
        my_tot_E_sus_points += num_sigmas[E_stuff][j] * chunks[i]->chiP[0]->ntot;
      }
    }
  }
  size_t my_tot_H_sus_points = 0;
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine() && chunks[i]->chiP[H_stuff]) {
      for (size_t j = 0; j < my_num_sus[H_stuff]; ++j) {
        my_tot_H_sus_points += num_sigmas[H_stuff][j] * chunks[i]->chiP[0]->ntot;
      }
    }
  }

  size_t my_E_sus_start = partial_sum_to_all(my_tot_E_sus_points) - my_tot_E_sus_points;
  size_t my_H_sus_start = partial_sum_to_all(my_tot_H_sus_points) - my_tot_H_sus_points;
  size_t E_sus_ntotal = sum_to_all(my_tot_E_sus_points);
  size_t H_sus_ntotal = sum_to_all(my_tot_H_sus_points);

  // Write sigma data
  file.create_data("E_sigma", 1, &E_sus_ntotal);
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine() && chunks[i]->chiP[E_stuff]) {
      susceptibility *sus = chunks[i]->chiP[E_stuff];
      while (sus) {
        sus->dump(&file, &my_E_sus_start);
        sus = sus->next;
      }
    }
  }

  file.create_data("H_sigma", 1, &H_sus_ntotal);
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine() && chunks[i]->chiP[H_stuff]) {
      susceptibility *sus = chunks[i]->chiP[H_stuff];
      while (sus) {
        sus->dump(&file, &my_H_sus_start);
        sus = sus->next;
      }
    }
  }

  for (int ft = 0; ft < 2; ++ft) {
    for (size_t i = 0; i < my_num_sus[ft]; ++i) {
      delete[] sigma_cd[ft][i];
    }
    delete[] sigma_cd[ft];
    delete[] num_sigmas[ft];
  }

  // Get number of susceptibility params
  size_t E_params_ntotal = 0;
  size_t H_params_ntotal = 0;
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine()) {
      for (int ft = 0; ft < NUM_FIELD_TYPES; ++ft) {
        if (chunks[i]->chiP[ft]) {
          susceptibility *sus = chunks[i]->chiP[ft];
          while (sus) {
            if (ft == 0)
              E_params_ntotal += sus->get_num_params() + 1;
            else if (ft == 1)
              H_params_ntotal += sus->get_num_params() + 1;
            sus = sus->next;
          }
        }
      }
    }
  }

  // Write params for E_stuff susceptibilities
  size_t params_start = 0;
  file.create_data("E_params", 1, &E_params_ntotal);
  if (am_master()) {
    for (int i = 0; i < num_chunks; ++i) {
      if (chunks[i]->chiP[E_stuff] && chunks[i]->is_mine()) {
        susceptibility *sus = chunks[i]->chiP[E_stuff];
        while (sus) {
          sus->dump_params(&file, &params_start);
          sus = sus->next;
        }
      }
    }
  }

  // Write params for H_stuff susceptibilities
  params_start = 0;
  file.create_data("H_params", 1, &H_params_ntotal);
  if (am_master()) {
    for (int i = 0; i < num_chunks; ++i) {
      if (chunks[i]->chiP[H_stuff] && chunks[i]->is_mine()) {
        susceptibility *sus = chunks[i]->chiP[H_stuff];
        while (sus) {
          sus->dump_params(&file, &params_start);
          sus = sus->next;
        }
      }
    }
  }
}

susceptibility *make_sus_list_from_params(h5file *file, int rank, size_t *start, size_t dims[3], size_t ntot) {
  susceptibility *sus = NULL;
  susceptibility *res = NULL;
  while (*start < dims[0] - 1) {
    size_t num_params_dims[3] = {1, 0, 0};
    realnum num_params;
    file->read_chunk(rank, start, num_params_dims, &num_params);
    *start += num_params_dims[0];

    if (num_params == 4) {
      // This is a lorentzian_susceptibility and the next 4 values in the dataset
      // are id, omega_0, gamma, and no_omega_0_denominator.
      size_t lorentzian_dims[3] = {4, 0, 0};
      realnum lorentzian_params[4];
      file->read_chunk(rank, start, lorentzian_dims, lorentzian_params);
      *start += lorentzian_dims[0];

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
      if (sus->next)
        sus = sus->next;
    }
    else if (num_params == 5) {
      // This is a noisy_lorentzian_susceptibility and the next 5 values in the dataset
      // are id, noise_amp, omega_0, gamma, and no_omega_0_denominator.
      size_t noisy_lorentzian_dims[3] = {5, 0, 0};
      realnum noisy_lorentzian_params[5];
      file->read_chunk(rank, start, noisy_lorentzian_dims, noisy_lorentzian_params);
      *start += noisy_lorentzian_dims[0];

      int id = (int)noisy_lorentzian_params[0];
      double noise_amp = noisy_lorentzian_params[1];
      double omega_0 = noisy_lorentzian_params[2];
      double gamma = noisy_lorentzian_params[3];
      bool no_omega_0_denominator = (bool)noisy_lorentzian_params[4];
      if (sus) {
        sus->next = new noisy_lorentzian_susceptibility(noise_amp, omega_0, gamma, no_omega_0_denominator);
        sus->next->ntot = ntot;
        sus->next->set_id(id);
      }
      else {
        sus = new noisy_lorentzian_susceptibility(noise_amp, omega_0, gamma, no_omega_0_denominator);
        sus->ntot = ntot;
        sus->set_id(id);
        res = sus;
      }
      if (sus->next)
        sus = sus->next;
    }
    else {
      abort("Invalid number of susceptibility parameters in structure::load");
    }
  }
  return res;
}

void structure::set_chiP_from_file(h5file *file, const char *dataset, field_type ft) {
  int rank = 0;
  size_t dims[3] = {0, 0, 0};
  size_t start = 0;

  file->read_size(dataset, &rank, dims, 1);
  if (rank != 1)
    abort("inconsistent data size in structure::load");

  if (dims[0] != 0) {
    for (int i = 0; i < num_chunks; ++i) {
      if (chunks[i]->is_mine()) {
        chunks[i]->chiP[ft] = make_sus_list_from_params(file, rank, &start, dims, chunks[i]->gv.ntot());
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

  // Create susceptibilites from params datasets
  set_chiP_from_file(&file, "E_params", E_stuff);
  set_chiP_from_file(&file, "H_params", H_stuff);

  // Get number of H and E susceptibilites on this processor
  size_t my_num_E_sus = 0;
  int num_E_sus_rank = 0;
  size_t num_E_sus_dims[] = {0, 0, 0};
  file.read_size("num_E_sus", &num_E_sus_rank, num_E_sus_dims, 1);
  if (num_E_sus_dims[0] > 0) {
    for (int i = 0; i < num_chunks; ++i) {
      if (chunks[i]->is_mine()) {
        size_t start = i;
        size_t count = 1;
        file.read_chunk(num_E_sus_rank, &start, &count, &my_num_E_sus);
      }
    }
  }
  size_t my_num_H_sus = 0;
  int num_H_sus_rank = 0;
  size_t num_H_sus_dims[] = {0, 0, 0};
  file.read_size("num_H_sus", &num_H_sus_rank, num_H_sus_dims, 1);
  if (num_H_sus_dims[0] > 0) {
    for (int i = 0; i < num_chunks; ++i) {
      if (chunks[i]->is_mine()) {
        size_t start = i;
        size_t count = 1;
        file.read_chunk(num_H_sus_rank, &start, &count, &my_num_H_sus);
      }
    }
  }
  file.prevent_deadlock();

  // Get non-null sigma entry data
  size_t *num_E_sigmas = new size_t[my_num_E_sus];
  size_t my_E_sigma_start = partial_sum_to_all(my_num_E_sus) - my_num_E_sus;

  int num_E_sigma_rank = 0;
  size_t num_E_sigma_dims[] = {0, 0, 0};
  file.read_size("num_E_sigmas", &num_E_sigma_rank, num_E_sigma_dims, 1);
  if (num_E_sigma_dims[0] != num_chunks * my_num_E_sus)
    abort("inconsistent data size in structure::load");

  if (my_num_E_sus > 0) {
    for (int i = 0; i < num_chunks; ++i) {
      if (chunks[i]->is_mine()) {
        file.read_chunk(num_E_sigma_rank, &my_E_sigma_start, &my_num_E_sus, num_E_sigmas);
      }
    }
  }

  file.prevent_deadlock();

  size_t *num_H_sigmas = new size_t[my_num_H_sus];
  size_t my_H_sigma_start = partial_sum_to_all(my_num_H_sus) - my_num_H_sus;

  int num_H_sigma_rank = 0;
  size_t num_H_sigma_dims[] = {0, 0, 0};
  file.read_size("num_H_sigmas", &num_H_sigma_rank, num_H_sigma_dims, 1);
  if (num_H_sigma_dims[0] != num_chunks * my_num_H_sus)
    abort("inconsistent data size in structure::load");

  if (my_num_H_sus > 0) {
    for (int i = 0; i < num_chunks; ++i) {
      if (chunks[i]->is_mine()) {
        file.read_chunk(num_H_sigma_rank, &my_H_sigma_start, &my_num_H_sus, num_H_sigmas);
      }
    }
  }
  file.prevent_deadlock();

  // Get component and direction data of the non-null susceptibilities
  size_t **sigma_cd_E = new size_t*[my_num_E_sus];
  for (size_t i = 0; i < my_num_E_sus; ++i) {
    sigma_cd_E[i] = new size_t[num_E_sigmas[i] * 2];
  }

  size_t my_num_E_cd_pairs = 0;
  for (size_t i = 0; i < my_num_E_sus; ++i) {
    my_num_E_cd_pairs += num_E_sigmas[i];
  }

  size_t my_sigma_cd_E_start = partial_sum_to_all(my_num_E_cd_pairs * 2) - my_num_E_cd_pairs * 2;

  int sigma_cd_E_rank = 0;
  size_t sigma_cd_E_dims[] = {0, 0, 0};
  file.read_size("sigma_cd_E", &sigma_cd_E_rank, sigma_cd_E_dims, 1);

  if (my_num_E_sus > 0) {
    for (int i = 0; i < num_chunks; ++i) {
      if (chunks[i]->is_mine()) {
        for (size_t j = 0; j < my_num_E_sus; ++j) {
          size_t count = num_E_sigmas[j] * 2;
          file.read_chunk(sigma_cd_E_rank, &my_sigma_cd_E_start, &count, sigma_cd_E[j]);
          my_sigma_cd_E_start += count;
        }
      }
    }
  }
  file.prevent_deadlock();

  size_t **sigma_cd_H = new size_t*[my_num_H_sus];
  for (size_t i = 0; i < my_num_H_sus; ++i) {
    sigma_cd_H[i] = new size_t[num_H_sigmas[i] * 2];
  }

  size_t my_num_H_cd_pairs = 0;
  for (size_t i = 0; i < my_num_H_sus; ++i) {
    my_num_H_cd_pairs += num_H_sigmas[i];
  }

  size_t my_sigma_cd_H_start = partial_sum_to_all(my_num_H_cd_pairs * 2) - my_num_H_cd_pairs * 2;

  int sigma_cd_H_rank = 0;
  size_t sigma_cd_H_dims[] = {0, 0, 0};
  file.read_size("sigma_cd_H", &sigma_cd_H_rank, sigma_cd_H_dims, 1);

  if (my_num_H_sus > 0) {
    for (int i = 0; i < num_chunks; ++i) {
      if (chunks[i]->is_mine()) {
        for (size_t j = 0; j < my_num_H_sus; ++j) {
          size_t count = num_H_sigmas[j] * 2;
          file.read_chunk(sigma_cd_H_rank, &my_sigma_cd_H_start, &count, sigma_cd_H[j]);
          my_sigma_cd_H_start += count;
        }
      }
    }
  }
  file.prevent_deadlock();

  // Get total sigma size on this processor
  size_t my_E_sigma_ntot = 0;
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine() && chunks[i]->chiP[E_stuff]) {
      for (size_t j = 0; j < my_num_E_sus; ++j) {
        my_E_sigma_ntot += num_E_sigmas[j] * chunks[i]->chiP[E_stuff]->ntot;
      }
    }
  }

  size_t my_H_sigma_ntot = 0;
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine() && chunks[i]->chiP[H_stuff]) {
      for (size_t j = 0; j < my_num_H_sus; ++j) {
        my_H_sigma_ntot += num_H_sigmas[j] * chunks[i]->chiP[H_stuff]->ntot;
      }
    }
  }

  size_t E_sigma_start = partial_sum_to_all(my_E_sigma_ntot) - my_E_sigma_ntot;
  size_t H_sigma_start = partial_sum_to_all(my_H_sigma_ntot) - my_H_sigma_ntot;

  // Load sigma data into susceptibilites
  int E_sigma_rank = 0;
  size_t E_sigma_dims[3] = {0, 0, 0};
  file.read_size("E_sigma", &E_sigma_rank, E_sigma_dims, 1);
  if (E_sigma_dims[0] > 0) {
    for (int i = 0; i < num_chunks; ++i) {
      if (chunks[i]->is_mine()) {
        susceptibility *sus = chunks[i]->chiP[E_stuff];
        size_t n = 0;
        while (sus) {
          sus->load(&file, &E_sigma_start, num_E_sigmas[n], sigma_cd_E[n]);
          sus = sus->next;
          n++;
        }
      }
    }
  }

  int H_sigma_rank = 0;
  size_t H_sigma_dims[3] = {0, 0, 0};
  file.read_size("H_sigma", &H_sigma_rank, H_sigma_dims, 1);
  if (H_sigma_dims[0] > 0) {
    for (int i = 0; i < num_chunks; ++i) {
      if (chunks[i]->is_mine()) {
        susceptibility *sus = chunks[i]->chiP[H_stuff];
        size_t n = 0;
        while (sus) {
          sus->load(&file, &H_sigma_start, num_H_sigmas[n], sigma_cd_H[n]);
          sus = sus->next;
          n++;
        }
      }
    }
  }

  delete[] num_E_sigmas;
  delete[] num_H_sigmas;
  for (size_t i = 0; i < my_num_E_sus; ++i) {
    delete[] sigma_cd_E[i];
  }
  delete[] sigma_cd_E;
  for (size_t i = 0; i < my_num_H_sus; ++i) {
    delete[] sigma_cd_H[i];
  }
  delete[] sigma_cd_H;
}
} // namespace meep
