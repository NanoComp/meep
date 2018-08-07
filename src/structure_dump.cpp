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

// Writes the number of susceptibilites to an hdf5 file. This is the length of the susceptibility
// lists pointed to by chunks[i]->chiP[E_stuff] and chunks[i]->chiP[H_stuff]
void structure::write_num_susceptibilites(h5file *file, const char *dname, size_t *source) {
  size_t len = *source == 0 ? 0 : num_chunks;
  file->create_data(dname, 1, &len);
  for (size_t i = 0; i < len; i++)
    if (chunks[i]->is_mine()) {
      size_t start = i;
      size_t ntot = 1;
      file->write_chunk(1, &start, &ntot, source);
    }
}

// The number of sigmas is the number of non-null elements of
// chunks[i]->chiP[E_stuff|H_stuff]->sigma
void structure::write_num_sigmas(h5file *file, const char *dname, size_t num_sus, size_t *num_sigmas) {
  size_t sus_ntot = sum_to_all(num_sus);
  size_t sigma_start = partial_sum_to_all(num_sus) - num_sus;

  file->create_data(dname, 1, &sus_ntot);
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine()) {
      file->write_chunk(1, &sigma_start, &num_sus, num_sigmas);
    }
  }
  file->prevent_deadlock();
}

// Write location (component and direction) of each non-null sigma (sigma[c][d])
void structure::write_component_direction_data(h5file *file, const char* dname, size_t *num_sus, size_t **num_sigmas,
                                               size_t ***sigma_cd, int EorH) {
  size_t my_sigma_cd_ntot = 0;
  for (size_t i = 0; i < num_sus[EorH]; ++i) {
    my_sigma_cd_ntot += num_sigmas[EorH][i] * 2;
  }

  size_t sigma_cd_ntot = sum_to_all(my_sigma_cd_ntot);
  size_t my_sigma_cd_start = partial_sum_to_all(my_sigma_cd_ntot) - my_sigma_cd_ntot;

  file->create_data(dname, 1, &sigma_cd_ntot);
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine() && chunks[i]->chiP[EorH]) {
      for (size_t j = 0; j < num_sus[EorH]; ++j) {
        size_t cd_count = num_sigmas[EorH][j] * 2;
        file->write_chunk(1, &my_sigma_cd_start, &cd_count, sigma_cd[EorH][j]);
        my_sigma_cd_start += cd_count;
      }
    }
  }
  file->prevent_deadlock();
}

// Write the actual data in a particular non-null sigma[c][d] for each susceptibility in this
// chunk's chiP lists.
void structure::write_sigma_data(h5file *file, const char *dname, size_t *num_sus, size_t **num_sigmas, int EorH) {
  size_t my_ntot = 0;
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine() && chunks[i]->chiP[EorH]) {
      for (size_t j = 0; j < num_sus[EorH]; ++j) {
        my_ntot += num_sigmas[EorH][j] * chunks[i]->chiP[0]->ntot;
      }
    }
  }

  size_t my_start = partial_sum_to_all(my_ntot ) - my_ntot;
  size_t ntotal = sum_to_all(my_ntot );

  // Write sigma data
  file->create_data(dname, 1, &ntotal);
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine() && chunks[i]->chiP[EorH]) {
      susceptibility *sus = chunks[i]->chiP[EorH];
      while (sus) {
        sus->dump(file, &my_start);
        sus = sus->next;
      }
    }
  }
  file->prevent_deadlock();
}

// Write the parameters required to reconstruct the susceptibility (id, noise_amp (for noisy), omega_0,
// gamma, no_omega_0_denominator)
void structure::write_susceptibility_params(h5file *file, const char *dname, int EorH) {
  // Get number of susceptibility params
  size_t params_ntotal = 0;
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine() && chunks[i]->chiP[EorH]) {
      susceptibility *sus = chunks[i]->chiP[EorH];
      while (sus) {
        params_ntotal += sus->get_num_params() + 1;
        sus = sus->next;
      }
    }
  }

  // Write params
  size_t params_start = 0;
  file->create_data(dname, 1, &params_ntotal);
  if (am_master()) {
    for (int i = 0; i < num_chunks; ++i) {
      if (chunks[i]->chiP[EorH] && chunks[i]->is_mine()) {
        susceptibility *sus = chunks[i]->chiP[EorH];
        while (sus) {
          sus->dump_params(file, &params_start);
          sus = sus->next;
        }
      }
    }
  }
}

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
  write_num_susceptibilites(&file, "num_E_sus", &my_num_sus[0]);
  write_num_susceptibilites(&file, "num_H_sus", &my_num_sus[1]);
  file.prevent_deadlock();

  // Get number of non-null sigma entries for each susceptibility of each chiP
  size_t *num_sigmas[2];
  for (int i = 0; i < 2; ++i) {
    num_sigmas[i] = new size_t[my_num_sus[i]];
    memset(num_sigmas[i], 0, sizeof(size_t) * my_num_sus[i]);
  }

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
  write_num_sigmas(&file, "num_E_sigmas", my_num_sus[E_stuff], num_sigmas[E_stuff]);
  write_num_sigmas(&file, "num_H_sigmas", my_num_sus[H_stuff], num_sigmas[H_stuff]);

  // Get component and direction of non-null sigmas
  size_t **sigma_cd[2] = {NULL, NULL};
  for (int ft = 0; ft < 2; ++ft) {
    sigma_cd[ft] = new size_t*[my_num_sus[ft]];
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

  write_component_direction_data(&file, "sigma_cd_E", my_num_sus, num_sigmas, sigma_cd, E_stuff);
  write_component_direction_data(&file, "sigma_cd_H", my_num_sus, num_sigmas, sigma_cd, H_stuff);

  write_sigma_data(&file, "E_sigma", my_num_sus, num_sigmas, E_stuff);
  write_sigma_data(&file, "H_sigma", my_num_sus, num_sigmas, H_stuff);

  write_susceptibility_params(&file, "E_params", E_stuff);
  write_susceptibility_params(&file, "H_params", H_stuff);

  for (int ft = 0; ft < 2; ++ft) {
    for (size_t i = 0; i < my_num_sus[ft]; ++i) {
      delete[] sigma_cd[ft][i];
    }
    delete[] sigma_cd[ft];
    delete[] num_sigmas[ft];
  }
}

// Reconstruct the chiP lists of susceptibilities from the params hdf5 data
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

// Read the number of susceptibilities in the chiP lists.
size_t structure::read_num_susceptibilities(h5file *file, const char *dname) {
  size_t result = 0;
  int rank = 0;
  size_t dims[] = {0, 0, 0};
  file->read_size(dname, &rank, dims, 1);
  if (dims[0] > 0) {
    for (int i = 0; i < num_chunks; ++i) {
      if (chunks[i]->is_mine()) {
        size_t start = i;
        size_t count = 1;
        file->read_chunk(rank, &start, &count, &result);
      }
    }
  }
  file->prevent_deadlock();

  return result;
}

// Read the num_sigma data that was created by write_num_sigmas
size_t *structure::read_num_sigmas(h5file *file, const char *dname, size_t num_sus) {
  size_t *result = new size_t[num_sus];
  size_t start = partial_sum_to_all(num_sus) - num_sus;

  int rank = 0;
  size_t dims[] = {0, 0, 0};
  file->read_size(dname, &rank, dims, 1);
  if (dims[0] != num_chunks * num_sus)
    abort("inconsistent data size in structure::load");

  if (num_sus > 0) {
    for (int i = 0; i < num_chunks; ++i) {
      if (chunks[i]->is_mine()) {
        file->read_chunk(rank, &start, &num_sus, result);
      }
    }
  }

  file->prevent_deadlock();
  return result;
}

// Read the cd data created by write_component_direction_data
size_t **structure::read_component_direction_data(h5file *file, const char *dname, size_t num_sus,
                                                  size_t *num_sigmas) {
  size_t **result = new size_t*[num_sus];
  for (size_t i = 0; i < num_sus; ++i) {
    result[i] = new size_t[num_sigmas[i] * 2];
  }

  size_t num_pairs = 0;
  for (size_t i = 0; i < num_sus; ++i) {
    num_pairs += num_sigmas[i];
  }

  size_t start = partial_sum_to_all(num_pairs * 2) - num_pairs * 2;

  int rank = 0;
  size_t dims[] = {0, 0, 0};
  file->read_size(dname, &rank, dims, 1);

  if (num_sus > 0) {
    for (int i = 0; i < num_chunks; ++i) {
      if (chunks[i]->is_mine()) {
        for (size_t j = 0; j < num_sus; ++j) {
          size_t count = num_sigmas[j] * 2;
          file->read_chunk(rank, &start, &count, result[j]);
          start += count;
        }
      }
    }
  }
  file->prevent_deadlock();
  return result;
}

// Read the sigma data into the non-null sigma[c][d] entries for each susceptibility
// in this chunk's chiP list.
void structure::read_sigma(h5file *file, const char *dname, size_t num_sus, size_t *num_sigmas, size_t **sigma_cd,
                           int EorH) {
  size_t my_ntot = 0;
  for (int i = 0; i < num_chunks; ++i) {
    if (chunks[i]->is_mine() && chunks[i]->chiP[EorH]) {
      for (size_t j = 0; j < num_sus; ++j) {
        my_ntot += num_sigmas[j] * chunks[i]->chiP[EorH]->ntot;
      }
    }
  }

  size_t start = partial_sum_to_all(my_ntot) - my_ntot;

  // Load sigma data into susceptibilites
  int rank = 0;
  size_t dims[3] = {0, 0, 0};
  file->read_size(dname, &rank, dims, 1);
  if (dims[0] > 0) {
    for (int i = 0; i < num_chunks; ++i) {
      if (chunks[i]->is_mine()) {
        susceptibility *sus = chunks[i]->chiP[EorH];
        size_t n = 0;
        while (sus) {
          sus->load(file, &start, num_sigmas[n], sigma_cd[n]);
          sus = sus->next;
          n++;
        }
      }
    }
  }
  file->prevent_deadlock();
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
  size_t my_num_E_sus = read_num_susceptibilities(&file, "num_E_sus");
  size_t my_num_H_sus = read_num_susceptibilities(&file, "num_H_sus");

  // Allocate and read non-null sigma entry data
  size_t *num_E_sigmas = read_num_sigmas(&file, "num_E_sigmas", my_num_E_sus);
  size_t *num_H_sigmas = read_num_sigmas(&file, "num_H_sigmas", my_num_H_sus);

  // Allocate and read component and direction data of the non-null susceptibilities
  size_t **sigma_cd_E = read_component_direction_data(&file, "sigma_cd_E", my_num_E_sus, num_E_sigmas);
  size_t **sigma_cd_H = read_component_direction_data(&file, "sigma_cd_H", my_num_H_sus, num_H_sigmas);

  // Load sigma data
  read_sigma(&file, "E_sigma", my_num_E_sus, num_E_sigmas, sigma_cd_E, E_stuff);
  read_sigma(&file, "H_sigma", my_num_H_sus, num_H_sigmas, sigma_cd_H, H_stuff);

  delete[] num_E_sigmas;
  delete[] num_H_sigmas;
  for (size_t i = 0; i < my_num_E_sus; ++i) {
    delete[] sigma_cd_E[i];
  }
  for (size_t i = 0; i < my_num_H_sus; ++i) {
    delete[] sigma_cd_H[i];
  }
  delete[] sigma_cd_E;
  delete[] sigma_cd_H;
}
} // namespace meep
