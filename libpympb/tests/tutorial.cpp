#include <cassert>
#include "pympb.hpp"

void test_create_maxwell_data() {
  int num_bands = 8;
  bool match_frequency = false;
  int parity = 0;
  double resolution = 32;
  double eigensolver_tol = 1.0e-7;
  vector3 lattice_size = {1, 1, 0};

  bool reset_fields = true;

  py_mpb::mode_solver ms(num_bands, match_frequency, parity, resolution,
                         lattice_size, eigensolver_tol);

  ms.init(parity, reset_fields);

  maxwell_data *md = ms.mdata;

  assert(md->nx == 32);
  assert(md->ny == 32);
  assert(md->nz == 1);
  assert(md->local_nx == 32);
  assert(md->local_ny == 32);
  assert(md->local_x_start == 0);
  assert(md->local_y_start == 0);
  assert(md->last_dim == 32);
  assert(md->last_dim_size == 32);
  assert(md->other_dims == 32);
  assert(md->num_bands == 8);
  assert(md->N == 1024);
  assert(md->local_N == 1024);
  assert(md->N_start == 0);
  assert(md->alloc_N == 1024);
  assert(md->fft_output_size == 1024);
  assert(md->max_fft_bands == 8);
  assert(md->num_fft_bands == 8);
  assert(current_k[0] == 0);
  assert(current_k[1] == 0);
  assert(current_k[2] == 0);
  assert(md->parity == 0);
  assert(md->nplans == 0);
  assert(md->zero_k == 0);
  assert(md->eps_inv->m00 == 0);
  assert(md->eps_inv->m01 == 0);
  assert(md->eps_inv->m02 == 0);
  assert(md->eps_inv->m11 == 0);
  assert(md->eps_inv->m12 == 0);
  assert(md->eps_inv->m22 == 0);
  assert(md->eps_inv_mean == 1);
  assert(md->mu_inv == NULL);
  assert(md->mu_inv_mean == 1);
}

int main() {

  test_create_maxwell_data();

  return 0;
}