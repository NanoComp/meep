#include "pympb.hpp"

bool test_create_maxwell_data() {
  int band_num = 8;
  vector3 kpoint = {0.0, 0.0, 0.0};
  bool match_frequency = false;
  int parity = 0;
  double resolution = 32;
  double eigensolver_tol = 1.0e-7;

  py_mpb::add_eigenmode_source(band_num, kpoint, match_frequency, parity, resolution, eigensolver_tol);

  // Expected values of maxwell_data
  //  nx = 32,
  //  ny = 32,
  //  nz = 1,
  //  local_nx = 32,
  //  local_ny = 32,
  //  local_x_start = 0,
  //  local_y_start = 0,
  //  last_dim = 32,
  //  last_dim_size = 32,
  //  other_dims = 32,
  //  num_bands = 8,
  //  N = 1024,
  //  local_N = 1024,
  //  N_start = 0,
  //  alloc_N = 1024,
  //  fft_output_size = 1024,
  //  max_fft_bands = 8,
  //  num_fft_bands = 8,
  //  current_k = {0, 0, 0},
  //  parity = 0,
  //  plans = {0x0 <repeats 32 times>},
  //  iplans = {0x0 <repeats 32 times>},
  //  nplans = 0,
  //  plans_howmany = {0 <repeats 32 times>},
  //  plans_stride = {0 <repeats 32 times>},
  //  plans_dist = {0 <repeats 32 times>},
  //  fft_data = 0x7ffff7e25040,
  //  fft_data2 = 0x7ffff7e25040,
  //  zero_k = 0,
  //  k_plus_G = 0x808020,
  //  k_plus_G_normsqr = 0x655280,
  //  eps_inv = 0x7fc010,
  //  eps_inv_mean = 1,
  //  mu_inv = 0x0,
  //  mu_inv_mean = 1

  return false;
}

int main() {

  assert(test_create_maxwell_data());

  return 0;
}