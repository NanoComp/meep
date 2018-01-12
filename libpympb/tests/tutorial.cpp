#include <cassert>
#include "pympb.hpp"

void test_create_maxwell_data() {
  static const int NUM_KPOINTS = 16;

  int num_bands = 8;
  int parity = 0;
  double resolution = 32;
  double tolerance = 1.0e-7;
  lattice lat = {};
  lat.size.x = 1;
  lat.size.y = 1;
  lat.basis_size.x = 1;
  lat.basis_size.y = 1;
  lat.basis_size.z = 1;
  lat.basis1.x = 1;
  lat.basis2.y = 1;
  lat.basis3.z = 1;
  lat.basis.c0 = lat.basis1;
  lat.basis.c1 = lat.basis2;
  lat.basis.c2 = lat.basis3;

  vector3 k_points[NUM_KPOINTS] = {
     {0, 0, 0},
     {0.1, 0, 0},
     {0.2, 0, 0},
     {0.3, 0, 0},
     {0.4, 0, 0},
     {0.5, 0, 0},
     {0.5, 0.1, 0},
     {0.5, 0.2, 0},
     {0.5, 0.3, 0},
     {0.5, 0.4, 0},
     {0.5, 0.5, 0},
     {0.4, 0.4, 0},
     {0.3, 0.3, 0},
     {0.2, 0.2, 0},
     {0.1, 0.1, 0},
     {0, 0, 0},
  };

  meep_geom::material_type mat = meep_geom::make_dielectric(1);
  geometric_object_list g;
  g.num_items = 0;
  g.items = NULL;

  bool reset_fields = true;

  py_mpb::mode_solver ms(num_bands, parity, resolution, lat, tolerance, mat, g);

  ms.init(parity, reset_fields);

  for (int i = 0; i < NUM_KPOINTS; ++i) {
    ms.solve_kpoint(k_points[i]);
  }

  // maxwell_data *md = ms.mdata;

  // assert(md->nx == 32);
  // assert(md->ny == 32);
  // assert(md->nz == 1);
  // assert(md->local_nx == 32);
  // assert(md->local_ny == 32);
  // assert(md->local_x_start == 0);
  // assert(md->local_y_start == 0);
  // assert(md->last_dim == 32);
  // assert(md->last_dim_size == 32);
  // assert(md->other_dims == 32);
  // assert(md->num_bands == 8);
  // assert(md->N == 1024);
  // assert(md->local_N == 1024);
  // assert(md->N_start == 0);
  // assert(md->alloc_N == 1024);
  // assert(md->fft_output_size == 1024);
  // assert(md->max_fft_bands == 8);
  // assert(md->num_fft_bands == 8);
  // assert(md->current_k[0] == 0);
  // assert(md->current_k[1] == 0);
  // assert(md->current_k[2] == 0);
  // assert(md->parity == 0);
  // assert(md->nplans == 0);
  // assert(md->zero_k == 0);
  // assert(md->eps_inv->m00 == 0);
  // assert(md->eps_inv->m01 == 0);
  // assert(md->eps_inv->m02 == 0);
  // assert(md->eps_inv->m11 == 0);
  // assert(md->eps_inv->m12 == 0);
  // assert(md->eps_inv->m22 == 0);
  // assert(md->eps_inv_mean == 1);
  // assert(md->mu_inv == NULL);
  // assert(md->mu_inv_mean == 1);

  free(mat);
}

int main() {

  test_create_maxwell_data();

  return 0;
}