#include <cassert>
#include <complex>
#include "pympb.hpp"

void test_mode_solver() {
  static const int NUM_KPOINTS = 16;

  int num_bands = 8;
  int parity = 1;
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

  py_mpb::mode_solver ms(num_bands, parity, resolution, lat, tolerance, mat, g,
                         reset_fields, true);

  ms.get_epsilon();
  ms.output_field_to_file(-1, (char *)"tutorial-");

  for (int i = 0; i < NUM_KPOINTS; ++i) {
    ms.solve_kpoint(k_points[i]);
  }

  int size = ms.mdata->fft_output_size * 3;
  std::complex<mpb_real> *h_field = new std::complex<mpb_real>[size];
  ms.get_h_field(h_field, size);

  for (int i = 0; i < ms.mdata->fft_output_size * 3; ++i) {
    printf("<%f, %fi>\n", real(h_field[i]), imag(h_field[i]));
  }

  free(mat);
}

int main(int argc, char** argv) {

  meep::initialize mpi(argc, argv);
  test_mode_solver();
  return 0;
}
