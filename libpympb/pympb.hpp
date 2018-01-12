#ifndef PYMPB_H
#define PYMPB_H

#include "ctlgeom.h"
#include "mpb.h"
#include "mpb/maxwell.h"
// #include "mpb/scalar.h"
#include "meepgeom.hpp"

namespace py_mpb {

struct mode_solver {
  static const int MAX_NWORK = 10;

  int num_bands;
  int parity;
  double resolution;
  lattice lat;
  double tolerance;

  int n[3];
  int local_N;
  int N_start;
  int alloc_N;
  int nwork_alloc;

  // TODO: Passed in from python
  int eigensolver_nwork;
  int eigensolver_block_size;

  geometric_object_list geometry;

  mpb_real R[3][3];
  mpb_real G[3][3];

  maxwell_data *mdata;

  evectmatrix H;
  evectmatrix Hblock;
  evectmatrix muinvH;
  evectmatrix W[MAX_NWORK];

  mode_solver(int num_bands, int parity, double resolution, lattice lat, double tolerance,
              meep_geom::material_data *_default_material, geometric_object_list geom);
  ~mode_solver();
  bool using_mup();
  void init(int p, bool reset_fields);
  void randomize_fields();
  void solve_kpoint(vector3 kpoint);
};
} // namespace py_mpb
#endif
