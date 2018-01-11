#ifndef PYMPB_H
#define PYMPB_H

#include "ctlgeom.h"
#include "mpb.h"
#include "mpb/maxwell.h"
// #include "mpb/scalar.h"
#include "meepgeom.hpp"

namespace py_mpb {

struct mode_solver {
  int num_bands;
  int parity;
  double resolution;
  lattice lat;
  double tolerance;

  int n[3];
  int local_N;
  int N_start;
  int alloc_N;

  geometric_object_list geometry;

  mpb_real R[3][3];
  mpb_real G[3][3];

  maxwell_data *mdata;

  mode_solver(int num_bands, int parity, double resolution, lattice lat, double tolerance,
              meep_geom::material_type _default_material, geometric_object_list geom);
  ~mode_solver();
  void init(int p, bool reset_fields);
  void solve_kpoint(vector3 kpoint);
};
} // namespace py_mpb
#endif
