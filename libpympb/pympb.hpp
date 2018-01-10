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
  bool match_frequency;
  int parity;
  double resolution;
  lattice lat;
  double tolerance;

  int n[3];
  int local_N;
  int N_start;
  int alloc_N;

  maxwell_data *mdata;

  mode_solver(int num_bands, bool match_frequency, int parity, double resolution,
              lattice lat, double tolerance);
  ~mode_solver();
  void init(int p, bool reset_fields);
  // void add_eigenmode_source(int band_num, const vector3 &kpoint, bool match_frequency, int parity,
  //                           double resolution, vector3 lattice_size, double tolerance);
};
} // namespace py_mpb
#endif
