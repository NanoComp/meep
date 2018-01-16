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
  static const char epsilon_CURFIELD_TYPE = 'n';

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

  // TODO: Get from python ?
  int eigensolver_nwork;
  int eigensolver_block_size;

  int last_parity;

  // Output variable
  bool negative_epsilon_ok;

  geometric_object_list geometry;

  geom_box_tree geometry_tree;
  geom_box_tree restricted_tree;

  mpb_real R[3][3];
  mpb_real G[3][3];

  maxwell_data *mdata;
  maxwell_target_data *mtdata;

  scalar_complex *curfield;
  int curfield_band;
  char curfield_type;

  matrix3x3 Gm;

  evectmatrix H;
  evectmatrix Hblock;
  evectmatrix muinvH;
  evectmatrix W[MAX_NWORK];

  meep_geom::material_data *default_md;

  mode_solver(int num_bands, int parity, double resolution, lattice lat, double tolerance,
              meep_geom::material_data *_default_material, geometric_object_list geom, bool reset_fields);
  ~mode_solver();

  bool using_mu();
  void init(int p, bool reset_fields);
  void set_parity(int p);
  void set_kpoint_index(int i);
  void get_epsilon();
  void get_material_pt(meep_geom::material_type &material, vector3 p);
  void material_epsmu(meep_geom::material_type material, symmetric_matrix *epsmu,
                      symmetric_matrix *epsmu_inv);
  void epsilon_file_material(meep_geom::material_data *md, vector3 p);
  void randomize_fields();
  void solve_kpoint(vector3 kpoint);

private:
  int kpoint_index;
};
} // namespace py_mpb
#endif
