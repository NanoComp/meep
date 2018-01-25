#ifndef PYMPB_H
#define PYMPB_H

#include <vector>

#include "ctlgeom.h"
#include "mpb.h"
#include "mpb/maxwell.h"
#include "meepgeom.hpp"

namespace py_mpb {

// TODO: Temporary matrixio stuff
#if defined(HAVE_HDF5)
/* don't use new HDF5 1.8 API (which isn't even fully documented yet, grrr) */
#  define H5_USE_16_API 1
#  include <hdf5.h>
typedef hid_t matrixio_id_;
/* HDF5 changed this datatype in their interfaces starting in version 1.6.4 */
#  if H5_VERS_MAJOR > 1 \
     || (H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6) \
     || (H5_VERS_MAJOR == 1 && H5_VERS_MINOR == 6 && H5_VERS_RELEASE > 3)
typedef hsize_t start_t;
#  else
typedef hssize_t start_t;
#  endif
#else /* no HDF */
typedef int matrixio_id_; /* dummy */
#endif

typedef struct {
     matrixio_id_ id;
     int parallel;
} matrixio_id;

#define TWOPI 6.2831853071795864769252867665590057683943388


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
  int mesh_size;

  int last_parity;

  // Output variable
  bool negative_epsilon_ok;
  int iterations;

  geometric_object_list geometry;

  geom_box_tree geometry_tree;
  geom_box_tree restricted_tree;

  mpb_real vol;

  mpb_real R[3][3];
  mpb_real G[3][3];

  maxwell_data *mdata;
  maxwell_target_data *mtdata;

  int curfield_band;

  vector3 cur_kvector;
  matrix3x3 Rm;
  matrix3x3 Gm;

  evectmatrix H;
  evectmatrix Hblock;
  evectmatrix muinvH;
  evectmatrix W[MAX_NWORK];

  meep_geom::material_data *default_md;

  std::vector<mpb_real> freqs;

  bool verbose;
  bool deterministic;

  mode_solver(int num_bands, int parity, double resolution, lattice lat, double tolerance,
              meep_geom::material_data *_default_material, geometric_object_list geom,
              bool reset_fields, bool deterministic);
  ~mode_solver();

  void init(int p, bool reset_fields);
  void solve_kpoint(vector3 kpoint);
  bool using_mu();
  void set_parity(int p);
  void set_kpoint_index(int i);
  void get_epsilon();
  void get_epsilon_tensor(int c1, int c2, int imag, int inv);
  void get_material_pt(meep_geom::material_type &material, vector3 p);
  void material_epsmu(meep_geom::material_type material, symmetric_matrix *epsmu,
                      symmetric_matrix *epsmu_inv);
  void randomize_fields();
  void init_epsilon();
  void reset_epsilon();
  void curfield_reset();
  void output_field_to_file(int which_component, char *filename_prefix);
  void output_scalarfield(mpb_real *vals,
                          const int dims[3],
                          const int local_dims[3],
                          const int start[3],
                          matrixio_id file_id,
                          const char *dataname,
                          int last_dim_index,
                          int last_dim_start,
                          int last_dim_size,
                          int first_dim_start,
                          int first_dim_size,
                          int write_start0_special);
  char *fix_fname(const char *fname, const char *prefix, maxwell_data *d, int parity_suffix);
  void load_eigenvectors(char *filename);

  size_t get_field_size();

  std::vector<mpb_real> get_freqs();
  // std::complex<mpb_real>* get_e_field();
  // std::complex<mpb_real>* get_d_field();
  void get_h_field(std::complex<mpb_real> *cdata, int size);

private:
  int kpoint_index;
  scalar_complex *curfield;
  char curfield_type;
};
} // namespace py_mpb
#endif
