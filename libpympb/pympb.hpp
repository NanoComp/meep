/* Copyright (C) 2005-2023 Massachusetts Institute of Technology
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
#ifndef PYMPB_H
#define PYMPB_H

#include <string>
#include <vector>

#include "ctlgeom.h"
#include "mpb.h"
#include "meepgeom.hpp"

namespace py_mpb {

// A global "verbosity" level.
//    0: minimal output
//    1: a little (default)
//    2: a lot
//    3: debugging
//
// This is redeclared here so SWIG will see it.
extern "C" int mpb_verbosity;

#define TWOPI 6.2831853071795864769252867665590057683943388

void map_data(mpb_real *d_in_re, int size_in_re, mpb_real *d_in_im, int size_in_im, int n_in[3],
              mpb_real *d_out_re, int size_out_re, mpb_real *d_out_im, int size_out_im,
              int n_out[3], matrix3x3 coord_map, mpb_real *kvector, bool pick_nearest, bool verbose,
              bool multiply_bloch_phase);

bool with_hermitian_epsilon();

typedef mpb_real (*field_integral_energy_func)(mpb_real, mpb_real, vector3, void *);
typedef cnumber (*field_integral_func)(cvector3, mpb_real, vector3, void *);

struct mode_solver {
  static constexpr int MAX_NWORK = 10;
  static constexpr char epsilon_CURFIELD_TYPE = 'n';
  static constexpr char mu_CURFIELD_TYPE = 'm';
  static constexpr int NUM_FFT_BANDS = 20;

  int num_bands;
  double resolution[3];
  double target_freq;
  double tolerance;
  int mesh_size;
  bool negative_epsilon_ok;
  std::string epsilon_input_file;
  std::string mu_input_file;
  bool force_mu;
  bool use_simple_preconditioner;
  vector3 grid_size;

  int n[3];
  int local_N;
  int N_start;
  int alloc_N;
  int nwork_alloc;

  int eigensolver_nwork;
  int eigensolver_block_size;

  int last_parity;

  // Output variable
  int iterations;
  double eigensolver_flops;

  geometric_object_list geometry_list;
  geom_box_tree geometry_tree;

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

  std::vector<mpb_real> freqs;

  bool verbose;
  bool deterministic;

  mode_solver(int num_bands, double resolution[3], lattice lat, double tolerance, int mesh_size,
              meep_geom::material_data *_default_material, bool deterministic, double target_freq,
              int dims, bool verbose, bool periodicity, double flops, bool negative_epsilon_ok,
              std::string epsilon_input_file, std::string mu_input_file, bool force_mu,
              bool use_simple_preconditioner, vector3 grid_size, int eigensolver_nwork,
              int eigensolver_block_size);

  ~mode_solver();

  // Initializes the mode_solver. Geometric objects are moved out of `geometry`
  // and deallocated upon destruction of the mode_solver.
  void init(int p, bool reset_fields, geometric_object_list *geometry,
            meep_geom::material_data *_default_material);
  void solve_kpoint(vector3 kpoint);
  bool using_mu();
  void set_parity(int p);
  void set_num_bands(int nb);
  int get_kpoint_index();
  void set_kpoint_index(int i);
  void get_epsilon();
  void get_mu();
  void get_epsilon_tensor(int c1, int c2, int imag, int inv);
  void get_material_pt(meep_geom::material_type &material, vector3 p);
  void material_epsmu(meep_geom::material_type material, symmetric_matrix *epsmu,
                      symmetric_matrix *epsmu_inv, bool eps = true);
  int mean_epsilon(symmetric_matrix *meps, symmetric_matrix *meps_inv, mpb_real n[3], mpb_real d1,
                   mpb_real d2, mpb_real d3, mpb_real tol, const mpb_real r[3]);

  void randomize_fields();
  void init_epsilon(geometric_object_list *geometry_in);
  void reset_epsilon(geometric_object_list *geometry);
  bool has_mu(geometric_object_list *geometry);
  bool material_has_mu(void *mt);
  void curfield_reset();

  size_t get_field_size();

  std::vector<mpb_real> get_freqs();
  double get_eigensolver_flops();
  int get_iterations();
  void get_efield(int band);
  void get_dfield(int band);
  void get_hfield(int band);
  void get_bfield(int band);
  void get_efield_from_dfield();

  void get_curfield(double *data, int size);
  void get_curfield_cmplx(std::complex<mpb_real> *cdata, int size);
  void set_curfield(double *data, int size);
  void set_curfield_cmplx(std::complex<mpb_real> *cdata, int size);

  void get_lattice(double data[3][3]);
  void get_eigenvectors(int p_start, int p, std::complex<mpb_real> *cdata, int size);
  std::vector<int> get_eigenvectors_slice_dims(int num_bands);
  void set_eigenvectors(int b_start, std::complex<mpb_real> *cdata, int size);

  std::vector<mpb_real> compute_field_energy();
  double compute_energy_in_objects(geometric_object_list objects);

  char get_curfield_type();
  void set_curfield_type(char t);
  std::string get_parity_string();
  std::vector<int> get_dims();
  void set_grid_size(vector3 gs);
  int get_libctl_dimensions();
  void set_libctl_dimensions(int val);
  bool get_libctl_ensure_periodicity();
  void set_libctl_ensure_periodicity(bool val);
  void set_libctl_geometry_lattice(lattice val);
  std::vector<mpb_real> get_output_k();

  mpb_real get_val(int ix, int iy, int iz, int nx, int ny, int nz, int last_dim_size,
                   mpb_real *data, int stride, int conjugate);
  mpb_real interp_val(vector3 p, int nx, int ny, int nz, int last_dim_size, mpb_real *data,
                      int stride, int conjugate);
  scalar_complex interp_cval(vector3 p, int nx, int ny, int nz, int last_dim_size, mpb_real *data,
                             int stride);
  symmetric_matrix interp_eps_inv(vector3 p);

  mpb_real get_epsilon_point(vector3 p);
  cmatrix3x3 get_epsilon_inverse_tensor_point(vector3 p);
  mpb_real get_energy_point(vector3 p);
  cvector3 get_field_point(vector3 p);
  void get_bloch_field_point_(scalar_complex field[3], vector3 p);
  cvector3 get_bloch_field_point(vector3 p);

  void multiply_bloch_phase(std::complex<double> *cdata = NULL);
  void fix_field_phase();
  void compute_field_divergence();
  std::vector<mpb_real> compute_zparities();
  std::vector<mpb_real> compute_yparities();
  std::vector<mpb_real> compute_group_velocity_component(vector3 d);
  mpb_real compute_1_group_velocity_component(vector3 d, int b);
  vector3 compute_1_group_velocity(int b);
  vector3 compute_1_group_velocity_reciprocal(int b);
  mpb_real compute_energy_in_dielectric(mpb_real eps_low, mpb_real eps_high);

  cnumber compute_field_integral(field_integral_func field_func,
                                 field_integral_energy_func energy_func, void *py_func);
  number compute_energy_integral(field_integral_func field_func,
                                 field_integral_energy_func energy_func, void *py_func);

  vector3 get_dominant_planewave(int band);

  cnumber transformed_overlap(matrix3x3 W, vector3 w);
  cnumber compute_symmetry(int which_band, matrix3x3 W, vector3 w);

private:
  int kpoint_index;
  scalar_complex *curfield;
  char curfield_type;
  bool eps;

  double compute_field_energy_internal(mpb_real comp_sum[6]);
  void clear_geometry_list();
};
} // namespace py_mpb
#endif
