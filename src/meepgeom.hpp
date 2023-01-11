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
#ifndef MEEP_GEOM_H
#define MEEP_GEOM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "material_data.hpp"

namespace meep_geom {

// constants from meep-ctl-const.hpp
const int CYLINDRICAL = -2;

/* should be the same as meep::direction enum */
const int X_DIR = 0;
const int Y_DIR = 1;
const int Z_DIR = 2;
const int R_DIR = 4;
const int PHI_DIR = 5;

// constant used in meep.scm
const int ALL_SIDES = -1;
const int ALL_DIRECTIONS = -1;

// FIXME: we should really purge these ancient hacks and just use Inf and 0.0.
// large (but not strictly inf!) floating-point number for
// effectively infinite lengths
const double ENORMOUS = 1e20;
// tiny floating-point number for effectively zero lengths
const double TINY = 1e-20;

struct dft_data {
  int num_freqs;
  int num_components;
  std::vector<meep::volume> vols;

  dft_data(int freqs, int components, std::vector<meep::volume> volumes);
};

struct fragment_stats {
  static double tol;
  static int maxeval;
  static int resolution;
  static meep::ndim dims;
  static geometric_object_list geom;
  static std::vector<dft_data> dft_data_list;
  static std::vector<meep::volume> pml_1d_vols;
  static std::vector<meep::volume> pml_2d_vols;
  static std::vector<meep::volume> pml_3d_vols;
  static std::vector<meep::volume> absorber_vols;
  static material_type_list extra_materials;
  static bool split_chunks_evenly;
  static bool eps_averaging;

  static bool has_non_medium_material();

  size_t num_anisotropic_eps_pixels;
  size_t num_anisotropic_mu_pixels;
  size_t num_nonlinear_pixels;
  size_t num_susceptibility_pixels;
  size_t num_nonzero_conductivity_pixels;

  // Pixels in single PML regions
  size_t num_1d_pml_pixels;

  // Pixels where 2 PML regions overlap
  size_t num_2d_pml_pixels;

  // Pixels where 3 PML regions overlap
  size_t num_3d_pml_pixels;
  size_t num_dft_pixels;
  size_t num_pixels_in_box;
  geom_box box;

  fragment_stats() {}
  fragment_stats(geom_box &bx);

  void compute();
  double cost() const;
  void print_stats() const;

private:
  void update_stats_from_material(material_type mat, size_t pixels,
                                  bool anisotropic_pixels_already_added = false);
  void compute_stats();
  bool count_anisotropic_pixels(medium_struct *med, size_t pixels);
  bool count_nonlinear_pixels(medium_struct *med, size_t pixels);
  bool count_susceptibility_pixels(medium_struct *med, size_t pixels);
  bool count_nonzero_conductivity_pixels(medium_struct *med, size_t pixels);
  void compute_dft_stats();
  void compute_pml_stats();
  void compute_absorber_stats();
};

fragment_stats
compute_fragment_stats(geometric_object_list geom, meep::grid_volume *gv, vector3 cell_size,
                       vector3 cell_center, material_type default_mat,
                       std::vector<dft_data> dft_data_list, std::vector<meep::volume> pml_1d_vols,
                       std::vector<meep::volume> pml_2d_vols, std::vector<meep::volume> pml_3d_vols,
                       std::vector<meep::volume> absorber_vols, material_type_list extra_materials,
                       double tol, int maxeval, bool ensure_per, bool eps_averaging);

/***************************************************************/
/* these routines create and append absorbing layers to an     */
/* optional list of absorbing layers which is added to the     */
/* material geometry by set_materials_from_geometry.           */
/***************************************************************/
typedef struct absorber {
  double thickness;
  int direction;
  int side;
  double R_asymptotic;
  double mean_stretch;
  meep::pml_profile_func pml_profile;
  void *pml_profile_data;
} absorber;

typedef std::vector<absorber> absorber_list_type;
typedef absorber_list_type *absorber_list;

absorber_list create_absorber_list();
void destroy_absorber_list(absorber_list alist);
void add_absorbing_layer(absorber_list alist, double thickness, int direction = ALL_DIRECTIONS,
                         int side = ALL_SIDES, double R_asymptotic = 1.0e-15,
                         double mean_stretch = 1.0,
                         meep::pml_profile_func func = meep::pml_quadratic_profile,
                         void *func_data = 0);

/***************************************************************/
/***************************************************************/
/***************************************************************/
inline vector3 make_vector3(double x = 0.0, double y = 0.0, double z = 0.0) {
  vector3 v = {x, y, z};
  return v;
}

typedef struct {
  double m00, m01, m02, m11, m12, m22;
} symm_matrix;

struct pol {
  meep_geom::susceptibility user_s;
  struct pol *next;
};

// structure to hold a conductivity profile (for scalar absorbing layers)
struct cond_profile {
  double L;     // thickness
  int N;        // number of points prof[n] from 0..N corresponding to 0..L
  double *prof; // (NULL if none)
};

class geom_epsilon : public meep::material_function {

public:
  double u_p = 0;
  geom_box_tree geometry_tree;
  geom_box_tree restricted_tree;
  geometric_object_list geometry;
  cond_profile cond[5][2]; // [direction][side]
  double tol = DEFAULT_SUBPIXEL_TOL;
  int maxeval = DEFAULT_SUBPIXEL_MAXEVAL;

  geom_epsilon(geometric_object_list g, material_type_list mlist, const meep::volume &v);
  geom_epsilon(const geom_epsilon &geps1); // copy constructor
  virtual ~geom_epsilon();

  virtual void set_cond_profile(meep::direction, meep::boundary_side, double L, double dx,
                                double (*prof)(int, double *, void *), void *, double R);

  virtual void set_volume(const meep::volume &v);
  virtual void unset_volume(void);

  bool has_chi(meep::component c, int p);
  virtual bool has_chi3(meep::component c);
  virtual bool has_chi2(meep::component c);

  double chi(meep::component c, const meep::vec &r, int p);
  virtual double chi3(meep::component c, const meep::vec &r);
  virtual double chi2(meep::component c, const meep::vec &r);

  virtual bool has_mu();

  virtual bool has_conductivity(meep::component c);
  virtual double conductivity(meep::component c, const meep::vec &r);

  virtual double chi1p1(meep::field_type ft, const meep::vec &r);
  virtual void eff_chi1inv_row(meep::component c, double chi1inv_row[3], const meep::volume &v,
                               double tol, int maxeval);

  void eff_chi1inv_matrix(meep::component c, symm_matrix *chi1inv_matrix, const meep::volume &v,
                          double tol, int maxeval, bool &fallback);

  void fallback_chi1inv_row(meep::component c, double chi1inv_row[3], const meep::volume &v,
                            double tol, int maxeval);

  virtual void sigma_row(meep::component c, double sigrow[3], const meep::vec &r);
  void add_susceptibilities(meep::structure *s);
  void add_susceptibilities(meep::field_type ft, meep::structure *s);

  void get_material_pt(material_type &material, const meep::vec &r);

private:
  material_type_list extra_materials;
  pol *current_pol;
};

void set_dimensions(int dims);
geom_epsilon *make_geom_epsilon(meep::structure *s, geometric_object_list *g,
                                vector3 center = make_vector3(), bool ensure_periodicity = false,
                                material_type _default_material = vacuum,
                                material_type_list extra_materials = material_type_list());
//, geometric_object_list g, material_type_list extra_materials
void set_materials_from_geometry(meep::structure *s, geometric_object_list g,
                                 vector3 center = make_vector3(),
                                 bool use_anisotropic_averaging = true,
                                 double tol = DEFAULT_SUBPIXEL_TOL,
                                 int maxeval = DEFAULT_SUBPIXEL_MAXEVAL,
                                 bool ensure_periodicity = false,
                                 material_type _default_material = vacuum, absorber_list alist = 0,
                                 material_type_list extra_materials = material_type_list());
void set_materials_from_geom_epsilon(meep::structure *s, geom_epsilon *geps,
                                     bool use_anisotropic_averaging = true,
                                     double tol = DEFAULT_SUBPIXEL_TOL,
                                     int maxeval = DEFAULT_SUBPIXEL_MAXEVAL,
                                     absorber_list alist = 0);

material_type make_dielectric(double epsilon);
material_type make_user_material(user_material_func user_func, void *user_data, bool do_averaging);
material_type make_file_material(const char *eps_input_file);
material_type make_material_grid(bool do_averaging, double beta, double eta, double damping);

vector3 vec_to_vector3(const meep::vec &pt);
meep::vec vector3_to_vec(const vector3 v3);

void set_default_material(material_type _default_material);
void unset_default_material(void);
void epsilon_material_grid(material_data *md, double u);
void epsilon_file_material(material_data *md, vector3 p);
bool susceptibility_equal(const susceptibility &s1, const susceptibility &s2);
bool susceptibility_list_equal(const susceptibility_list &s1, const susceptibility_list &s2);
bool medium_struct_equal(const medium_struct *m1, const medium_struct *m2);
void material_gc(material_type m);
void material_free(material_type m);
bool material_type_equal(const material_type m1, const material_type m2);
bool is_material_grid(material_type mt);
bool is_material_grid(void *md);
bool is_variable(material_type mt);
bool is_variable(void *md);
bool is_file(material_type md);
bool is_file(void *md);
bool is_medium(material_type md, medium_struct **m);
bool is_medium(void *md, medium_struct **m);
bool is_metal(meep::field_type ft, const material_type *material);
geom_box gv2box(const meep::volume &v);
void get_epsilon_grid(geometric_object_list gobj_list, material_type_list mlist,
                      material_type _default_material, bool _ensure_periodicity,
                      meep::grid_volume gv, vector3 cell_size, vector3 cell_center, int nx,
                      const double *x, int ny, const double *y, int nz, const double *z,
                      std::complex<double> *grid_vals, double frequency = 0);
void init_libctl(material_type default_mat, bool ensure_per, meep::grid_volume *gv,
                 vector3 cell_size, vector3 cell_center, geometric_object_list *geom_list);

/***************************************************************/
// material grid functions
/***************************************************************/
void update_weights(material_type matgrid, double *weights);
meep::vec matgrid_grad(vector3 p, geom_box_tree tp, int oi, material_data *md);
meep::vec material_grid_grad(vector3 p, material_data *md, const geometric_object *o);
double matgrid_val(vector3 p, geom_box_tree tp, int oi, material_data *md);
double material_grid_val(vector3 p, material_data *md);
geom_box_tree calculate_tree(const meep::volume &v, geometric_object_list g);
void material_grids_addgradient(double *v, size_t ng, size_t nf,
                                std::vector<meep::dft_fields *> fields_a,
                                std::vector<meep::dft_fields *> fields_f, double *frequencies,
                                double scalegrad, meep::grid_volume &gv, geom_epsilon *geps,
                                double du = 1e-6);

/***************************************************************/
/* routines in GDSIIgeom.cc ************************************/
/***************************************************************/
bool with_libGDSII();
meep::grid_volume set_geometry_from_GDSII(double resolution, const char *GDSIIFile,
                                          const char *Text, int Layer = -1, double zsize = 0.0);
meep::grid_volume set_geometry_from_GDSII(double resolution, const char *GDSIIFile, int Layer,
                                          double zsize = 0.0);
geometric_object_list get_GDSII_prisms(material_type material, const char *GDSIIFile,
                                       int Layer = -1, double zmin = 0.0, double zmax = 0.0);
geometric_object get_GDSII_prism(material_type material, const char *GDSIIFile, const char *Text,
                                 int Layer = -1, double zmin = 0.0, double zmax = 0.0);
geometric_object get_GDSII_prism(material_type material, const char *GDSIIFile, int Layer,
                                 double zmin = 0.0, double zmax = 0.0);
meep::volume get_GDSII_volume(const char *GDSIIFile, const char *Text, int Layer = -1,
                              double zmin = 0.0, double zmax = 0.0);
meep::volume get_GDSII_volume(const char *GDSIIFile, int Layer, double zmin = 0.0,
                              double zmax = 0.0);
std::vector<int> get_GDSII_layers(const char *GDSIIFile);

}; // namespace meep_geom

#endif // #ifndef MEEP_GEOM_H
