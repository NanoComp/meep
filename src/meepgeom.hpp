/* Copyright (C) 2005-2019 Massachusetts Institute of Technology
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

#include "meep.hpp"
#include "material_data.hpp"

namespace meep_geom {

#ifndef cdouble
typedef std::complex<double> cdouble;
#endif

// constants from meep-ctl-const.hpp
#define CYLINDRICAL -2

/* should be the same as meep::direction enum */
#define X_DIR 0
#define Y_DIR 1
#define Z_DIR 2
#define R_DIR 4
#define PHI_DIR 5

// constant used in meep.scm
#define ALL_SIDES -1
#define ALL_DIRECTIONS -1

// large (but not strictly inf!) floating-point number for
// effectively infinite lengths
#define ENORMOUS 1e20

// tiny floating-point number for effectively zero lengths
#define TINY 1e-20

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
  static bool split_chunks_evenly;

  static bool has_non_medium_material();
  static void init_libctl(meep_geom::material_type default_mat, bool ensure_per,
                          meep::grid_volume *gv, vector3 cell_size, vector3 cell_center,
                          geometric_object_list *geom_list);

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
                                  bool anisotropic_pixels_already_added=false);
  void compute_stats();
  void count_anisotropic_pixels(medium_struct *med, size_t pixels);
  void count_nonlinear_pixels(medium_struct *med, size_t pixels);
  void count_susceptibility_pixels(medium_struct *med, size_t pixels);
  void count_nonzero_conductivity_pixels(medium_struct *med, size_t pixels);
  void compute_dft_stats();
  void compute_pml_stats();
  void compute_absorber_stats();
};

fragment_stats
compute_fragment_stats(geometric_object_list geom, meep::grid_volume *gv, vector3 cell_size,
                       vector3 cell_center, material_type default_mat,
                       std::vector<dft_data> dft_data_list, std::vector<meep::volume> pml_1d_vols,
                       std::vector<meep::volume> pml_2d_vols, std::vector<meep::volume> pml_3d_vols,
                       std::vector<meep::volume> absorber_vols, double tol, int maxeval,
                       bool ensure_per);

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

void set_dimensions(int dims);
void set_materials_from_geometry(meep::structure *s, geometric_object_list g,
                                 vector3 center = make_vector3(),
                                 bool use_anisotropic_averaging = true,
                                 double tol = DEFAULT_SUBPIXEL_TOL,
                                 int maxeval = DEFAULT_SUBPIXEL_MAXEVAL,
                                 bool ensure_periodicity = false, bool verbose = false,
                                 material_type _default_material = vacuum, absorber_list alist = 0,
                                 material_type_list extra_materials = material_type_list());

material_type make_dielectric(double epsilon);
material_type make_user_material(user_material_func user_func, void *user_data, bool do_averaging);
material_type make_file_material(const char *eps_input_file);

vector3 vec_to_vector3(const meep::vec &pt);
meep::vec vector3_to_vec(const vector3 v3);

void epsilon_file_material(material_data *md, vector3 p);
bool susceptibility_equal(const susceptibility &s1, const susceptibility &s2);
bool susceptibility_list_equal(const susceptibility_list &s1, const susceptibility_list &s2);
bool medium_struct_equal(const medium_struct *m1, const medium_struct *m2);
void material_gc(material_type m);
bool material_type_equal(const material_type m1, const material_type m2);
bool is_variable(material_type mt);
bool is_variable(void *md);
bool is_file(material_type md);
bool is_file(void *md);
bool is_medium(material_type md, medium_struct **m);
bool is_medium(void *md, medium_struct **m);
bool is_metal(meep::field_type ft, const material_type *material);
void check_offdiag(medium_struct *m);
geom_box gv2box(const meep::volume &v);

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
