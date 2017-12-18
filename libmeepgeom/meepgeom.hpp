/* Copyright (C) 2005-2015 Massachusetts Institute of Technology
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

#include <meep.hpp>
#include <ctlgeom.h>

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

// tiny floating-point number for effectively zero lengthsh
#define TINY 1e-20

/***************************************************************/
/* these routines create and append absorbing layers to an     */
/* optional list of absorbing layers which is added to the     */
/* material geometry by set_materials_from_geometry.           */
/***************************************************************/
typedef struct absorber {
  double thickness;
  int direction;
  int side;
  double strength;
  double R_asymptotic;
  double mean_stretch;
  meep::pml_profile_func pml_profile;
  void *pml_profile_data;
} absorber;

typedef std::vector<absorber> absorber_list_type;
typedef absorber_list_type *absorber_list;

absorber_list create_absorber_list();
void destroy_absorber_list(absorber_list alist);
void add_absorbing_layer(absorber_list alist,
                         double thickness,
                         int direction=ALL_DIRECTIONS, int side=ALL_SIDES,
                         double strength=1.0, double R_asymptotic=1.0e-15, double mean_stretch=1.0,
                         meep::pml_profile_func func=meep::pml_quadratic_profile, void *func_data=0);

/***************************************************************/
/***************************************************************/
/***************************************************************/
void set_dimensions(int dims);
void set_materials_from_geometry(meep::structure *s,
                                 geometric_object_list g,
                                 bool use_anisotropic_averaging=true,
                                 double tol=DEFAULT_SUBPIXEL_TOL,
                                 int maxeval=DEFAULT_SUBPIXEL_MAXEVAL,
                                 bool ensure_periodicity=false,
                                 bool verbose=false,
                                 material_type _default_material=vacuum,
                                 absorber_list alist=0,
                                 material_type_list extra_materials=material_type_list());

material_type make_dielectric(double epsilon);
material_type make_user_material(user_material_func user_func,
                                 void *user_data);
material_type make_file_material(const char *eps_input_file);

vector3 vec_to_vector3(const meep::vec &pt);
meep::vec vector3_to_vec(const vector3 v3);

}; // namespace meep_geom

#endif // #ifndef MEEP_GEOM_H
