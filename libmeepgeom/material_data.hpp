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
#ifndef MATERIAL_DATA_H
#define MATERIAL_DATA_H

/***************************************************************/
/* the void *data field in geometric_object_struct points to   */
/* a material_data structure, defined below.                   */
/***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <meep.hpp>
#include <ctlgeom.h>

#include "meep.hpp"

namespace meep_geom {

/* FIXME: we don't have to emulate the Scheme/libctl code here, which was
  limited to C functionality.  These types, especially the material types,
  should be proper C++ classes */

typedef struct susceptibility_struct {
  vector3 sigma_offdiag;
  vector3 sigma_diag;
  double frequency;
  double gamma;
  double noise_amp;
  bool drude;
  bool is_file;
} susceptibility;

struct susceptibility_list {
  int num_items;
  susceptibility *items;

  susceptibility_list(): num_items(0), items(NULL) {}
};

struct medium_struct {
  vector3 epsilon_diag;
  cvector3 epsilon_offdiag;
  vector3 mu_diag;
  cvector3 mu_offdiag;
  susceptibility_list E_susceptibilities;
  susceptibility_list H_susceptibilities;
  vector3 E_chi2_diag;
  vector3 E_chi3_diag;
  vector3 H_chi2_diag;
  vector3 H_chi3_diag;
  vector3 D_conductivity_diag;
  vector3 B_conductivity_diag;

  medium_struct(double epsilon=1): E_susceptibilities(), H_susceptibilities() {
    epsilon_diag.x = epsilon;
    epsilon_diag.y = epsilon;
    epsilon_diag.z = epsilon;

    mu_diag.x = 1;
    mu_diag.y = 1;
    mu_diag.z = 1;

    epsilon_offdiag.x.re = 0;
    epsilon_offdiag.x.im = 0;
    epsilon_offdiag.y.re = 0;
    epsilon_offdiag.y.im = 0;
    epsilon_offdiag.z.re = 0;
    epsilon_offdiag.z.im = 0;

    mu_offdiag.x.re = 0;
    mu_offdiag.x.im = 0;
    mu_offdiag.y.re = 0;
    mu_offdiag.y.im = 0;
    mu_offdiag.z.re = 0;
    mu_offdiag.z.im = 0;

    E_chi2_diag.x = 0;
    E_chi2_diag.y = 0;
    E_chi2_diag.z = 0;

    E_chi3_diag.x = 0;
    E_chi3_diag.y = 0;
    E_chi3_diag.z = 0;

    H_chi2_diag.x = 0;
    H_chi2_diag.y = 0;
    H_chi2_diag.z = 0;

    H_chi3_diag.x = 0;
    H_chi3_diag.y = 0;
    H_chi3_diag.z = 0;

    D_conductivity_diag.x = 0;
    D_conductivity_diag.y = 0;
    D_conductivity_diag.z = 0;

    B_conductivity_diag.x = 0;
    B_conductivity_diag.y = 0;
    B_conductivity_diag.z = 0;
  }
};

// prototype for user-defined material function,
// which should fill in medium as appropriate to
// describe the material properties at point x
typedef void (*user_material_func)(vector3 x, void *user_data,
                                   medium_struct *medium);

// the various types of materials are as follows:
//  MEDIUM:        material properties independent of position. In
//                 this case the 'medium' field below is
//                 initialized once and doesn't change.
//  MATERIAL_FILE: material properties position-dependent, described
//                 by user-supplied data file. In this case the
//                 'medium' field is filled in appropriately at
//                 each evaluation point by interpolating file data.
//  MATERIAL_USER: material properties position-dependent, described
//                 by user-supplied function. In this case the
//                 'medium' field is filled in appropriately at
//                 each evaluation point by calling the user's
//                 routine.
//  PERFECT_METAL: the 'medium' field is never referenced in this case.
struct material_data
 {
   enum { MEDIUM,
          MATERIAL_FILE,      // formerly MATERIAL_TYPE_SELF
          MATERIAL_USER,      // formerly MATERIAL_FUNCTION
          PERFECT_METAL
        } which_subclass;

   // this field is used for all material types except PERFECT_METAL
   medium_struct medium;

   // these fields used only if which_subclass==MATERIAL_USER
   user_material_func user_func;
   void *             user_data;

   // these fields used only if which_subclass==MATERIAL_FILE
   meep::realnum *epsilon_data;
   int epsilon_dims[3];

   material_data(): which_subclass(MEDIUM), medium(), user_data(NULL), epsilon_data(NULL) {
     epsilon_dims[0] = 0;
     epsilon_dims[1] = 0;
     epsilon_dims[2] = 0;
   }
};

typedef material_data *material_type;

struct material_type_list {
  material_type *items;
  int num_items;

  material_type_list(): items(NULL), num_items(0) {}
};

// global variables
extern material_type vacuum;

// exported functions for creating particular material types
material_type make_dielectric(double epsilon);
material_type make_user_material(user_material_func user_func,
                                 void *user_data);
material_type make_file_material(char *epsilon_input_file);
void epsilon_file_material(vector3 p);
void read_epsilon_file(const char *eps_input_file);




}; // namespace meep_geom

#endif // #ifndef MATERIAL_DATA_H
