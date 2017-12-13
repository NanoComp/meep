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
  bool is_self;
} susceptibility;

typedef struct {
  int num_items;
  susceptibility *items;
} susceptibility_list;

// prototype for user-supplied material function
// TODO: this prototype only allows user-defined scalar dielectric permeabilities
//       (which seems to be all that is possible via the libctl interface as well).
//       extend to allow more complicated user-specified materials?
//
typedef std::complex<double> (*user_material_func)(vector3 x, void *user_data);

typedef struct medium_struct {
  vector3 epsilon_diag;
  vector3 epsilon_offdiag;
  vector3 mu_diag;
  vector3 mu_offdiag;
  susceptibility_list E_susceptibilities;
  susceptibility_list H_susceptibilities;
  vector3 E_chi2_diag;
  vector3 E_chi3_diag;
  vector3 H_chi2_diag;
  vector3 H_chi3_diag;
  vector3 D_conductivity_diag;
  vector3 B_conductivity_diag;
} medium_struct;

typedef struct material_data_struct
 {
   enum { MATERIAL_TYPE_SELF, MATERIAL_FUNCTION, PERFECT_METAL, MEDIUM } which_subclass;

   // these fields used only if which_subclass==MATERIAL_FUNCTION
   user_material_func user_func;
   void *             user_data;

   // used only if which_subclass==MEDIUM
   medium_struct *medium;

 } material_data;

typedef material_data *material_type;

struct material_type_list {
  material_type *items;
  int num_items;

  material_type_list(): items(NULL), num_items(0) {}
};

// global variables and exported functions
extern material_type vacuum;
material_type make_dielectric(double epsilon);

void read_epsilon_file(const char *eps_input_file);


}; // namespace meep_geom

#endif // #ifndef MATERIAL_DATA_H
