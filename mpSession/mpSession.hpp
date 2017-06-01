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

/*
 * mpSession.h -- header for libmpSession, a library that sits
 *             -- atop libmeep to provide to C++/python codes
 *             -- the same higher-level functionality that is
 *             -- provided to scheme codes by libctl
 */
#ifndef MPSESSION_H
#define MPSESSION_H

#include <stdio.h>
#include <math.h>

#include "config.h"
#include "meep.hpp"
#include "ctl-noscheme.h"
#include "ctlgeom-types.h"
#include "ctlgeom.h"

#include "dataStructures.hpp"

using namespace meep;

namespace meepSession {

#define ALL_DIRECTIONS -1

/***************************************************************/
/* dataStructures.cpp ******************************************/
/***************************************************************/
typedef amplitude_function std::complex<double> A(const vec &);

material_type *make_dielectric(double epsilon);
void init_pml(ctl_pml *p, double thickness, int direction=ALL_DIRECTIONS);
ctl_pml *make_pml(double thickness);

/***************************************************************/
/* some convenient global variables ****************************/
/***************************************************************/
extern vector3 v3_zeroes;
extern vector3 v3_xaxis;
extern vector3 v3_yaxis;
extern vector3 v3_zaxis;
extern vector3 v3_ones;

/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
/* an mpSession ("meep-python session") stores roughly the     */
/* state that is present globally in a libctl-driven meep      */
/* session.                                                    */
/*                                                             */
/* the names of public methods in mpSession are chosen in      */
/* analogy with the names of corresponding statements one      */
/* might make in a .ctl file: for example,                     */
/*                                                             */
/*  .ctl file statement         |  mpSession equivalent        */
/* ----------------------------------------------------------- */
/*                              |                              */
/* (set! geometry_lattice ... ) |  S.set_geometry_lattice(...) */
/*                              |                              */
/* (set! pml-layers  ... )      |  S.add_pml_layer(...)        */
/*                              |  S.add_pml_layer(...)        */
/*                              |                              */
/* (set! sources                |  S.add_source( source1 )     */
/*       (list                  |  S.add_source( source2 )     */
/*         source1 source2...)) |                              */
/*                              |                              */
/* (run-until T ...)            |  S.run_until(T, ...)         */
/***************************************************************/
class mpSession
 { 
   public:
     mpSession(); 

     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     void set_dimensions(int dims);
     void set_resolution(int res);

     void set_geometry_lattice(double L1, double L2=0.0, double L3=0.0);

     void add_pml_layer(double thickness, int direction=ALL_DIRECTIONS );

     void add_continuous_src(double fcen, double df, component cp,
                             vector3 
     void add_object(char *format, ...);

     //step_func_list step_funcs;
     void run_until(double T);

     // static fields used as global variables
     static material_type *vacuum, *air;

   private:

     // data fields used to construct the_structure
     int dimensions;
     int        resolution;
     bool       eps_averaging;
     double     subpixel_tol;
     int        subpixel_maxeval;
     bool       ensure_periodicity_p;
     material_type_list extra_materials;
     material_type default_mat;
     const char *eps_input_file;
     ctl_pml_list pml_layers;
     ctl_symmetry_list symmetries;
     int num_chunks;
     double Courant;
     double global_D_conductivity;
     double global_B_conductivity;

     // data fields used to construct the_fields
     vector3 k_point;
     int m;
     bool accurate_fields_near_cylorigin;
     bool special_kz;

     geometric_object_list objects;

     // 
     structure  *the_structure;
     fields     *the_fields;

     // private methods
     void initStructure(double *kPoint=0);
     void initFields(double *kPoint=0);

 };
} // namespace meepSession

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
meep::structure *make_structure(int dims, vector3 size, vector3 center,
				double resolution, bool enable_averaging,
				double subpixel_tol, int subpixel_maxeval,
				bool ensure_periodicity_p,
				geometric_object_list geometry,
				material_type_list extra_materials,
				material_type default_mat,
				const char *eps_input_file,
				ctl_pml_list pml_layers,
				ctl_symmetry_list symmetries,
				int num_chunks, double Courant,
				double global_D_conductivity_,
				double global_B_conductivity_);

#endif //MPSESSION_H
