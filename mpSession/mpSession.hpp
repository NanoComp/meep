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
 * mpSession.h -- header for libmpSession.
 */
#ifndef MPSESSION_H
#define MPSESSION_H

#include <stdio.h>
#include <math.h>

#include "config.h"
#include "meep.hpp"
#include "ctl-io.h"

using namespace ctlio;
using namespace meep;

namespace meepSession {

#define ALL_DIRECTIONS -1

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
/* (set! pml-layers  ... )      |  S.set_pml_layers(...)       */
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

     void set_geometry_lattice(double L1, double L2=0.0, double L3=0.0);
     void add_pml_layer(double thickness, int direction=ALL_DIRECTIONS );

     //step_func_list step_funcs;
     void run_until(double T);

   private:

     // data fields used to construct the_structure
     vector3    k_point;
     double     resolution;
     bool       eps_averaging;
     double     subpixel_tol;
     int        subpixel_maxeval;
     bool       ensure_periodicity_p;
     ctlio::material_type_list extra_materials;
     ctlio::material_type default_mat;
     const char *eps_input_file;
     ctlio::pml_list pml_layers;
     ctlio::symmetry_list symmetries;
     int num_chunks;
     double Courant;
     double global_D_conductivity;
     double global_B_conductivity;

     // data fields used to construct the_fields
     int m; 
     bool accurate_fields_near_cylorigin;
     bool special_kz;

     // 
     structure  *the_structure;
     fields     *the_fields;

     // private methods
     void initStructure(double *kPoint=0);
     void initFields(double *kPoint=0);

 };

/*--------------------------------------------------------------*/
/*- standalone constructors for data structures initialized     */
/*- in meep.scm                                                 */
/*--------------------------------------------------------------*/
ctlio::pml *make_pml(number thickness);

} // namespace meepSession

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/*--------------------------------------------------------------*/
/*- function prototypes and global variables in libmeep-ctl    -*/
/*- that are not exposed in any of the libctl headers          -*/
/*--------------------------------------------------------------*/
namespace ctlio {
extern lattice geometry_lattice;
}

extern "C" {
void geom_fix_lattice0(lattice *L);
}

meep::structure *make_structure(int dims, vector3 size, vector3 center,
				double resolution, bool enable_averaging,
				double subpixel_tol, int subpixel_maxeval,
				bool ensure_periodicity_p,
				geometric_object_list geometry,
				material_type_list extra_materials,
				material_type default_mat,
				const char *eps_input_file,
				pml_list pml_layers,
				symmetry_list symmetries,
				int num_chunks, double Courant,
				double global_D_conductivity_,
				double global_B_conductivity_);

#endif //MPSESSION_H
