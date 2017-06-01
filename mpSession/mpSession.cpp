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
%  GNU General Public License for more details.  %
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include "mpSession.hpp"

using namespace meep;
using namespace meepSession;

material_type *mpSession::vacuum;
material_type *mpSession::air;

/*--------------------------------------------------------------*/
/*- mpSession constructor --------------------------------------*/
/*--------------------------------------------------------------*/
mpSession::mpSession()
{
  // initialize some static class variables used as globals
  mpSession::vacuum=make_dielectric(1.0);
  mpSession::air=mpSession::vacuum;

  // initialize global variables in libctlgeom
  geometry_lattice.size=v3_zeroes;

  // default values for data fields used to construct the_structure
  resolution=10;
  eps_averaging=true;
  subpixel_tol=1.0e-4;
  subpixel_maxeval=100000;
  ensure_periodicity=true;
  material_type_copy(vacuum, &default_mat);
  eps_input_file=0;
  num_chunks=0;
  Courant=0.5;
  global_D_conductivity=0.0;
  global_B_conductivity=0.0;

  // default values for data fields used to construct the_fields
  m=0;
  accurate_fields_near_cylorigin=false;
  special_kz=false;

}

void mpSession::set_dimensions(int dims)
{ 
  dimensions=dims; // global variable in libctlgeom-mpSession
}

void mpSession::set_resolution(int res)
{ 
  resolution=res;  // class field in mpSession
}

// this routine modeled on geom_cartesian_lattice0(lattice *L)
// in geom.cpp
void initLattice(lattice *L, double L1, double L2, double L3)
{
  L->basis1.x = L1; L->basis1.y = 0;  L->basis1.z = 0;
  L->basis2.x = 0;  L->basis2.y = L2; L->basis2.z = 0;
  L->basis3.x = 0;  L->basis3.y = 0;  L->basis3.z = L3;
  L->basis_size.x = L1;
  L->basis_size.y = L2;
  L->basis_size.z = L3;
  geom_fix_lattice0(L);
}

void mpSession::set_geometry_lattice(double L1, double L2, double L3)
{
  initLattice(&(geometry_lattice), L1, L2, L3);
}

void mpSession::add_pml_layer(double thickness, int direction)
{
  int n = pml_layers.num_items++;
  pml_layers.items= (ctl_pml *)realloc(pml_layers.items, (n+1)*sizeof(ctl_pml));
  init_pml( pml_layers.items + n, thickness, direction);
}

int infer_dimensions(double *k)
{
  if (!k) return 3;
  if (k[1]==HUGE_VAL) return 1;
  if (k[2]==HUGE_VAL) return 2;
  return 3;
}
/***************************************************************/

/***************************************************************/
/* based on                                                    */
/*  (define (init-structure . k_)                              */
/* in meep.scm                                                 */
/***************************************************************/
void mpSession::initStructure(double *kPoint)
{
  the_structure
   = make_structure( infer_dimensions(kPoint),
                     geometry_lattice.size, geometry_center,
                     resolution,
                     eps_averaging, subpixel_tol, subpixel_maxeval,
                     (kPoint!=0 && ensure_periodicity),
                     geometry, extra_materials,
                     default_mat, eps_input_file,
                     pml_layers, symmetries, num_chunks, Courant,
                     global_D_conductivity, global_B_conductivity);
}

/***************************************************************/
/* based on                                                    */
/*  (define (init-fields) ...                                  */
/* in meep.scm                                                 */
/***************************************************************/
void mpSession::initFields(double *kPoint)
{
  if (!the_structure)
   initStructure(kPoint);

  the_fields = new fields(the_structure,
                          (dimensions==CYLINDRICAL) ? m : 0.0,
                          (kPoint && special_kz) ? kPoint[2] : 0.0,
                          !accurate_fields_near_cylorigin
                         );

  if (verbose) the_fields->verbose();
}

/***************************************************************/
/* run_until is based on                                       */
/*   (define (run-until cond? . step-funcs) ...                */
/* in libctl/meep.scm                                          */
/***************************************************************/
void mpSession::run_until(double T)
{
  if (!the_fields)
   initFields();
}

//} // namespace meepSession
