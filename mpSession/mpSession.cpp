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

#include "mpSession.hpp"
#include "meep-ctl-const.hpp"

using namespace meep;
using namespace meepSession;

/*--------------------------------------------------------------*/
/* the next several routines create and initialize instances of */
/* structures defined in ctl-io.h; they replicate the behavior  */
/* of (define-class ...) statements in libctl/xx.scm files, and */
/* calling them from C++/python codes is the equivalent of      */
/* saying (make ...) in .ctl scripts.                           */
/*--------------------------------------------------------------*/

/* (define-class medium material-type ...                       */
ctlio::material_type *make_dielectric(double epsilon)
 { 
   ctlio::material_type *mt = (ctlio::material_type *)malloc(sizeof(ctlio::material_type));
   mt->which_subclass=ctlio::material_type_struct::MEDIUM;

   mt->subclass.medium_data = (ctlio::medium *)malloc(sizeof(ctlio::medium));

   vector3 ones={1.0, 1.0, 1.0};
   vector3 zeroes={0.0, 0.0, 0.0};

   mt->subclass.medium_data->epsilon_diag=vector3_scale(epsilon,ones);
   mt->subclass.medium_data->mu_diag=ones;
   mt->subclass.medium_data->E_susceptibilities.num_items=0;
   mt->subclass.medium_data->H_susceptibilities.num_items=0;
   mt->subclass.medium_data->E_chi2_diag=zeroes;
   mt->subclass.medium_data->E_chi3_diag=zeroes;
   mt->subclass.medium_data->H_chi2_diag=zeroes;
   mt->subclass.medium_data->H_chi3_diag=zeroes;
   mt->subclass.medium_data->D_conductivity_diag=zeroes;
   mt->subclass.medium_data->B_conductivity_diag=zeroes;

   return mt;
 }

// (define vacuum (make dielectric (epsilon 1.0)))
ctlio::material_type *make_vacuum()
 { return make_dielectric(1.0); }

// (define-class perfect-metal material-type)
ctlio::material_type *make_perfect_metal()
{ 
  ctlio::material_type *mt=(ctlio::material_type *)malloc(sizeof(ctlio::material_type));
  mt->which_subclass=ctlio::material_type::PERFECT_METAL;

  // looks like perfect_metal is an empty structure, so just
  // set this field to NULL?
  mt->subclass.perfect_metal_data=0; // (ctlio::perfect_metal *)malloc(sizeof(ctlio::perfect_metal);

  return mt;
}

/* (define-class pml ... */
SCM _wrap_meep_pml_quadratic_profile (SCM s_0, SCM s_1);
void init_pml(ctlio::pml *p, number thickness, integer direction)
{ 
  p->thickness=thickness;
  p->direction=direction;
  p->side=ALL_DIRECTIONS;
  p->strength=1.0;
  p->R_asymptotic=1.0e-15;
  p->mean_stretch=1.0;
  //p->pml_profile=_wrap_meep_pml_quadratic_profile;
  p->which_subclass=ctlio::pml::ABSORBER;
  p->subclass.absorber_data=0;
}

ctlio::pml *make_pml(number thickness, integer direction)
{ 
  ctlio::pml *p = (ctlio::pml *)malloc(sizeof(ctlio::pml));
  init_pml(p, thickness, direction);
  return p;
}

/* (define-class symmetry... */
#define NO_SYMMETRY      0
#define MIRROR_SYMMETRY  1
#define ROTATE4_SYMMETRY 2
#define ROTATE2_SYMMETRY 3
ctlio::symmetry *make_symmetry(integer direction,
                               int WhichSymmetry=NO_SYMMETRY)
{
  ctlio::symmetry *s = (ctlio::symmetry*)malloc(sizeof(*s));
  s->direction=direction;
  s->phase=make_cnumber(1.0, 0.0);
  if (WhichSymmetry==MIRROR_SYMMETRY)
   s->which_subclass = ctlio::symmetry::MIRROR_SYM;
  else if (WhichSymmetry==ROTATE4_SYMMETRY)
   s->which_subclass = ctlio::symmetry::ROTATE4_SYM;
  else if (WhichSymmetry==ROTATE2_SYMMETRY)
   s->which_subclass = ctlio::symmetry::ROTATE2_SYM;
  else
   s->which_subclass = ctlio::symmetry::SYMMETRY_SELF;
  return s;
}

ctlio::symmetry *make_mirror_symmetry(integer direction)
{ return make_symmetry(direction, MIRROR_SYMMETRY); }

/*--------------------------------------------------------------*/
/*- mpSession class methods ------------------------------------*/
/*--------------------------------------------------------------*/
mpSession::mpSession()
{
  // initialize global variables in libmeep-ctl
  epsilon_input_file=0;
  dimensions=3;
  ctlio::material_type *vacuum=make_vacuum();
  ctlio::material_type_copy(vacuum, &default_material);
  geometry_center.x=geometry_center.y=geometry_center.z=0.0;
  set_geometry_lattice(1.0, 1.0, 1.0);
  geometry.num_items=0;
  ensure_periodicity=true;

  // default values for data fields used to construct the_structure
  resolution=0.0;
  eps_averaging=true;
  subpixel_tol=1.0e-4;
  subpixel_maxeval=100000;
  ensure_periodicity=true;
  ctlio::material_type_copy(&default_material, &default_mat);
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
  initLattice(&(ctlio::geometry_lattice), L1, L2, L3);
}

void mpSession::add_pml_layer(double thickness, int direction)
{
  int n = pml_layers.num_items++;
  pml_layers.items= (ctlio::pml *)realloc(pml_layers.items, (n+1)*sizeof(ctlio::pml));
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
                     default_material, epsilon_input_file,
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
