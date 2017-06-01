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

/*
  dataStructures.cpp -- reimplementation of short-term
                     -- data structures defined by scheme code in
                     -- the libctl-based meep binary
                     -- 
                     -- i call these "short-term" data structures
                     -- because they are only used to store
                     -- information between the launch of a
                     -- C++ program / python script and the
                     -- start of an actual meep calculation
                     -- (i.e. a call to run_until or similar)
                     -- at which point these data structures 
                     -- are used to initialize low-level data 
                     -  structures in libmeep and then discarded.
                     --  
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpSession.hpp"

using namespace meep;
namespace meepSession 
{

/***************************************************************/
/* some convenient global variables ****************************/
/***************************************************************/
vector3 v3_zeroes = {0.0,0.0,0.0};
vector3 v3_xaxis  = {1.0,0.0,0.0};
vector3 v3_yaxis  = {0.0,1.0,0.0};
vector3 v3_zaxis  = {0.0,0.0,1.0};
vector3 v3_ones   = {1.0,1.0,1.0};

/***************************************************************/
/* constructors ************************************************/
/***************************************************************/

/* (define-class medium material-type ...                       */
material_type *make_dielectric(double epsilon)
 { 
   material_type *mt = (material_type *)malloc(sizeof(material_type));
   mt->which_subclass=material_type_struct::MEDIUM;

   mt->subclass.medium_data = (medium *)malloc(sizeof(medium));

   mt->subclass.medium_data->epsilon_diag=vector3_scale(epsilon,v3_ones);
   mt->subclass.medium_data->mu_diag=v3_ones;
   mt->subclass.medium_data->E_susceptibilities.num_items=0;
   mt->subclass.medium_data->H_susceptibilities.num_items=0;
   mt->subclass.medium_data->E_chi2_diag=v3_zeroes;
   mt->subclass.medium_data->E_chi3_diag=v3_zeroes;
   mt->subclass.medium_data->H_chi2_diag=v3_zeroes;
   mt->subclass.medium_data->H_chi3_diag=v3_zeroes;
   mt->subclass.medium_data->D_conductivity_diag=v3_zeroes;
   mt->subclass.medium_data->B_conductivity_diag=v3_zeroes;

   return mt;
 }

// (define-class perfect-metal material-type)
material_type *make_perfect_metal()
{ 
  material_type *mt=(material_type *)malloc(sizeof(material_type));
  mt->which_subclass=material_type::PERFECT_METAL;

  // looks like perfect_metal is an empty structure, so just
  // set this field to NULL?
  mt->subclass.perfect_metal_data=0; // (perfect_metal *)malloc(sizeof(perfect_metal);

  return mt;
}

/* (define-class pml ... */
void init_pml(ctl_pml *p, number thickness, integer direction)
{ 
  p->thickness=thickness;
  p->direction=direction;
  p->side=ALL_DIRECTIONS;
  p->strength=1.0;
  p->R_asymptotic=1.0e-15;
  p->mean_stretch=1.0;
  //p->pml_profile=_wrap_meep_pml_quadratic_profile;
  p->which_subclass=ctl_pml::ABSORBER;
  p->subclass.absorber_data=0;
}

ctl_pml *make_pml(number thickness, integer direction)
{ 
  ctl_pml *p = (ctl_pml *)malloc(sizeof(ctl_pml));
  init_pml(p, thickness, direction);
  return p;
}

/* (define-class symmetry... */
#define NO_SYMMETRY      0
#define MIRROR_SYMMETRY  1
#define ROTATE4_SYMMETRY 2
#define ROTATE2_SYMMETRY 3
ctl_symmetry *make_symmetry(integer direction,
                               int WhichSymmetry=NO_SYMMETRY)
{
  ctl_symmetry *s = (ctl_symmetry*)malloc(sizeof(*s));
  s->direction=direction;
  s->phase=make_cnumber(1.0, 0.0);
  if (WhichSymmetry==MIRROR_SYMMETRY)
   s->which_subclass = ctl_symmetry::MIRROR_SYM;
  else if (WhichSymmetry==ROTATE4_SYMMETRY)
   s->which_subclass = ctl_symmetry::ROTATE4_SYM;
  else if (WhichSymmetry==ROTATE2_SYMMETRY)
   s->which_subclass = ctl_symmetry::ROTATE2_SYM;
  else
   s->which_subclass = ctl_symmetry::SYMMETRY_SELF;
  return s;
}

ctl_symmetry *make_mirror_symmetry(integer direction)
{ return make_symmetry(direction, MIRROR_SYMMETRY); }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
boolean perfect_metal_equal(const perfect_metal *o0, const perfect_metal *o)
{
;
return 1;
}

boolean medium_equal(const medium *o0, const medium *o)
{
if (!vector3_equal(o->epsilon_diag, o0->epsilon_diag)) return 0;
if (!vector3_equal(o->epsilon_offdiag, o0->epsilon_offdiag)) return 0;
if (!vector3_equal(o->mu_diag, o0->mu_diag)) return 0;
if (!vector3_equal(o->mu_offdiag, o0->mu_offdiag)) return 0;
{
int i_t;
if (o->E_susceptibilities.num_items != o0->E_susceptibilities.num_items) return 0;
for (i_t = 0; i_t < o->E_susceptibilities.num_items; i_t++) {
if (!ctl_susceptibility_equal(&o0->E_susceptibilities.items[i_t], &o->E_susceptibilities.items[i_t])) return 0;
}
}
{
int i_t;
if (o->H_susceptibilities.num_items != o0->H_susceptibilities.num_items) return 0;
for (i_t = 0; i_t < o->H_susceptibilities.num_items; i_t++) {
if (!ctl_susceptibility_equal(&o0->H_susceptibilities.items[i_t], &o->H_susceptibilities.items[i_t])) return 0;
}
}
if (!vector3_equal(o->E_chi2_diag, o0->E_chi2_diag)) return 0;
if (!vector3_equal(o->E_chi3_diag, o0->E_chi3_diag)) return 0;
if (!vector3_equal(o->H_chi2_diag, o0->H_chi2_diag)) return 0;
if (!vector3_equal(o->H_chi3_diag, o0->H_chi3_diag)) return 0;
if (!vector3_equal(o->D_conductivity_diag, o0->D_conductivity_diag)) return 0;
if (!vector3_equal(o->B_conductivity_diag, o0->B_conductivity_diag)) return 0;
;
return 1;
}

boolean ctl_noisy_drude_susceptibility_equal(const ctl_noisy_drude_susceptibility *o0, const ctl_noisy_drude_susceptibility *o)
{
if (o->noise_amp != o0->noise_amp) return 0;
;
return 1;
}

boolean ctl_noisy_lorentzian_susceptibility_equal(const ctl_noisy_lorentzian_susceptibility *o0, const ctl_noisy_lorentzian_susceptibility *o)
{
if (o->noise_amp != o0->noise_amp) return 0;
;
return 1;
}

boolean ctl_drude_susceptibility_equal(const ctl_drude_susceptibility *o0, const ctl_drude_susceptibility *o)
{
if (o->frequency != o0->frequency) return 0;
if (o->gamma != o0->gamma) return 0;
if (o0->which_subclass != o->which_subclass) return 0;
if (o0->which_subclass == ctl_drude_susceptibility::NOISY_DRUDE_SUSCEPTIBILITY) {
if (!ctl_noisy_drude_susceptibility_equal(o0->subclass.noisy_drude_susceptibility_data, o->subclass.noisy_drude_susceptibility_data)) return 0;
}
else ;
return 1;
}

boolean ctl_lorentzian_susceptibility_equal(const ctl_lorentzian_susceptibility *o0, const ctl_lorentzian_susceptibility *o)
{
if (o->frequency != o0->frequency) return 0;
if (o->gamma != o0->gamma) return 0;
if (o0->which_subclass != o->which_subclass) return 0;
if (o0->which_subclass == ctl_lorentzian_susceptibility::NOISY_LORENTZIAN_SUSCEPTIBILITY) {
if (!ctl_noisy_lorentzian_susceptibility_equal(o0->subclass.noisy_lorentzian_susceptibility_data, o->subclass.noisy_lorentzian_susceptibility_data)) return 0;
}
else ;
return 1;
}

boolean ctl_susceptibility_equal(const ctl_susceptibility *o0, const ctl_susceptibility *o)
{
if (!vector3_equal(o->sigma_offdiag, o0->sigma_offdiag)) return 0;
if (!vector3_equal(o->sigma_diag, o0->sigma_diag)) return 0;
if (o0->which_subclass != o->which_subclass) return 0;
if (o0->which_subclass == ctl_susceptibility::DRUDE_SUSCEPTIBILITY) {
if (!ctl_drude_susceptibility_equal(o0->subclass.drude_susceptibility_data, o->subclass.drude_susceptibility_data)) return 0;
}
else if (o0->which_subclass == ctl_susceptibility::LORENTZIAN_SUSCEPTIBILITY) {
if (!ctl_lorentzian_susceptibility_equal(o0->subclass.lorentzian_susceptibility_data, o->subclass.lorentzian_susceptibility_data)) return 0;
}
else ;
return 1;
}

boolean material_type_equal(const material_type *o0, const material_type *o)
{
if (o0->which_subclass != o->which_subclass) return 0;
if (o0->which_subclass == material_type::MATERIAL_FUNCTION) {
//if (!material_function_equal(o0->subclass.material_function_data, o->subclass.material_function_data)) return 0;
}
else if (o0->which_subclass == material_type::PERFECT_METAL) {
if (!perfect_metal_equal(o0->subclass.perfect_metal_data, o->subclass.perfect_metal_data)) return 0;
}
else if (o0->which_subclass == material_type::MEDIUM) {
if (!medium_equal(o0->subclass.medium_data, o->subclass.medium_data)) return 0;
}
else ;
return 1;
}

void material_type_destroy(material_type o)
{
  // FIXME
}

void ctl_susceptibility_copy(const ctl_susceptibility *src, ctl_susceptibility *dest)
{ memcpy(dest,src,sizeof(ctl_susceptibility));
}

void ctl_susceptibility_destroy(ctl_susceptibility o)
{ 
  // FIXME
}
