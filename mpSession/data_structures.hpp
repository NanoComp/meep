/*
 * Copyright (C) 1998-2014 Massachusetts Institute of Technology and Steven G. Johnson
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA  02111-1307, USA.
 *
 * Steven G. Johnson can be contacted at stevenj@alum.mit.edu.
 */

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- code for manipulating short-term data structures            */
/*- taken from meep/libctl/ctl_io.h and modified                */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

#ifndef DATASTRUCTURES_H
#define DATASTRUCTURES_H

#ifdef __cpluscplus
 extern "C" {
#endif

/******* Type declarations *******/
typedef double (*pml_profile_func)(double u, void *func_data);
typedef struct ctl_pml_profile_data_struct
 { enum { QUADRATIC, USER } type;
   pml_profile_func user_func;
   void *user_data;
 } ctl_pml_profile_data;

typedef struct ctl_pml_struct {
number thickness;
integer direction;
integer side;
number strength;
number R_asymptotic;
number mean_stretch;
ctl_pml_profile_data ppdata;
enum { PML_SELF, ABSORBER } which_subclass;
union {
struct absorber_struct *absorber_data;
} subclass;
} ctl_pml;

typedef struct {
int num_items;
ctl_pml *items;
} ctl_pml_list;

typedef struct ctl_symmetry_struct {
integer direction;
cnumber phase;
enum { SYMMETRY_SELF, MIRROR_SYM, ROTATE4_SYM, ROTATE2_SYM } which_subclass;
union {
struct mirror_sym_struct *mirror_sym_data;
struct rotate4_sym_struct *rotate4_sym_data;
struct rotate2_sym_struct *rotate2_sym_data;
} subclass;
} ctl_symmetry;

typedef struct {
int num_items;
ctl_symmetry *items;
} ctl_symmetry_list;

typedef struct ctl_susceptibility_struct {
vector3 sigma_offdiag;
vector3 sigma_diag;
enum { SUSCEPTIBILITY_SELF, DRUDE_SUSCEPTIBILITY, LORENTZIAN_SUSCEPTIBILITY } which_subclass;
union {
struct ctl_drude_susceptibility_struct *drude_susceptibility_data;
struct ctl_lorentzian_susceptibility_struct *lorentzian_susceptibility_data;
} subclass;
} ctl_susceptibility;

typedef struct ctl_lorentzian_susceptibility_struct {
number frequency;
number gamma;
enum { LORENTZIAN_SUSCEPTIBILITY_SELF, NOISY_LORENTZIAN_SUSCEPTIBILITY } which_subclass;
union {
struct ctl_noisy_lorentzian_susceptibility_struct *noisy_lorentzian_susceptibility_data;
} subclass;
} ctl_lorentzian_susceptibility;

typedef struct ctl_drude_susceptibility_struct {
number frequency;
number gamma;
enum { DRUDE_SUSCEPTIBILITY_SELF, NOISY_DRUDE_SUSCEPTIBILITY } which_subclass;
union {
struct ctl_noisy_drude_susceptibility_struct *noisy_drude_susceptibility_data;
} subclass;
} ctl_drude_susceptibility;

typedef struct ctl_noisy_lorentzian_susceptibility_struct {
number noise_amp;
} ctl_noisy_lorentzian_susceptibility;

typedef struct ctl_noisy_drude_susceptibility_struct {
number noise_amp;
} ctl_noisy_drude_susceptibility;

typedef struct {
int num_items;
ctl_susceptibility *items;
} susceptibility_list;

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
} medium;

typedef struct material_type_struct {
enum { MATERIAL_TYPE_SELF, MATERIAL_FUNCTION, PERFECT_METAL, MEDIUM } which_subclass;
union {
struct material_function_struct *material_function_data;
struct perfect_metal_struct *perfect_metal_data;
struct medium_struct *medium_data;
} subclass;
} material_type;

typedef struct perfect_metal_struct {
} perfect_metal;

typedef material_type (*material_func)(vector3 p);
typedef struct material_function_struct {
material_func func;
} material_function;

typedef struct {
int num_items;
material_type *items;
} material_type_list;

extern boolean perfect_metal_equal(const perfect_metal *o0, const perfect_metal *o);
extern boolean medium_equal(const medium *o0, const medium *o);
extern boolean ctl_noisy_drude_susceptibility_equal(const ctl_noisy_drude_susceptibility *o0, const ctl_noisy_drude_susceptibility *o);
extern boolean ctl_noisy_lorentzian_susceptibility_equal(const ctl_noisy_lorentzian_susceptibility *o0, const ctl_noisy_lorentzian_susceptibility *o);
extern boolean ctl_drude_susceptibility_equal(const ctl_drude_susceptibility *o0, const ctl_drude_susceptibility *o);
extern boolean ctl_lorentzian_susceptibility_equal(const ctl_lorentzian_susceptibility *o0, const ctl_lorentzian_susceptibility *o);
extern boolean ctl_susceptibility_equal(const ctl_susceptibility *o0, const ctl_susceptibility *o);

extern void ctl_susceptibility_copy(const ctl_susceptibility *o0, ctl_susceptibility *o);
extern void ctl_susceptibility_destroy(ctl_susceptibility o);

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- section (D): content from meep/libctl/meep-ctl-const.hpp   -*/
/*-                       and meep/libctl/meep-ctl-const.hpp   -*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
#define CYLINDRICAL -2

/* should be the same as meep::direction enum */
#define X_DIR 0
#define Y_DIR 1
#define Z_DIR 2
#define R_DIR 4
#define PHI_DIR 5

#define CK(ex, msg) \
    (void)((ex) || (meep::abort(msg), 0))

#ifdef __cpluscplus
 } // extern "C" 
#endif

#endif // DATASTRUCTURES_H
