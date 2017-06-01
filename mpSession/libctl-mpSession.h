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

/*
 * HR 20170521
 * libctl-mpSession is (for now) a hand-selected hodgepodge
 * of content from (a) libctl itself, (b) the libctl subdirectory
 * of the meep repository.
 * 
 * More specifically, it includes the following bits and pieces:
 * 
 * (A) from the standalone libctl package: 
 *       code for manipulating vector3/matrix3x3 objects
 *
 * (B) from generated C++ code (ctl_io.cpp) in meep/libctl:
 *       code for manipulating short-term data structures:
 *        pml, symmetry, susceptibility, material_type, etc.
 *
 * (C) from native C++ code (geom.cpp) in meep/libctl:
 *       code for manipulating geometric_objects
 *
 * (D) from meep/libctl/meep-ctl-const.hpp: 
 *       constant definitions 
 */
#ifndef LIBCTL_MPSESSION_H
#define LIBCTL_MPSESSION_H

#ifdef __cpluscplus
 extern "C" {
#endif

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- SECTION (A) vector3 / matrix3x3 code from ctl.h in the      */
/*- standalone libctl package                                   */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
typedef int integer;
typedef double number;
typedef struct { number re, im; } cnumber; /* complex numbers! */
typedef short boolean;
typedef char *string;

/* define vector3 as a structure, not an array, so that it can
   be a function return value and so that simple assignment works. */
typedef struct { number x,y,z; } vector3;

  /* similarly for matrix3x3 */
typedef struct { vector3 c0, c1, c2; /* the columns */ } matrix3x3;

/* define complex equivalents: */
typedef struct { cnumber x,y,z; } cvector3;
typedef struct { cvector3 c0, c1, c2; /* the columns */ } cmatrix3x3;

/* vector3 and matrix3x3 utilities: */

extern number vector3_dot(vector3 v1,vector3 v2);
extern number vector3_norm(vector3 v);
extern vector3 vector3_scale(number s, vector3 v);
extern vector3 unit_vector3(vector3 v);
extern vector3 vector3_cross(vector3 v1,vector3 v2);
extern vector3 vector3_plus(vector3 v1,vector3 v2);
extern vector3 vector3_minus(vector3 v1,vector3 v2);
extern int vector3_equal(vector3 v1, vector3 v2);

extern vector3 matrix3x3_vector3_mult(matrix3x3 m, vector3 v);
extern vector3 matrix3x3_transpose_vector3_mult(matrix3x3 m, vector3 v);
extern matrix3x3 matrix3x3_mult(matrix3x3 m1, matrix3x3 m2);
extern matrix3x3 matrix3x3_transpose(matrix3x3 m);
extern number matrix3x3_determinant(matrix3x3 m);
extern matrix3x3 matrix3x3_inverse(matrix3x3 m);
extern int matrix3x3_equal(matrix3x3 m1, matrix3x3 m2);

extern vector3 matrix3x3_row1(matrix3x3 m);
extern vector3 matrix3x3_row2(matrix3x3 m);
extern vector3 matrix3x3_row3(matrix3x3 m);

/* complex number utilities */

extern cnumber make_cnumber(number r, number i);
extern cnumber cnumber_conj(cnumber c);
extern int cnumber_equal(cnumber c1, cnumber c2);
#define cnumber_re(c) ((c).re)
#define cnumber_im(c) ((c).im)

extern vector3 cvector3_re(cvector3 cv);
extern vector3 cvector3_im(cvector3 cv);
extern cvector3 make_cvector3(vector3 vr, vector3 vi);
extern int cvector3_equal(cvector3 v1, cvector3 v2);

extern matrix3x3 cmatrix3x3_re(cmatrix3x3 cm);
extern matrix3x3 cmatrix3x3_im(cmatrix3x3 cm);
extern cmatrix3x3 make_cmatrix3x3(matrix3x3 mr, matrix3x3 mi);
cmatrix3x3 make_hermitian_cmatrix3x3(number m00, number m11, number m22,
                                     cnumber m01, cnumber m02, cnumber m12);
extern int cmatrix3x3_equal(cmatrix3x3 m1, cmatrix3x3 m2);

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- SECTION (B) code for manipulating short-term data structures*/
/*-          taken from meep/libctl/ctl_io.h and modified       */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

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

typedef material_type (*ctl_material_func)(vector3 p);
typedef struct material_function_struct {
ctl_material_func material_func;
} material_function;

typedef struct {
int num_items;
material_type *items;
} material_type_list;

typedef struct geometric_object_struct {
material_type material;
vector3 center;
enum { GEOMETRIC_OBJECT_SELF, BLOCK, SPHERE, CYLINDER, COMPOUND_GEOMETRIC_OBJECT } which_subclass;
union {
struct block_struct *block_data;
struct sphere_struct *sphere_data;
struct cylinder_struct *cylinder_data;
struct compound_geometric_object_struct *compound_geometric_object_data;
} subclass;
} geometric_object;

typedef struct {
int num_items;
geometric_object *items;
} geometric_object_list;

typedef struct compound_geometric_object_struct {
geometric_object_list component_objects;
} compound_geometric_object;

typedef struct cylinder_struct {
vector3 axis;
number radius;
number height;
enum { CYLINDER_SELF, WEDGE, CONE } which_subclass;
union {
struct wedge_struct *wedge_data;
struct cone_struct *cone_data;
} subclass;
} cylinder;

typedef struct cone_struct {
number radius2;
} cone;

typedef struct wedge_struct {
number wedge_angle;
vector3 wedge_start;
vector3 e1;
vector3 e2;
} wedge;

typedef struct sphere_struct {
number radius;
} sphere;

typedef struct block_struct {
vector3 e1;
vector3 e2;
vector3 e3;
vector3 size;
matrix3x3 projection_matrix;
enum { BLOCK_SELF, ELLIPSOID } which_subclass;
union {
struct ellipsoid_struct *ellipsoid_data;
} subclass;
} block;

typedef struct ellipsoid_struct {
vector3 inverse_semi_axes;
} ellipsoid;

typedef struct lattice_struct {
vector3 basis1;
vector3 basis2;
vector3 basis3;
vector3 size;
vector3 basis_size;
vector3 b1;
vector3 b2;
vector3 b3;
matrix3x3 basis;
matrix3x3 metric;
} lattice;

/******* Input variables *******/
extern integer dimensions;
extern material_type default_material;
extern vector3 geometry_center;
extern lattice geometry_lattice;
extern geometric_object_list geometry;
extern boolean ensure_periodicity;

extern vector3 v3_zeroes;
extern vector3 v3_xaxis;
extern vector3 v3_yaxis;
extern vector3 v3_zaxis;
extern vector3 v3_ones;

extern matrix3x3 square_basis(matrix3x3, vector3);
extern number range_overlap_with_object(vector3, vector3, geometric_object, number, integer);
extern void display_geometric_object_info(integer, geometric_object);
extern boolean point_in_periodic_objectp(vector3, geometric_object);
extern vector3 normal_to_object(vector3, geometric_object);
extern boolean point_in_objectp(vector3, geometric_object);
extern void export_external_functions(void);

/******* class copy function prototypes *******/

extern void lattice_copy(const lattice *o0,lattice *o);
extern void ellipsoid_copy(const ellipsoid *o0,ellipsoid *o);
extern void block_copy(const block *o0,block *o);
extern void sphere_copy(const sphere *o0,sphere *o);
extern void wedge_copy(const wedge *o0,wedge *o);
extern void cone_copy(const cone *o0,cone *o);
extern void cylinder_copy(const cylinder *o0,cylinder *o);
extern void compound_geometric_object_copy(const compound_geometric_object *o0,compound_geometric_object *o);
extern void geometric_object_copy(const geometric_object *o0,geometric_object *o);
extern void material_type_copy(const material_type *o0,material_type *o);

/******* class equal function prototypes *******/

extern boolean lattice_equal(const lattice *o0, const lattice *o);
extern boolean ellipsoid_equal(const ellipsoid *o0, const ellipsoid *o);
extern boolean block_equal(const block *o0, const block *o);
extern boolean sphere_equal(const sphere *o0, const sphere *o);
extern boolean wedge_equal(const wedge *o0, const wedge *o);
extern boolean cone_equal(const cone *o0, const cone *o);
extern boolean cylinder_equal(const cylinder *o0, const cylinder *o);
extern boolean compound_geometric_object_equal(const compound_geometric_object *o0, const compound_geometric_object *o);
extern boolean geometric_object_equal(const geometric_object *o0, const geometric_object *o);
extern boolean material_type_equal(const material_type *o0, const material_type *o);

/******* class destruction function prototypes *******/

extern void lattice_destroy(lattice o);
extern void ellipsoid_destroy(ellipsoid o);
extern void block_destroy(block o);
extern void sphere_destroy(sphere o);
extern void wedge_destroy(wedge o);
extern void cone_destroy(cone o);
extern void cylinder_destroy(cylinder o);
extern void compound_geometric_object_destroy(compound_geometric_object o);
extern void geometric_object_destroy(geometric_object o);
extern void material_type_destroy(material_type o);

extern boolean perfect_metal_equal(const perfect_metal *o0, const perfect_metal *o);
extern boolean medium_equal(const medium *o0, const medium *o);
extern boolean ctl_noisy_drude_susceptibility_equal(const ctl_noisy_drude_susceptibility *o0, const ctl_noisy_drude_susceptibility *o);
extern boolean ctl_noisy_lorentzian_susceptibility_equal(const ctl_noisy_lorentzian_susceptibility *o0, const ctl_noisy_lorentzian_susceptibility *o);
extern boolean ctl_drude_susceptibility_equal(const ctl_drude_susceptibility *o0, const ctl_drude_susceptibility *o);
extern boolean ctl_lorentzian_susceptibility_equal(const ctl_lorentzian_susceptibility *o0, const ctl_lorentzian_susceptibility *o);
extern boolean ctl_susceptibility_equal(const ctl_susceptibility *o0, const ctl_susceptibility *o);
extern boolean material_type_equal(const material_type *o0, const material_type *o);

extern void material_type_copy(const material_type *o0,material_type *o);
extern void ctl_susceptibility_copy(const ctl_susceptibility *o0, ctl_susceptibility *o);

extern void ctl_susceptibility_destroy(ctl_susceptibility o);
extern void material_type_destroy(material_type o);


/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- section (C): geometric primitives from libctlgeom ----------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
#define MATERIAL_TYPE material_type
#define GEOMETRIC_OBJECT geometric_object
#define GEOMETRIC_OBJECT_LIST geometric_object_list
#define LATTICE lattice

extern void geom_initialize(void);
extern void geom_fix_object(GEOMETRIC_OBJECT o);
extern void geom_fix_objects(void);
extern void geom_fix_objects0(GEOMETRIC_OBJECT_LIST geometry);
extern void geom_fix_lattice(void);
extern void geom_fix_lattice0(LATTICE *L);
extern void geom_cartesian_lattice(void);
extern void geom_cartesian_lattice0(LATTICE *L);
extern boolean point_in_objectp(vector3 p, GEOMETRIC_OBJECT o);
extern boolean point_in_periodic_objectp(vector3 p, GEOMETRIC_OBJECT o);
extern boolean point_in_fixed_objectp(vector3 p, GEOMETRIC_OBJECT o);
extern boolean point_in_fixed_pobjectp(vector3 p, GEOMETRIC_OBJECT *o);
extern boolean point_in_periodic_fixed_objectp(vector3 p, GEOMETRIC_OBJECT o);
extern vector3 to_geom_object_coords(vector3 p, GEOMETRIC_OBJECT o);
extern vector3 from_geom_object_coords(vector3 p, GEOMETRIC_OBJECT o);
extern vector3 normal_to_object(vector3 p, GEOMETRIC_OBJECT o);
extern vector3 normal_to_fixed_object(vector3 p, GEOMETRIC_OBJECT o);
extern int intersect_line_with_object(vector3 p, vector3 d, GEOMETRIC_OBJECT o,
				      double s[2]);
extern MATERIAL_TYPE material_of_point_inobject(vector3 p, boolean *inobject);
extern MATERIAL_TYPE material_of_point_inobject0(
     GEOMETRIC_OBJECT_LIST geometry, vector3 p, boolean *inobject);
extern MATERIAL_TYPE material_of_point(vector3 p);
extern MATERIAL_TYPE material_of_point0(GEOMETRIC_OBJECT_LIST geometry,
					vector3 p);
GEOMETRIC_OBJECT object_of_point0(GEOMETRIC_OBJECT_LIST geometry, vector3 p,
				   vector3 *shiftby);
GEOMETRIC_OBJECT object_of_point(vector3 p, vector3 *shiftby);
vector3 shift_to_unit_cell(vector3 p);
extern matrix3x3 square_basis(matrix3x3 lattice_basis, vector3 size);

typedef struct {
     vector3 low, high;
} geom_box;

typedef struct {
     geom_box box;
     const GEOMETRIC_OBJECT *o;
     vector3 shiftby;
     int precedence;
} geom_box_object;

typedef struct geom_box_tree_struct {
     geom_box b, b1, b2;
     struct geom_box_tree_struct *t1, *t2;
     int nobjects;
     geom_box_object *objects;
} *geom_box_tree;

extern void destroy_geom_box_tree(geom_box_tree t);
extern geom_box_tree create_geom_box_tree(void);
extern geom_box_tree create_geom_box_tree0(GEOMETRIC_OBJECT_LIST geometry,
					   geom_box b0);
extern geom_box_tree restrict_geom_box_tree(geom_box_tree, const geom_box *);
extern geom_box_tree geom_tree_search(vector3 p, geom_box_tree t, int *oindex);
extern geom_box_tree geom_tree_search_next(vector3 p, geom_box_tree t, int *oindex);
extern MATERIAL_TYPE material_of_point_in_tree_inobject(vector3 p, geom_box_tree t, boolean *inobject);
extern MATERIAL_TYPE material_of_point_in_tree(vector3 p, geom_box_tree t);
extern MATERIAL_TYPE material_of_unshifted_point_in_tree_inobject(vector3 p, geom_box_tree t, boolean *inobject);
const GEOMETRIC_OBJECT *object_of_point_in_tree(vector3 p, geom_box_tree t,
						vector3 *shiftby,
						int *precedence);
extern vector3 to_geom_box_coords(vector3 p, geom_box_object *gbo);
extern void display_geom_box_tree(int indentby, geom_box_tree t);
extern void geom_box_tree_stats(geom_box_tree t, int *depth, int *nobjects);

extern void geom_get_bounding_box(GEOMETRIC_OBJECT o, geom_box *box);
extern number box_overlap_with_object(geom_box b, GEOMETRIC_OBJECT o, number tol, integer maxeval);
extern number ellipsoid_overlap_with_object(geom_box b, GEOMETRIC_OBJECT o, number tol, integer maxeval);
extern number range_overlap_with_object(vector3 low, vector3 high,
					GEOMETRIC_OBJECT o, number tol,
					integer maxeval);

extern vector3 get_grid_size(void);
extern vector3 get_resolution(void);
extern void get_grid_size_n(int *nx, int *ny, int *nz);

GEOMETRIC_OBJECT make_geometric_object(MATERIAL_TYPE material, vector3 center);
GEOMETRIC_OBJECT make_cylinder(MATERIAL_TYPE material, vector3 center,
			       number radius, number height, vector3 axis);
GEOMETRIC_OBJECT make_cone(MATERIAL_TYPE material, vector3 center,
			   number radius, number height, vector3 axis,
			   number radius2);
GEOMETRIC_OBJECT make_sphere(MATERIAL_TYPE material, vector3 center,
			     number radius);
GEOMETRIC_OBJECT make_block(MATERIAL_TYPE material, vector3 center,
			    vector3 e1, vector3 e2, vector3 e3,
			    vector3 size);
GEOMETRIC_OBJECT make_ellipsoid(MATERIAL_TYPE material, vector3 center,
				vector3 e1, vector3 e2, vector3 e3,
				vector3 size);

typedef double (*integrand) (unsigned ndim, const double *x, void *);
number adaptive_integration(integrand f, number *xmin, number *xmax,
			    integer n, void *fdata,
			    number abstol, number reltol, integer maxnfe,
			    number *esterr, integer *errflag);

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

extern int verbose;

#ifdef __cpluscplus
 } // extern "C" 
#endif

#endif // LIBCTL_MPSESSION_H
