/* THIS FILE WAS AUTOMATICALLY GENERATED.  DO NOT MODIFY! */
/* generated from the file: ./geom.scm */

#ifndef CTL_IO_H
#define CTL_IO_H

#include <ctl-noscheme.h>
#include "data_structures.hpp"

#ifdef __cplusplus
extern "C" {
#endif                          /* __cplusplus */

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

/******* Output variables *******/

extern int num_read_input_vars;
extern int num_write_output_vars;

/******* external-functions *******/

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

#ifdef __cplusplus
}                               /* extern "C" */
#endif                          /* __cplusplus */

#endif                          /* CTL_IO_H */

