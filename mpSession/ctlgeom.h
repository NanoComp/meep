/* libctl: flexible Guile-based control files for scientific software 
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

#ifndef GEOM_H
#define GEOM_H

#include <ctlgeom-types.h>

/***************************************************************/
/* global variables in libctl-geom, there defined by scheme code*/
/***************************************************************/
extern integer dimensions;
extern material_type default_material;
extern vector3 geometry_center;
extern lattice geometry_lattice;
extern geometric_object_list geometry;
extern boolean ensure_periodicity;
extern int verbose;

#define MATERIAL_TYPE material_type
#define GEOMETRIC_OBJECT geometric_object
#define GEOMETRIC_OBJECT_LIST geometric_object_list
#define LATTICE lattice

#ifdef __cplusplus
extern "C" {
#endif                          /* __cplusplus */

/**************************************************************************/

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


/**************************************************************************/

#ifdef __cplusplus
}                               /* extern "C" */
#endif                          /* __cplusplus */

#endif /* GEOM_H */
