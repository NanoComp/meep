/* Copyright (C) 2005-2017 Massachusetts Institute of Technology  
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software Foundation,
 *  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

%module meep

%{
#define SWIG_FILE_WITH_INIT

#include <string>

#include "meep/vec.hpp"
#include "meep.hpp"
#include "ctl-math.h"
#include "ctlgeom.h"
#include "meepgeom.hpp"

using namespace meep;
using namespace meep_geom;

extern boolean point_in_objectp(vector3 p, GEOMETRIC_OBJECT o);

%}

%include "numpy.i"

%init %{
  import_array();
%}

%{
PyObject *py_callback = NULL;

static PyObject* vec2py(const meep::vec &v);
static double py_callback_wrap(const meep::vec &v);
static int pyv3_to_v3(PyObject *po, vector3 *v);

static int get_attr_v3(PyObject *py_obj, vector3 *v, const char *name);
static int get_attr_dbl(PyObject *py_obj, double *result, const char *name);
// static material_type get_attr_material(PyObject *po);
// static material_type pymaterial_to_material(PyObject *po);

// static susceptibility_struct py_susceptibility_to_susceptibility(PyObject *po);
// static susceptibility_list py_list_to_susceptibility_list(PyObject *po);

// static geometric_object pysphere_to_sphere(PyObject *py_sphere);
static int pycylinder_to_cylinder(PyObject *py_cyl, geometric_object *o);
// static geometric_object pywedge_to_wedge(PyObject *py_wedge);
// static geometric_object pycone_to_cone(PyObject *py_cone);
// static geometric_object pyblock_to_block(PyObject *py_blk);
// static geometric_object pyellipsoid_to_ellipsoid(PyObject *py_ell);
// static geometric_object py_gobj_to_gobj(PyObject *po);
// static geometric_object_list py_list_to_gobj_list(PyObject *po);
// static geometric_object pycgo_to_cgo(PyObject *py_cgo);

#include "typemap_utils.cpp"

%}

// Typemap suite for double func(meep::vec &)

%typemap(in) double (*)(const meep::vec &) {
  $1 = py_callback_wrap;
  py_callback = $input;
  Py_INCREF(py_callback);
}
%typemap(freearg) double (*)(const meep::vec &) {
  Py_XDECREF(py_callback);
}
%typecheck(SWIG_TYPECHECK_POINTER) double (*)(const meep::vec &) {
  $1 = PyCallable_Check($input);
}

// Typemap suite for vector3

%typemap(in) vector3 {
    // TODO(chogan): Accept inputs of tuple, np.array, list?
    if(!pyv3_to_v3($input, &$1)) {
        SWIG_fail;
    }
}

// Typemap suite for GEOMETRIC_OBJECT

%typecheck(SWIG_TYPECHECK_POINTER) GEOMETRIC_OBJECT {
    PyObject *geom_mod = PyImport_ImportModule("geom");
    PyObject *geom_dict = PyModule_GetDict(geom_mod);
    PyObject *py_sphere = PyDict_GetItemString(geom_dict, "GeometricObject");
    $1 = PyObject_IsInstance($input, py_sphere);
    Py_XDECREF(geom_mod);
}

%typemap(in) GEOMETRIC_OBJECT {
    if(!py_gobj_to_gobj($input, &$1)) {
        SWIG_fail;
    }
}

%typemap(freearg) GEOMETRIC_OBJECT {
    geometric_object_destroy($1);
}

%typemap(out) boolean {
    long b = $1 == 0 ? 0 : 1;
    $result = PyBool_FromLong(b);
}

// Typemap suite for geometric_object_list

// %typecheck(SWIG_TYPECHECK_POINTER) geometric_object_list {
//     $1 = PyList_Check($input);
// }

// %typemap(in) geometric_object_list {
//     $1 = py_list_to_gobj_list($input);
// }

// %typemap(freearg) geometric_object_list {
//     for(int i = 0; i < $1.num_items; i++) {
//         geometric_object_destroy($1.items[i]);
//     }
//     delete[] $1.items;
// }

// Rename python builtins
%rename(br_apply) meep::boundary_region::apply;
%rename(_is) meep::dft_chunk::is;
%rename(Meep_None) meep::None;

// Operator renaming
%rename(boundary_region_assign) meep::boundary_region::operator=;

// TODO:  Fix these with a typemap when necessary
%feature("immutable") meep::fields_chunk::connections;
%feature("immutable") meep::fields_chunk::num_connections;

%include "vec.i"
%include "meep.hpp"
%include "meepgeom.hpp"

extern boolean point_in_objectp(vector3 p, GEOMETRIC_OBJECT o);

%ignore eps_func;
%ignore inveps_func;
