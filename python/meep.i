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
#include "meep/mympi.hpp"
#include "ctl-math.h"
#include "ctlgeom.h"
#include "meepgeom.hpp"

using namespace meep;
using namespace meep_geom;

extern boolean point_in_objectp(vector3 p, GEOMETRIC_OBJECT o);
extern number adaptive_integration(multivar_func f, number *xmin, number *xmax,
                                   integer n, void *fdata, number abstol,
                                   number reltol, integer maxnfe, number *esterr,
                                   integer *errflag);

extern double py_pml_profile(double u, void *f);
extern double py_pml_profile2(int dim, double *u, void *f);
extern double py_pml_profile2u(int dim, double *u, void *f);
%}

%inline %{
double py_pml_helper(PyObject *f, PyObject *d) {
    if(!PyCallable_Check(f)) {
        PyErr_SetString(PyExc_TypeError, "py_pml_profile: Object is not callable");
        // TODO(chogan): Fix this error handling.
        throw;
    }
    PyObject *pyret = PyObject_CallFunctionObjArgs(f, d, NULL);
    double ret = PyFloat_AsDouble(pyret);
    Py_XDECREF(pyret);
    Py_XDECREF(d);
    return ret;
}

// Wrapper for Python PML profile function
double py_pml_profile(double u, void *f) {
    PyObject *func = (PyObject *)f;
    PyObject *d = PyFloat_FromDouble(u);

    return py_pml_helper(func, d);
}

// For passing to multidimensional integration routine
double py_pml_profile2(int dim, double *u, void *f) {
    PyObject *func = (PyObject *)f;
    PyObject *d = PyFloat_FromDouble(*u);
    (void)dim;

    return py_pml_helper(func, d);
}

// For integrating profile(u) * u
double py_pml_profile2u(int dim, double *u, void *f) {
    PyObject *func = (PyObject *)f;
    PyObject *d = PyFloat_FromDouble(*u);
    (void)dim;

    return py_pml_helper(func, d) * (*u);
}

number f_py_wrapper(integer n, number *x, void *f_py_p) {
    PyObject *f_py = (PyObject *)f_py_p;

    PyObject *py_list = PyList_New(n);

    for(int i = 0; i < n; i++) {
        PyObject *num = PyFloat_FromDouble(x[i]);
        Py_INCREF(num);
        PyList_SetItem(py_list, i, num);
    }

    // TODO(chogan): Double check the ref counting for py_list and its items.
    PyObject *py_ret = PyObject_CallFunctionObjArgs(f_py, py_list, NULL);

    number ret = PyFloat_AsDouble(py_ret);
    Py_XDECREF(py_ret);
    return ret;
}

// Python wrapper for adaptive_integration
PyObject *py_adaptive_integration(PyObject *py_func, PyObject *py_xmin,
                                  PyObject *py_xmax, PyObject *py_abstol,
                                  PyObject *py_reltol, PyObject *py_maxnfe) {
    int n, maxnfe, errflag;
    number abstol, reltol, integral;
    number *xmin, *xmax;

    if(!PyList_Check(py_xmin) || !PyList_Check(py_xmax)) {
        PyErr_SetString(PyExc_TypeError, "adaptive_integration: Expected a list");
        throw;
    }

    n = (int)PyList_Size(py_xmin);
    abstol = fabs(PyFloat_AsDouble(py_abstol));
    reltol = fabs(PyFloat_AsDouble(py_reltol));
    maxnfe = PyInt_AsLong(py_maxnfe);

    // if(PyErr_Occurred()) {
    //     PyErr_Print();
    //     PyErr_SetString(PyExc_RuntimeError, "adaptive_integration: Couldn't convert Python maxnfe to an integer");
    //     throw;
    // }

    if((int)PyList_Size(py_xmax) != n) {
        PyErr_SetString(PyExc_ValueError, "adaptive_integration: xmin/xmax dimension mismatch");
        throw;
    }

    xmin = (number *)malloc(sizeof(number) * n);
    xmax = (number *)malloc(sizeof(number) * n);

    if(!xmin || !xmax) {
        PyErr_SetString(PyExc_RuntimeError, "adaptive_integration: error, out of memory!");
        throw;
    }

    for(int i = 0; i < n; i++) {
        xmin[i] = PyFloat_AsDouble(PyList_GetItem(py_xmin, i));
        xmax[i] = PyFloat_AsDouble(PyList_GetItem(py_xmax, i));
    }

    integral = adaptive_integration(f_py_wrapper, xmin, xmax, n, (void*)py_func, abstol,
                                    reltol, maxnfe, &abstol, &errflag);

    free(xmax);
    free(xmin);

    switch(errflag) {
        case 3:
            PyErr_SetString(PyExc_RuntimeError, "adaptive_integration: invalid inputs");
            throw;
        case 1:
            PyErr_SetString(PyExc_RuntimeError, "adaptive_integration: maxnfe too small");
            throw;
        case 2:
            PyErr_SetString(PyExc_RuntimeError, "adaptive_integration: lenwork too small");
            throw;
        default:
            break;
    }

    return Py_BuildValue("[dd]", integral, abstol);
}
%}

%include "numpy.i"

%init %{
  import_array();
%}

%{
PyObject *py_callback = NULL;

static PyObject *py_geometric_object();
static PyObject *py_source_time_object();
static PyObject* vec2py(const meep::vec &v);
static double py_callback_wrap(const meep::vec &v);
static int pyv3_to_v3(PyObject *po, vector3 *v);

static int get_attr_v3(PyObject *py_obj, vector3 *v, const char *name);
static int get_attr_dbl(PyObject *py_obj, double *result, const char *name);
static int get_attr_material(PyObject *po, material_type *m);
static int pymaterial_to_material(PyObject *po, material_type *mt);

static int py_susceptibility_to_susceptibility(PyObject *po, susceptibility_struct *s);
static int py_list_to_susceptibility_list(PyObject *po, susceptibility_list *sl);

static int pysphere_to_sphere(PyObject *py_sphere, geometric_object *go);
static int pycylinder_to_cylinder(PyObject *py_cyl, geometric_object *o);
static int pywedge_to_wedge(PyObject *py_wedge, geometric_object *w);
static int pycone_to_cone(PyObject *py_cone, geometric_object *cone);
static int pyblock_to_block(PyObject *py_blk, geometric_object *blk);
static int pyellipsoid_to_ellipsoid(PyObject *py_ell, geometric_object *e);
static std::string py_class_name_as_string(PyObject *po);
static int py_gobj_to_gobj(PyObject *po, geometric_object *o);
static int py_list_to_gobj_list(PyObject *po, geometric_object_list *l);

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
    if(!pyv3_to_v3($input, &$1)) {
        SWIG_fail;
    }
}

// Typemap suite for GEOMETRIC_OBJECT

%typecheck(SWIG_TYPECHECK_POINTER) GEOMETRIC_OBJECT {
    $1 = PyObject_IsInstance($input, py_geometric_object());
}

%typemap(in) GEOMETRIC_OBJECT {
    if(!py_gobj_to_gobj($input, &$1)) {
        SWIG_fail;
    }
}

%typemap(freearg) GEOMETRIC_OBJECT {
    if($1.subclass.sphere_data || $1.subclass.cylinder_data || $1.subclass.block_data) {
        geometric_object_destroy($1);
    }
}

// Typemap suite for boolean

%typemap(out) boolean {
    $result = PyBool_FromLong($1);
}

// Typemap suite for geometric_object_list

%typecheck(SWIG_TYPECHECK_POINTER) geometric_object_list {
    $1 = PyList_Check($input);
}

%typemap(in) geometric_object_list {
    if(!py_list_to_gobj_list($input, &$1)) {
        SWIG_fail;
    }
}

%typemap(freearg) geometric_object_list {
    for(int i = 0; i < $1.num_items; i++) {
        geometric_object_destroy($1.items[i]);
    }
    delete[] $1.items;
}

// Typemap suite for susceptibility_list

%typecheck(SWIG_TYPECHECK_POINTER) susceptibility_list {
    $1 = PyList_Check($input);
}

%typemap(in) susceptibility_list {
    if(!py_list_to_susceptibility_list($input, &$1)) {
        SWIG_fail;
    }
}

%typemap(freearg) susceptibility_list {
    delete[] $1.items;
}

// Typemap suite for sources

%typecheck(SWIG_TYPECHECK_POINTER) const meep::src_time & {
    int py_source_time = PyObject_IsInstance($input, py_source_time_object());
    int swig_src_time = PyObject_IsInstance($input, py_meep_src_time_object());

    $1 = py_source_time || swig_src_time;
}

%typemap(in) const meep::src_time & {
    PyObject *swig_obj = NULL;
    void *tmp_ptr = 0;
    int tmp_res = 0;

    if(PyObject_IsInstance($input, py_source_time_object())) {
        swig_obj = PyObject_GetAttrString($input, "swigobj");
    } else if(PyObject_IsInstance($input, py_meep_src_time_object())) {
        swig_obj = $input;
        Py_XINCREF(swig_obj);
    } else {
        PyErr_SetString(PyExc_TypeError, "Expected a meep.source.SourceTime or a meep.src_time\n");
        SWIG_fail;
    }

    tmp_res = SWIG_ConvertPtr(swig_obj, &tmp_ptr, $1_descriptor, 0);
    Py_XDECREF(swig_obj);

    if(!SWIG_IsOK(tmp_res)) {
        SWIG_exception_fail(SWIG_ArgError(tmp_res), "Couldn't convert Python object to meep::src_time");
    }
    $1 = reinterpret_cast<meep::src_time *>(tmp_ptr);
}

// Typemap suite for boundary_region
%typecheck(SWIG_TYPECHECK_POINTER) double (*)(double u, void *func_data) {
    $1 = PyCallable_Check($input);
}

%typemap(in) double (*)(double u, void *func_data) {
    $1 = (pml_profile_func)$input;
}

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
%include "meep/mympi.hpp"
%include "meepgeom.hpp"

extern boolean point_in_objectp(vector3 p, GEOMETRIC_OBJECT o);
extern number adaptive_integration(multivar_func f, number *xmin, number *xmax,
                                   integer n, void *fdata, number abstol,
                                   number reltol, integer maxnfe, number *esterr,
                                   integer *errflag);
%ignore eps_func;
%ignore inveps_func;
