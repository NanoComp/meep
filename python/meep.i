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

using namespace meep;

extern boolean point_in_objectp(vector3 p, GEOMETRIC_OBJECT o);

%}

%include "numpy.i"

%init %{
  import_array();
%}

%{
PyObject *py_callback = NULL;

static PyObject* vec2py(const meep::vec &v) {

    double x = 0, y = 0, z = 0;

    switch (v.dim) {
      case meep::D1:
        z = v.z();
        break;
      case meep::D2:
        x = v.x();
        y = v.y();
        break;
      case meep::D3:
        x = v.x();
        y = v.y();
        z = v.z();
        break;
      case meep::Dcyl:
        return Py_BuildValue("(ddd)", v.r(), v.z(), 0);
    }

    return Py_BuildValue("(ddd)", x, y, z);
}

static double py_callback_wrap(const meep::vec &v) {
    PyObject *pyv = vec2py(v);
    PyObject *pyret = PyObject_CallFunctionObjArgs(py_callback, pyv, NULL);
    Py_XDECREF(pyv);
    double ret = PyFloat_AsDouble(pyret);
    Py_XDECREF(pyret);
    return ret;
}


static vector3 pyv3_to_v3(PyObject *po) {
    PyObject *py_x = PyObject_GetAttrString(po, "x");
    PyObject *py_y = PyObject_GetAttrString(po, "y");
    PyObject *py_z = PyObject_GetAttrString(po, "z");
    double x = PyFloat_AsDouble(py_x);
    double y = PyFloat_AsDouble(py_y);
    double z = PyFloat_AsDouble(py_z);
    Py_DECREF(py_x);
    Py_DECREF(py_y);
    Py_DECREF(py_z);

    vector3 v = {x, y, z};
    return v;
}

static material_type pymaterial_to_material(PyObject *po) {
    material_type m;
    // TODO(chogan): I don't think this will work unless the original python
    // data was created from PyLong_FromVoidPtr(). Need to test.
    m.data = PyLong_AsVoidPtr(PyObject_GetAttrString(po, "data"));
    return m;
}

static geometric_object pysphere_to_sphere(PyObject *py_sphere) {
    material_type material = pymaterial_to_material(PyObject_GetAttrString(py_sphere, "material"));

    PyObject *py_center = PyObject_GetAttrString(py_sphere, "center");
    vector3 center = pyv3_to_v3(py_center);

    PyObject *py_radius = PyObject_GetAttrString(py_sphere, "radius");
    double radius = PyFloat_AsDouble(py_radius);

    Py_XDECREF(py_center);
    Py_XDECREF(py_radius);

    return make_sphere(material, center, radius);
}
%}

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

%typemap(in) vector3 {
    // TODO(chogan): Accept inputs of tuple, np.array, list?
    $1 = pyv3_to_v3($input);
}

%typemap(in) GEOMETRIC_OBJECT {
    PyObject *py_type = PyObject_Type($input);
    PyObject *name = PyObject_GetAttrString(py_type, "__name__");

    std::string go_type(PyString_AsString(name));

    if(go_type == "Sphere") {
        $1 = pysphere_to_sphere($input);
    }
    else if(go_type == "Cylinder") {
        // TODO(chogan)
        // $1 = pycylinder_to_cylinder($input);
    }
    else if(go_type == "Wedge") {
        // TODO(chogan)
        // $1 = pywedge_to_wedge($input);
    }
    else if(go_type == "Cone") {
        // TODO(chogan)
        // $1 = pycone_to_cone($input);
    }
    else if(go_type == "Block") {
        // TODO(chogan)
        // $1 = pyblock_to_block($input);
    }
    else if(go_type == "Ellipsoid") {
        // TODO(chogan)
        // $1 = pyellipsoid_to_ellipsoid($input);
    }
    else if(go_type == "CompoundGeometricObject") {
        // TODO(chogan)
        // $1 = pycgo_to_cgo($input);
    }
    else {
        // TODO(chogan): Fail
    }
    Py_XDECREF(py_type);
    Py_XDECREF(name);
}

%typemap(out) boolean {
    int b = $1 == 0 ? 0 : 1;
    $result = PyBool_FromLong(b);
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

extern boolean point_in_objectp(vector3 p, GEOMETRIC_OBJECT o);

%ignore eps_func;
%ignore inveps_func;
