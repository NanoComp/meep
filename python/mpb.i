/* Copyright (C) 2005-2018 Massachusetts Institute of Technology
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

%module(package="meep.mpb") mpb

%{
#include "pympb.hpp"
#include "meepgeom.hpp"

using namespace py_mpb;
%}

// TODO: Don't repeat this code here and in meep.i
%{
using namespace meep;
using namespace meep_geom;

typedef struct {
    PyObject *func;
    int num_components;
} py_field_func_data;

PyObject *py_callback = NULL;
PyObject *py_callback_v3 = NULL;
PyObject *py_amp_func = NULL;

static PyObject *py_source_time_object();
static PyObject *py_material_object();
static PyObject *vec2py(const meep::vec &v);
static double py_callback_wrap(const meep::vec &v);
static std::complex<double> py_amp_func_wrap(const meep::vec &v);
static std::complex<double> py_field_func_wrap(const std::complex<double> *fields,
                                               const meep::vec &loc,
                                               void *data_);
static void py_user_material_func_wrap(vector3 x, void *user_data, medium_struct *medium);
static void py_epsilon_func_wrap(vector3 x, void *user_data, medium_struct *medium);
static int pyv3_to_v3(PyObject *po, vector3 *v);

static int get_attr_v3(PyObject *py_obj, vector3 *v, const char *name);
static int get_attr_dbl(PyObject *py_obj, double *result, const char *name);
static int get_attr_int(PyObject *py_obj, int *result, const char *name);
static int get_attr_material(PyObject *po, material_type *m);
static int pymaterial_to_material(PyObject *po, material_type *mt);
static int pymedium_to_medium(PyObject *po, medium_struct *m);
static int pyabsorber_to_absorber(PyObject *py_absorber, meep_geom::absorber *a);
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

%include "numpy.i"
%import "meep.i"

%apply material_type { meep_geom::material_data *};

%typemap(in) lattice {
    if (!pylattice_to_lattice($input, &$1)) {
        PyErr_PrintEx(0);
        SWIG_fail;
    }
}

%typemap(out) std::vector<mpb_real> py_mpb::mode_solver::get_freqs {
    Py_ssize_t n = $1.size();

    $result = PyList_New(n);

    for (Py_ssize_t i = 0; i < n; ++i) {
        PyObject *freq = PyFloat_FromDouble($1[i]);
        PyList_SetItem($result, i, freq);
    }
}

%include "pympb.hpp"

%pythoncode %{
    from .solver import (
        ModeSolver,
    )
%}
