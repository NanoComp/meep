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
#define SWIG_FILE_WITH_INIT

#include "pympb.hpp"
#include "meepgeom.hpp"

using namespace py_mpb;
%}

%{
using namespace meep;
using namespace meep_geom;

#include "typemap_utils.cpp"

static int pymatrix_to_matrix(PyObject *po, matrix3x3 *m) {
    vector3 c1, c2, c3;

    PyObject *py_c1 = PyObject_GetAttrString(po, "c1");
    PyObject *py_c2 = PyObject_GetAttrString(po, "c2");
    PyObject *py_c3 = PyObject_GetAttrString(po, "c3");

    if (!pyv3_to_v3(py_c1, &c1) ||
        !pyv3_to_v3(py_c2, &c2) ||
        !pyv3_to_v3(py_c3, &c3)) {

        return 0;
    }

    m->c0 = c1;
    m->c1 = c2;
    m->c2 = c3;

    Py_DECREF(py_c1);
    Py_DECREF(py_c2);
    Py_DECREF(py_c3);

    return 1;
}

static int get_attr_matrix(PyObject *py_obj, matrix3x3 *m, const char *name) {
    PyObject *py_attr = PyObject_GetAttrString(py_obj, name);

    if (!py_attr) {
        PyErr_Format(PyExc_ValueError, "Class attribute '%s' is None\n", name);
        return 0;
    }

    if (!pymatrix_to_matrix(py_attr, m)) {
        return 0;
    }

    Py_XDECREF(py_attr);
    return 1;
}

static int pylattice_to_lattice(PyObject *py_lat, lattice *l) {
    vector3 basis1, basis2, basis3, size, basis_size, b1, b2, b3;
    matrix3x3 basis, metric;

    if (!get_attr_v3(py_lat, &basis1, "basis1") ||
        !get_attr_v3(py_lat, &basis2, "basis2") ||
        !get_attr_v3(py_lat, &basis3, "basis3") ||
        !get_attr_v3(py_lat, &size, "size") ||
        !get_attr_v3(py_lat, &basis_size, "basis_size") ||
        !get_attr_v3(py_lat, &b1, "b1") ||
        !get_attr_v3(py_lat, &b2, "b2") ||
        !get_attr_v3(py_lat, &b3, "b3") ||
        !get_attr_matrix(py_lat, &basis, "basis") ||
        !get_attr_matrix(py_lat, &metric, "metric")) {

        return 0;
    }

    l->basis1 = basis1;
    l->basis2 = basis2;
    l->basis3 = basis3;
    l->size = size;
    l->basis_size = basis_size;
    l->b1 = b1;
    l->b2 = b2;
    l->b3 = b3;
    l->basis = basis;
    l->metric = metric;

    return 1;
}
%}

%include "std_string.i"
%include "numpy.i"
%init %{
    import_array();
%}

%import "meep.i"

%numpy_typemaps(std::complex<mpb_real>, NPY_CDOUBLE, int);

%apply (std::complex<mpb_real>* INPLACE_ARRAY1, int DIM1) {
    (std::complex<mpb_real>* cdata, int size)
};

%apply (double* INPLACE_ARRAY1, int DIM1) {
    (double* data, int size)
};

%apply material_type {
    meep_geom::material_data*
};

%typemap(in) lattice {
    if (!pylattice_to_lattice($input, &$1)) {
        PyErr_PrintEx(0);
        SWIG_fail;
    }
}

%typemap(out) std::vector<mpb_real> {
    Py_ssize_t n = $1.size();

    $result = PyList_New(n);

    for (Py_ssize_t i = 0; i < n; ++i) {
        PyObject *freq = PyFloat_FromDouble($1[i]);
        PyList_SetItem($result, i, freq);
    }
}

%typemap(out) std::vector<int> {
    Py_ssize_t n = $1.size();

    $result = PyList_New(n);

    for (Py_ssize_t i = 0; i < n; ++i) {
        PyObject *dim = PyInteger_FromLong($1[i]);
        PyList_SetItem($result, i, dim);
    }
}

%include "pympb.hpp"

%pythoncode %{
    from .solver import (
        ModeSolver,
        output_hfield,
        output_hfield_x,
        output_hfield_y,
        output_hfield_z,
        output_bfield,
        output_bfield_x,
        output_bfield_y,
        output_bfield_z,
        output_dfield,
        output_dfield_x,
        output_dfield_y,
        output_dfield_z,
        output_efield,
        output_efield_x,
        output_efield_y,
        output_efield_z,
        output_charge_density,
        output_bpwr,
        output_dpwr,
        output_at_kpoint,
        fix_hfield_phase,
        fix_bfield_phase,
        fix_dfield_phase,
        fix_efield_phase,
    )
%}
