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

static PyObject* v3_to_pyv3(vector3 *v) {
    PyObject *geom_mod = PyImport_ImportModule("meep.geom");
    PyObject *v3_class = PyObject_GetAttrString(geom_mod, "Vector3");
    PyObject *args = Py_BuildValue("(ddd)", v->x, v->y, v->z);
    PyObject *py_v = PyObject_Call(v3_class, args, NULL);

    Py_DECREF(geom_mod);
    Py_DECREF(args);
    Py_DECREF(v3_class);

    return py_v;
}

static PyObject* cv3_to_pyv3(cvector3 *cv) {
    PyObject *geom_mod = PyImport_ImportModule("meep.geom");
    PyObject *v3_class = PyObject_GetAttrString(geom_mod, "Vector3");

    vector3 r = cvector3_re(*cv);
    vector3 i = cvector3_im(*cv);

    Py_complex x, y, z;
    x.real = r.x;
    x.imag = i.x;
    y.real = r.y;
    y.imag = i.y;
    z.real = r.z;
    z.imag = i.z;

    PyObject *args = Py_BuildValue("(DDD)", &x, &y, &z);
    PyObject *py_v = PyObject_Call(v3_class, args, NULL);

    Py_DECREF(geom_mod);
    Py_DECREF(args);
    Py_DECREF(v3_class);

    return py_v;
}

static PyObject* cmatrix3x3_to_pymatrix(cmatrix3x3 *m) {
    PyObject *c1 = cv3_to_pyv3(&m->c0);
    PyObject *c2 = cv3_to_pyv3(&m->c1);
    PyObject *c3 = cv3_to_pyv3(&m->c2);

    PyObject *geom_mod = PyImport_ImportModule("meep.geom");
    PyObject *matrix_class = PyObject_GetAttrString(geom_mod, "Matrix");

    PyObject *args = Py_BuildValue("(OOO)", c1, c2, c3);
    PyObject *res = PyObject_Call(matrix_class, args, NULL);

    Py_DECREF(c1);
    Py_DECREF(c2);
    Py_DECREF(c3);
    Py_DECREF(geom_mod);
    Py_DECREF(matrix_class);
    Py_DECREF(args);

    return res;
}
static int pyv3_to_cv3(PyObject *po, cvector3 *v) {
    PyObject *py_x = PyObject_GetAttrString(po, "x");
    PyObject *py_y = PyObject_GetAttrString(po, "y");
    PyObject *py_z = PyObject_GetAttrString(po, "z");

    if (!py_x || !py_y || !py_z) {
        PyErr_SetString(PyExc_ValueError, "Vector3 is not initialized");
        return 0;
    }

    std::complex<double> x = std::complex<double>(PyComplex_RealAsDouble(py_x),
                                                  PyComplex_ImagAsDouble(py_x));
    std::complex<double> y = std::complex<double>(PyComplex_RealAsDouble(py_y),
                                                  PyComplex_ImagAsDouble(py_y));
    std::complex<double> z = std::complex<double>(PyComplex_RealAsDouble(py_z),
                                                  PyComplex_ImagAsDouble(py_z));
    Py_DECREF(py_x);
    Py_DECREF(py_y);
    Py_DECREF(py_z);

    v->x.re = x.real();
    v->x.im = x.imag();
    v->y.re = y.real();
    v->y.im = y.imag();
    v->z.re = z.real();
    v->z.im = z.imag();

    return 1;
}

static int get_attr_v3_cmplx(PyObject *py_obj, cvector3 *v, const char *name) {

    PyObject *py_attr = PyObject_GetAttrString(py_obj, name);

    if (!py_attr) {
        PyErr_Format(PyExc_ValueError, "Class attribute '%s' is None\n", name);
        return 0;
    }

    if (!pyv3_to_cv3(py_attr, v)) {
        return 0;
    }

    Py_XDECREF(py_attr);
    return 1;
}

static int pyhermmedium_to_medium(PyObject *po, hermitian_medium *m) {

    if (!get_attr_v3(po, &m->epsilon_diag, "epsilon_diag") ||
        !get_attr_v3(po, &m->mu_diag, "mu_diag")) {

        return 0;
    }
    if (!get_attr_v3_cmplx(po, &m->mu_offdiag, "mu_offdiag") ||
        !get_attr_v3_cmplx(po, &m->epsilon_offdiag, "epsilon_offdiag")) {

        return 0;
    }

    // TODO: Reuse pymedium_to_medium
    PyObject *py_e_susceptibilities = PyObject_GetAttrString(po, "E_susceptibilities");
    PyObject *py_h_susceptibilities = PyObject_GetAttrString(po, "H_susceptibilities");

    if (!py_e_susceptibilities || !py_h_susceptibilities) {
        return 0;
    }

    if (!py_list_to_susceptibility_list(py_e_susceptibilities, &m->E_susceptibilities) ||
       !py_list_to_susceptibility_list(py_h_susceptibilities, &m->H_susceptibilities)) {

        return 0;
    }

    Py_XDECREF(py_e_susceptibilities);
    Py_XDECREF(py_h_susceptibilities);

    if (!get_attr_v3(po, &m->E_chi2_diag, "E_chi2_diag") ||
        !get_attr_v3(po, &m->E_chi3_diag, "E_chi3_diag") ||
        !get_attr_v3(po, &m->H_chi2_diag, "H_chi2_diag") ||
        !get_attr_v3(po, &m->H_chi3_diag, "H_chi3_diag") ||
        !get_attr_v3(po, &m->D_conductivity_diag, "D_conductivity_diag") ||
        !get_attr_v3(po, &m->B_conductivity_diag, "B_conductivity_diag")) {

        return 0;
    }

    return 1;
}

static int pyhermmaterial_to_material(PyObject *po, hermitian_material_data **mdp) {
    hermitian_material_data *md;

    if (PyObject_IsInstance(po, py_material_object())) {
        md = new hermitian_material_data();
        if (!pyhermmedium_to_medium(po, &md->medium)) {
            return 0;
        }
    }
    // TODO: user_material and eps_input_file

    *mdp = md;

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

%apply double INPLACE_ARRAY2[ANY][ANY] {
    double data[3][3]
};

%apply material_type {
    meep_geom::material_data*
};

%apply double { mpb_real };

// Typemap suite for hermitian_material_data

/* %typecheck(SWIG_TYPECHECK_POINTER) hermitian_material_data* = material_type; */

%typemap(in) hermitian_material_data* {
    if (!pyhermmaterial_to_material($input, &$1)) {
        SWIG_fail;
    }
}

%typemap(freearg) hermitian_material_data* {

    if ($1->medium.E_susceptibilities.items) {
        delete[] $1->medium.E_susceptibilities.items;
    }
    if ($1->medium.H_susceptibilities.items) {
        delete[] $1->medium.H_susceptibilities.items;
    }
    free($1);
}

%typemap(in) lattice {
    if (!pylattice_to_lattice($input, &$1)) {
        PyErr_PrintEx(0);
        SWIG_fail;
    }
}

%typemap(in) double resolution[3] (double temp[3]) {

    for (Py_ssize_t i = 0; i < 3; ++i) {
        temp[i] = PyFloat_AsDouble(PyList_GetItem($input, i));
    }
    $1 = &temp[0];
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

%typemap(out) vector3 {
    $result = v3_to_pyv3(&$1);

    if (!$result) {
        SWIG_fail;
    }
}

%typemap(out) cvector3 {
    $result = cv3_to_pyv3(&$1);

    if (!$result) {
        SWIG_fail;
    }
}

%typemap(out) cmatrix3x3 {
    $result = cmatrix3x3_to_pymatrix(&$1);

    if (!$result) {
        SWIG_fail;
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
        output_tot_pwr,
        output_dpwr_in_objects,
        output_poynting,
        output_poynting_x,
        output_poynting_y,
        output_poynting_z,
        output_at_kpoint,
        display_group_velocities,
        display_yparities,
        display_zparities,
        fix_hfield_phase,
        fix_bfield_phase,
        fix_dfield_phase,
        fix_efield_phase,
    )
%}
