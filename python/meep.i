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

#include <complex>
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
%}

%include "numpy.i"

%init %{
  import_array();
%}

%{
typedef struct {
    PyObject *func;
    int num_components;
} py_field_func_data;

PyObject *py_callback = NULL;
PyObject *py_callback_v3 = NULL;
PyObject *py_amp_func = NULL;

static PyObject *py_geometric_object();
static PyObject *py_source_time_object();
static PyObject *py_material_object();
static PyObject* vec2py(const meep::vec &v);
static double py_callback_wrap(const meep::vec &v);
static std::complex<double> py_amp_func_wrap(const meep::vec &v);
static std::complex<double> py_field_func_wrap(const std::complex<double> *fields,
                                               const meep::vec &loc,
                                               void *data_);
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

// Wrapper for Python PML profile function
double py_pml_profile(double u, void *f) {
    PyObject *func = (PyObject *)f;
    PyObject *d = PyFloat_FromDouble(u);

    if(!PyCallable_Check(func)) {
        PyErr_SetString(PyExc_TypeError, "py_pml_profile: Object is not callable");
        // TODO(chogan): Fix this error handling.
        throw;
    }

    PyObject *pyret = PyObject_CallFunctionObjArgs(func, d, NULL);
    double ret = PyFloat_AsDouble(pyret);
    Py_XDECREF(pyret);
    Py_XDECREF(d);
    return ret;
}

PyObject *py_do_harminv(PyObject *vals, double dt, double f_min, double f_max, int maxbands,
                     double spectral_density, double Q_thresh, double rel_err_thresh,
                     double err_thresh, double rel_amp_thresh, double amp_thresh) {

    std::complex<double> *amp = new std::complex<double>[maxbands];
    double *freq_re = new double[maxbands];
    double *freq_im = new double[maxbands];
    double *freq_err = new double[maxbands];

    Py_ssize_t n = PyList_Size(vals);
    std::complex<double> *items = new std::complex<double>[n];

    for(int i = 0; i < n; i++) {
        Py_complex py_c = PyComplex_AsCComplex(PyList_GetItem(vals, i));
        std::complex<double> c(py_c.real, py_c.imag);
        items[i] = c;
    }

    maxbands = do_harminv(items, n, dt, f_min, f_max, maxbands, amp,
                          freq_re, freq_im, freq_err, spectral_density, Q_thresh,
                          rel_err_thresh, err_thresh, rel_amp_thresh, amp_thresh);

    PyObject *res = PyList_New(maxbands);

    for(int i = 0; i < maxbands; i++) {
        Py_complex pyfreq = {freq_re[i], freq_im[i]};
        Py_complex pyamp = {amp[i].real(), amp[i].imag()};
        Py_complex pyfreq_err = {freq_err[i], 0};

        PyObject *pyobj = Py_BuildValue("(DDD)", &pyfreq, &pyamp, &pyfreq_err);
        PyList_SetItem(res, i, pyobj);
    }

    delete[] freq_err;
    delete[] freq_im;
    delete[] freq_re;
    delete[] amp;
    delete[] items;

    return res;
}

// Wrapper around meep::dft_near2far::farfield
PyObject *_get_farfield(meep::dft_near2far *f, const meep::vec & v) {
    Py_ssize_t len = f->Nfreq * 6;
    PyObject *res = PyList_New(len);

    std::complex<double> *ff_arr = f->farfield(v);

    for (Py_ssize_t i = 0; i < len; i++) {
        PyList_SetItem(res, i, PyComplex_FromDoubles(ff_arr[i].real(), ff_arr[i].imag()));
    }

    delete ff_arr;

    return res;
}

// Wrapper around meep::dft_ldos::ldos
PyObject *_dft_ldos_ldos(meep::dft_ldos *f) {
    Py_ssize_t len = f->Nomega;
    PyObject *res = PyList_New(len);

    double *tmp = f->ldos();

    for (Py_ssize_t i = 0; i < len; i++) {
        PyList_SetItem(res, i, PyFloat_FromDouble(tmp[i]));
    }

    delete tmp;

    return res;
}

// Wrapper around meep::dft_ldos_F
PyObject *_dft_ldos_F(meep::dft_ldos *f) {
    Py_ssize_t len = f->Nomega;
    PyObject *res = PyList_New(len);

    std::complex<double> *tmp = f->F();

    for (Py_ssize_t i = 0; i < len; i++) {
        PyList_SetItem(res, i, PyComplex_FromDoubles(tmp[i].real(), tmp[i].imag()));
    }

    delete tmp;

    return res;
}

// Wrapper arond meep::dft_ldos_J
PyObject *_dft_ldos_J(meep::dft_ldos *f) {
    Py_ssize_t len = f->Nomega;
    PyObject *res = PyList_New(len);

    std::complex<double> *tmp = f->J();

    for (Py_ssize_t i = 0; i < len; i++) {
        PyList_SetItem(res, i, PyComplex_FromDoubles(tmp[i].real(), tmp[i].imag()));
    }

    delete tmp;

    return res;
}

/* This is a wrapper function to fool SWIG...since our list constructor
   takes ownership of the next pointer, we have to make sure that SWIG
   does not garbage-collect volume_list objects.  We do
   this by wrapping a "helper" function around the constructor which
   does not have the %newobject SWIG attribute.   Note that we then
   need to deallocate the list explicitly in Python. */
meep::volume_list *make_volume_list(const meep::volume &v, int c,
                                    std::complex<double> weight,
                                    meep::volume_list *next) {

    return new meep::volume_list(v, c, weight, next);
}

%}

// This is necessary so that SWIG wraps py_pml_profile as a SWIG function
// pointer object instead of as a built-in function
%constant double py_pml_profile(double u, void *f);
%ignore py_pml_profile;
double py_pml_profile(double u, void *f);

PyObject *py_do_harminv(PyObject *vals, double dt, double f_min, double f_max, int maxbands,
                     double spectral_density, double Q_thresh, double rel_err_thresh,
                     double err_thresh, double rel_amp_thresh, double amp_thresh);

PyObject *_get_farfield(meep::dft_near2far *f, const meep::vec & v);
PyObject *_dft_ldos_ldos(meep::dft_ldos *f);
PyObject *_dft_ldos_F(meep::dft_ldos *f);
PyObject *_dft_ldos_J(meep::dft_ldos *f);
meep::volume_list *make_volume_list(const meep::volume &v, int c,
                                    std::complex<double> weight,
                                    meep::volume_list *next);

// Typemap suite for do_harminv

%typecheck(SWIG_TYPECHECK_POINTER) PyObject *vals {
    $1 = PyList_Check($input);
}

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

// Typemap suite for amplitude function

%typecheck(SWIG_TYPECHECK_POINTER) std::complex<double> (*)(const meep::vec &) {
  $1 = PyCallable_Check($input);
}

%typemap(in) std::complex<double> (*)(const meep::vec &) {
    $1 = py_amp_func_wrap;
    py_amp_func = $input;
    Py_INCREF(py_amp_func);
}

%typemap(freearg) std::complex<double> (*)(const meep::vec &) {
    Py_XDECREF(py_amp_func);
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
        delete[] ((material_data *)$1.material)->medium->E_susceptibilities.items;
        delete[] ((material_data *)$1.material)->medium->H_susceptibilities.items;
        delete ((material_data *)$1.material)->medium;
        delete (material_data *)$1.material;
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
        delete[] ((material_data *)$1.items[i].material)->medium->E_susceptibilities.items;
        delete[] ((material_data *)$1.items[i].material)->medium->H_susceptibilities.items;
        delete ((material_data *)$1.items[i].material)->medium;
        delete (material_data *)$1.items[i].material;
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

%typecheck(SWIG_TYPECHECK_POINTER) void *pml_profile_data {
    $1 = PyCallable_Check($input);
}

%typemap(in) void *pml_profile_data {
    $1 = (void*)$input;
}

// Typemap suite for dft_flux

%typemap(out) double* flux {
    int size = arg1->Nfreq;
    $result = PyList_New(size);
    for(int i = 0; i < size; i++) {
        PyList_SetItem($result, i, PyFloat_FromDouble($1[i]));
    }

    delete $1;
}

// Typemap suite for dft_force

%typemap(out) double* force {
    int size = arg1->Nfreq;
    $result = PyList_New(size);
    for(int i = 0; i < size; i++) {
        PyList_SetItem($result, i, PyFloat_FromDouble($1[i]));
    }

    delete $1;
}

// Typemap suite for material_type

%typecheck(SWIG_TYPECHECK_POINTER) material_type {
    $1 = PyObject_IsInstance($input, py_material_object());
}

%typemap(in) material_type {
    if(!pymaterial_to_material($input, &$1)) {
        SWIG_fail;
    }
}

%typemap(freearg) material_type {
    delete[] $1->medium->E_susceptibilities.items;
    delete[] $1->medium->H_susceptibilities.items;
    delete $1->medium;
    delete $1;
}

// Typemap suite for array_slice

%typecheck(SWIG_TYPECHECK_POINTER, fragment="NumPy_Fragments") double* slice {
    $1 = is_array($input);
}

%typemap(in, fragment="NumPy_Macros") double* slice {
    $1 = (double *)array_data($input);
}

%typecheck(SWIG_TYPECHECK_POINTER, fragment="NumPy_Fragments") std::complex<double>* slice {
    $1 = is_array($input);
}

%typemap(in) std::complex<double>* slice {
    $1 = (std::complex<double> *)array_data($input);
}

%typecheck(SWIG_TYPECHECK_POINTER) meep::component {
    $1 = PyInteger_Check($input) && PyInteger_AsLong($input) < 100;
}

%typemap(in) meep::component {
    $1 = static_cast<meep::component>(PyInteger_AsLong($input));
}

%typecheck(SWIG_TYPECHECK_POINTER) meep::derived_component {
    $1 = PyInteger_Check($input) && PyInteger_AsLong($input) >= 100;
}

%typemap(in) meep::derived_component {
    $1 = static_cast<meep::derived_component>(PyInteger_AsLong($input));
}


%typemap(freearg) std::complex<double> (*)(const meep::vec &) {
    Py_XDECREF(py_amp_func);
}

%apply int INPLACE_ARRAY1[ANY] { int [3] };

// typemap suite for field functions

%typecheck(SWIG_TYPECHECK_POINTER) (int num_fields, const meep::component *components,
                                    meep::field_function fun, void *fun_data_) {
    $1 = PySequence_Check($input) &&
         PySequence_Check(PyList_GetItem($input, 0)) &&
         PyCallable_Check(PyList_GetItem($input, 1));
}
%typemap(in) (int num_fields, const meep::component *components, meep::field_function fun, void *fun_data_)
    (py_field_func_data tmp_data) {

    if (!PySequence_Check($input)) {
        PyErr_SetString(PyExc_ValueError, "Expected a sequence");
        SWIG_fail;
    }

    PyObject *cs = PyList_GetItem($input, 0);

    if (!PySequence_Check(cs)) {
        PyErr_SetString(PyExc_ValueError, "Expected first item in list to be a list");
        SWIG_fail;
    }

    PyObject *func = PyList_GetItem($input, 1);

    if (!PyCallable_Check(func)) {
        PyErr_SetString(PyExc_ValueError, "Expected a function");
        SWIG_fail;
    }

    $1 = PyList_Size(cs);
    $2 = new meep::component[$1];

    for (Py_ssize_t i = 0; i < $1; i++) {
        $2[i] = (meep::component)PyInteger_AsLong(PyList_GetItem(cs, i));
    }

    $3 = py_field_func_wrap;

    tmp_data.num_components = $1;
    tmp_data.func = func;
    Py_INCREF(tmp_data.func);
    $4 = &tmp_data;
}

%typemap(freearg) (int num_fields, const meep::component *components, meep::field_function fun, void *fun_data_) {
    delete[] $2;
    Py_XDECREF(tmp_data$argnum.func);
}

%rename(_dft_ldos) meep::dft_ldos::dft_ldos;

// Rename python builtins
%rename(br_apply) meep::boundary_region::apply;
%rename(_is) meep::dft_chunk::is;
%rename(Meep_None) meep::None;

// Operator renaming
%rename(boundary_region_assign) meep::boundary_region::operator=;

%rename(get_field_from_comp) meep::fields::get_field(component, const vec &) const;

%feature("python:cdefaultargs") meep::fields::add_eigenmode_source;

%feature("immutable") meep::fields_chunk::connections;
%feature("immutable") meep::fields_chunk::num_connections;

%include "vec.i"
%include "meep.hpp"
%include "meep/mympi.hpp"
%include "meepgeom.hpp"

extern boolean point_in_objectp(vector3 p, GEOMETRIC_OBJECT o);

%ignore eps_func;
%ignore inveps_func;

%pythoncode %{
    from .geom import (
        Block,
        Cone,
        Cylinder,
        DrudeSusceptibility,
        Ellipsoid,
        GeometricObject,
        LorentzianSusceptibility,
        Medium,
        NoisyDrudeSusceptibility,
        NoisyLorentzianSusceptibility,
        Sphere,
        Susceptibility,
        Vector3,
        Wedge,
        check_nonnegative,
    )
    from .simulation import (
        NO_PARITY,
        EVEN_Z,
        ODD_Z,
        EVEN_Y,
        ODD_Y,
        TE,
        TM,
        Absorber,
        FluxRegion,
        ForceRegion,
        Harminv,
        Identity,
        Mirror,
        Near2FarRegion,
        PML,
        Rotate2,
        Rotate4,
        Simulation,
        Symmetry,
        Volume,
        after_sources,
        after_sources_and_time,
        after_time,
        at_beginning,
        at_end,
        at_every,
        at_time,
        dft_ldos,
        display_progress,
        during_sources,
        get_flux_freqs,
        get_fluxes,
        get_force_freqs,
        get_forces,
        get_near2far_freqs,
        get_ldos_freqs,
        in_point,
        in_volume,
        inf,
        interpolate,
        output_epsilon,
        output_mu,
        output_hpwr,
        output_dpwr,
        output_tot_pwr,
        output_bfield,
        output_bfield_x,
        output_bfield_y,
        output_bfield_z,
        output_bfield_r,
        output_bfield_p,
        output_dfield,
        output_dfield_x,
        output_dfield_y,
        output_dfield_z,
        output_dfield_r,
        output_dfield_p,
        output_efield,
        output_efield_x,
        output_efield_y,
        output_efield_z,
        output_efield_r,
        output_efield_p,
        output_hfield,
        output_hfield_x,
        output_hfield_y,
        output_hfield_z,
        output_hfield_r,
        output_hfield_p,
        output_png,
        output_poynting,
        output_poynting_x,
        output_poynting_y,
        output_poynting_z,
        output_poynting_r,
        output_poynting_p,
        output_sfield,
        output_sfield_x,
        output_sfield_y,
        output_sfield_z,
        output_sfield_r,
        output_sfield_p,
        py_v3_to_vec,
        scale_flux_fields,
        scale_force_fields,
        scale_near2far_fields,
        stop_when_fields_decayed,
        synchronized_magnetic,
        to_appended,
        when_true,
        when_false,
        with_prefix
    )
    from .source import (
        ALL_COMPONENTS,
        ContinuousSource,
        CustomSource,
        EigenModeSource,
        GaussianSource,
        Source,
        SourceTime,
        check_positive,
    )

    if with_mpi():
        try:
            from mpi4py import MPI
        except ImportError:
            print('\n**\n** failed to load python MPI module (mpi4py)\n**\n')
            pass
        else:
            # this variable reference is needed for lazy initialization of MPI
            comm = MPI.COMM_WORLD
            master_printf('\n**\n** successfully loaded python MPI module (mpi4py)\n**\n')

            if not am_master():
                import os
                import sys
                saved_stdout = sys.stdout
                sys.stdout = open(os.devnull, 'w')

    vacuum = Medium(epsilon=1)
    air = Medium(epsilon=1)
    metal = Medium(epsilon=-inf)
    perfect_electric_conductor = Medium(epsilon=-inf)
    perfect_magnetic_conductor = Medium(mu=-inf)
%}
