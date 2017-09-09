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
PyObject *py_callback = NULL;

static PyObject *py_geometric_object();
static PyObject *py_source_time_object();
static PyObject *py_material_object();
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
%}

// This is necessary so that SWIG wraps py_pml_profile as a SWIG function
// pointer object instead of as a built-in function
%constant double py_pml_profile(double u, void *f);
%ignore py_pml_profile;
double py_pml_profile(double u, void *f);

PyObject *py_do_harminv(PyObject *vals, double dt, double f_min, double f_max, int maxbands,
                     double spectral_density, double Q_thresh, double rel_err_thresh,
                     double err_thresh, double rel_amp_thresh, double amp_thresh);

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

%typecheck(SWIG_TYPECHECK_POINTER) void *pml_profile_data {
    $1 = PyCallable_Check($input);
}

%typemap(in) void *pml_profile_data {
    $1 = (void*)$input;
}

// Typemap suite for dtf_flux

%typemap(out) double* flux {
    int size = arg1->Nfreq;
    $result = PyList_New(size);
    for(int i = 0; i < size; i++) {
        PyList_SetItem($result, i, PyFloat_FromDouble($1[i]));
    }

  delete $1;
}

%typecheck(SWIG_TYPECHECK_POINTER) material_type _default_material {
    $1 = PyObject_IsInstance($input, py_material_object());
}

%typemap(in) material_type _default_material {
    if(!pymaterial_to_material($input, &$1)) {
        SWIG_fail;
    }
}

// Rename python builtins
%rename(br_apply) meep::boundary_region::apply;
%rename(_is) meep::dft_chunk::is;
%rename(Meep_None) meep::None;

// Operator renaming
%rename(boundary_region_assign) meep::boundary_region::operator=;

%rename(get_field_from_comp) meep::fields::get_field(component, const vec &) const;

// TODO:  Fix these with a typemap when necessary
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
        Absorber,
        FluxRegion,
        Harminv,
        Identity,
        Mirror,
        Pml,
        Rotate2,
        Rotate4,
        Simulation,
        Symmetry,
        Volume,
        after_sources,
        after_time,
        at_beginning,
        at_end,
        at_every,
        during_sources,
        display_progress,
        get_flux_freqs,
        get_fluxes,
        in_volume,
        interpolate,
        output_epsilon,
        output_hfield_z,
        output_efield_z,
        py_v3_to_vec,
        stop_when_fields_decayed,
        to_appended
    )
    from .source import (
        ContinuousSource,
        CustomSource,
        EigenModeSource,
        GaussianSource,
        Source,
        SourceTime,
        check_positive,
    )
%}
