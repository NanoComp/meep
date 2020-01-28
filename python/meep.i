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

%module meep

%import "config.h"

%{
#define SWIG_FILE_WITH_INIT
#define SWIG_PYTHON_2_UNICODE

#include <complex>
#include <string>

#include "config.h"
#include "meep/vec.hpp"
#include "meep.hpp"
#include "meep/mympi.hpp"
#include "ctl-math.h"
#include "ctlgeom.h"
#include "meepgeom.hpp"

namespace meep {
    size_t dft_chunks_Ntotal(dft_chunk *dft_chunks, size_t *my_start);
    typedef std::complex<double> (*amplitude_function)(const vec &);
    void (*master_printf_callback)(const char *s);
}

#ifdef HAVE_MPB
#include "mpb.h"

namespace meep {
    struct eigenmode_data {
        maxwell_data *mdata;
        scalar_complex *fft_data_H, *fft_data_E;
        evectmatrix H;
        int n[3];
        double s[3];
        double Gk[3];
        vec center;
        amplitude_function amp_func;
        int band_num;
        double omega;
        double group_velocity;
    };
}
#else
namespace meep {
    struct eigenmode_data {};
}
#endif

using namespace meep;
using namespace meep_geom;

extern boolean point_in_objectp(vector3 p, GEOMETRIC_OBJECT o);
extern boolean point_in_periodic_objectp(vector3 p, GEOMETRIC_OBJECT o);
void display_geometric_object_info(int indentby, GEOMETRIC_OBJECT o);

%}

%include "numpy.i"
%include "std_vector.i"

%init %{
  import_array();
%}

%{
typedef struct {
    PyObject *func;
    int num_components;
} py_field_func_data;


#include "typemap_utils.cpp"

static PyObject *py_source_time_object() {
    static PyObject *source_time_object = NULL;
    if (source_time_object == NULL) {
        PyObject *source_mod = PyImport_ImportModule("meep.source");
        source_time_object = PyObject_GetAttrString(source_mod, "SourceTime");
        Py_XDECREF(source_mod);
    }
    return source_time_object;
}

static PyObject *py_meep_src_time_object() {
    static PyObject *src_time = NULL;
    if (src_time == NULL) {
        PyObject *meep_mod = PyImport_ImportModule("meep");
        src_time = PyObject_GetAttrString(meep_mod, "src_time");
        Py_XDECREF(meep_mod);
    }
    return src_time;
}

static double py_callback_wrap(const meep::vec &v) {
    PyObject *pyv = vec2py(v);
    PyObject *pyret = PyObject_CallFunctionObjArgs(py_callback, pyv, NULL);
    double ret = PyFloat_AsDouble(pyret);
    Py_XDECREF(pyret);
    return ret;
}

static std::complex<double> py_amp_func_wrap(const meep::vec &v) {
    PyObject *pyv = vec2py(v);
    PyObject *pyret = PyObject_CallFunctionObjArgs(py_amp_func, pyv, NULL);
    double real = PyComplex_RealAsDouble(pyret);
    double imag = PyComplex_ImagAsDouble(pyret);
    std::complex<double> ret(real, imag);
    Py_DECREF(pyret);
    return ret;
}

static std::complex<double> py_field_func_wrap(const std::complex<double> *fields,
                                               const meep::vec &loc,
                                               void *data_) {
    PyObject *pyv = vec2py(loc);

    py_field_func_data *data = (py_field_func_data *)data_;
    int len = data->num_components;

    PyObject *py_args = PyTuple_New(len + 1);
    // Increment here because PyTuple_SetItem steals a reference
    Py_INCREF(pyv);
    PyTuple_SetItem(py_args, 0, pyv);

    for (Py_ssize_t i = 1; i < len + 1; i++) {
        PyObject *cmplx = PyComplex_FromDoubles(fields[i - 1].real(), fields[i - 1].imag());
        PyTuple_SetItem(py_args, i, cmplx);
    }

    PyObject *pyret = PyObject_CallObject(data->func, py_args);

    if (!pyret) { abort_with_stack_trace(); }

    double real = PyComplex_RealAsDouble(pyret);
    double imag = PyComplex_ImagAsDouble(pyret);
    std::complex<double> ret(real, imag);
    Py_DECREF(pyret);
    Py_DECREF(py_args);
    return ret;
}

static std::complex<double> py_src_func_wrap(double t, void *f) {
    PyObject *py_t = PyFloat_FromDouble(t);
    PyObject *pyres = PyObject_CallFunctionObjArgs((PyObject *)f, py_t, NULL);
    double real = PyComplex_RealAsDouble(pyres);
    double imag = PyComplex_ImagAsDouble(pyres);
    std::complex<double> ret(real, imag);
    Py_DECREF(py_t);
    Py_DECREF(pyres);

    return ret;
}

static meep::vec py_kpoint_func_wrap(double freq, int mode, void *user_data) {
    PyObject *py_freq = PyFloat_FromDouble(freq);
    PyObject *py_mode = PyInteger_FromLong(mode);

    PyObject *py_result = PyObject_CallFunctionObjArgs((PyObject*)user_data, py_freq, py_mode, NULL);

    if (!py_result) {
        PyErr_PrintEx(0);
        Py_DECREF(py_freq);
        Py_DECREF(py_mode);
        return meep::vec(0, 0, 0);
    }

    vector3 v3;
    if (!pyv3_to_v3(py_result, &v3)) {
        PyErr_PrintEx(0);
        Py_DECREF(py_freq);
        Py_DECREF(py_mode);
        Py_XDECREF(py_result);
        return meep::vec(0, 0, 0);
    }

    meep::vec result(v3.x, v3.y, v3.z);

    Py_DECREF(py_freq);
    Py_DECREF(py_mode);
    Py_DECREF(py_result);

    return result;
}

void py_master_printf_wrap(const char *s) {
    PySys_WriteStdout("%s", s);
    static PyObject *py_stdout = NULL;
    if (py_stdout == NULL) {
      py_stdout = PySys_GetObject("stdout");
    }
    PyObject *result = PyObject_CallMethod(py_stdout, "flush", NULL);
    Py_XDECREF(result);
}

void set_ctl_printf_callback(void (*callback)(const char *s)) {
#if HAVE_CTL_PRINTF_CALLBACK
  ctl_printf_callback = callback;
#else
  (void)callback;
#endif
}

void set_mpb_printf_callback(void (*callback)(const char *s)) {
#if HAVE_MPB_PRINTF_CALLBACK
  mpb_printf_callback = callback;
#else
  (void)callback;
#endif
}

static int pyabsorber_to_absorber(PyObject *py_absorber, meep_geom::absorber *a) {

    if (!get_attr_dbl(py_absorber, &a->thickness, "thickness") ||
        !get_attr_int(py_absorber, &a->direction, "direction") ||
        !get_attr_int(py_absorber, &a->side, "side") ||
        !get_attr_dbl(py_absorber, &a->R_asymptotic, "R_asymptotic") ||
        !get_attr_dbl(py_absorber, &a->mean_stretch, "mean_stretch")) {

        return 0;
    }

    PyObject *py_pml_profile_func = PyObject_GetAttrString(py_absorber, "pml_profile");

    if (!py_pml_profile_func) {
      meep::abort("Class attribute 'pml_profile' is None\n");
    }

    a->pml_profile_data = py_pml_profile_func;

    return 1;
}

// Wrapper for Python PML profile function
double py_pml_profile(double u, void *f) {
    PyObject *func = (PyObject *)f;
    PyObject *d = PyFloat_FromDouble(u);

    PyObject *pyret = PyObject_CallFunctionObjArgs(func, d, NULL);

    if (!pyret) { abort_with_stack_trace(); }

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

    delete[] ff_arr;

    return res;
}

// Wrapper around meep::dft_near2far::get_farfields_array
 PyObject *_get_farfields_array(meep::dft_near2far *n2f, const meep::volume &where,
                                double resolution) {
    size_t dims[4] = {1, 1, 1, 1};
    int rank = 0;
    size_t N = 1;

    // TODO: Support single precision?
    if (sizeof(realnum) == sizeof(float)) abort("Single precision not supported for get_farfields");

    meep::realnum *EH = n2f->get_farfields_array(where, rank, dims, N, resolution);

    if (!EH) return PyArray_SimpleNew(0, 0, NPY_CDOUBLE);

    // frequencies are the last dimension
    if (n2f->Nfreq > 1) dims[rank++] = n2f->Nfreq;

    // Additional rank to store all 12 E/H x/y/z r/i arrays.
    rank++;
    npy_intp *arr_dims = new npy_intp[rank];
    arr_dims[0] = 12;
    for (int i = 1; i < rank; ++i) {
        arr_dims[i] = dims[i - 1];
    }

    PyObject *py_arr = PyArray_SimpleNew(rank, arr_dims, NPY_DOUBLE);
    memcpy(PyArray_DATA((PyArrayObject*)py_arr), EH, sizeof(meep::realnum) * 2 * N * 6 * n2f->Nfreq);

    delete[] arr_dims;
    delete[] EH;
    return py_arr;

}

// Wrapper around meep::dft_ldos::ldos
PyObject *_dft_ldos_ldos(meep::dft_ldos *f) {
    Py_ssize_t len = f->Nomega;
    PyObject *res = PyList_New(len);

    double *tmp = f->ldos();

    for (Py_ssize_t i = 0; i < len; i++) {
        PyList_SetItem(res, i, PyFloat_FromDouble(tmp[i]));
    }

    delete[] tmp;

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

    delete[] tmp;

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

    delete[] tmp;

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

template<typename dft_type>
PyObject *_get_dft_array(meep::fields *f, dft_type dft, meep::component c, int num_freq) {
    int rank;
    size_t dims[3];
    std::complex<double> *dft_arr = f->get_dft_array(dft, c, num_freq, &rank, dims);

    if (rank==0 || dft_arr==NULL) // this can happen e.g. if component c vanishes by symmetry
     return PyArray_SimpleNew(0, 0, NPY_CDOUBLE);

    size_t length = 1;
    npy_intp *arr_dims = new npy_intp[rank];
    for (int i = 0; i < rank; ++i) {
        arr_dims[i] = dims[i];       // implicit size_t -> int cast, presumed safe for individual array dimensions
        length *= dims[i];
    }

    PyObject *py_arr = PyArray_SimpleNew(rank, arr_dims, NPY_CDOUBLE);
    memcpy(PyArray_DATA((PyArrayObject*)py_arr), dft_arr, sizeof(std::complex<double>) * length);
    delete[] dft_arr;
    if (arr_dims) delete[] arr_dims;

    return py_arr;
}

size_t _get_dft_data_size(meep::dft_chunk *dc) {
    size_t istart;
    return meep::dft_chunks_Ntotal(dc, &istart) / 2;
}

void _get_dft_data(meep::dft_chunk *dc, std::complex<meep::realnum> *cdata, int size) {
    size_t istart;
    size_t n = meep::dft_chunks_Ntotal(dc, &istart) / 2;
    istart /= 2;

    if (n != (size_t)size) {
        meep::abort("Total dft_chunks size does not agree with size allocated for output array.\n");
    }

    for (meep::dft_chunk *cur = dc; cur; cur = cur->next_in_dft) {
        size_t Nchunk = cur->N * cur->Nomega;
        for (size_t i = 0; i < Nchunk; ++i) {
            cdata[i + istart] = cur->dft[i];
        }
        istart += Nchunk;
    }
}

void _load_dft_data(meep::dft_chunk *dc, std::complex<meep::realnum> *cdata, int size) {
    size_t istart;
    size_t n = meep::dft_chunks_Ntotal(dc, &istart) / 2;
    istart /= 2;

    if (n != (size_t)size) {
        meep::abort("Total dft_chunks size does not agree with size allocated for output array.\n");
    }

    for (meep::dft_chunk *cur = dc; cur; cur = cur->next_in_dft) {
        size_t Nchunk = cur->N * cur->Nomega;
        for (size_t i = 0; i < Nchunk; ++i) {
            cur->dft[i] = cdata[i + istart];
        }
        istart += Nchunk;
    }
}

struct kpoint_list {
    meep::vec *kpoints;
    size_t n;
    meep::vec *kdom;
    size_t num_bands;
};

kpoint_list get_eigenmode_coefficients_and_kpoints(meep::fields *f, meep::dft_flux flux, const meep::volume &eig_vol,
                                                   int *bands, int num_bands, int parity, double eig_resolution,
                                                   double eigensolver_tol, std::complex<double> *coeffs,
                                                   double *vgrp, meep::kpoint_func user_kpoint_func,
                                                   void *user_kpoint_data, meep::direction d) {

    size_t num_kpoints = num_bands * flux.Nfreq;
    meep::vec *kpoints = new meep::vec[num_kpoints];
    meep::vec *kdom = new meep::vec[num_kpoints];

    f->get_eigenmode_coefficients(flux, eig_vol, bands, num_bands, parity, eig_resolution, eigensolver_tol,
                                  coeffs, vgrp, user_kpoint_func, user_kpoint_data, kpoints, kdom, d);

    kpoint_list res = {kpoints, num_kpoints, kdom, num_kpoints};

    return res;
}

PyObject *_get_array_slice_dimensions(meep::fields *f, const meep::volume &where, size_t dims[3],
                                      bool collapse_empty_dimensions, bool snap_empty_dimensions) {
    meep::direction dirs[3] = {meep::X, meep::X, meep::X};
    int rank = f->get_array_slice_dimensions(where, dims, dirs, collapse_empty_dimensions, snap_empty_dimensions);

    PyObject *py_dirs = PyList_New(3);
    for (Py_ssize_t i = 0; i < 3; ++i) {
        PyList_SetItem(py_dirs, i, PyInteger_FromLong(static_cast<int>(dirs[i])));
    }

    return Py_BuildValue("(iO)", rank, py_dirs);
}

#ifdef HAVE_MPB
meep::eigenmode_data *_get_eigenmode(meep::fields *f, double omega_src, meep::direction d, const meep::volume where,
                                     const meep::volume eig_vol, int band_num, const meep::vec &_kpoint,
                                     bool match_frequency, int parity, double resolution, double eigensolver_tol,
                                     double kdom[3]) {

    void *data = f->get_eigenmode(omega_src, d, where, eig_vol, band_num, _kpoint, match_frequency,
                                  parity, resolution, eigensolver_tol, kdom);
    return (meep::eigenmode_data *)data;
}

PyObject *_get_eigenmode_Gk(meep::eigenmode_data *emdata) {
    PyObject *v3_class = py_vector3_object();
    PyObject *args = Py_BuildValue("(ddd)", emdata->Gk[0], emdata->Gk[1], emdata->Gk[2]);
    PyObject *result = PyObject_Call(v3_class, args, NULL);
    Py_DECREF(args);
    return result;
}

#else
void _get_eigenmode(meep::fields *f, double omega_src, meep::direction d, const meep::volume where,
                    const meep::volume eig_vol, int band_num, const meep::vec &_kpoint,
                    bool match_frequency, int parity, double resolution, double eigensolver_tol,
                    double kdom[3]) {
    (void) f; (void) omega_src; (void) d; (void) where; (void) eig_vol; (void) band_num; (void) _kpoint;
    (void) match_frequency; (void) parity; (void) resolution; (void) eigensolver_tol;
    (void) kdom;
    meep::abort("Must compile Meep with MPB for get_eigenmode");
}
#endif
%}

%numpy_typemaps(std::complex<meep::realnum>, NPY_CDOUBLE, int);
%numpy_typemaps(std::complex<double>, NPY_CDOUBLE, size_t);

%apply (std::complex<meep::realnum> *INPLACE_ARRAY1, int DIM1) {(std::complex<meep::realnum> *cdata, int size)};

// add_volume_source
%apply (std::complex<double> *INPLACE_ARRAY3, size_t DIM1, size_t DIM2, size_t DIM3) {
    (std::complex<double> *arr, size_t dim1, size_t dim2, size_t dim3)
};

// This is necessary so that SWIG wraps py_pml_profile as a SWIG function
// pointer object instead of as a built-in function
%constant double py_pml_profile(double u, void *f);
%ignore py_pml_profile;
double py_pml_profile(double u, void *f);

%constant void py_master_printf_wrap(const char *s);
namespace meep {
    void (*master_printf_callback)(const char *s);
}
void set_ctl_printf_callback(void (*callback)(const char *s));
void set_mpb_printf_callback(void (*callback)(const char *s));

PyObject *py_do_harminv(PyObject *vals, double dt, double f_min, double f_max, int maxbands,
                     double spectral_density, double Q_thresh, double rel_err_thresh,
                     double err_thresh, double rel_amp_thresh, double amp_thresh);

PyObject *_get_farfield(meep::dft_near2far *f, const meep::vec & v);
PyObject *_get_farfields_array(meep::dft_near2far *n2f, const meep::volume &where, double resolution);
PyObject *_dft_ldos_ldos(meep::dft_ldos *f);
PyObject *_dft_ldos_F(meep::dft_ldos *f);
PyObject *_dft_ldos_J(meep::dft_ldos *f);
template<typename dft_type>
PyObject *_get_dft_array(meep::fields *f, dft_type dft, meep::component c, int num_freq);
size_t _get_dft_data_size(meep::dft_chunk *dc);
void _get_dft_data(meep::dft_chunk *dc, std::complex<meep::realnum> *cdata, int size);
void _load_dft_data(meep::dft_chunk *dc, std::complex<meep::realnum> *cdata, int size);
meep::volume_list *make_volume_list(const meep::volume &v, int c,
                                    std::complex<double> weight,
                                    meep::volume_list *next);

// Typemap suite for get_eigenmode_coefficients_and_kpoints

%typemap(out) kpoint_list {

    PyObject *py_kpoints = PyList_New($1.n);
    PyObject *py_kdom = PyList_New($1.num_bands);

    for (size_t i = 0; i < $1.n; ++i) {
        PyList_SetItem(py_kpoints, i, vec2py($1.kpoints[i], true));
    }
    for (size_t i = 0; i < $1.num_bands; ++i) {
        PyList_SetItem(py_kdom, i, vec2py($1.kdom[i], true));
    }

    $result = Py_BuildValue("(O,O)", py_kpoints, py_kdom);

    delete[] $1.kpoints;
    delete[] $1.kdom;
}

// Typemap suite for do_harminv

%typecheck(SWIG_TYPECHECK_POINTER) PyObject *vals {
    $1 = PyList_Check($input);
}

// Typemap suite for double func(meep::vec &)

%typemap(in) double (*)(const meep::vec &) {
  if ($input == Py_None) {
    $1 = NULL;
    py_callback = NULL;
  } else {
    $1 = py_callback_wrap;
    py_callback = $input;
    Py_INCREF(py_callback);
  }
}

%typemap(freearg) double (*)(const meep::vec &) {
  Py_XDECREF(py_callback);
}

%typecheck(SWIG_TYPECHECK_POINTER) double (*)(const meep::vec &) {
  $1 = PyCallable_Check($input) || $input == Py_None;
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

%typecheck (SWIG_TYPECHECK_POINTER) vector3 {
    $1 = PyObject_IsInstance($input, py_vector3_object());
}

%typemap(in) vector3 {
    if(!pyv3_to_v3($input, &$1)) {
        SWIG_fail;
    }
}

// Typemap suite for GEOMETRIC_OBJECT

%typemap(in) GEOMETRIC_OBJECT {
    if(!py_gobj_to_gobj($input, &$1)) {
        SWIG_fail;
    }
}

%typemap(freearg) GEOMETRIC_OBJECT {
    if($1.subclass.sphere_data || $1.subclass.cylinder_data || $1.subclass.block_data) {
        if (((material_data *)$1.material)->medium.E_susceptibilities.items) {
            delete[] ((material_data *)$1.material)->medium.E_susceptibilities.items;
        }
        if (((material_data *)$1.material)->medium.H_susceptibilities.items) {
            delete[] ((material_data *)$1.material)->medium.H_susceptibilities.items;
        }
        delete[] ((material_data *)$1.material)->epsilon_data;
        delete (material_data *)$1.material;
        geometric_object_destroy($1);
    }
}

%typemap(out) geometric_object {
    $result = gobj_to_py_obj(&$1);

    if (!$result) {
        SWIG_fail;
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
        if (((material_data *)$1.items[i].material)->medium.E_susceptibilities.items) {
            delete[] ((material_data *)$1.items[i].material)->medium.E_susceptibilities.items;
        }
        if (((material_data *)$1.items[i].material)->medium.H_susceptibilities.items) {
            delete[] ((material_data *)$1.items[i].material)->medium.H_susceptibilities.items;
        }
        delete[] ((material_data *)$1.items[i].material)->epsilon_data;
        delete (material_data *)$1.items[i].material;
        geometric_object_destroy($1.items[i]);
    }
    delete[] $1.items;
}

%typemap(out) geometric_object_list {
    $result = gobj_list_to_py_list(&$1);

    if (!$result) {
        SWIG_fail;
    }
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
      meep::abort("Expected a meep.source.SourceTime or a meep.src_time\n");
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

    delete[] $1;
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
    int py_material = PyObject_IsInstance($input, py_material_object());
    int user_material = PyFunction_Check($input);
    int file_material = IsPyString($input);
    int numpy_material = PyArray_Check($input);

    $1 = py_material || user_material || file_material || numpy_material;
}

%typemap(in) material_type {
    if(!pymaterial_to_material($input, &$1)) {
        SWIG_fail;
    }
}

%typemap(freearg) material_type {
    if ($1->medium.E_susceptibilities.items) {
        delete[] $1->medium.E_susceptibilities.items;
    }
    if ($1->medium.H_susceptibilities.items) {
        delete[] $1->medium.H_susceptibilities.items;
    }
    delete[] $1->epsilon_data;
    delete $1;
}

// For some reason SWIG needs the namespaced version too
%apply material_type { meep_geom::material_type };

// Typemap suite for get_array_metadata

%typecheck(SWIG_TYPECHECK_POINTER, fragment="NumPy_Fragments") double* xtics {
    $1 = is_array($input);
}
%typemap(in, fragment="NumPy_Macros") double* xtics {
    $1 = (double *)array_data($input);
}
%typecheck(SWIG_TYPECHECK_POINTER, fragment="NumPy_Fragments") double* ytics {
    $1 = is_array($input);
}
%typemap(in, fragment="NumPy_Macros") double* ytics {
    $1 = (double *)array_data($input);
}
%typecheck(SWIG_TYPECHECK_POINTER, fragment="NumPy_Fragments") double* ztics {
    $1 = is_array($input);
}
%typemap(in, fragment="NumPy_Macros") double* ztics {
    $1 = (double *)array_data($input);
}
%typecheck(SWIG_TYPECHECK_POINTER, fragment="NumPy_Fragments") double* weights {
    $1 = is_array($input);
}
%typemap(in, fragment="NumPy_Macros") double* weights {
    $1 = (double *)array_data($input);
}

// Typemap suite for array_slice

%typecheck(SWIG_TYPECHECK_POINTER, fragment="NumPy_Fragments") size_t dims[3] {
    $1 = is_array($input);
}

%typemap(in, fragment="NumPy_Macros") size_t dims[3] {
    $1 = (size_t *)array_data($input);
}

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
%apply double INPLACE_ARRAY1[ANY] { double [3] };

//--------------------------------------------------
// typemaps needed for get_eigenmode_coefficients
//--------------------------------------------------
%apply (int *IN_ARRAY1, int DIM1) {(int *bands, int num_bands)};

%typecheck(SWIG_TYPECHECK_POINTER, fragment="NumPy_Fragments") std::complex<double>* coeffs {
    $1 = is_array($input);
}

%typemap(in, fragment="NumPy_Macros") std::complex<double>* coeffs {
    $1 = (std::complex<double> *)array_data($input);
}

%typecheck(SWIG_TYPECHECK_POINTER, fragment="NumPy_Fragments") double* vgrp {
    $1 = is_array($input);
}

%typemap(in, fragment="NumPy_Macros") double* vgrp {
    $1 = (double *)array_data($input);
}

//--------------------------------------------------
// end typemaps for get_eigenmode_coefficients
//--------------------------------------------------

//--------------------------------------------------
// typemaps needed for add_dft_fields
//--------------------------------------------------

%typecheck(SWIG_TYPECHECK_POINTER) const volume where {
    int py_material = PyObject_IsInstance($input, py_volume_object());
    $1 = py_material;
}

%typecheck(SWIG_TYPECHECK_POINTER) meep::component *components {
    $1 = PyList_Check($input);
}

%typemap(in) (meep::component *components, int num_components) {
    if (!PyList_Check($input)) {
        meep::abort("Expected a list");
    }
    $2 = PyList_Size($input);
    $1 = new meep::component[$2];

    for (Py_ssize_t i = 0; i < $2; i++) {
        $1[i] = (meep::component)PyInteger_AsLong(PyList_GetItem($input, i));
    }
}

%typemap(freearg) (meep::component *components, int num_components) {
    delete[] $1;
}
//--------------------------------------------------
// end typemaps for add_dft_fields
//--------------------------------------------------

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
        meep::abort("Expected a sequence");
    }

    PyObject *cs = PyList_GetItem($input, 0);

    if (!PySequence_Check(cs)) {
        meep::abort("Expected first item in list to be a list");
    }

    PyObject *func = PyList_GetItem($input, 1);

    if (!PyCallable_Check(func)) {
        meep::abort("Expected a function");
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

// integrate2
%typecheck(SWIG_TYPECHECK_POINTER) (int num_fields1, const meep::component *components1, int num_fields2,
                                    const meep::component *components2, meep::field_function integrand,
                                    void *integrand_data_) {
    $1 = PySequence_Check($input) &&
         PySequence_Check(PyList_GetItem($input, 0)) &&
         PySequence_Check(PyList_GetItem($input, 1)) &&
         PyCallable_Check(PyList_GetItem($input, 2));
}

%typemap(in) (int num_fields1, const meep::component *components1, int num_fields2,
              const meep::component *components2, meep::field_function integrand,
              void *integrand_data_) (py_field_func_data data) {

    if (!PySequence_Check($input)) {
        meep::abort("Expected a sequence");
    }

    PyObject *cs1 = PyList_GetItem($input, 0);

    if (!PySequence_Check(cs1)) {
        meep::abort("Expected 1st item in list to be a sequence");
    }

    PyObject *cs2 = PyList_GetItem($input, 1);

    if (!PySequence_Check(cs2)) {
        meep::abort("Expected 2nd item in list to be a sequence");
    }

    PyObject *func = PyList_GetItem($input, 2);

    if (!PyCallable_Check(func)) {
        meep::abort("Expected 3rd item in list to be a function");
    }

    $1 = PyList_Size(cs1);
    $3 = PyList_Size(cs2);

    $2 = new meep::component[$1];
    $4 = new meep::component[$3];

    for (Py_ssize_t i = 0; i < $1; i++) {
        $2[i] = (meep::component)PyInteger_AsLong(PyList_GetItem(cs1, i));
    }

    for (Py_ssize_t i = 0; i < $3; i++) {
        $4[i] = (meep::component)PyInteger_AsLong(PyList_GetItem(cs2, i));
    }

    $5 = py_field_func_wrap;

    data.num_components = $1 + $3;
    data.func = func;
    Py_INCREF(func);
    $6 = &data;
}

%typemap(freearg) (int num_fields1, const meep::component *components1, int num_fields2,
                   const meep::component *components2, meep::field_function integrand, void *integrand_data_) {
    if ($2) {
        delete[] $2;
    }
    if ($4) {
        delete[] $4;
    }
    Py_XDECREF(data$argnum.func);
}

// Typemap suite for absorber_list

%typecheck(SWIG_TYPECHECK_POINTER) meep_geom::absorber_list {
    $1 = PySequence_Check($input);
}

%typemap(in) meep_geom::absorber_list {

    Py_ssize_t len = PyList_Size($input);

    if (len == 0) {
        $1 = 0;
    } else {
        $1 = create_absorber_list();

        for (Py_ssize_t i = 0; i < len; i++) {
            absorber a;
            PyObject *py_absorber = PyList_GetItem($input, i);

            if (!pyabsorber_to_absorber(py_absorber, &a)) {
                SWIG_fail;
            }

            add_absorbing_layer($1, a.thickness, a.direction, a.side,
                                a.R_asymptotic, a.mean_stretch, py_pml_profile,
                                a.pml_profile_data);
            Py_DECREF((PyObject *)a.pml_profile_data);
        }
    }
}

%typemap(freearg) meep_geom::absorber_list {
    if ($1) {
        destroy_absorber_list($1);
    }
}

// Typemap suite for material_type_list

%typecheck(SWIG_TYPECHECK_POINTER) material_type_list {
    $1 = PySequence_Check($input);
}

%typemap(in) material_type_list {
    Py_ssize_t len = PyList_Size($input);

    if (len == 0) {
        $1 = material_type_list();
    } else {
        material_type_list mtl;
        mtl.num_items = len;
        mtl.items = new material_type[len];
        for (Py_ssize_t i = 0; i < len; i++) {
            PyObject *py_material = PyList_GetItem($input, i);
            if (!pymaterial_to_material(py_material, &mtl.items[i])) {
                SWIG_fail;
            }
        }
        $1 = mtl;
    }
}

%typemap(freearg) material_type_list {
    if ($1.num_items != 0) {
        for (int i = 0; i < $1.num_items; i++) {
            if ($1.items[i]->medium.E_susceptibilities.items) {
                delete[] $1.items[i]->medium.E_susceptibilities.items;
            }
            if ($1.items[i]->medium.H_susceptibilities.items) {
                delete[] $1.items[i]->medium.H_susceptibilities.items;
            }
            delete[] $1.items[i]->epsilon_data;
        }
        delete[] $1.items;
    }
}

// For some reason SWIG needs the namespaced version too
%apply material_type_list { meep_geom::material_type_list };

// Typemap suite for custom_src_time

%typecheck(SWIG_TYPECHECK_POINTER) (std::complex<double> (*func)(double t, void *), void *data) {
    $1 = PyFunction_Check($input);
}

%typemap(in) (std::complex<double> (*func)(double t, void *), void *data) {
  $1 = py_src_func_wrap;
  $2 = (void *)$input;
}

// Typemap suite for kpoint_func

%typecheck(SWIG_TYPECHECK_POINTER) (meep::kpoint_func user_kpoint_func, void *user_kpoint_data) {
    $1 = PyFunction_Check($input) || $input == Py_None;
}

%typemap(in) (meep::kpoint_func user_kpoint_func, void *user_kpoint_data) {
    if ($input == Py_None) {
        $1 = NULL;
        $2 = NULL;
    }
    else {
        $1 = py_kpoint_func_wrap;
        $2 = (void*)$input;
    }
}

%apply double *flux {
    double *electric,
    double *magnetic,
    double *total
};

%exception {
  try {
    $action
  } catch (std::runtime_error &e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    SWIG_fail;
  }
}

// Tells Python to take ownership of the h5file* this function returns so that
// it gets garbage collected and the file gets closed.
%newobject meep::fields::open_h5file;

%newobject _get_eigenmode;

%rename(_vec) meep::vec::vec;
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

%ignore susceptibility_equal;
%ignore susceptibility_list_equal;
%ignore medium_struct_equal;
%ignore material_gc;
%ignore material_type_equal;
%ignore is_variable;
%ignore is_variable;
%ignore is_file;
%ignore is_file;
%ignore is_medium;
%ignore is_medium;
%ignore is_metal;
%ignore meep::infinity;

%ignore std::vector<meep::volume>::vector(size_type);
%ignore std::vector<meep::volume>::resize;
%ignore std::vector<meep_geom::dft_data>::vector(size_type);
%ignore std::vector<meep_geom::dft_data>::resize;

// template instantiations
%template(get_dft_flux_array) _get_dft_array<meep::dft_flux>;
%template(get_dft_fields_array) _get_dft_array<meep::dft_fields>;
%template(get_dft_force_array) _get_dft_array<meep::dft_force>;
%template(get_dft_near2far_array) _get_dft_array<meep::dft_near2far>;
%template(FragmentStatsVector) std::vector<meep_geom::fragment_stats>;
%template(DftDataVector) std::vector<meep_geom::dft_data>;
%template(VolumeVector) std::vector<meep::volume>;
%template(GridVolumeVector) std::vector<meep::grid_volume>;
%template(IntVector) std::vector<int>;
%template(DoubleVector) std::vector<double>;

%include "vec.i"
%include "meep.hpp"
%include "meep/mympi.hpp"
%include "meepgeom.hpp"

struct vector3 {
    double x;
    double y;
    double z;
};

struct geom_box {
    vector3 low;
    vector3 high;
};

%rename(is_point_in_object) point_in_objectp(vector3 p, GEOMETRIC_OBJECT o);
%rename(is_point_in_periodic_object) point_in_periodic_objectp(vector3 p, GEOMETRIC_OBJECT o);

#ifdef HAVE_MPB
namespace meep {
struct eigenmode_data {
    maxwell_data *mdata;
    scalar_complex *fft_data_H, *fft_data_E;
    evectmatrix H;
    int n[3];
    double s[3];
    double Gk[3];
    vec center;
    amplitude_function amp_func;
    int band_num;
    double omega;
    double group_velocity;
};
}

meep::eigenmode_data *_get_eigenmode(meep::fields *f, double omega_src, meep::direction d, const meep::volume where,
                                     const meep::volume eig_vol, int band_num, const meep::vec &_kpoint,
                                     bool match_frequency, int parity, double resolution, double eigensolver_tol,
                                     double kdom[3]);
PyObject *_get_eigenmode_Gk(meep::eigenmode_data *emdata);

%extend meep::eigenmode_data {
    ~eigenmode_data() {
        meep::destroy_eigenmode_data($self);
    }
}

#else
void _get_eigenmode(meep::fields *f, double omega_src, meep::direction d, const meep::volume where,
                    const meep::volume eig_vol, int band_num, const meep::vec &_kpoint,
                    bool match_frequency, int parity, double resolution, double eigensolver_tol,
                    double kdom[3]);
#endif // HAVE_MPB

// Make omega members of meep::dft_ldos available as 'freq' in python
%extend meep::dft_ldos {

    double get_omega_min() {
        return $self->omega_min;
    }
    double get_domega() {
        return $self->domega;
    }
    int get_Nomega() {
        return $self->Nomega;
    }
    %pythoncode %{
        def freqs(self):
            import math
            import numpy as np
            start = self.omega_min / (2 * math.pi)
            stop = start + (self.domega / (2 * math.pi)) * self.Nomega
            return np.linspace(start, stop, num=self.Nomega, endpoint=False).tolist()

        __swig_getmethods__["freq_min"] = get_omega_min
        __swig_getmethods__["nfreq"] = get_Nomega
        __swig_getmethods__["dfreq"] = get_domega
        if _newclass: freq_min = property(get_omega_min)
        if _newclass: nfreq = property(get_Nomega)
        if _newclass: dfreq = property(get_domega)
    %}
}

%extend meep::fields {
  bool is_periodic(boundary_side side, direction dir) {
    return $self->boundaries[side][dir] == meep::Periodic;
  }
}

extern boolean point_in_objectp(vector3 p, GEOMETRIC_OBJECT o);
extern boolean point_in_periodic_objectp(vector3 p, GEOMETRIC_OBJECT o);
void display_geometric_object_info(int indentby, GEOMETRIC_OBJECT o);
kpoint_list get_eigenmode_coefficients_and_kpoints(meep::fields *f, meep::dft_flux flux,
                                                   const meep::volume &eig_vol, int *bands, int num_bands,
                                                   int parity, double eig_resolution, double eigensolver_tol,
                                                   std::complex<double> *coeffs, double *vgrp,
                                                   meep::kpoint_func user_kpoint_func, void *user_kpoint_data,
                                                   meep::direction d);
PyObject *_get_array_slice_dimensions(meep::fields *f, const meep::volume &where, size_t dims[3],
                                      bool collapse_empty_dimensions, bool snap_empty_dimensions);

%ignore eps_func;
%ignore inveps_func;

%pythoncode %{
    AUTOMATIC = -1
    CYLINDRICAL = -2
    ALL = -1
    ALL_COMPONENTS = Dielectric

    # MPB definitions
    NO_PARITY = 0
    EVEN_Z = 1
    ODD_Z = 2
    EVEN_Y = 4
    ODD_Y = 8
    TE = EVEN_Z
    TM = ODD_Z
    PREV_PARITY = -1

    inf = 1.0e20

    from .geom import (
        Block,
        Cone,
        Cylinder,
        DrudeSusceptibility,
        Ellipsoid,
        FreqRange,
        GeometricObject,
        GyrotropicDrudeSusceptibility,
        GyrotropicLorentzianSusceptibility,
        GyrotropicSaturatedSusceptibility,
        Lattice,
        LorentzianSusceptibility,
        Matrix,
        Medium,
        MultilevelAtom,
        NoisyDrudeSusceptibility,
        NoisyLorentzianSusceptibility,
        Prism,
        Sphere,
        Susceptibility,
        Transition,
        Vector3,
        Wedge,
        check_nonnegative,
        geometric_object_duplicates,
        geometric_objects_duplicates,
        geometric_objects_lattice_duplicates,
        cartesian_to_lattice,
        lattice_to_cartesian,
        lattice_to_reciprocal,
        reciprocal_to_lattice,
        cartesian_to_reciprocal,
        reciprocal_to_cartesian,
        find_root_deriv,
        get_rotation_matrix,
    )
    from .simulation import (
        Absorber,
        Ldos,
        EnergyRegion,
        FluxRegion,
        ForceRegion,
        Harminv,
        Identity,
        Mirror,
        ModeRegion,
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
        before_time,
        complexarray,
        dft_ldos,
        display_progress,
        during_sources,
        GDSII_layers,
        GDSII_vol,
        get_center_and_size,
        get_eigenmode_freqs,
        get_electric_energy,
        get_energy_freqs,
        get_flux_freqs,
        get_fluxes,
        get_force_freqs,
        get_forces,
        get_magnetic_energy,
        get_near2far_freqs,
        get_total_energy,
        in_point,
        in_volume,
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
        quiet,
        scale_energy_fields,
        scale_flux_fields,
        scale_force_fields,
        scale_near2far_fields,
        stop_after_walltime,
        stop_on_interrupt,
        stop_when_fields_decayed,
        synchronized_magnetic,
        to_appended,
        vec,
        verbosity,
        when_true,
        when_false,
        with_prefix
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
    from .visualization import (
        plot2D,
        plot3D,
        plot_fields,
        Animate2D
    )

    if with_mpi():
        try:
            from mpi4py import MPI
        except ImportError as e:
            print('\n**\n** failed to load python MPI module (mpi4py)\n**', e, '\n**\n')
            pass
        else:
            # this variable reference is needed for lazy initialization of MPI
            comm = MPI.COMM_WORLD
            if am_master():
                Procs=comm.Get_size()
                (Major,Minor)=MPI.Get_version();
                print('Using MPI version {}.{}, {} processes'.format(Major, Minor, Procs));

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
    _t_start = wall_time()

    def report_elapsed_time():
        print("\nElapsed run time = {:.4f} s".format(wall_time() - _t_start))

    import atexit
    atexit.register(report_elapsed_time)
%}

%newobject create_structure_and_set_materials;
%inline %{

size_t get_realnum_size() {
  return sizeof(meep::realnum);
}

meep::structure *create_structure_and_set_materials(vector3 cell_size,
                                                    std::vector<meep_geom::dft_data> dft_data_list_,
                                                    std::vector<meep::volume> pml_1d_vols_,
                                                    std::vector<meep::volume> pml_2d_vols_,
                                                    std::vector<meep::volume> pml_3d_vols_,
                                                    std::vector<meep::volume> absorber_vols_,
                                                    meep::grid_volume &gv,
                                                    const meep::boundary_region &br,
                                                    const meep::symmetry &sym,
                                                    int num_chunks,
                                                    double Courant,
                                                    bool use_anisotropic_averaging,
                                                    double tol,
                                                    int maxeval,
                                                    geometric_object_list gobj_list,
                                                    vector3 center,
                                                    bool _ensure_periodicity,
                                                    meep_geom::material_type _default_material,
                                                    meep_geom::absorber_list alist,
                                                    meep_geom::material_type_list extra_materials,
                                                    bool split_chunks_evenly,
                                                    bool set_materials,
                                                    meep::structure *existing_s) {
    // Initialize fragment_stats static members (used for creating chunks in choose_chunkdivision)
    meep_geom::fragment_stats::geom = gobj_list;
    meep_geom::fragment_stats::dft_data_list = dft_data_list_;
    meep_geom::fragment_stats::pml_1d_vols = pml_1d_vols_;
    meep_geom::fragment_stats::pml_2d_vols = pml_2d_vols_;
    meep_geom::fragment_stats::pml_3d_vols = pml_3d_vols_;
    meep_geom::fragment_stats::absorber_vols = absorber_vols_;
    meep_geom::fragment_stats::tol = tol;
    meep_geom::fragment_stats::maxeval = maxeval;
    meep_geom::fragment_stats::resolution = gv.a;
    meep_geom::fragment_stats::dims = gv.dim;
    meep_geom::fragment_stats::split_chunks_evenly = split_chunks_evenly;
    meep_geom::fragment_stats::init_libctl(_default_material, _ensure_periodicity,
                                           &gv, cell_size, center, &gobj_list);
    meep::structure *s;
    if (existing_s) {
      s = existing_s;
    }
    else {
      s = new meep::structure(gv, NULL, br, sym, num_chunks, Courant,
                              use_anisotropic_averaging, tol, maxeval);
    }
    s->shared_chunks = true;

    if (set_materials) {
      meep_geom::set_materials_from_geometry(s, gobj_list, center, use_anisotropic_averaging, tol,
                                             maxeval, _ensure_periodicity, _default_material,
                                             alist, extra_materials);
    }

    if (meep::verbosity > 1 && !split_chunks_evenly && set_materials) {
      int num_procs = meep::count_processors();
      double *costs = new double[num_procs];
      for (int i = 0; i < num_procs; i++)
        costs[i] = 0;
      for (int i = 0; i < s->num_chunks; i++)
        costs[s->chunks[i]->n_proc()] += s->chunks[i]->gv.get_cost();
      double sum = 0, sumsq = 0;
      master_printf("estimated costs per process: ");
      for (int i = 0; i < num_procs; i++) {
        double cost = costs[i];
        sum += cost;
        sumsq += cost*cost;
        master_printf("%g%s", cost, i == num_procs - 1 ? "\n" : ", ");
      }
      delete[] costs;
      double mean = sum / num_procs;
      double stddev = sumsq - num_procs * mean * mean;
      stddev = num_procs == 1 || stddev <= 0 ? 0.0 : sqrt(stddev / (num_procs - 1));
      master_printf("estimated cost mean = %g, stddev = %g\n", mean, stddev);
    }

    // Return params to default state
    meep_geom::fragment_stats::resolution = 0;
    meep_geom::fragment_stats::split_chunks_evenly = false;

    return s;
}
%}
