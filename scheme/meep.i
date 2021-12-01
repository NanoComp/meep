// -*- C++ -*-
%module meep
%{
#include "meep-ctl.hpp"

#if SCM_MAJOR_VERSION >= 2
#  define SCM_VECTORP(o) scm_vector_p(o)
#  define SCM_VECTOR_LENGTH(o) scm_c_vector_length(o)
#endif

static inline int SwigComplex_Check(SCM o) {
  return SCM_REALP(o) || SCM_COMPLEXP(o);
}

static inline int SwigVector3_Check(SCM o) {
  return SCM_VECTORP(o) && SCM_VECTOR_LENGTH(o) == 3;
}

/* Unfortunately, this is not re-entrant.  Damn dynamic scoping.
   Hopefully, it should be good enough for our purposes. */
static SCM my_complex_func_scm;
static inline std::complex<double> my_complex_func(meep::vec const &v) {
  SCM ret = gh_call1(my_complex_func_scm,
		     ctl_convert_vector3_to_scm(vec_to_vector3(v)));
  cnumber cret = ctl_convert_cnumber_to_c(ret);
  return std::complex<double>(cret.re, cret.im);
}

static inline std::complex<double> my_complex_func2(double t, void *f) {
  SCM ret = gh_call1((SCM) f, ctl_convert_number_to_scm(t));
  cnumber cret = ctl_convert_cnumber_to_c(ret);
  return std::complex<double>(cret.re, cret.im);
}

typedef struct { SCM func; int nf; } my_field_func_data;
static inline std::complex<double> my_field_func(const std::complex<meep::realnum> *fields,
                                                 const meep::vec &loc,
                                                 void *data_) {
  my_field_func_data *data = (my_field_func_data *) data_;
  int num_items = data->nf;
  cnumber *items = new cnumber[num_items];
  for (int i = 0; i < num_items; ++i)
    items[i] = make_cnumber(real(fields[i]), imag(fields[i]));
  SCM ret = scm_apply_0(data->func,
		     scm_cons(ctl_convert_vector3_to_scm(vec_to_vector3(loc)),
			      make_cnumber_list(num_items, items)));
  delete[] items;
  cnumber cret = ctl_convert_cnumber_to_c(ret);
  return std::complex<double>(cret.re, cret.im);
}

/* Unfortunately, this is not re-entrant.  Damn dynamic scoping.
   Hopefully, it should be good enough for our purposes. */
static SCM my_complex_func3_scm;
static inline std::complex<double> my_complex_func3(std::complex<double> x) {
  cnumber cx;
  cx.re = real(x); cx.im = imag(x);
  SCM ret = gh_call1(my_complex_func3_scm, ctl_convert_cnumber_to_scm(cx));
  cnumber cret = ctl_convert_cnumber_to_c(ret);
  return std::complex<double>(cret.re, cret.im);
}

static meep::vec my_kpoint_func(double freq, int mode, void *user_data) {
  SCM scm_res = gh_call2((SCM)user_data, ctl_convert_number_to_scm(freq),
                         ctl_convert_number_to_scm(mode));
  vector3 v3 = ctl_convert_vector3_to_c(scm_res);
  meep::vec result(v3.x, v3.y, v3.z);
  return result;
}

%}

%typecheck(SWIG_TYPECHECK_COMPLEX) std::complex<double> {
  $1 = SwigComplex_Check($input);
}

%typemap(out) complex, std::complex<double>, std::complex<double> {
  $result = scm_make_rectangular(ctl_convert_number_to_scm($1.real()),
				 ctl_convert_number_to_scm($1.imag()));
}
%typemap(in) complex, std::complex<double>, std::complex<double> {
  cnumber cnum = ctl_convert_cnumber_to_c($input);
  $1 = std::complex<double>(cnum.re, cnum.im);
}

%typemap(in) std::complex<double>(*)(meep::vec const &) {
  my_complex_func_scm = $input;
  $1 = my_complex_func;
}
%typecheck(SWIG_TYPECHECK_POINTER) std::complex<double>(*)(meep::vec const &) {
  $1 = SCM_NFALSEP(scm_procedure_p($input));
}

%typemap(in) std::complex<double>(*)(std::complex<double>) {
  my_complex_func3_scm = $input;
  $1 = my_complex_func3;
}
%typecheck(SWIG_TYPECHECK_POINTER) std::complex<double>(*)(std::complex<double>) {
  $1 = SCM_NFALSEP(scm_procedure_p($input));
}

%typemap(in) (std::complex<double> (*func)(double t, void *), void *data) {
  $1 = my_complex_func2;
  $2 = (void *) $input; // input is SCM pointer to Scheme function
}
%typecheck(SWIG_TYPECHECK_POINTER) (std::complex<double> (*func)(double t, void *), void *data) {
  $1 = SCM_NFALSEP(scm_procedure_p($input));
}

%typemap(in) meep::vec {
  $1 = vector3_to_vec(ctl_convert_vector3_to_c($input));
}
%typemap(out) meep::vec {
  $result = ctl_convert_vector3_to_scm(vec_to_vector3($1));
}
%typemap(in) meep::vec const & %{
  meep::vec vec__$1 = vector3_to_vec(ctl_convert_vector3_to_c($input));
  $1 = &vec__$1;
%}
%typecheck(SWIG_TYPECHECK_COMPLEX) meep::vec, meep::vec const & {
  $1 = SwigVector3_Check($input);
}

/* field_function arguments are passed as a cons pair of (components . func)
   in order to set all four arguments at once. */
%typemap(in) (int num_fields, const meep::component *components, meep::field_function fun, void *fun_data_) (my_field_func_data data) {
  $1 = list_length(gh_car($input));
  $2 = new meep::component[$1];
  for (int i = 0; i < $1; ++i)
    $2[i] = meep::component(integer_list_ref(gh_car($input), i));
  data.nf = $1;
  data.func = gh_cdr($input);
  $3 = my_field_func;
  $4 = &data;
}
%typemap(freearg) (int num_fields, const meep::component *components, meep::field_function fun, void *fun_data_) {
  if ($2) delete[] $2;
}
%typecheck(SWIG_TYPECHECK_POINTER) (int num_fields, const meep::component *components, meep::field_function fun, void *fun_data_) {
  $1 = SCM_NFALSEP(scm_pair_p($input)) &&
       SCM_NFALSEP(scm_list_p(gh_car($input))) &&
       SCM_NFALSEP(scm_procedure_p(gh_cdr($input)));
}

/* integrate2 arguments are passed as a cons pair of
   ((components1 . components2) . func)
   in order to set all six arguments at once. */
%typemap(in) (int num_fields1, const meep::component *components1, int num_fields2, const meep::component *components2, meep::field_function integrand, void *integrand_data_) (my_field_func_data data) {
  $1 = list_length(gh_car(gh_car($input)));
  $2 = new meep::component[$1];
  for (int i = 0; i < $1; ++i)
    $2[i] = meep::component(integer_list_ref(gh_car(gh_car($input)), i));
  $3 = list_length(gh_cdr(gh_car($input)));
  $4 = new meep::component[$3];
  for (int i = 0; i < $3; ++i)
    $4[i] = meep::component(integer_list_ref(gh_cdr(gh_car($input)), i));
  data.nf = $1 + $3;
  data.func = gh_cdr($input);
  $5 = my_field_func;
  $6 = &data;
}
%typemap(freearg) (int num_fields1, const meep::component *components1, int num_fields2, const meep::component *components2, meep::field_function integrand, void *integrand_data_) (my_field_func_data data) {
  if ($2) delete[] $2;
  if ($4) delete[] $4;
}
%typecheck(SWIG_TYPECHECK_POINTER) (int num_fields1, const meep::component *components1, int num_fields2, const meep::component *components2, meep::field_function integrand, void *integrand_data_) (my_field_func_data data) {
  $1 = SCM_NFALSEP(scm_pair_p($input)) &&
       SCM_NFALSEP(scm_pair_p(gh_car($input))) &&
       SCM_NFALSEP(scm_list_p(gh_car(gh_car($input)))) &&
       SCM_NFALSEP(scm_list_p(gh_cdr(gh_car($input)))) &&
       SCM_NFALSEP(scm_procedure_p(gh_cdr($input)));
}

%typecheck(SWIG_TYPECHECK_POINTER) (meep::component *components, int num_components) {
    $1 = SCM_NFALSEP(scm_list_p($input));
}

%typemap(in) (meep::component *components, int num_components) {
  $2 = list_length($input);
  $1 = new meep::component[$2];

  for (int i = 0; i < $2; ++i) {
    $1[i] = meep::component(integer_list_ref($input, i));
  }
}

%typemap(freearg) (meep::component *components, int num_components) {
  if ($1) {
    delete[] $1;
  }
}

// add_dft_fields

%typemap(in) (component *components, int num_components) {
  $2 = list_length($input);
  $1 = new meep::component[$2];

  for (int i = 0; i < $2; ++i) {
    $1[i] = integer_list_ref($input, i);
  }
}

%typemap(freearg) (component *components, int num_components) {
  if ($1) {
    delete[] $1;
  }
}

// do_get_eigenmode_coefficients

%typemap(in) (int *bands, int num_bands) {
  $2 = list_length($input);
  $1 = new int[$2];

  for (int i = 0; i < $2; ++i) {
    $1[i] = integer_list_ref($input, i);
  }
}

%typemap(freearg) (int *bands, int num_bands) {
  if ($1) {
    delete[] $1;
  }
}

%typemap(in) (meep::kpoint_func user_kpoint_func, void *user_kpoint_data) {
  if (SCM_NFALSEP(scm_procedure_p($input))) {
    $1 = my_kpoint_func;
    $2 = (void*)$input;
  }
  else {
    $1 = NULL;
    $2 = NULL;
  }
}

%typemap(in, noblock=1) std::complex<double> *coeffs {
  scm_t_array_handle coeffs_handle;
  scm_array_get_handle($input, &coeffs_handle);
  $1 = (std::complex<double>*)scm_array_handle_uniform_writable_elements(&coeffs_handle);
}

%typemap(in, noblock=1) double *vgrp {
  scm_t_array_handle vgrp_handle;
  scm_array_get_handle($input, &vgrp_handle);
  $1 = (double*)scm_array_handle_uniform_writable_elements(&vgrp_handle);
}

%typemap(freearg) std::complex<double> *coeffs {
  scm_array_handle_release(&coeffs_handle);
}

%typemap(freearg) double *vgrp {
  scm_array_handle_release(&vgrp_handle);
}

%typemap(out) kpoint_list {
  int i;
  list kpoint_list = SCM_EOL;
  list kdom_list = SCM_EOL;
  for (i = $1.n - 1; i >= 0; --i) {
    kpoint_list = gh_cons(vector32scm(vec_to_vector3($1.kpoints[i])), kpoint_list);
  }
  for (i = $1.num_bands - 1; i >= 0; --i) {
    kdom_list = gh_cons(vector32scm(vec_to_vector3($1.kdom[i])), kdom_list);
  }

  $result = scm_list_2(kpoint_list, kdom_list);

  delete[] $1.kpoints;
  delete[] $1.kdom;
}

// Need to tell SWIG about any method that returns a new object
// which needs to be garbage-collected.
%newobject meep::fields::open_h5file;
%newobject *::clone;
%newobject meep::dft_flux::flux;

%warnfilter(302,325,451,503,509);

%ignore meep::all_in_or_out;
%ignore meep::all_connect_phases;
%ignore meep::choose_chunkdivision;
%ignore meep::comms_key;
%ignore meep::comms_key_hash_fn;
%ignore meep::comms_manager;
%ignore meep::comms_operation;
%ignore meep::comms_sequence;
%ignore meep::create_comms_manager;
%ignore meep::fields_chunk;
%ignore meep_geom::fragment_stats;

%include "meep_renames.i"
%include "meep_enum_renames.i"
%include "meep_op_renames.i"
%include "meep_swig_bug_workaround.i"

%include "meep/vec.hpp"
%include "meep/mympi.hpp"
%include "meep.hpp"

%include "ctl-io.i"

%{
#include "meep-ctl-swig.hpp"
%}

%include "meep-ctl-swig.hpp"
