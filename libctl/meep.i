// -*- C++ -*-   
%module meep
%{
#include "meep-ctl.hpp"

static inline int SwigComplex_Check(SCM o) {
  return SCM_REALP(o) || SCM_COMPLEXP(o);
}

static inline int SwigVector3_Check(SCM o) {
  return SCM_VECTORP(o) && SCM_VECTOR_LENGTH(o) == 3;
}

/* Unfortunately, this is not re-entrant.  Damn dynamic scoping. 
   Hopefully, it should be good enough for our purposes. */
static SCM my_complex_func_scm;
static inline complex<double> my_complex_func(meep::vec const &v) {
  SCM ret = gh_call1(my_complex_func_scm, 
		     ctl_convert_vector3_to_scm(vec_to_vector3(v)));
  cnumber cret = ctl_convert_cnumber_to_c(ret);
  return std::complex<double>(cret.re, cret.im);
}

static inline complex<double> my_complex_func2(double t, void *f) {
  SCM ret = gh_call1((SCM) f, gh_double2scm(t));
  cnumber cret = ctl_convert_cnumber_to_c(ret);
  return std::complex<double>(cret.re, cret.im);
}

typedef struct { SCM func; int nf; } my_field_func_data;
static inline complex<double> my_field_func(const complex<double> *fields,
					    const meep::vec &loc,
					    void *data_) {
  my_field_func_data *data = (my_field_func_data *) data_;
  int num_items = data->nf;
  cnumber *items = new cnumber[num_items];
  for (int i = 0; i < num_items; ++i)
    items[i] = make_cnumber(real(fields[i]), imag(fields[i]));
  SCM ret = gh_apply(data->func,
		     scm_cons(ctl_convert_vector3_to_scm(vec_to_vector3(loc)),
			      make_cnumber_list(num_items, items)));
  delete[] items;
  cnumber cret = ctl_convert_cnumber_to_c(ret);
  return std::complex<double>(cret.re, cret.im);
}

%}

%typecheck(SWIG_TYPECHECK_COMPLEX) complex<double> {
  $1 = SwigComplex_Check($input);
}

%typemap(guile,out) complex, complex<double>, std::complex<double> {
  $result = scm_make_rectangular(gh_double2scm($1.real()),
				 gh_double2scm($1.imag()));
}
%typemap(guile,in) complex, complex<double>, std::complex<double> {
  cnumber cnum = ctl_convert_cnumber_to_c($input);
  $1 = std::complex<double>(cnum.re, cnum.im);
}

%typemap(guile,in) complex<double>(*)(meep::vec const &) {
  my_complex_func_scm = $input;
  $1 = my_complex_func;
}
%typecheck(SWIG_TYPECHECK_POINTER) complex<double>(*)(meep::vec const &) {
  $1 = SCM_NFALSEP(scm_procedure_p($input));
}

%typemap(guile,in) complex<double>(*)(double, void*) {
  $1 = my_complex_func2; // the actual function had better be the next arg
}
%typemap(guile,in) void* { $1 = (void*) ($input); }

%typemap(guile,in) meep::vec { 
  $1 = vector3_to_vec(ctl_convert_vector3_to_c($input));
}
%typemap(guile,out) meep::vec { 
  $result = ctl_convert_vector3_to_scm(vec_to_vector3($1));
}
%typemap(guile,in) meep::vec const & %{ 
  meep::vec vec__$1 = vector3_to_vec(ctl_convert_vector3_to_c($input));
  $1 = &vec__$1;
%}
%typecheck(SWIG_TYPECHECK_COMPLEX) meep::vec, meep::vec const & {
  $1 = SwigVector3_Check($input);
}

/* field_function arguments are passed as a cons pair of (components . func)
   in order to set all four arguments at once. */
%typemap(guile,in) (int num_fields, const meep::component *components, meep::field_function fun, void *fun_data_) (my_field_func_data data) {
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

// Need to tell SWIG about any method that returns a new object
// which needs to be garbage-collected.
%newobject meep::fields::open_h5file;
%newobject *::clone;
%newobject meep::dft_flux::flux;

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
