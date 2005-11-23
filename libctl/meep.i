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
