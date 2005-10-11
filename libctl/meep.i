// -*- C++ -*-   
%module meep
%{
#include "meep-ctl.hpp"

static inline int SwigComplex_Check(SCM o) {
  return SCM_COMPLEXP(o);
}

/* Unfortunately, this is not re-entrant.  Damn dynamic scoping. 
   Hopefully, it should be good enough for our purposes. */
static SCM my_complex_func_scm;
static inline complex<double> my_complex_func(meep::vec const &v) {
  SCM ret = gh_call1(my_complex_func_scm, 
		     ctl_convert_vector3_to_scm(vec2vector3(v)));
  return std::complex<double>(gh_scm2double(scm_real_part(ret)),
			      gh_scm2double(scm_imag_part(ret)));
}

static inline complex<double> my_complex_func2(double t, void *f) {
  SCM ret = gh_call1((SCM) f, gh_double2scm(t));
  return std::complex<double>(gh_scm2double(scm_real_part(ret)),
			      gh_scm2double(scm_imag_part(ret)));
}
%}

%typecheck(SWIG_TYPECHECK_COMPLEX) complex<double> {
  $1 = SwigComplex_Check($input);
}

%typemap(guile,out) complex, complex<double>, std::complex<double> {
  $result = scm_make_rectangular(gh_double2scm($1.real()),
				 gh_double2scm($1.imag()) );
}
%typemap(guile,in) complex, complex<double>, std::complex<double> {
  $1 = std::complex<double>(gh_scm2double(scm_real_part($input)),
			    gh_scm2double(scm_imag_part($input)));
}

%typemap(guile,in) complex<double>(*)(meep::vec const &) {
  my_complex_func_scm = $input;
  $1 = my_complex_func;
}

%typemap(guile,in) complex<double>(*)(double, void*) {
  $1 = my_complex_func2; // the actual function had better be the next arg
}
%typemap(guile,in) void* { $1 = (void*) ($input); }

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
