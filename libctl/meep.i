%module meep
%{
#include "meep.hpp"
using namespace meep;

static inline int
SwigComplex_Check(SCM o)
{
  return SCM_COMPLEXP(o);
}
%}

%typecheck(SWIG_TYPECHECK_COMPLEX) complex<double> {
  $1 = SwigComplex_Check($input);
}

%typemap(guile,out) complex, complex<double>, std::complex<double> {
  $result = scm_make_rectangular( gh_double2scm ($1.real ()),
           gh_double2scm ($1.imag ()) );
}
%typemap(guile,in) complex, complex<double>, std::complex<double> {
  $1 = std::complex<double>( gh_scm2double (scm_real_part ($input)),
           gh_scm2double (scm_imag_part ($input)) );
}

%include "meep_renames.i"

%include "meep/vec.hpp"
%include "meep/mympi.hpp"
%include "meep.hpp"

%include "ctl-io.i"

%{
#include "meep-ctl-swig.hpp"
%}

%include "meep-ctl-swig.hpp"
