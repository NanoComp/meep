# Process this file with autoconf to produce a configure script.
AC_INIT(meep, 0.2, meep@ab-initio.mit.edu)
AC_CONFIG_SRCDIR(src/step.cpp)

# Shared-library version number; indicates api compatibility, and is
# not the same as the "public" version number.  (Don't worry about this
# except for public releases.)
SHARED_VERSION_INFO="0:0:0"

AM_INIT_AUTOMAKE(1.6)
AM_CONFIG_HEADER(config.h)
AM_MAINTAINER_MODE
AC_SUBST(SHARED_VERSION_INFO)

DEFAULTFLAGS="no"
if test "x" == "x$CXXFLAGS"; then
   DEFAULTFLAGS=""
fi

# Checks for programs.
AC_PROG_CXX
AC_PROG_LIBTOOL

AC_CHECK_PROG(LATEX2HTML, latex2html, latex2html)
if test -z "$LATEX2HTML"; then
    AC_MSG_WARN([Cannot find latex2html in your path!])
    # FIXME: use standard 'missing' script from automake
    LATEX2HTML='echo not running latex2html...'
fi
AC_SUBST(LATEX2HTML)

AC_CHECK_LIB(m, sqrt)

AC_CHECK_LIB(fftw3, fftw_plan_dft_1d, [],
	[AC_CHECK_LIB(fftw, fftw_create_plan)])

###########################################################################

# Checks for BLAS/LAPACK libraries:

AC_F77_WRAPPERS

sinclude(m4/acx_blas.m4)
sinclude(m4/acx_lapack.m4)

ACX_BLAS
ACX_LAPACK

# just use lapack/blas in default libs since we always need them
LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

##############################################################################
# GNU Scientific Library

AC_CHECK_FUNC(cblas_cgemm, [], [AC_CHECK_LIB(gslcblas, cblas, cblas_cgemm)])

AC_CHECK_LIB(gsl, gsl_sf_bessel_Jn)

##############################################################################
# Check for MPI library

sinclude(m4/acx_mpi.m4)

AC_ARG_WITH(mpi, [AC_HELP_STRING([--with-mpi],[enable MPI parallelization])],
                  with_mpi=$withval, with_mpi=no)

if test "x$with_mpi" = "xyes"; then
	AC_LANG_PUSH([C++])
	ACX_MPI([],[AC_MSG_ERROR([could not find mpi library for --with-mpi])])
	CXX=$MPICXX
	AC_LANG_POP([C++])
fi

##############################################################################
# Figure out the number of processors so we can quickly build the docs and
# run the tests...  (FIXME: use more standard autoconf macros)

RUNCODE=""

if test "x$with_mpi" = "xyes"; then
        cat > test_num_processors.c <<MYEOF
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
  long nprocs;
  nprocs = sysconf(_SC_NPROCESSORS_ONLN);
  if (nprocs < 1)
    nprocs = 1;
  printf ("%ld\n",nprocs);
  exit (EXIT_SUCCESS);
}
MYEOF
        $CC -o test_num_processors test_num_processors.c
        NUMPROCS=`./test_num_processors`
        rm -f test_num_processors test_num_processors.c
        echo Looks like we have got "$NUMPROCS" processors

        RUNCODE="mpirun -np $NUMPROCS"
fi

AC_SUBST(RUNCODE)

##############################################################################

# FIXME: see how FFTW does it

# Try to guess some good optimization flags if CXXFLAGS is empty...
if test "x" == "x$DEFAULTFLAGS"; then
   CXXFLAGS_my_saved="-g"
   CFLAGS_my_saved="$CFLAGS"
   CFLAGS="$CXXFLAGS_my_saved -O3"
   AC_COMPILE_IFELSE([int main() { return 0; }],
       [echo Looks like -O3 is good && CXXFLAGS_my_saved="$CFLAGS"],
       [echo Looks like -O3 is bad  && CFLAGS="$CXXFLAGS_my_saved"])
   CFLAGS="$CFLAGS -funroll-loops"
   AC_COMPILE_IFELSE([int main() { return 0; }],
       [echo Looks like -funroll-loops is good && CXXFLAGS_my_saved="$CFLAGS"],
       [echo Looks like -funroll-loops is bad  && CFLAGS="$CXXFLAGS_my_saved"])
   CFLAGS="$CFLAGS -msse"
   AC_COMPILE_IFELSE([int main() { return 0; }],
       [echo Looks like -msse is good && CXXFLAGS_my_saved="$CFLAGS"],
       [echo Looks like -msse is bad  && CFLAGS="$CXXFLAGS_my_saved"])
   CFLAGS="$CFLAGS -mfpmath=sse"
   AC_COMPILE_IFELSE([int main() { return 0; }],
       [echo Looks like -mfpmath=sse is good && CXXFLAGS_my_saved="$CFLAGS"],
       [echo Looks like -mfpmath=sse is bad  && CFLAGS="$CXXFLAGS_my_saved"])
   CFLAGS="$CFLAGS -mcpu=pentium4"
   AC_COMPILE_IFELSE([int main() { return 0; }],
       [echo Looks like -mcpu=pentium4 is good && CXXFLAGS_my_saved="$CFLAGS"],
       [echo Looks like -mcpu=pentium4 is bad  && CFLAGS="$CXXFLAGS_my_saved"])

   CXXFLAGS="$CFLAGS"
   CFLAGS="$CFLAGS_my_saved"
fi

##############################################################################

AC_ARG_WITH(timings, [AC_HELP_STRING([--with-timings],[calculate timings of various operations])],with_timings=$withval,with_timings=no)

AC_MSG_CHECKING([calculate timings of various operations])
AC_MSG_RESULT($with_timings)
if test "x$with_timings" = "xyes"; then
    CPPFLAGS="$CPPFLAGS -DWITH_TIMINGS"
fi

##############################################################################

AC_CHECK_PROG(HC, ghc, ghc)

##############################################################################
# Libraries and flags (other than -lmeep) required to link Meep:

MEEPLIBS="$LDFLAGS $LIBS"
AC_SUBST(MEEPLIBS)

##############################################################################

AC_CONFIG_FILES([
	Makefile
	meep.pc
	src/Makefile
	hsrc/Makefile
	doc/Makefile
	tests/Makefile
])

AC_OUTPUT
