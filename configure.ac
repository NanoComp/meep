# Process this file with autoconf to produce a configure script.
AC_INIT(meep, 0.9, meep@ab-initio.mit.edu)
AC_CONFIG_SRCDIR(src/step.cpp)

# Shared-library version number; indicates api compatibility, and is
# not the same as the "public" version number.  (Don't worry about this
# except for public releases.)
SHARED_VERSION_INFO="0:0:0"

AM_INIT_AUTOMAKE(1.6)
AM_CONFIG_HEADER(config.h)
AM_MAINTAINER_MODE
AC_SUBST(SHARED_VERSION_INFO)
AM_ENABLE_SHARED(no) dnl shared libs cause too many headaches to be default

MEEP_SUFFIX=""

##############################################################################

AC_ARG_WITH(doc,
	[AC_HELP_STRING([--with-doc], [include old documentation in dist])],
	with_doc=$withval,with_doc=no)

AM_CONDITIONAL(WITH_DOC, test "x$with_doc" = xyes)

##############################################################################

# Check for mpiCC immediately after getting C++ compiler...
AC_PROG_CXX

# Check for MPI library

AC_ARG_WITH(mpi, [AC_HELP_STRING([--with-mpi],[enable MPI parallelization])],
                  with_mpi=$withval, with_mpi=no)

if test "x$with_mpi" = "xyes"; then
	AC_LANG_PUSH([C++])
	ACX_MPI([],[AC_MSG_ERROR([could not find mpi library for --with-mpi])])
	CXX=$MPICXX
	AC_LANG_POP([C++])
	MEEP_SUFFIX="${MEEP_SUFFIX}_mpi"
fi

##############################################################################
# More checks

AC_PROG_LIBTOOL

AC_CHECK_PROG(LATEX2HTML, latex2html, latex2html)
if test -z "$LATEX2HTML"; then
    AC_MSG_WARN([Cannot find latex2html in your path!])
    # FIXME: use standard 'missing' script from automake
    LATEX2HTML='echo not running latex2html...'
fi
AC_SUBST(LATEX2HTML)

AC_CHECK_LIB(m, sin)

AC_CHECK_LIB(fftw3, fftw_plan_dft_1d, [],
  [AC_CHECK_LIB(dfftw, fftw_create_plan, [],
    [AC_CHECK_LIB(fftw, fftw_create_plan, [],
       [AC_MSG_WARN([Missing FFTW library...slow DFT will be used for monitor-point Fourier transform])])])])

###########################################################################
# Harminv library

# First, try pkg-config, if it is installed.  Note that, annoyingly,
# pkg-config doesn't look in /usr/local by default, so we have to
# add this to the path for that common case.  It also doesn't give
# us a way to print its default path, grr, so we have to assume that this
# is /usr/lib/pkgconfig.
PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/usr/lib/pkgconfig:/usr/local/lib/pkgconfig
export PKG_CONFIG_PATH
PKG_CHECK_MODULES(HARMINV, harminv >= 1.1, [CFLAGS="$CFLAGS $HARMINV_CFLAGS"; LIBS="$HARMINV_LIBS $LIBS"; have_harminv=yes], [have_harminv=no])

# As a fallback, check manually for BLAS/LAPACK and harminv libraries:
if test $have_harminv = no; then
	AC_F77_WRAPPERS
	ACX_BLAS
	ACX_LAPACK([], [AC_MSG_WARN([BLAS/LAPACK needed for harminv])])

        save_CC=$CC
	CC=$CXX # harminv test must use $CXX since harminv may be a C++ lib

	# just use lapack/blas in default libs since we always need them
	AC_CHECK_LIB(harminv, harminv_get_freq_error, [have_harminv=yes; LIBS="-lharminv $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"], [AC_MSG_WARN([harminv support is disabled])], [$LAPACK_LIBS $BLAS_LIBS $FLIBS])

	CC=$save_CC
fi

if test $have_harminv = yes; then
	AC_DEFINE([HAVE_HARMINV],[1],[Define if you have libharminv])
fi

##############################################################################
# GNU Scientific Library

AC_CHECK_FUNC(cblas_cgemm, [], [AC_CHECK_LIB(gslcblas, cblas_cgemm)])

AC_CHECK_LIB(gsl, gsl_sf_bessel_Jn, [],
	[AC_MSG_WARN([Missing GNU GSL library...Bessel-function field initialization will not be supported.])])

##############################################################################
# Check for HDF5 library
AC_ARG_WITH(hdf5, [AC_HELP_STRING([--without-hdf5],[do not allow HDF5 output])], ok=$withval, ok=yes)
if test "$ok" = "yes"; then
	AC_CHECK_LIB(z,deflate, [],[AC_MSG_WARN([zlib is required for HDF5!])])

	save_CC=$CC
	if test "x$with_mpi" = "xyes"; then
	   CC=$MPICXX
	fi
	save_LIBS_0="$LIBS" # save, to check later if we found any library
	AC_CHECK_LIB(hdf5, H5Pcreate, [
		AC_CHECK_HEADERS(hdf5.h, [LIBS="-lhdf5 $LIBS"
                     AC_DEFINE(HAVE_HDF5,1,[Define if we have & link HDF5])])])
	if test x"$save_LIBS_0" = x"$LIBS"; then
		AC_MSG_WARN([Couldn't find the HDF5 library!!  Switching to --without-hdf5.])
	fi
	CC=$save_CC

        if test "x$with_mpi" = "xyes"; then
		AC_CHECK_FUNCS(H5Pset_mpi H5Pset_fapl_mpio)
	fi
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
# Compiler flags

AC_ARG_ENABLE(debug, [AC_HELP_STRING([--enable-debug,compile for debugging])],
	      enable_debug=$enableval, enable_debug=no)
if test "$enable_debug" = "yes"; then
	CFLAGS="-g"
	CXXFLAGS="-g"
	FFLAGS="-g"
	AC_DEFINE(DEBUG,1,[define to enable debugging code])
else
  # Pick "good" compiler flags(?)
  AX_CXX_MAXOPT
fi

# Add lots of compiler warnings in maintainer mode if we are using gcc:
# (The variable $GXX is set to "yes" by AC_PROG_CXX if we are using g++.)
if test "$GXX" = "yes" && test "$USE_MAINTAINER_MODE" = yes; then
        CXXFLAGS="$CXXFLAGS -Wall -W"
fi

##############################################################################

AC_CHECK_PROG(HC, ghc, ghc)
if test x"$HC" = x && test "$USE_MAINTAINER_MODE" = yes; then
	AC_MSG_ERROR([the ghc Haskell compiler is required in maintainer mode])
fi

##############################################################################
# Libraries and flags (other than -lmeep) required to link Meep:

MEEPLIBS="$LDFLAGS $LIBS"
AC_SUBST(MEEPLIBS)

##############################################################################

# subsequent libraries are only for libctl front end
meep_save_LIBS=$LIBS

AC_ARG_WITH(libctl, [AC_HELP_STRING([--with-libctl=<dir>],[specify libctl directory])], with_libctl=$withval, with_libctl=yes)
if test "x$with_libctl" = "xno"; then
        :
elif test "x$with_libctl" != "xyes"; then
	LIBCTL_DIR="$with_libctl"
fi

AM_CONDITIONAL(WITH_LIBCTL, test "x$with_libctl" != xno)

##############################################################################
# Check for Guile library and its behavior (for libctl front end):

if test "x$with_libctl" != xno; then

AC_CHECK_PROG(HAVE_GUILE_CONFIG, guile-config, yes, no)
if test "$HAVE_GUILE_CONFIG" = "yes"; then
        CPPFLAGS="$CPPFLAGS `guile-config compile`"
        LIBS="`guile-config link` $LIBS"
	AC_MSG_CHECKING([if linking to guile works])
	AC_TRY_LINK_FUNC(gh_enter, AC_MSG_RESULT(yes), [
	AC_MSG_RESULT(no)
	AC_MSG_ERROR(guile-config is broken)
	])
else
        AC_CHECK_LIB(readline, readline)
        AC_CHECK_LIB(dl, dlopen)
        AC_CHECK_LIB(guile, gh_eval_str)
fi

# Check how smob types work in this Guile version:

AC_CHECK_FUNCS(scm_make_smob_type)
AC_MSG_CHECKING([for SCM_SMOB_PREDICATE])
AC_TRY_LINK([#include <guile/gh.h>
#include <libguile/smob.h>],
                    [int x; SCM_SMOB_PREDICATE(1,x);],
                    ok=yes, ok=no)
AC_MSG_RESULT($ok)
if test $ok = yes; then
	AC_DEFINE(HAVE_SCM_SMOB_PREDICATE, 1, [define if we have SCM_SMOB_PREDICATE])
fi
AC_MSG_CHECKING([for SCM_SMOB_DATA])
AC_TRY_LINK([#include <guile/gh.h>
#include <libguile/smob.h>],
                    [int x; SCM_SMOB_DATA(x);],
                    ok=yes, ok=no)
AC_MSG_RESULT($ok)
if test $ok = yes; then
	AC_DEFINE(HAVE_SCM_SMOB_DATA, 1, [define if we have SCM_SMOB_DATA])
fi
AC_MSG_CHECKING([for SCM_NEWSMOB])
AC_TRY_LINK([#include <guile/gh.h>
#include <libguile/smob.h>],
                    [int x; SCM_NEWSMOB(x,1,0);],
                    ok=yes, ok=no)
AC_MSG_RESULT($ok)
if test $ok = yes; then
	AC_DEFINE(HAVE_SCM_NEWSMOB, 1, [define if we have SCM_NEWSMOB])
fi

# Check how to activate Guile readline support:

ACTIVATE_READLINE=""
AC_MSG_CHECKING(how to activate readline in Guile)
ractivate="(use-modules (readline-activator)) (activate-readline)"
if guile -c "$ractivate" > /dev/null 2>&1; then
	AC_MSG_RESULT(readline-activator)
	ACTIVATE_READLINE="$ractivate"
else
	ractivate="(use-modules (ice-9 readline)) (activate-readline)"
	if guile -c "$ractivate" >& /dev/null; then
	        AC_MSG_RESULT(ice-9 readline)
	        ACTIVATE_READLINE="$ractivate"
	else
		AC_MSG_RESULT(cannot)
		ACTIVATE_READLINE=""
	fi
fi
AC_SUBST(ACTIVATE_READLINE)

fi # if with_libctl

##############################################################################
# Check for libctl library and files

if test "x$with_libctl" != xno; then

AC_MSG_CHECKING([for libctl dir])
if test x != x"$LIBCTL_DIR" -a ! -r "$LIBCTL_DIR/base/ctl.scm"; then
	LIBCTL_DIR=""
fi
if test x = x"$LIBCTL_DIR" -a -r /usr/local/share/libctl/base/ctl.scm; then
	LIBCTL_DIR="/usr/local/share/libctl"
fi
if test x = x"$LIBCTL_DIR" -a -r /usr/share/libctl/base/ctl.scm; then
        LIBCTL_DIR="/usr/share/libctl"
fi
if test x = x"$LIBCTL_DIR"; then
	AC_MSG_ERROR([could not find libctl files; use --with-libctl=<dir>])
fi
AC_MSG_RESULT($LIBCTL_DIR)
AC_SUBST(LIBCTL_DIR)

AC_CHECK_PROGS(GEN_CTL_IO, gen-ctl-io, echo)
if test x"$GEN_CTL_IO" = xecho; then
	AC_MSG_ERROR([could not find gen-ctl-io program; check your PATH])
fi

# check for -lctl:
AC_CHECK_LIB(ctl, ctl_get_vector3, [], [AC_MSG_ERROR([Couldn't find the required libctl library.])])

AC_CHECK_HEADER(ctl.h, [], [AC_MSG_ERROR([Couldn't find the libctl header.])])

# Check libctl version >= LIBCTL_MAJOR.LIBCTL_MINOR.LIBCTL_BUGFIX
LIBCTL_MAJOR=3; LIBCTL_MINOR=0; LIBCTL_BUGFIX=0
AC_MSG_CHECKING([whether libctl version is at least ${LIBCTL_MAJOR}.${LIBCTL_MINOR}.${LIBCTL_BUGFIX}])
AC_EGREP_CPP(yes, [[
#include <ctl.h>
#if LIBCTL_MAJOR_VERSION > $LIBCTL_MAJOR || (LIBCTL_MAJOR_VERSION == $LIBCTL_MAJOR && (LIBCTL_MINOR_VERSION > $LIBCTL_MINOR || (LIBCTL_MINOR_VERSION == $LIBCTL_MINOR && LIBCTL_BUGFIX_VERSION >= $LIBCTL_BUGFIX)))
  yes
#endif
]], [AC_MSG_RESULT(ok)], [AC_MSG_ERROR([libctl version ${LIBCTL_MAJOR}.${LIBCTL_MINOR}.${LIBCTL_BUGFIX} or later is required])])

fi # if with_libctl

##############################################################################

LIBCTL_LIBS=$LIBS
AC_SUBST(LIBCTL_LIBS)
LIBS=$meep_save_LIBS

##############################################################################
# The following function is used only for debugging.  Note that
# we must test for it *after* setting the compiler flags (which
# affect whether it is declared, as it is a GNU extension).
# We need to #include <stdio.h> because that somehow affects whether
# the function is declared with older gcc versions.  We need
# to use AC_TRY_COMPILE because the test in AC_HAVE_DECL seems
# to be optimized out.

AC_CHECK_FUNCS(feenableexcept) dnl GNU libc fp exception control function
AC_LANG_PUSH([C++])
AC_MSG_CHECKING([whether feenableexcept declaration is usable])
feenableexcept_decl_ok=yes
AC_TRY_COMPILE([#include <stdio.h>
#define _GNU_SOURCE 1
#include <fenv.h>], [feenableexcept(0);],
[AC_DEFINE([HAVE_DECL_FEENABLEEXCEPT],[1],[Define if fenv.h declares this.])],
[feenableexcept_decl_ok=no])
AC_MSG_RESULT($feenableexcept_decl_ok)
AC_LANG_POP([C++])

##############################################################################
# See if we need to catch SIGFPE to avoid crashing on underflow exceptions
# etcetera, since some stupid operating systems (e.g. Alpha/Tru64) crash on 
# SIGFPE by default.

AC_MSG_CHECKING([whether to catch and ignore SIGFPE signals])
ignore_sigfpe=yes
AC_RUN_IFELSE([AC_LANG_PROGRAM([], [volatile double x=1,y=0; return x/y;])], 
              [ignore_sigfpe=no])
AC_MSG_RESULT($ignore_sigfpe)
if test $ignore_sigfpe = yes; then
    AC_DEFINE(IGNORE_SIGFPE, [1], [Define to catch and ignore SIGFPE signals])
fi

##############################################################################
# Miscellaneous function and header checks

AC_HEADER_TIME
AC_CHECK_HEADERS([sys/time.h])
AC_CHECK_FUNCS([BSDgettimeofday gettimeofday cblas_ddot cblas_daxpy])

##############################################################################

AC_SUBST(MEEP_SUFFIX)
program_transform_name="s,_,-,g;$program_transform_name"

##############################################################################

AC_CONFIG_FILES([
	Makefile
	meep.pc
	src/Makefile
	hsrc/Makefile
	doc/Makefile
	tests/Makefile
	examples/Makefile
	libctl/Makefile
	libctl/meep.scm
])

AC_OUTPUT
