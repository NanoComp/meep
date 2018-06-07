# Process this file with autoconf to produce a configure script.
AC_INIT(meep, 1.5, meep@ab-initio.mit.edu)
AC_CONFIG_SRCDIR(src/step.cpp)

# Shared-library version number; indicates api compatibility, and is
# not the same as the "public" version number.  (Don't worry about this
# except for public releases.) Note that any change to a C++ class
# definition (in the .hpp file) generally breaks binary compatibility.
SHARED_VERSION_INFO="10:0:0" # CURRENT:REVISION:AGE

AM_INIT_AUTOMAKE([foreign color-tests parallel-tests silent-rules 1.11])
AM_SILENT_RULES(yes)
AM_CONFIG_HEADER(config.h)
AC_CONFIG_MACRO_DIR([m4])
AM_MAINTAINER_MODE
AC_SUBST(SHARED_VERSION_INFO)
# need shared libraries for python interfaces
AC_ENABLE_SHARED

##############################################################################

# Check for mpiCC immediately after getting C++ compiler...
AC_PROG_CXX
AC_LANG([C++])

# Check for MPI library

AC_ARG_WITH(mpi, [AC_HELP_STRING([--with-mpi],[enable MPI parallelization])],
                  with_mpi=$withval, with_mpi=no)

if test "x$with_mpi" = "xyes"; then
  ACX_MPI([],[AC_MSG_ERROR([could not find mpi library for --with-mpi])])
  CXX=$MPICXX

  # Detect stupid MPICH problem with SEEK_SET conflicting between
  # the mpi.h header file and the stdio.h header file in C++

  AC_MSG_CHECKING([for extra flag needed to combine stdio.h and mpi.h])
  for flag in none -DMPICH_IGNORE_CXX_SEEK unknown; do
    if test "x$flag" = xunknown; then break; fi
    save_CPPFLAGS=$CPPFLAGS
    if test "x$flag" != xnone; then CPPFLAGS="$CPPFLAGS $flag"; fi
    AC_TRY_COMPILE([#include <stdio.h>
#include <mpi.h>], [], [break], [])
    CPPFLAGS=$save_CPPFLAGS
  done
  AC_MSG_RESULT($flag)
  if test "x$flag" = xunknown; then
    AC_MSG_CHECKING([if we can combine stdio.h and mpi.h by undef'ing SEEK_*])
    AC_TRY_COMPILE([#include <stdio.h>
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#include <mpi.h>], [],
	 [AC_MSG_RESULT(yes)
          AC_DEFINE([NEED_UNDEF_SEEK_FOR_MPI], 1, [Define if mpi.h needs SEEK macros to be undefined])],
	 [AC_MSG_RESULT(no)
          AC_MSG_ERROR([can't include stdio.h and mpi.h in same file])])
  fi
fi

##############################################################################
# Compiler flags

# Pick "good" compiler flags(?)
AX_CXX_MAXOPT

AC_ARG_ENABLE(debug,
              [AC_HELP_STRING([--enable-debug],[compile for debugging])],
	      enable_debug=$enableval, enable_debug=no)
if test "$enable_debug" = "yes"; then
	CFLAGS="-g"
	CXXFLAGS="-g"
	FFLAGS="-g"
	AC_DEFINE(DEBUG,1,[define to enable debugging code])
fi

# Add lots of compiler warnings in maintainer mode if we are using gcc:
# (The variable $GXX is set to "yes" by AC_PROG_CXX if we are using g++.)
if test "$GXX" = "yes" && test "$USE_MAINTAINER_MODE" = yes; then
        CXXFLAGS="$CXXFLAGS -Wall -W"
fi

# For some annoying reason, g++ requires you to compile
# all code with -march if you compile any code with -march,
# otherwise segfaults can occur (observed in g++ 3.3.5).
# We also need to be careful about -malign-double.
ARCHFLAG=""
for flg in $CXXFLAGS; do
    case $flg in
      -march=*|-mcpu=*|-malign*) ARCHFLAG="$ARCHFLAG $flg" ;;
    esac
done
AC_SUBST(ARCHFLAG)

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
       [AC_MSG_WARN([FFTW needed for MPB])])])])

###########################################################################

AC_PROG_F77
AC_F77_WRAPPERS
if test "x$with_lapack" != xno -o "x$with_blas" != xno; then
   ACX_BLAS
   ACX_LAPACK
   LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
fi

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
        if test "x$acx_lapack_ok" = x; then
	   ACX_BLAS
	   ACX_LAPACK([], [AC_MSG_WARN([BLAS/LAPACK needed for harminv])])
	   LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
	elif test "x$acx_lapack_ok" = xno; then
	   AC_MSG_WARN([BLAS/LAPACK needed for harminv])
	fi

        save_CC=$CC
	CC=$CXX # harminv test must use $CXX since harminv may be a C++ lib

	# just use lapack/blas in default libs since we always need them
	AC_CHECK_LIB(harminv, harminv_get_freq_error, [have_harminv=yes; LIBS="-lharminv $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"], [AC_MSG_WARN([harminv support is disabled])], [$LAPACK_LIBS $BLAS_LIBS $FLIBS])

	CC=$save_CC
fi

if test $have_harminv = yes; then
	AC_DEFINE([HAVE_HARMINV],[1],[Define if you have libharminv])
fi

###########################################################################
# MPB library

AC_CHECK_HEADER(mpb.h, [have_mpb=maybe], [have_mpb=no])

if test $have_mpb = maybe; then
  if test "x$acx_lapack_ok" = x -a "x$have_harminv" = xno; then
	ACX_BLAS
	ACX_LAPACK([], [AC_MSG_WARN([BLAS/LAPACK needed for MPB])])
        LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
  elif test "x$acx_lapack_ok" = xno; then
       AC_MSG_WARN([BLAS/LAPACK needed for MPB])
  fi
  AC_CHECK_LIB(mpb, maxwell_target_operator1, [have_mpb=yes; LIBS="-lmpb $LIBS"],
               [AC_CHECK_LIB(mpb, maxwell_target_operator1,
                   [have_mpb=yes; LIBS="-lmpb $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"], [AC_MSG_WARN([MPB support is disabled])], [$LAPACK_LIBS $BLAS_LIBS $FLIBS])])
  if test $have_mpb = yes; then
     AC_DEFINE([HAVE_MPB], [1], [Define if you have libmpb])
  fi
fi

AM_CONDITIONAL(WITH_MPB, test "x$have_mpb" = "xyes")

##############################################################################
# GNU Scientific Library

AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_sf_bessel_Jn], [],
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

RUNCODE=""
if test "x$with_mpi" = "xyes"; then
   AC_CHECK_PROG(MPIRUN, mpirun, mpirun)
   # always use 2 processes for 'make check'
   RUNCODE="$MPIRUN -np 2"
fi
AC_SUBST(RUNCODE)
AM_CONDITIONAL(WITH_MPI, test "x$with_mpi" = "xyes")

##############################################################################
# Libraries and flags (other than -lmeep) required to link Meep:

MEEPLIBS="$LDFLAGS $LIBS"
AC_SUBST(MEEPLIBS)

##############################################################################

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

save_LIBS=$LIBS
save_CPPFLAGS=$CPPFLAGS

AC_CHECK_PROG(GUILE_CONFIG, guile-config, guile-config, unknown)
if test "x$GUILE_CONFIG" = "xunknown"; then
  AC_CHECK_LIB(readline, readline)
  AC_CHECK_LIB(dl, dlopen)
  AC_CHECK_LIB(guile-ltdl, scm_lt_dlinit, [], [AC_CHECK_LIB(ltdl, dlopen)])
  AC_CHECK_LIB(guile, gh_eval_str, [],
    [AC_CHECK_LIB(guile, scm_eval_string, [],
    [AC_MSG_ERROR([linking to guile failed])])])
  GUILE_LIBS="-lguile -lguile-ltdl -lreadline"
  GUILE_CPPFLAGS=""
else
  GUILE_LIBS=`$GUILE_CONFIG link`
  GUILE_CPPFLAGS=`$GUILE_CONFIG compile`
  CPPFLAGS="$CPPFLAGS $GUILE_CPPFLAGS"
  LIBS="$GUILE_LIBS $LIBS"
  AC_MSG_CHECKING([if linking to guile works])
  AC_TRY_LINK_FUNC(gh_enter, AC_MSG_RESULT(yes),
    [AC_TRY_LINK_FUNC(scm_boot_guile, AC_MSG_RESULT(yes),
    [AC_MSG_RESULT(no)
     AC_MSG_ERROR([guile-config is broken])])])
fi
AC_CHECK_HEADERS([libguile.h guile/gh.h])

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
if test "x$cross_compiling" = xno; then
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
fi
AC_SUBST(ACTIVATE_READLINE)

AC_SUBST(GUILE_LIBS)
AC_SUBST(GUILE_CPPFLAGS)
LIBS=$save_LIBS
CPPFLAGS=$save_CPPFLAGS

fi # if with_libctl

##############################################################################
# Check for libctl library and files

if test "x$with_libctl" != xno; then

save_LIBS=$LIBS
save_CPPFLAGS=$CPPFLAGS
LIBS="$GUILE_LIBS $LIBS"
CPPFLAGS="$CPPFLAGS $GUILE_CPPFLAGS"

AC_MSG_CHECKING([for libctl dir])
if test x != x"$LIBCTL_DIR" -a ! -r "$LIBCTL_DIR/base/ctl.scm"; then
	LIBCTL_DIR=""
fi
if test x = x"$LIBCTL_DIR" -a -r /usr/local/share/libctl/base/ctl.scm; then
	LIBCTL_DIR="/usr/local/share/libctl"
fi
if test x = x"$LIBCTL_DIR" -a -r /usr/share/libctl3/base/ctl.scm; then
        LIBCTL_DIR="/usr/share/libctl3"
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

# check for ctl.h (Fedora puts it in /usr/include/ctl, grrr)
CTL_H_CPPFLAG=""
AC_CHECK_HEADER(ctl.h, [CTL_H_CPPFLAG=" "], [
  save_CPPFLAGS=$CPPFLAGS
  AC_MSG_CHECKING([for alternate ctl.h location])
  for ctl_h_path in /usr/include/ctl /usr/local/include/ctl; do
     CPPFLAGS="$save_CPPFLAGS -I$ctl_h_path"
     AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <ctl.h>]],[])], [CTL_H_CPPFLAG=-I$ctl_h_path; break])
  done
  AC_MSG_RESULT([$CTL_H_CPPFLAG])
  CPPFLAGS=$save_CPPFLAGS])
if test "x$CTL_H_CPPFLAG" = "x"; then
  AC_MSG_ERROR([Couldn't find the ctl.h header file for libctl.])
fi
AC_SUBST(CTL_H_CPPFLAG)

save_CPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS $CTL_H_CPPFLAG"
# Check libctl version >= LIBCTL_MAJOR.LIBCTL_MINOR.LIBCTL_BUGFIX
LIBCTL_MAJOR=4; LIBCTL_MINOR=0; LIBCTL_BUGFIX=0
AC_MSG_CHECKING([whether libctl version is at least ${LIBCTL_MAJOR}.${LIBCTL_MINOR}.${LIBCTL_BUGFIX}])
AC_EGREP_CPP(yes, [[
#include <ctl.h>
#if LIBCTL_MAJOR_VERSION > $LIBCTL_MAJOR || (LIBCTL_MAJOR_VERSION == $LIBCTL_MAJOR && (LIBCTL_MINOR_VERSION > $LIBCTL_MINOR || (LIBCTL_MINOR_VERSION == $LIBCTL_MINOR && LIBCTL_BUGFIX_VERSION >= $LIBCTL_BUGFIX)))
  yes
#endif
]], [AC_MSG_RESULT(ok)], [AC_MSG_ERROR([libctl version ${LIBCTL_MAJOR}.${LIBCTL_MINOR}.${LIBCTL_BUGFIX} or later is required])])
CPPFLAGS=$save_CPPFLAGS

# Check for libctl_quiet feature (in libctl 3.1.1 or later):
AC_MSG_CHECKING([for libctl_quiet feature])
if grep libctl_quiet "$LIBCTL_DIR/base/main.c" 2>&1 > /dev/null; then
   AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_LIBCTL_QUIET], [1], [If we have the libctl_quiet variable])
else
   AC_MSG_RESULT(no)
fi

LIBS=$save_LIBS
CPPFLAGS=$save_CPPFLAGS

fi # if with_libctl

##############################################################################

LIBCTL_LIBS="-lctl $GUILE_LIBS"
LIBCTL_CPPFLAGS="$GUILE_CPPFLAGS"
AC_SUBST(LIBCTL_LIBS)
AC_SUBST(LIBCTL_CPPFLAGS)

##############################################################################

# check for -lctlgeom:
AC_CHECK_LIB(ctlgeom, vector3_dot, [have_libctlgeom=yes], [have_libctlgeom=no])
AM_CONDITIONAL(WITH_LIBCTLGEOM, test x"$have_libctlgeom" = "xyes")

LIBCTLGEOM_LIBS="-lctlgeom"
AC_SUBST(LIBCTLGEOM_LIBS)

##############################################################################
# check for libGDSII
AC_CHECK_HEADER(libGDSII.h, [have_gdsii=maybe], [have_gdsii=no])
if test $have_gdsii = maybe; then
  AC_CHECK_LIB(GDSII, libGDSIIExists)
fi

##############################################################################
# The following function is used only for debugging.  Note that
# we must test for it *after* setting the compiler flags (which
# affect whether it is declared, as it is a GNU extension).
# We need to #include <stdio.h> because that somehow affects whether
# the function is declared with older gcc versions.  We need
# to use AC_TRY_COMPILE because the test in AC_HAVE_DECL seems
# to be optimized out.

AC_CHECK_FUNCS(feenableexcept) dnl GNU libc fp exception control function
AC_MSG_CHECKING([whether feenableexcept declaration is usable])
feenableexcept_decl_ok=yes
AC_TRY_COMPILE([#include <stdio.h>
#define _GNU_SOURCE 1
#include <fenv.h>], [feenableexcept(0);],
[AC_DEFINE([HAVE_DECL_FEENABLEEXCEPT],[1],[Define if fenv.h declares this.])],
[feenableexcept_decl_ok=no])
AC_MSG_RESULT($feenableexcept_decl_ok)

##############################################################################
# See if we need to catch SIGFPE to avoid crashing on underflow exceptions
# etcetera, since some stupid operating systems (e.g. Alpha/Tru64) crash on
# SIGFPE by default.

AC_MSG_CHECKING([whether to catch and ignore SIGFPE signals])
AC_RUN_IFELSE([AC_LANG_PROGRAM([], [volatile double x=1,y=0; return x/y;])],
              [ignore_sigfpe=no], [ignore_sigfpe=yes], [ignore_sigfpe=no])
AC_MSG_RESULT($ignore_sigfpe)
if test $ignore_sigfpe = yes; then
    AC_DEFINE(IGNORE_SIGFPE, [1], [Define to catch and ignore SIGFPE signals])
fi

##############################################################################
# Miscellaneous function and header checks

AC_HEADER_TIME
AC_CHECK_HEADERS([sys/time.h])
AC_CHECK_FUNCS([BSDgettimeofday gettimeofday cblas_ddot cblas_daxpy jn])

##############################################################################
# check for restrict keyword in C++
AC_PREREQ(2.58) dnl for ac_c_restrict
AC_C_RESTRICT

# check for C99 _Pragma support, so that we can emit pragmas in macros
AC_TRY_COMPILE([], [_Pragma("ivdep")], [], [AC_DEFINE([_Pragma],[],[define to nothing if C99 _Pragma is not supported])])

##############################################################################

# checks for python
##############################################################################
AC_ARG_WITH(python,
	[AC_HELP_STRING([--without-python], [compile without Python interface])],
        with_python=$withval,with_python=yes)

if test "x$with_python" = xno; then
  have_python=no
else
  if test "x$have_libctlgeom" = xno; then
    AC_MSG_ERROR([libctlgeom was not found, is required for Python interface])
  fi

  if test "$enable_shared" = no; then
     AC_MSG_WARN([Python interface requires --enable-shared; disabling])
     have_python=no
  else

    AM_PATH_PYTHON([],[have_python=yes],[have_python=no])
    if test $have_python = yes; then

      AC_MSG_CHECKING([for Python include directory])
      pinc=`echo "import distutils.sysconfig; print (distutils.sysconfig.get_python_inc())" | $PYTHON - 2>/dev/null`
      AC_MSG_RESULT([${pinc:-unknown}])
      test -n "$pinc" && PYTHON_INCLUDES="-I$pinc"
      save_CPPFLAGS=$CPPFLAGS
      CPPFLAGS="$CPPFLAGS $PYTHON_INCLUDES"
      AC_CHECK_HEADER([Python.h], [], [AC_MSG_WARN([disabling Python wrappers])
                                       have_python=no])

      if test $have_python = yes; then
        AC_MSG_CHECKING([for Numpy include directory])
        pinc=`echo "import numpy; print (numpy.get_include())" | $PYTHON - 2>/dev/null`
        AC_MSG_RESULT([${pinc:-unknown}])
        test -n "$pinc" && PYTHON_INCLUDES="$PYTHON_INCLUDES -I$pinc"
        CPPFLAGS="$save_CPPFLAGS $PYTHON_INCLUDES"
        AC_CHECK_HEADER([numpy/arrayobject.h],[],[
          AC_MSG_WARN([disabling Python wrappers])
          have_python=no],[#include <Python.h>])
      fi

      CPPFLAGS=$save_CPPFLAGS
    fi # have_python

  fi # enabled_shared
fi # with_python

AC_SUBST(PYTHON_INCLUDES)
AM_CONDITIONAL(WITH_PYTHON, test x"$have_python" = "xyes")

# Copy/symlink casimir.scm to builddir for out-of-tree builds
AC_CONFIG_LINKS(scheme/casimir.scm:scheme/casimir.scm)

##############################################################################

AC_CONFIG_FILES([
	Makefile
	meep-pkgconfig
	src/Makefile
	tests/Makefile
	scheme/Makefile
	scheme/meep.scm
	libmeepgeom/Makefile
	libpympb/Makefile
	python/Makefile
])

AC_OUTPUT
