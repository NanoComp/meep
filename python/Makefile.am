HPPFILES=                                    \
 $(top_srcdir)/src/meep.hpp                  \
 $(top_srcdir)/src/meep_internals.hpp        \
 $(top_srcdir)/src/bicgstab.hpp              \
 $(top_srcdir)/src/meep/vec.hpp              \
 $(top_srcdir)/src/meep/mympi.hpp            \
 $(top_srcdir)/src/meepgeom.hpp              \
 $(top_srcdir)/src/material_data.hpp

BUILT_SOURCES = meep-python.cxx __init__.py
EXTRA_DIST = $(BUILT_SOURCES) typemap_utils.cpp materials.py examples tests

CLEANFILES = $(BUILT_SOURCES) meep.py .coverage

if WITH_MPB
  PYMPBINCLUDE=-I$(top_srcdir)/libpympb
  LIBPYMPB=$(top_builddir)/libpympb/libpympb.la
  BUILT_SOURCES += mpb-python.cxx mpb.py
  EXTRA_DIST += solver.py mpb_data.py
endif

LIBMEEP = $(top_builddir)/src/libmeep.la

AM_CPPFLAGS = -I$(top_srcdir)/src                   \
              $(PYMPBINCLUDE)                       \
              -I$(top_builddir) # for config.h

_meep_la_SOURCES = meep-python.cxx
_meep_la_LIBADD = $(LIBMEEP) $(PYTHON_LIBS) @LIBCTLGEOM_LIBS@
_meep_la_LDFLAGS = -module -version-info @SHARED_VERSION_INFO@
_meep_la_CPPFLAGS = $(PYTHON_INCLUDES) $(AM_CPPFLAGS)

if WITH_MPB
_mpb_la_SOURCES = mpb-python.cxx
_mpb_la_LIBADD = $(PYTHON_LIBS) $(LIBPYMPB) @LIBCTLGEOM_LIBS@
_mpb_la_LDFLAGS = -module -version-info @SHARED_VERSION_INFO@
_mpb_la_CPPFLAGS = $(PYTHON_INCLUDES) $(AM_CPPFLAGS)
MPB_LA = _mpb.la
endif # WITH_MPB

TEST_DIR = tests

# material_dispersion.py test must be excluded from test suite for MPI build
if WITH_MPI
  MDPYTEST=
  MPBPYTEST=
else
  MDPYTEST=$(TEST_DIR)/material_dispersion.py
if WITH_MPB
  MPBPYTEST=$(TEST_DIR)/mpb.py
else
  MPBPYTEST=
endif # WITH_MPB
endif # WITH_MPI

if WITH_MPB
  BINARY_GRATING_TEST = $(TEST_DIR)/binary_grating.py
  KDOM_TEST = $(TEST_DIR)/kdom.py
  MODE_COEFFS_TEST = $(TEST_DIR)/mode_coeffs.py
  MODE_DECOMPOSITION_TEST = $(TEST_DIR)/mode_decomposition.py
  WVG_SRC_TEST = $(TEST_DIR)/wvg_src.py
else
  BINARY_GRATING_TEST =
  KDOM_TEST =
  MODE_COEFFS_TEST =
  MODE_DECOMPOSITION_TEST =
  WVG_SRC_TEST =
endif

TESTS =                                   \
    $(TEST_DIR)/3rd_harm_1d.py            \
    $(TEST_DIR)/absorber_1d.py            \
    $(TEST_DIR)/antenna_radiation.py      \
    $(TEST_DIR)/bend_flux.py              \
    $(BINARY_GRATING_TEST)                \
    $(TEST_DIR)/cavity_arrayslice.py      \
    $(TEST_DIR)/chunks.py                 \
    $(TEST_DIR)/cyl_ellipsoid.py          \
    $(TEST_DIR)/dft_fields.py             \
    $(TEST_DIR)/field_functions.py        \
    $(TEST_DIR)/force.py                  \
    $(TEST_DIR)/fragment_stats.py         \
    $(TEST_DIR)/geom.py                   \
    $(TEST_DIR)/holey_wvg_bands.py        \
    $(TEST_DIR)/holey_wvg_cavity.py       \
    $(KDOM_TEST)                          \
    $(TEST_DIR)/ldos.py                   \
    $(MDPYTEST)                           \
    $(MPBPYTEST)                          \
    $(MODE_COEFFS_TEST)                   \
    $(MODE_DECOMPOSITION_TEST)            \
    $(TEST_DIR)/multilevel_atom.py        \
    $(TEST_DIR)/oblique_source.py        \
    $(TEST_DIR)/physical.py               \
    $(TEST_DIR)/pw_source.py              \
    $(TEST_DIR)/refl_angular.py           \
    $(TEST_DIR)/ring.py                   \
    $(TEST_DIR)/ring_cyl.py               \
    $(TEST_DIR)/simulation.py             \
    $(TEST_DIR)/source.py                 \
    $(TEST_DIR)/user_defined_material.py  \
    $(WVG_SRC_TEST)

if WITH_COVERAGE
  PY_LOG_COMPILER = coverage run -a --omit=$(top_srcdir)/python/tests/*,${HOME}/virtualenv/*,$(top_srcdir)/python/examples/*
else
  PY_LOG_COMPILER = $(RUNCODE) $(PYTHON)
endif

TEST_EXTENSIONS = .py
TESTS_ENVIRONMENT = export PYTHONPATH=$(abs_top_builddir)/python:$$PYTHONPATH;

if WITH_PYTHON
  pkgpython_PYTHON = geom.py __init__.py source.py simulation.py materials.py
  pkgpyexec_LTLIBRARIES = _meep.la
if WITH_MPB
  pympbdir = $(pkgpythondir)/mpb
  pympb_DATA = mpb.py solver.py mpb_data.py
  pympb_LTLIBRARIES = _mpb.la
endif # WITH_MPB
endif # WITH_PYTHON

if WITH_MPB
install-data-hook:
	mv $(DESTDIR)$(pkgpythondir)/mpb/mpb.py $(DESTDIR)$(pkgpythondir)/mpb/__init__.py

uninstall-hook:
	rm -f $(DESTDIR)$(pkgpythondir)/mpb/__init__.py
endif

if MAINTAINER_MODE

SWIG ?= swig
SWIG_VERSION = $(shell $(SWIG) -version | grep Version | awk '{print $$3}')
MEEP_SWIG_SRC = meep.i numpy.i vec.i

meep-python.cxx: $(MEEP_SWIG_SRC) $(HPPFILES)
	$(SWIG) -Wextra $(AM_CPPFLAGS) -outdir $(builddir) -c++ -python -o $@ $(srcdir)/meep.i

if WITH_MPB
MPB_SWIG_SRC = mpb.i
mpb-python.cxx: $(MPB_SWIG_SRC) $(top_srcdir)/libpympb/pympb.hpp
	$(SWIG) -Wextra $(AM_CPPFLAGS) $(PYMPBINCLUDE) -outdir $(builddir) -c++ -python -o $@ $(srcdir)/mpb.i

mpb.py: mpb-python.cxx
MPB_PY = mpb.py
endif # WITH_MPB

meep.py: meep-python.cxx

__init__.py: meep.py $(MPB_PY)
	cp $< $@
	echo "__version__ = '$(shell git describe --tags | sed 's/^v//')'" >> $@
	if [[ "${SWIG_VERSION}" = 3.0.12 ]]; then \
		sed -i.bak '/^if _swig_python_version_info >= (2, 7, 0):/,/^else:/d' $@; \
		sed -i.bak 's/    import _meep/from . import _meep/' $@; \
	fi
if WITH_MPB
	if [[ "${SWIG_VERSION}" == 3.0.12 ]]; then \
		sed -i.bak '/^if _swig_python_version_info >= (2, 7, 0):/,/^else:/d' mpb.py; \
		sed -i.bak 's/    import _mpb/from . import _mpb/' mpb.py; \
	fi
endif

INIT_PY = __init__.py
MPB_PY_PATH = $(builddir)/mpb.py
else

INIT_PY = $(srcdir)/__init__.py
MPB_PY_PATH = $(srcdir)/mpb.py

endif # MAINTAINER_MODE

HL_IFACE =                  \
    $(srcdir)/geom.py       \
    $(srcdir)/simulation.py \
    $(srcdir)/source.py     \
    $(srcdir)/materials.py

PY_PKG_FILES =      \
    $(INIT_PY)      \
    $(HL_IFACE)     \
    .libs/_meep.so

meep: _meep.la $(MPB_LA) __init__.py $(HL_IFACE)
	mkdir -p meep
	cp $(PY_PKG_FILES) meep
if WITH_MPB
	mkdir -p meep/mpb
	cp .libs/_mpb.so meep/mpb
	cp $(MPB_PY_PATH) meep/mpb/__init__.py
	cp $(srcdir)/solver.py $(srcdir)/mpb_data.py meep/mpb
endif

all-local: meep

clean-local:
	rm -rf meep __init__.py.bak

distclean-local:
	rm -f *.h5
