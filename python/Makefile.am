HPPFILES=                                    \
 $(top_srcdir)/src/meep_internals.hpp        \
 $(top_srcdir)/src/bicgstab.hpp              \
 $(top_srcdir)/src/meep/vec.hpp              \
 $(top_srcdir)/src/meep/mympi.hpp            \
 $(top_srcdir)/libmeepgeom/meepgeom.hpp      \
 $(top_srcdir)/libmeepgeom/material_data.hpp

BUILT_SOURCES = meep-python.cpp mpb-python.cpp __init__.py mpb.py
EXTRA_DIST = $(BUILT_SOURCES) typemap_utils.cpp solver.py mpb_data.py examples tests

CLEANFILES = $(BUILT_SOURCES) meep.py

if WITH_MPB
  PYMPBINCLUDE=-I$(top_srcdir)/libpympb
  LIBPYMPB=$(top_builddir)/libpympb/libpympb.la
endif

LIBMEEPGEOM = $(top_builddir)/libmeepgeom/libmeepgeom.la
LIBMEEP = $(top_builddir)/src/libmeep.la

AM_CPPFLAGS = -I$(top_srcdir)/src                   \
              -I$(top_srcdir)/libmeepgeom           \
              $(PYMPBINCLUDE)                       \
              -I$(top_builddir) # for config.h

_meep_la_SOURCES = meep-python.cpp
_meep_la_LIBADD = $(LIBMEEP) $(LIBMEEPGEOM) $(PYTHON_LIBS) @LIBCTLGEOM_LIBS@
_meep_la_LDFLAGS = -module -version-info @SHARED_VERSION_INFO@
_meep_la_CPPFLAGS = $(PYTHON_INCLUDES) $(AM_CPPFLAGS)

_mpb_la_SOURCES = mpb-python.cpp
_mpb_la_LIBADD = $(LIBMEEPGEOM) $(PYTHON_LIBS) $(LIBPYMPB) @LIBCTLGEOM_LIBS@
_mpb_la_LDFLAGS = -module -version-info @SHARED_VERSION_INFO@
_mpb_la_CPPFLAGS = $(PYTHON_INCLUDES) $(AM_CPPFLAGS)

# material_dispersion.py test must be excluded from test suite for MPI build
if WITH_MPI
  MDPYTEST=
else
  MDPYTEST=$(TEST_DIR)/material_dispersion.py
endif

TEST_DIR = tests
TESTS =                                   \
    $(TEST_DIR)/3rd_harm_1d.py            \
    $(TEST_DIR)/absorber_1d.py            \
    $(TEST_DIR)/antenna_radiation.py      \
    $(TEST_DIR)/bend_flux.py              \
    $(TEST_DIR)/cavity_arrayslice.py      \
    $(TEST_DIR)/cyl_ellipsoid.py          \
    $(TEST_DIR)/field_functions.py        \
    $(TEST_DIR)/force.py                  \
    $(TEST_DIR)/geom.py                   \
    $(TEST_DIR)/holey_wvg_bands.py        \
    $(TEST_DIR)/holey_wvg_cavity.py       \
    $(TEST_DIR)/ldos.py                   \
    $(MDPYTEST)                           \
    $(TEST_DIR)/mpb.py                    \
    $(TEST_DIR)/physical.py               \
    $(TEST_DIR)/pw_source.py              \
    $(TEST_DIR)/ring.py                   \
    $(TEST_DIR)/ring_cyl.py               \
    $(TEST_DIR)/simulation.py             \
    $(TEST_DIR)/source.py                 \
    $(TEST_DIR)/user_defined_material.py  \
    $(TEST_DIR)/wvg_src.py

TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(RUNCODE) $(PYTHON)
TESTS_ENVIRONMENT = export PYTHONPATH=$(abs_top_builddir)/python:$$PYTHONPATH;

if WITH_PYTHON
  pympbdir = $(pkgpythondir)/mpb
  pkgpython_PYTHON = geom.py __init__.py source.py simulation.py
  pkgpyexec_LTLIBRARIES = _meep.la
  pympb_DATA = mpb.py solver.py mpb_data.py
  pympb_LTLIBRARIES = _mpb.la
endif

install-data-hook:
	mv $(DESTDIR)$(pkgpythondir)/mpb/mpb.py $(DESTDIR)$(pkgpythondir)/mpb/__init__.py

uninstall-hook:
	rm -f $(DESTDIR)$(pkgpythondir)/mpb/__init__.py

if MAINTAINER_MODE

SWIG ?= swig
SWIG_VERSION = $(shell $(SWIG) -version | grep Version | awk '{print $$3}')
MEEP_SWIG_SRC = meep.i numpy.i vec.i

meep-python.cpp: $(MEEP_SWIG_SRC) $(HPPFILES)
	$(SWIG) -Wextra $(AM_CPPFLAGS) -outdir $(builddir) -c++ -python -o $@ $(srcdir)/meep.i

MPB_SWIG_SRC = mpb.i

mpb-python.cpp: $(MPB_SWIG_SRC) $(top_srcdir)/libpympb/pympb.hpp
	$(SWIG) -Wextra $(AM_CPPFLAGS) $(PYMPBINCLUDE) -outdir $(builddir) -c++ -python -o $@ $(srcdir)/mpb.i

meep.py: meep-python.cpp

mpb.py: mpb-python.cpp

__init__.py: meep.py
	cp $< $@
	if [[ "${SWIG_VERSION}" = 3.0.12 ]]; then \
		sed -i.bak '/^if _swig_python_version_info >= (2, 7, 0):/,/^else:/d' $@; \
		sed -i.bak 's/    import _meep/from . import _meep/' $@; \
	fi


INIT_PY = __init__.py
MPB_PY = $(builddir)/mpb.py
else

INIT_PY = $(srcdir)/__init__.py
MPB_PY = $(srcdir)/mpb.py

endif # MAINTAINER_MODE

PY_PKG_FILES =               \
    $(INIT_PY)               \
    $(srcdir)/geom.py        \
    $(srcdir)/simulation.py  \
    $(srcdir)/source.py      \
    .libs/_meep.so

meep: _meep.la _mpb.la __init__.py
	mkdir -p meep/mpb
	cp $(PY_PKG_FILES) meep
	cp .libs/_mpb.so meep/mpb
	cp $(MPB_PY) meep/mpb/__init__.py
	cp $(srcdir)/solver.py $(srcdir)/mpb_data.py meep/mpb

all-local: meep

clean-local:
	rm -rf meep __init__.py.bak

distclean-local:
	rm -f *.h5
