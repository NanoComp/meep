SWIG_SRC = meep.i numpy.i vec.i
EXTRA_DIST = $(SWIG_SRC)
HPPFILES=                                    \
 $(top_srcdir)/src/meep_internals.hpp        \
 $(top_srcdir)/src/bicgstab.hpp              \
 $(top_srcdir)/src/meep/vec.hpp              \
 $(top_srcdir)/src/meep/mympi.hpp            \
 $(top_srcdir)/libmeepgeom/meepgeom.hpp      \
 $(top_srcdir)/libmeepgeom/material_data.hpp 

BUILT_SOURCES = meep-python.cpp meep.py __init__.py

CLEANFILES = $(BUILT_SOURCES) meep-python.h

AM_CPPFLAGS = -I$(top_srcdir)/src                   \
              -I$(top_srcdir)/libmeepgeom           \
              -I$(top_builddir) # for config.h

_meep_la_SOURCES = meep-python.cpp
_meep_la_LIBADD = $(top_builddir)/src/libmeep@MEEP_SUFFIX@.la $(top_builddir)/libmeepgeom/libmeepgeom.la $(PYTHON_LIBS) -lctlgeom
_meep_la_LDFLAGS = -module -version-info @SHARED_VERSION_INFO@
_meep_la_CPPFLAGS = $(PYTHON_INCLUDES) $(AM_CPPFLAGS)

TEST_DIR = tests
TESTS =                                    \
    $(TEST_DIR)/cyl_ellipsoid.py      \
    $(TEST_DIR)/geom.py               \
    $(TEST_DIR)/holy_wvg_bands.py     \
    $(TEST_DIR)/physical.py           \
    $(TEST_DIR)/ring.py               \
    $(TEST_DIR)/simulation.py         \
    $(TEST_DIR)/source.py

TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(PYTHON)
TESTS_ENVIRONMENT = export PYTHONPATH=$(abs_top_builddir)/python:$$PYTHONPATH;

if WITH_PYTHON
  pkgpython_PYTHON = geom.py __init__.py source.py
  pkgpyexec_LTLIBRARIES = _meep.la
endif

if MAINTAINER_MODE

PY_PKG_FILES =               \
    __init__.py              \
    $(srcdir)/geom.py        \
    $(srcdir)/simulation.py  \
    $(srcdir)/source.py      \
    .libs/_meep.so

meep-python.cpp: $(SWIG_SRC) $(HPPFILES)
	swig -Wextra $(AM_CPPFLAGS) -outdir $(builddir) -c++ -python -o $@ $(srcdir)/meep.i

meep.py: meep-python.cpp

__init__.py: meep.py
	cp $< $@

pymeep: .libs/_meep.so __init__.py
	test -d meep || mkdir meep
	cp $(PY_PKG_FILES) meep

all-local: pymeep

clean-local:
	rm -rf meep

endif # MAINTAINER_MODE

