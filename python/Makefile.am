SWIG_SRC = meep.i numpy.i vec.i
EXTRA_DIST = $(SWIG_SRC)
HPPFILES=                                    \
 $(top_srcdir)/src/meep_internals.hpp        \
 $(top_srcdir)/src/bicgstab.hpp              \
 $(top_srcdir)/src/meep/vec.hpp              \
 $(top_srcdir)/src/meep/mympi.hpp            \
 $(top_srcdir)/libmeepgeom/meepgeom.hpp      \
 $(top_srcdir)/libmeepgeom/material_data.hpp 

BUILT_SOURCES = meep-python.cpp meep.py

CLEANFILES = $(BUILT_SOURCES) meep-python.h

AM_CPPFLAGS = -I$(top_srcdir)/src                   \
              -I$(top_srcdir)/libmeepgeom           \
              -I$(top_builddir) # for config.h

_meep_la_SOURCES = meep-python.cpp
_meep_la_LIBADD = $(top_builddir)/src/libmeep@MEEP_SUFFIX@.la -L$(top_srcdir)/libmeepgeom/.libs $(PYTHON_LIBS) -lctlgeom -lmeepgeom
_meep_la_LDFLAGS = -module -version-info @SHARED_VERSION_INFO@
_meep_la_CPPFLAGS = $(PYTHON_INCLUDES) $(AM_CPPFLAGS)

TEST_DIR = tests
TESTS =                               \
	$(TEST_DIR)/test_cyl_ellipsoid.py \
    $(TEST_DIR)/test_geom.py          \
    $(TEST_DIR)/test_physical.py      \
    $(TEST_DIR)/test_ring.py

TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(PYTHON)
AM_TESTS_ENVIRONMENT = PYTHONPATH='$(top_srcdir)/python/.libs'; export PYTHONPATH;

if WITH_PYTHON
  python_PYTHON = meep.py
  pyexec_LTLIBRARIES = _meep.la
endif

if MAINTAINER_MODE

meep-python.cpp:	$(SWIG_SRC) $(HPPFILES)
	swig -Wextra $(AM_CPPFLAGS) -outdir $(builddir) -c++ -python -o $@ $(srcdir)/meep.i

meep.py:       		meep-python.cpp

endif # MAINTAINER_MODE

