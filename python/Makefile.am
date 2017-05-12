SWIG_SRC = meep.i numpy.i vec.i
EXTRA_DIST = $(SWIG_SRC)
HPPFILES=	                          \
 $(top_srcdir)/src/meep.hpp               \
 $(top_srcdir)/src/meep_internals.hpp     \
 $(top_srcdir)/src/bicgstab.hpp           \
 $(top_srcdir)/src/meep/vec.hpp           \
 $(top_srcdir)/src/meep/mympi.hpp

BUILT_SOURCES = meep-python.cpp meep.py

AM_CPPFLAGS = -I$(top_srcdir)/src                   \
              -I$(top_srcdir)/src/meep              \
              -I$(top_builddir) # for config.h

_meep_la_SOURCES = meep-python.cpp
_meep_la_LIBADD = $(top_builddir)/src/libmeep.la $(PYTHON_LIBS)
_meep_la_LDFLAGS = -module -version-info @SHARED_VERSION_INFO@
_meep_la_CPPFLAGS = $(PYTHON_INCLUDES) $(AM_CPPFLAGS)

if WITH_PYTHON
  python_PYTHON = meep.py
  pyexec_LTLIBRARIES = _meep.la
endif

if MAINTAINER_MODE

meep-python.cpp:	$(SWIG_SRC) $(HPPFILES)
	swig -Wall $(AM_CPPFLAGS) -c++ -python -o $@ $(srcdir)/meep.i

meep.py:       		meep-python.cpp

endif # MAINTAINER_MODE

