LIBCTLGEOM   = -lctlgeom
LIBMEEP      = $(top_builddir)/src/libmeep.la
MEEPLIBS     = $(LIBCTLGEOM) $(LIBMEEP)

lib_LTLIBRARIES = libmeep_geom.la
pkginclude_HEADERS = meep_geom.hpp material_data.hpp

AM_CPPFLAGS = -I$(top_srcdir)/src -std=c++11

libmeep_geom_la_SOURCES = \
 meep_geom.cpp meep_geom.hpp material_data.hpp

libmeep_geom_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@

noinst_PROGRAMS      = pw-source-ll bend-flux-ll

pw_source_ll_SOURCES = pw-source-ll.cpp
pw_source_ll_LDADD   = libmeep_geom.la $(MEEPLIBS)

bend_flux_ll_SOURCES = bend-flux-ll.cpp
bend_flux_ll_LDADD   = libmeep_geom.la $(MEEPLIBS)
