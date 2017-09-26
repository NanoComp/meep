LIBCTLGEOM   = -lctlgeom
LIBMEEP      = $(top_builddir)/src/libmeep@MEEP_SUFFIX@.la
MEEPLIBS     = $(LIBCTLGEOM) $(LIBMEEP)

lib_LTLIBRARIES = libmeepgeom.la
pkginclude_HEADERS = meepgeom.hpp material_data.hpp

#AM_CPPFLAGS = -I$(top_srcdir)/src -std=c++11
AM_CPPFLAGS = -I$(top_srcdir)/src -DDATADIR=\"$(srcdir)/\"

libmeepgeom_la_SOURCES = \
 meepgeom.cpp meepgeom.hpp material_data.hpp

libmeepgeom_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@ 

check_PROGRAMS      = pw-source-ll ring-ll cyl-ellipsoid-ll absorber-1d-ll

absorber_1d_ll_SOURCES = absorber-1d-ll.cpp
absorber_1d_ll_LDADD   = libmeepgeom.la $(MEEPLIBS)

cyl_ellipsoid_ll_SOURCES = cyl-ellipsoid-ll.cpp
cyl_ellipsoid_ll_LDADD   = libmeepgeom.la $(MEEPLIBS)

pw_source_ll_SOURCES = pw-source-ll.cpp
pw_source_ll_LDADD   = libmeepgeom.la $(MEEPLIBS)

ring_ll_SOURCES = ring-ll.cpp
ring_ll_LDADD   = libmeepgeom.la $(MEEPLIBS)

TESTS           = cyl-ellipsoid-ll

noinst_PROGRAMS = bend-flux-ll cavity_arrayslice_ll
bend_flux_ll_SOURCES = bend-flux-ll.cpp
bend_flux_ll_LDADD   = libmeepgeom.la $(MEEPLIBS)

cavity_arrayslice_ll_SOURCES = cavity-arrayslice-ll.cpp
cavity_arrayslice_ll_LDADD = libmeepgeom.la $(MEEPLIBS)

noinst_DATA = cyl-ellipsoid-eps-ref.h5
