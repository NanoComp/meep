LIBCTLGEOM   = -lctlgeom
LIBMEEP      = $(top_builddir)/src/libmeep.la
MEEPLIBS     = $(LIBCTLGEOM) $(LIBMEEP)

lib_LTLIBRARIES = libmeepgeom.la
pkginclude_HEADERS = meepgeom.hpp material_data.hpp

#AM_CPPFLAGS = -I$(top_srcdir)/src -std=c++11
AM_CPPFLAGS = -I$(top_srcdir)/src -DDATADIR=\"$(srcdir)/\"

libmeepgeom_la_SOURCES = \
 meepgeom.cpp meepgeom.hpp material_data.hpp

libmeepgeom_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@

check_PROGRAMS      = pw-source-ll ring-ll cyl-ellipsoid-ll absorber-1d-ll array-slice-ll user-defined-material dft-fields

absorber_1d_ll_SOURCES = absorber-1d-ll.cpp
absorber_1d_ll_LDADD   = libmeepgeom.la $(MEEPLIBS)

cyl_ellipsoid_ll_SOURCES = cyl-ellipsoid-ll.cpp
cyl_ellipsoid_ll_LDADD   = libmeepgeom.la $(MEEPLIBS)

pw_source_ll_SOURCES = pw-source-ll.cpp
pw_source_ll_LDADD   = libmeepgeom.la $(MEEPLIBS)

ring_ll_SOURCES = ring-ll.cpp
ring_ll_LDADD   = libmeepgeom.la $(MEEPLIBS)

array_slice_ll_SOURCES = array-slice-ll.cpp
array_slice_ll_LDADD   = libmeepgeom.la $(MEEPLIBS)

dft_fields_SOURCES = dft-fields.cpp
dft_fields_LDADD   = libmeepgeom.la $(MEEPLIBS)

TESTS           = cyl-ellipsoid-ll array-slice-ll
TESTS_ENVIRONMENT = export LD_PRELOAD=$(abs_top_builddir)/src/.libs/libmeep.so;

#LOG_COMPILER = $(RUNCODE)

noinst_PROGRAMS = bend-flux-ll

bend_flux_ll_SOURCES = bend-flux-ll.cpp
bend_flux_ll_LDADD   = libmeepgeom.la $(MEEPLIBS)

user_defined_material_SOURCES = user-defined-material.cpp
user_defined_material_LDADD   = libmeepgeom.la $(MEEPLIBS)

dist_noinst_DATA = cyl-ellipsoid-eps-ref.h5 array-slice-ll.h5

distclean-local:
	rm -f eps-000000000.h5
