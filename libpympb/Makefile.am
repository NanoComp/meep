# LIBCTLGEOM   = -lctlgeom
# LIBMEEP      = $(top_builddir)/src/libmeep.la
MEEPLIBS     = $(LIBCTLGEOM) $(LIBMEEP)

AM_CPPFLAGS = -Wfatal-errors -I$(top_srcdir)/libmeepgeom
lib_LTLIBRARIES = libpympb.la
# pkginclude_HEADERS = meepgeom.hpp material_data.hpp

# AM_CPPFLAGS = -I$(top_srcdir)/src -DDATADIR=\"$(srcdir)/\"

libpympb_la_SOURCES = pympb.cpp
libpympb_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@

# check_PROGRAMS =
# TESTS =
