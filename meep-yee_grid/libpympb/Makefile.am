LIBMEEP = $(top_builddir)/src/libmeep.la

AM_CPPFLAGS =                   \
    -Wfatal-errors              \
    -I$(top_srcdir)/src         \
    -I$(top_builddir)

lib_LTLIBRARIES = libpympb.la
pkginclude_HEADERS = pympb.hpp

libpympb_la_SOURCES = pympb.cpp pympb.hpp
libpympb_la_LIBADD = $(LIBMEEP)
libpympb_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@
libpympb_la_CPPFLAGS = $(AM_CPPFLAGS)
