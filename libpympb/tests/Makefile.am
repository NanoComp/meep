LIBCTLGEOM = @LIBCTLGEOM_LIBS@
LIBMEEPGEOM = $(top_builddir)/libmeepgeom/libmeepgeom.la
LIBPYMPB = $(top_builddir)/libpympb/libpympb.la

AM_CPPFLAGS = -Wfatal-errors -I$(top_srcdir)/libmeepgeom -I$(top_srcdir)/libpympb
bin_PROGRAMS = tutorial

tutorial_SOURCES = tutorial.cpp
tutorial_LDADD = $(LIBPYMPB) $(LIBMEEPGEOM) $(LIBCTLGEOM)

TESTS = tutorial
