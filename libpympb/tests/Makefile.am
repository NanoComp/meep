LIBCTLGEOM = @LIBCTLGEOM_LIBS@
LIBMEEPGEOM = $(top_builddir)/libmeepgeom/libmeepgeom.la
LIBPYMPB = $(top_builddir)/libpympb/libpympb.la
LIBMEEP = $(top_builddir)/src/libmeep.la

AM_CPPFLAGS =                   \
    -Wfatal-errors              \
    -I$(top_srcdir)/src         \
    -I$(top_srcdir)/libmeepgeom \
    -I$(top_srcdir)/libpympb

bin_PROGRAMS = tutorial

tutorial_SOURCES = tutorial.cpp
tutorial_LDADD = $(LIBPYMPB) $(LIBMEEPGEOM) $(LIBCTLGEOM) $(LIBMEEP)
tutorial_CPPFLAGS = $(AM_CPPFLAGS)

TESTS = tutorial

clean-local:
	rm -f *.h5
