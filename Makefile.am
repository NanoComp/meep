if WITH_LIBCTL
LIBCTL = libctl
else
LIBCTL =
endif

ACLOCAL_AMFLAGS=-I m4

SUBDIRS = src $(LIBCTL) tests examples libmeepgeom

if WITH_PYTHON
  SUBDIRS += python python_examples
endif

EXTRA_DIST = COPYRIGHT m4 meep-pkgconfig.in

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = meep@MEEP_SUFFIX@.pc

meep@MEEP_SUFFIX@.pc: meep-pkgconfig
	cp -f $(top_builddir)/meep-pkgconfig $@

libmeep:
	cd src && $(MAKE)

DISTCLEANFILES = meep@MEEP_SUFFIX@.pc
