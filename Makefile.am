ACLOCAL_AMFLAGS=-I m4

SUBDIRS = src tests examples

if WITH_LIBCTL
  SUBDIRS += libctl
endif
if WITH_LIBCTLGEOM
  SUBDIRS += libmeepgeom
endif
if WITH_PYTHON
  SUBDIRS += python
endif

EXTRA_DIST = COPYRIGHT m4 meep-pkgconfig.in

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = meep@MEEP_SUFFIX@.pc

meep@MEEP_SUFFIX@.pc: meep-pkgconfig
	cp -f $(top_builddir)/meep-pkgconfig $@

libmeep:
	cd src && $(MAKE)

DISTCLEANFILES = meep@MEEP_SUFFIX@.pc
