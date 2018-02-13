ACLOCAL_AMFLAGS=-I m4

SUBDIRS = src tests examples

if WITH_LIBCTL
  SUBDIRS += scheme
endif
if WITH_LIBCTLGEOM
  SUBDIRS += libmeepgeom
endif
if WITH_PYTHON
  SUBDIRS += python
endif

EXTRA_DIST = NEWS.md LICENSE COPYRIGHT m4 meep-pkgconfig.in

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = meep.pc

meep.pc: meep-pkgconfig
	cp -f $(top_builddir)/meep-pkgconfig $@

libmeep:
	cd src && $(MAKE)

DISTCLEANFILES = meep.pc
