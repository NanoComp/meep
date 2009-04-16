if WITH_LIBCTL
LIBCTL = libctl
else
LIBCTL =
endif

ACLOCAL_AMFLAGS=-I m4

if WITH_DOC
DOC = doc
else
DOC = 
endif

SUBDIRS = src $(LIBCTL) tests examples $(DOC)

EXTRA_DIST = COPYRIGHT m4 meep-pkgconfig.in

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = meep@MEEP_SUFFIX@.pc

meep@MEEP_SUFFIX@.pc: meep-pkgconfig
	cp -f $(top_builddir)/meep-pkgconfig $@

libmeep:
	cd src && $(MAKE)

# Somewhat hackish.  The "right" way to do this is by a dist-hook target,
# but then darcs check will fail because it doesn't run in the darcs
# repository.
darcs-dist: distdir
	darcs changes --summary > $(distdir)/ChangeLog
	tardir=$(distdir) && $(am__tar) | GZIP=$(GZIP_ENV) gzip -c >$(distdir).tar.gz
	$(am__remove_distdir)

