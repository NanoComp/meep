# recompile Haskell source only in maintainer mode
if MAINTAINER_MODE
HSRC = hsrc
else
HSRC =
endif

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

SUBDIRS = $(HSRC) src $(LIBCTL) tests examples $(DOC)

EXTRA_DIST = COPYRIGHT m4 meep-pkgconfig.in

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = meep@MEEP_SUFFIX@.pc

meep@MEEP_SUFFIX@.pc: meep-pkgconfig
	cp -f $(top_builddir)/meep-pkgconfig $@

if MAINTAINER_MODE
generator:
	cd hsrc && $(MAKE) depend && $(MAKE)
else
generator:
endif

libmeep: generator
	cd src && $(MAKE)

# add a fake dependency to make sure stuff is built before 'make dist'
if WITH_DOC
COPYRIGHT: generator libmeep
else
COPYRIGHT: generator
endif

# Somewhat hackish.  The "right" way to do this is by a dist-hook target,
# but then darcs check will fail because it doesn't run in the darcs
# repository.
darcs-dist: distdir
	darcs changes --summary > $(distdir)/ChangeLog
	tardir=$(distdir) && $(am__tar) | GZIP=$(GZIP_ENV) gzip -c >$(distdir).tar.gz
	$(am__remove_distdir)

