# recompile Haskell source only in maintainer mode
if MAINTAINER_MODE
HSRC = hsrc
else
HSRC =
endif

SUBDIRS = $(HSRC) src doc tests

EXTRA_DIST = COPYRIGHT m4 meep.pc.in

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = meep.pc

if MAINTAINER_MODE
generator:
	cd hsrc && $(MAKE) depend && $(MAKE)
else
generator:
endif

libmeep: generator
	cd src && $(MAKE)

# add a fake dependency to make sure stuff is built before 'make dist'
COPYRIGHT: generator libmeep

