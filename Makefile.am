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

dist-hook: $(HSRC) src

