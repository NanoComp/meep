bin_PROGRAMS = meep
EXTRA_DIST = meep.scm.in

LIBMEEP = $(top_builddir)/src/libmeep.la

HDRS = meep-ctl.h meep-ctl-const.h

LIBHDRS = $(top_srcdir)/src/meep.h $(top_srcdir)/src/meep/vec.h $(top_srcdir)/src/meep/mympi.h

meep_SOURCES = meep.cpp structure.cpp meep_wrap.cxx $(HDRS) meep.i
nodist_meep_SOURCES = main.cpp geom.cpp ctl-io.cpp ctl-io.h meep_renames.i
meep_LDADD = $(LIBMEEP) @LIBCTL_LIBS@

BUILT_SOURCES = $(nodist_meep_SOURCES) meep_wrap.cxx

if WITH_LIBCTL
meep_wrap.cxx: meep.i meep_renames.i ctl-io.i meep-ctl-swig.h $(LIBHDRS)
	swig -I$(top_srcdir)/src -c++ -guile -o $@ meep.i
else
meep_wrap.cxx:
	echo "/* need --with-libctl to generate this file */" > $@
endif

meep_renames.i: $(LIBHDRS)
	(echo "// AUTOMATICALLY GENERATED -- DO NOT EDIT"; sed 's/^ *class *\([A-Za-z_0-9:]*\) *{ *$$/%rename(meep_\1) meep::\1;/' $(LIBHDRS) | grep "%rename") > $@

##############################################################################
#      stuff below is based on example Makefile.in included with libctl
##############################################################################

# what is printed out when invoking your program with --version:
VERSION_STRING = "Meep @VERSION@, Copyright (C) 2003 Massachusetts Insitute of Technology."

MY_DEFS = -DHAVE_CTL_HOOKS=1 -DHAVE_CTL_EXPORT_HOOK=1

PROGRAM_NAME = meep
SPECIFICATION_FILE = $(PROGRAM_NAME).scm

##############################################################################

AM_CPPFLAGS = -I$(top_srcdir)/src $(MY_DEFS) $(CTL_DEFS)

# libctl install. dir., e.g. /usr/local/share/libctl
LIBCTL_DIR = @LIBCTL_DIR@

# gen-ctl-io program
GEN_CTL_IO = @GEN_CTL_IO@

##############################################################################

CTL_DEFS = -DCTL_SCM='"'$(LIBCTL_DIR)/base/ctl.scm'"' \
       -DINCLUDE_SCM='"'$(LIBCTL_DIR)/base/include.scm'"' \
          -DSPEC_SCM='"'$(prefix)/share/libctl/specs/$(SPECIFICATION_FILE)'"' \
    -DVERSION_STRING='"'$(VERSION_STRING)'"'

##############################################################################

ctl-io.cpp ctl-io.h ctl-io.i: $(SPECIFICATION_FILE)
	$(GEN_CTL_IO) --cxx --swig $(srcdir)/$(SPECIFICATION_FILE) $(LIBCTL_DIR)

main.cpp: $(LIBCTL_DIR)/base/main.c ctl-io.h
	cp -f $(LIBCTL_DIR)/base/main.c $@

geom.cpp: $(LIBCTL_DIR)/utils/geom.c ctl-io.h
	cp -f $(LIBCTL_DIR)/utils/geom.c $@

clean-local:
	rm -f ctl-io.* main.* geom.*
