lib_LTLIBRARIES = libmeep@MEEP_SUFFIX@-ctl.la
bin_PROGRAMS = meep@MEEP_SUFFIX@
EXTRA_DIST = meep.scm.in meep-enums.scm casimir.scm meep_wrap.patch

LIBMEEP = $(top_builddir)/src/libmeep@MEEP_SUFFIX@.la

HDRS = meep-ctl.hpp meep-ctl-const.hpp meep-ctl-swig.hpp

LIBHDRS = $(top_srcdir)/src/meep.hpp $(top_srcdir)/src/meep/vec.hpp $(top_srcdir)/src/meep/mympi.hpp
CTLHDRS = $(LIBHDRS) $(srcdir)/meep-ctl.hpp $(top_builddir)/config.h $(srcdir)/meep-ctl-const.hpp ctl-io.h

libmeep@MEEP_SUFFIX@_ctl_la_SOURCES = ctl-io.cpp geom.cpp meep.cpp structure.cpp ctl-io.h ctl-io.i $(HDRS) meep.i meep_op_renames.i meep_renames.i meep_enum_renames.i meep_swig_bug_workaround.i

nodist_meep@MEEP_SUFFIX@_SOURCES = main.cpp
meep@MEEP_SUFFIX@_LDADD = libmeep@MEEP_SUFFIX@-ctl.la $(LIBMEEP) @LIBCTL_LIBS@

BUILT_SOURCES = $(nodist_meep_SOURCES) meep_renames.i meep_enum_renames.i meep_swig_bug_workaround.i meep-enums.scm meep_wrap.cxx

# manual dependency hack to force the build order in certain cases
# ... for some reason automake's automatic dependencies are not working here
meep.o: meep.cpp meep-enums.scm $(CTLHDRS)
structure.o: structure.cpp $(CTLHDRS)
meep_wrap.o: meep_wrap.cxx $(srcdir)/meep-ctl-swig.hpp $(CTLHDRS)

if WITH_LIBCTL
if MAINTAINER_MODE
meep_wrap.cxx: meep.i meep_op_renames.i meep_enum_renames.i meep_renames.i ctl-io.i meep-ctl-swig.hpp meep_swig_bug_workaround.i $(LIBHDRS)
	swig -I$(top_srcdir)/src -c++ -guile -o $@ meep.i
	patch -p0 $@ < $(srcdir)/meep_wrap.patch
else
meep_wrap.cxx:
	echo "#error need --with-maintainer-mode to generate this file" 1>&2
	exit 1
endif
else
meep_wrap.cxx:
	echo "#error need --with-libctl to generate this file" > $@
endif

# workaround missing namespace prefix in swig
meep_renames.i: $(LIBHDRS)
	(echo "// AUTOMATICALLY GENERATED -- DO NOT EDIT"; sed 's/^ *class \+\([A-Za-z_0-9:]*\)\( *\| *:[^{]*\){.*$$/%rename(meep_\1) meep::\1;/' $(LIBHDRS) | grep "%rename" | sort -u; echo; grep -hv typedef $(LIBHDRS) | sed 's/\(inline\|const\|extern\|static\) \+//g' | sed 's/^[A-Za-z_0-9:<>]\+[* ]\+\([A-Za-z_0-9:]*\) *(.*$$/%rename(meep_\1) meep::\1;/' | grep "%rename" | sort -u; ) > $@

# work around bug in swig, where it doesn't prepend namespace to friend funcs
meep_swig_bug_workaround.i: $(LIBHDRS)
	(echo "// AUTOMATICALLY GENERATED -- DO NOT EDIT"; grep -h friend $(LIBHDRS) | sed 's/^ *friend \+[A-Za-z_0-9:<>]\+[* ]\+\([A-Za-z_0-9:]*\) *(.*$$/%ignore \1;/' | grep "%ignore" | sort -u;) > $@

meep_enum_renames.i: $(LIBHDRS)
	(echo "// AUTOMATICALLY GENERATED -- DO NOT EDIT"; for f in $(LIBHDRS); do egrep "^enum" $$f | sed 's/enum \+\([A-Za-z_0-9:]\+\).*$$/\1/g' | while read enum; do cat $$f | tr -d '\n' | sed 's/.*enum \+'$${enum}' *{\([^}]*\)}.*/\1/g' | sed 's/= *[0-9]\+//g' |tr -d ' \t' | tr ',' '\n' | sed 's/^.*$$/'"%rename(meep_$${enum}_\0) meep::\0;/g"; echo; done; done;) > $@

meep-enums.scm: meep_enum_renames.i
	(sed 's/%rename(\([A-Za-z0-9_]*\)) *\([A-Za-z0-9:_]*\);$$/(define \2 (\1))/' meep_enum_renames.i | sed 's/[A-Za-z0-9:_]*:://g' | sed 's/_/-/g' | sed 's,//,;,'; echo "(define Centered Dielectric)") > $@

##############################################################################
#      stuff below is based on example Makefile.in included with libctl
##############################################################################

# what is printed out when invoking your program with --version:
VERSION_STRING = "Meep @VERSION@, Copyright (C) 2005-2015 Massachusetts Insitute of Technology."

MY_DEFS = -DHAVE_CTL_HOOKS=1 -DHAVE_CTL_EXPORT_HOOK=1

PROGRAM_NAME = meep
SPECIFICATION_FILE = $(PROGRAM_NAME).scm
nodist_pkgdata_DATA = $(SPECIFICATION_FILE) meep-enums.scm casimir.scm

##############################################################################

AM_CPPFLAGS = -I$(top_srcdir)/src $(MY_DEFS) $(CTL_DEFS) $(CTL_H_CPPFLAG)

# libctl install. dir., e.g. /usr/local/share/libctl
LIBCTL_DIR = @LIBCTL_DIR@

# gen-ctl-io program
GEN_CTL_IO = @GEN_CTL_IO@

##############################################################################

CTL_DEFS = -DCTL_SCM='"'$(LIBCTL_DIR)/base/ctl.scm'"' \
       -DINCLUDE_SCM='"'$(LIBCTL_DIR)/base/include.scm'"' \
          -DSPEC_SCM='"'$(pkgdatadir)/$(SPECIFICATION_FILE)'"' \
    -DVERSION_STRING='"'$(VERSION_STRING)'"'

##############################################################################

ctl-io.cpp: $(SPECIFICATION_FILE) $(LIBCTL_DIR)/utils/geom.scm ctl-io.h
	$(GEN_CTL_IO) --cxx --code -o $@ $(SPECIFICATION_FILE) $(LIBCTL_DIR)

ctl-io.h: $(SPECIFICATION_FILE) $(LIBCTL_DIR)/utils/geom.scm
	$(GEN_CTL_IO) --cxx --header -o $@ $(SPECIFICATION_FILE) $(LIBCTL_DIR)

ctl-io.i: $(SPECIFICATION_FILE) $(LIBCTL_DIR)/utils/geom.scm
	$(GEN_CTL_IO) --cxx --swig -o $@ $(SPECIFICATION_FILE) $(LIBCTL_DIR)

main.cpp: $(LIBCTL_DIR)/base/main.c ctl-io.h
	cp -f $(LIBCTL_DIR)/base/main.c $@

geom.cpp: $(LIBCTL_DIR)/utils/geom.c ctl-io.h
	cp -f $(LIBCTL_DIR)/utils/geom.c $@

clean-local:
	rm -f $(nodist_meep_SOURCES) meep-enums.scm ctl-io.* main.* geom.*

maintainer-clean-local:
	rm -f $(BUILT_SOURCES)
