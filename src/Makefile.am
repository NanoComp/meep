lib_LTLIBRARIES = libmeep@MEEP_SUFFIX@.la
include_HEADERS = meep.hpp
pkginclude_HEADERS = meep/mympi.hpp meep/vec.hpp

AM_CPPFLAGS = -I$(top_srcdir)/src

BUILT_SOURCES = update_from_e.hpp sphere-quad.h

HDRS = meep.hpp meep_internals.hpp meep/mympi.hpp meep/vec.hpp	\
bicgstab.hpp

libmeep@MEEP_SUFFIX@_la_SOURCES = anisotropic_averaging.cpp bands.cpp	\
boundaries.cpp bicgstab.cpp control_c.cpp cw_fields.cpp dft.cpp		\
energy_and_flux.cpp fields.cpp loop_in_chunks.cpp grace.cpp		\
h5fields.cpp h5file.cpp initialize.cpp integrate.cpp monitor.cpp	\
mympi.cpp output_directory.cpp polarization.cpp sources.cpp step.cpp	\
step_db.cpp structure.cpp time.cpp update_e_from_d.cpp	\
update_from_e.cpp vec.cpp step_generic.cpp update_h_from_b.cpp $(HDRS) $(BUILT_SOURCES)

libmeep@MEEP_SUFFIX@_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@

noinst_PROGRAMS = sphere_quad
sphere_quad_SOURCES = sphere-quad.cpp

PRELUDE = "/* This file was automatically generated --- DO NOT EDIT */"

SPHERE_QUAD = $(top_builddir)/src/sphere_quad
sphere-quad.h: sphere_quad$(EXEEXT)
	(echo $(PRELUDE); echo; $(SPHERE_QUAD)) > $@

# Generate built sources only in maintainer mode
if MAINTAINER_MODE

HSRC=$(top_builddir)/hsrc

update_from_e.hpp: $(HSRC)/update_from_e_gen
	(echo $(PRELUDE); echo; $(HSRC)/update_from_e_gen) > $@

maintainer-clean-local:
	rm -f $(BUILT_SOURCES)

endif # MAINTAINER_MODE
