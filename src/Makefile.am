lib_LTLIBRARIES = libmeep@MEEP_SUFFIX@.la
include_HEADERS = meep.hpp
pkginclude_HEADERS = meep/mympi.hpp meep/vec.hpp

AM_CPPFLAGS = -I$(top_srcdir)/src

BUILT_SOURCES = step_d.hpp step_h.hpp update_from_e.hpp	\
update_e_from_d_prepare.hpp update_e_from_d_update.hpp sphere-quad.h

HDRS = meep.hpp meep_internals.hpp meep/mympi.hpp ran.hpp meep/vec.hpp	\
bicgstab.hpp

libmeep@MEEP_SUFFIX@_la_SOURCES = anisotropic_averaging.cpp bands.cpp	\
boundaries.cpp bicgstab.cpp control_c.cpp cw_fields.cpp dft.cpp		\
energy_and_flux.cpp fields.cpp loop_in_chunks.cpp grace.cpp		\
h5fields.cpp h5file.cpp initialize.cpp integrate.cpp monitor.cpp	\
mympi.cpp output_directory.cpp polarization.cpp ran.cpp                 \
sources.cpp step.cpp step_d.cpp step_h.cpp structure.cpp time.cpp	\
update_e_from_d.cpp update_from_e.cpp vec.cpp $(HDRS) $(BUILT_SOURCES)	\
update_e_from_d_prepare.cpp update_e_from_d_sources.cpp			\
update_e_from_d_update.cpp

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

step_d.hpp: $(HSRC)/step_d_gen
	(echo $(PRELUDE); echo; $(HSRC)/step_d_gen) > $@

step_h.hpp: $(HSRC)/step_h_gen
	(echo $(PRELUDE); echo; $(HSRC)/step_h_gen) > $@

update_e_from_d_prepare.hpp: $(HSRC)/update_e_from_d_gen
	(echo $(PRELUDE); echo; $(HSRC)/update_e_from_d_gen prepare) > $@

update_e_from_d_update.hpp: $(HSRC)/update_e_from_d_gen
	(echo $(PRELUDE); echo; $(HSRC)/update_e_from_d_gen update) > $@

update_from_e.hpp: $(HSRC)/update_from_e_gen
	(echo $(PRELUDE); echo; $(HSRC)/update_from_e_gen) > $@

maintainer-clean-local:
	rm -f $(BUILT_SOURCES)

endif # MAINTAINER_MODE
