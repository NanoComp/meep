lib_LTLIBRARIES = libmeep.la

BUILT_SOURCES = step_d.h step_h.h update_from_e.h	\
	sphere-quad.h \
	update_e_from_d_prepare.h update_e_from_d_sources.h update_e_from_d_update.h

HDRS = meep.h meep_internals.h mympi.h ran.h threevec.h vec.h

libmeep_la_SOURCES = anisotropic_averaging.cpp bands.cpp		\
boundaries.cpp control_c.cpp dft.cpp energy_and_flux.cpp fields.cpp	\
fields_integrate.cpp grace.cpp h5fields.cpp h5file.cpp initialize.cpp	\
monitor.cpp mympi.cpp output_directory.cpp polarization.cpp ran.cpp	\
slices.cpp sources.cpp step.cpp step_d.cpp step_h.cpp structure.cpp	\
time.cpp update_e_from_d.cpp update_from_e.cpp vec.cpp $(HDRS)		\
$(BUILT_SOURCES) update_e_from_d_prepare.cpp				\
update_e_from_d_sources.cpp update_e_from_d_update.cpp

libmeep_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@

noinst_PROGRAMS = sphere_quad
sphere_quad_SOURCES = sphere-quad.cpp

PRELUDE = "/* This file was automatically generated --- DO NOT EDIT */"

SPHERE_QUAD = $(top_builddir)/src/sphere_quad
sphere-quad.h: sphere_quad
	(echo $(PRELUDE); echo; $(SPHERE_QUAD)) > $@

# Generate built sources only in maintainer mode
if MAINTAINER_MODE

HSRC=$(top_builddir)/hsrc

step_d.h: $(HSRC)/step_d_gen
	(echo $(PRELUDE); echo; $(HSRC)/step_d_gen) > $@

step_h.h: $(HSRC)/step_h_gen
	(echo $(PRELUDE); echo; $(HSRC)/step_h_gen) > $@

update_e_from_d_prepare.h: $(HSRC)/update_e_from_d_gen
	(echo $(PRELUDE); echo; $(HSRC)/update_e_from_d_gen prepare) > $@

update_e_from_d_sources.h: $(HSRC)/update_e_from_d_gen
	(echo $(PRELUDE); echo; $(HSRC)/update_e_from_d_gen sources) > $@

update_e_from_d_update.h: $(HSRC)/update_e_from_d_gen
	(echo $(PRELUDE); echo; $(HSRC)/update_e_from_d_gen update) > $@

update_from_e.h: $(HSRC)/update_from_e_gen
	(echo $(PRELUDE); echo; $(HSRC)/update_from_e_gen) > $@

maintainer-clean-local:
	rm -f $(BUILT_SOURCES)

endif # MAINTAINER_MODE
