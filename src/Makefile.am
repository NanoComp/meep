lib_LTLIBRARIES = libmeep.la

BUILT_SOURCES = step_d.h step_h.h update_e_from_d.h update_from_e.h

HDRS = harminv.h meep.h meep_internals.h mympi.h ran.h threevec.h	\
vec.h

libmeep_la_SOURCES = anisotropic_averaging.cpp bands.cpp		\
boundaries.cpp control_c.cpp energy_and_flux.cpp fields.cpp fluxes.cpp	\
grace.cpp harminv.cpp initialize.cpp material.cpp monitor.cpp		\
mympi.cpp output_directory.cpp polarization.cpp ran.cpp slices.cpp	\
sources.cpp step.cpp step_d.cpp step_h.cpp time.cpp			\
update_e_from_d.cpp update_from_e.cpp vec.cpp $(HDRS) $(BUILT_SOURCES)

libmeep_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@

# Generate built sources only in maintainer mode
if MAINTAINER_MODE

HSRC=$(top_builddir)/hsrc

PRELUDE = "/* This file was automatically generated --- DO NOT EDIT */"

step_d.h: $(HSRC)/step_d_gen
	(echo $(PRELUDE); echo; $(HSRC)/step_d_gen) > $@

step_h.h: $(HSRC)/step_h_gen
	(echo $(PRELUDE); echo; $(HSRC)/step_h_gen) > $@

update_e_from_d.h: $(HSRC)/update_e_from_d_gen
	(echo $(PRELUDE); echo; $(HSRC)/update_e_from_d_gen) > $@

update_from_e.h: $(HSRC)/update_from_e_gen
	(echo $(PRELUDE); echo; $(HSRC)/update_from_e_gen) > $@

maintainer-clean-local:
	rm -f $(BUILT_SOURCES)

endif # MAINTAINER_MODE
