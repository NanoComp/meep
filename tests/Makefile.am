SRC = aniso_disp.cpp bench.cpp bragg_transmission.cpp			\
convergence_cyl_waveguide.cpp cylindrical.cpp dump_load.cpp flux.cpp    \
harmonics.cpp integrate.cpp known_results.cpp near2far.cpp              \
one_dimensional.cpp physical.cpp stress_tensor.cpp symmetry.cpp 	\
three_d.cpp two_dimensional.cpp 2D_convergence.cpp h5test.cpp pml.cpp

EXTRA_DIST = $(SRC)

OBJ = $(SRC:.cpp=.o)
DAC = $(SRC:.cpp=.dac)
DONE = $(SRC:.cpp=.done)

LIBMEEP = $(top_builddir)/src/libmeep.la
MEEPLIBS = $(LIBMEEP)
AM_CPPFLAGS = -I$(top_srcdir)/src -DDATADIR=\"$(srcdir)/\"

.SUFFIXES = .dac .done

check_PROGRAMS = aniso_disp bench bragg_transmission convergence_cyl_waveguide cylindrical dump_load flux harmonics integrate known_results near2far one_dimensional physical stress_tensor symmetry three_d two_dimensional 2D_convergence h5test pml pw-source-ll ring-ll cyl-ellipsoid-ll absorber-1d-ll array-slice-ll user-defined-material dft-fields gdsII-3d bend-flux-ll array-metadata

array_metadata_SOURCES = array-metadata.cpp
array_metadata_LDADD   = $(MEEPLIBS)

aniso_disp_SOURCES = aniso_disp.cpp
aniso_disp_LDADD = $(MEEPLIBS)

bench_SOURCES = bench.cpp
bench_LDADD = $(MEEPLIBS)

bragg_transmission_SOURCES = bragg_transmission.cpp
bragg_transmission_LDADD = $(MEEPLIBS)

convergence_cyl_waveguide_SOURCES = convergence_cyl_waveguide.cpp
convergence_cyl_waveguide_LDADD = $(MEEPLIBS)

cylindrical_SOURCES = cylindrical.cpp
cylindrical_LDADD = $(MEEPLIBS)

dump_load_SOURCES = dump_load.cpp
dump_load_LDADD = $(MEEPLIBS)

flux_SOURCES = flux.cpp
flux_LDADD = $(MEEPLIBS)

harmonics_SOURCES = harmonics.cpp
harmonics_LDADD = $(MEEPLIBS)

integrate_SOURCES = integrate.cpp
integrate_LDADD = $(MEEPLIBS)

known_results_SOURCES = known_results.cpp
known_results_LDADD = $(MEEPLIBS)

near2far_SOURCES = near2far.cpp
near2far_LDADD = $(MEEPLIBS)

one_dimensional_SOURCES = one_dimensional.cpp
one_dimensional_LDADD = $(MEEPLIBS)

physical_SOURCES = physical.cpp
physical_LDADD = $(MEEPLIBS)

stress_tensor_SOURCES = stress_tensor.cpp
stress_tensor_LDADD = $(MEEPLIBS)

symmetry_SOURCES = symmetry.cpp
symmetry_LDADD = $(MEEPLIBS)

three_d_SOURCES = three_d.cpp
three_d_LDADD = $(MEEPLIBS)

two_dimensional_SOURCES = two_dimensional.cpp
two_dimensional_LDADD = $(MEEPLIBS)

2D_convergence_SOURCES = 2D_convergence.cpp
2D_convergence_LDADD = $(MEEPLIBS)

h5test_SOURCES = h5test.cpp
h5test_LDADD = $(MEEPLIBS)

pml_SOURCES = pml.cpp
pml_LDADD = $(MEEPLIBS)

absorber_1d_ll_SOURCES = absorber-1d-ll.cpp
absorber_1d_ll_LDADD   = $(MEEPLIBS)

cyl_ellipsoid_ll_SOURCES = cyl-ellipsoid-ll.cpp
cyl_ellipsoid_ll_LDADD   = $(MEEPLIBS)

pw_source_ll_SOURCES = pw-source-ll.cpp
pw_source_ll_LDADD   = $(MEEPLIBS)

ring_ll_SOURCES = ring-ll.cpp
ring_ll_LDADD   = $(MEEPLIBS)

array_slice_ll_SOURCES = array-slice-ll.cpp
array_slice_ll_LDADD   = $(MEEPLIBS)

dft_fields_SOURCES = dft-fields.cpp
dft_fields_LDADD   = $(MEEPLIBS)

gdsII_3d_SOURCES = gdsII-3d.cpp
gdsII_3d_LDADD   = $(MEEPLIBS)

bend_flux_ll_SOURCES = bend-flux-ll.cpp
bend_flux_ll_LDADD   = $(MEEPLIBS)

user_defined_material_SOURCES = user-defined-material.cpp
user_defined_material_LDADD   = $(MEEPLIBS)

dist_noinst_DATA = cyl-ellipsoid-eps-ref.h5 array-slice-ll-ref.h5 gdsII-3d.gds

TESTS = aniso_disp bench bragg_transmission convergence_cyl_waveguide cylindrical dump_load flux harmonics integrate known_results near2far one_dimensional physical stress_tensor symmetry three_d two_dimensional 2D_convergence h5test pml

if WITH_MPI
  LOG_COMPILER = $(RUNCODE)
else
  TESTS += cyl-ellipsoid-ll array-slice-ll ring-ll
endif

# Note: this requires GNU make

.o.dac:
	$(CXXLINK) $< $(LIBMEEP) $(LIBS)

$(DAC): $(LIBMEEP)

.dac.done:
	$(RUNCODE) ./$<
	touch $@

benchmark: bench
	$(RUNCODE) ./bench

dac: $(DAC)

clean-local::
	rm -f *.o *.dac debug_out_* *.done

distclean-local:
	rm -f $(shell ls *.h5 | sed '/.*ref.*/d')
