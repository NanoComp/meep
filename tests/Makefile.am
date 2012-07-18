SRC = aniso_disp.cpp bench.cpp bragg_transmission.cpp			\
convergence_cyl_waveguide.cpp cylindrical.cpp flux.cpp harmonics.cpp	\
integrate.cpp known_results.cpp one_dimensional.cpp physical.cpp	\
symmetry.cpp three_d.cpp two_dimensional.cpp 2D_convergence.cpp		\
h5test.cpp pml.cpp

EXTRA_DIST = $(SRC)

OBJ = $(SRC:.cpp=.o)
DAC = $(SRC:.cpp=.dac)
DONE = $(SRC:.cpp=.done)

LIBMEEP = $(top_builddir)/src/libmeep@MEEP_SUFFIX@.la
AM_CPPFLAGS = -I$(top_srcdir)/src

.SUFFIXES = .dac .done

check_PROGRAMS = aniso_disp bench bragg_transmission convergence_cyl_waveguide cylindrical flux harmonics integrate known_results one_dimensional physical symmetry three_d two_dimensional 2D_convergence h5test pml

aniso_disp_SOURCES = aniso_disp.cpp
aniso_disp_LDADD = $(LIBMEEP)

bench_SOURCES = bench.cpp
bench_LDADD = $(LIBMEEP)

bragg_transmission_SOURCES = bragg_transmission.cpp
bragg_transmission_LDADD = $(LIBMEEP)

convergence_cyl_waveguide_SOURCES = convergence_cyl_waveguide.cpp
convergence_cyl_waveguide_LDADD = $(LIBMEEP)

cylindrical_SOURCES = cylindrical.cpp
cylindrical_LDADD = $(LIBMEEP)

flux_SOURCES = flux.cpp
flux_LDADD = $(LIBMEEP)

harmonics_SOURCES = harmonics.cpp
harmonics_LDADD = $(LIBMEEP)

integrate_SOURCES = integrate.cpp
integrate_LDADD = $(LIBMEEP)

known_results_SOURCES = known_results.cpp
known_results_LDADD = $(LIBMEEP)

one_dimensional_SOURCES = one_dimensional.cpp
one_dimensional_LDADD = $(LIBMEEP)

physical_SOURCES = physical.cpp
physical_LDADD = $(LIBMEEP)

symmetry_SOURCES = symmetry.cpp
symmetry_LDADD = $(LIBMEEP)

three_d_SOURCES = three_d.cpp
three_d_LDADD = $(LIBMEEP)

two_dimensional_SOURCES = two_dimensional.cpp
two_dimensional_LDADD = $(LIBMEEP)

2D_convergence_SOURCES = 2D_convergence.cpp
2D_convergence_LDADD = $(LIBMEEP)

h5test_SOURCES = h5test.cpp
h5test_LDADD = $(LIBMEEP)

pml_SOURCES = pml.cpp
pml_LDADD = $(LIBMEEP)

TESTS = aniso_disp bench bragg_transmission convergence_cyl_waveguide cylindrical flux harmonics integrate known_results one_dimensional physical symmetry three_d two_dimensional 2D_convergence h5test pml

LOG_COMPILER = $(RUNCODE)

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
