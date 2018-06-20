SRC = aniso_disp.cpp bench.cpp bragg_transmission.cpp			\
convergence_cyl_waveguide.cpp cylindrical.cpp flux.cpp harmonics.cpp	\
integrate.cpp known_results.cpp near2far.cpp one_dimensional.cpp	\
physical.cpp stress_tensor.cpp symmetry.cpp three_d.cpp			\
two_dimensional.cpp 2D_convergence.cpp h5test.cpp pml.cpp

EXTRA_DIST = $(SRC)

OBJ = $(SRC:.cpp=.o)
DAC = $(SRC:.cpp=.dac)
DONE = $(SRC:.cpp=.done)

LIBMEEP = $(top_builddir)/src/libmeep.la
AM_CPPFLAGS = -I$(top_srcdir)/src
LDADD = $(LIBMEEP)

.SUFFIXES = .dac .done

check_PROGRAMS = aniso_disp bench bragg_transmission convergence_cyl_waveguide cylindrical flux harmonics integrate known_results near2far one_dimensional physical stress_tensor symmetry three_d two_dimensional 2D_convergence h5test pml

aniso_disp_SOURCES = aniso_disp.cpp

bench_SOURCES = bench.cpp

bragg_transmission_SOURCES = bragg_transmission.cpp

convergence_cyl_waveguide_SOURCES = convergence_cyl_waveguide.cpp

cylindrical_SOURCES = cylindrical.cpp

flux_SOURCES = flux.cpp

harmonics_SOURCES = harmonics.cpp

integrate_SOURCES = integrate.cpp

known_results_SOURCES = known_results.cpp

near2far_SOURCES = near2far.cpp

one_dimensional_SOURCES = one_dimensional.cpp

physical_SOURCES = physical.cpp

stress_tensor_SOURCES = stress_tensor.cpp

symmetry_SOURCES = symmetry.cpp

three_d_SOURCES = three_d.cpp

two_dimensional_SOURCES = two_dimensional.cpp

2D_convergence_SOURCES = 2D_convergence.cpp

h5test_SOURCES = h5test.cpp

pml_SOURCES = pml.cpp

TESTS = aniso_disp bench bragg_transmission convergence_cyl_waveguide cylindrical flux harmonics integrate known_results near2far one_dimensional physical stress_tensor symmetry three_d two_dimensional 2D_convergence h5test pml 

if WITH_MPI
 LOG_COMPILER = $(RUNCODE)
endif
#LOG_COMPILER = $(RUNCODE)

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
	rm -f *.h5
