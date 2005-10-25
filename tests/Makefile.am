SRC = bench.cpp bragg_transmission.cpp convergence_cyl_waveguide.cpp	\
cylindrical.cpp flux.cpp integrate.cpp known_results.cpp		\
one_dimensional.cpp physical.cpp symmetry.cpp three_d.cpp		\
two_dimensional.cpp 2D_convergence.cpp h5test.cpp

EXTRA_DIST = $(SRC)

OBJ = $(SRC:.cpp=.o)
DAC = $(SRC:.cpp=.dac)
DONE = $(SRC:.cpp=.done)

LIBMEEP = $(top_builddir)/src/libmeep@MEEP_SUFFIX@.la
AM_CPPFLAGS = -I$(top_srcdir)/src

.SUFFIXES = .dac .done

noinst_PROGRAMS = bench

bench_SOURCES = bench.cpp
bench_LDADD = $(LIBMEEP)

compile-tests: $(DAC)

check-local: $(DONE)

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
