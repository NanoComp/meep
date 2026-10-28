##################################################
# description of unit tests
##################################################
TEST_DIR = tests
EXAMPLE_DIR = examples

# material_dispersion.py test must be excluded from test suite for MPI build
MATERIAL_DISPERSION_TEST = $(TEST_DIR)/test_material_dispersion.py

if WITH_MPI
  MDPYTEST=
else
  MDPYTEST=$(MATERIAL_DISPERSION_TEST)
endif # WITH_MPI

MPB_TESTS =                                    \
    $(TEST_DIR)/test_binary_grating.py         \
    $(TEST_DIR)/test_diffracted_planewave.py   \
    $(TEST_DIR)/test_dispersive_eigenmode.py   \
    $(TEST_DIR)/test_kdom.py                   \
    $(TEST_DIR)/test_mode_coeffs.py            \
    $(TEST_DIR)/test_mode_decomposition.py     \
    $(TEST_DIR)/test_mpb.py                    \
    $(TEST_DIR)/test_wvg_src.py

if WITH_MPB
  MPBPYTEST= $(MPB_TESTS)
else
  MPBPYTEST =
endif

ADJOINT_TESTS =                                \
    $(TEST_DIR)/test_adjoint_solver.py         \
    $(TEST_DIR)/test_adjoint_utils.py          \
    $(TEST_DIR)/test_adjoint_cyl.py            \
    $(TEST_DIR)/test_adjoint_jax.py

TESTS =                                        \
    $(TEST_DIR)/test_3rd_harm_1d.py            \
    $(TEST_DIR)/test_absorber_1d.py            \
    $(TEST_DIR)/test_adjoint_chunks.py         \
    $(TEST_DIR)/test_antenna_radiation.py      \
    $(TEST_DIR)/test_array_metadata.py         \
    $(TEST_DIR)/test_bend_flux.py              \
    $(TEST_DIR)/test_binary_partition_utils.py \
    $(TEST_DIR)/test_boundaries_1D.py          \
    $(TEST_DIR)/test_cavity_arrayslice.py      \
    $(TEST_DIR)/test_cavity_farfield.py        \
    $(TEST_DIR)/test_chunk_balancer.py         \
    $(TEST_DIR)/test_chunk_layout.py           \
    $(TEST_DIR)/test_chunks.py                 \
    $(TEST_DIR)/test_conductivity.py           \
    $(TEST_DIR)/test_cyl_ellipsoid.py          \
    $(TEST_DIR)/test_dft_energy.py             \
    $(TEST_DIR)/test_dft_fields.py             \
    $(TEST_DIR)/test_divide_mpi_processes.py   \
    $(TEST_DIR)/test_dump_load.py              \
    $(TEST_DIR)/test_eigfreq.py                \
    $(TEST_DIR)/test_faraday_rotation.py       \
    $(TEST_DIR)/test_field_functions.py        \
    $(TEST_DIR)/test_force.py                  \
    $(TEST_DIR)/test_fragment_stats.py         \
    $(TEST_DIR)/test_gaussianbeam.py           \
    $(TEST_DIR)/test_geom.py                   \
    $(TEST_DIR)/test_get_point.py              \
    $(TEST_DIR)/test_get_epsilon_grid.py       \
    $(TEST_DIR)/test_holey_wvg_bands.py        \
    $(TEST_DIR)/test_holey_wvg_cavity.py       \
    $(TEST_DIR)/test_integrated_source.py      \
    $(TEST_DIR)/test_ldos.py                   \
    $(MDPYTEST)                                \
    $(TEST_DIR)/test_material_grid.py          \
    $(TEST_DIR)/test_materials_library.py      \
    $(TEST_DIR)/test_medium_evaluations.py     \
    $(TEST_DIR)/test_mesh.py                   \
    $(TEST_DIR)/test_multilevel_atom.py        \
    $(TEST_DIR)/test_n2f_periodic.py           \
    $(TEST_DIR)/test_oblique_source.py         \
    $(TEST_DIR)/test_pml_cyl.py                \
    $(TEST_DIR)/test_physical.py               \
    $(TEST_DIR)/test_planewave_1D.py           \
    $(TEST_DIR)/test_prism.py                  \
    $(TEST_DIR)/test_pw_source.py              \
    $(TEST_DIR)/test_refl_angular.py           \
    $(TEST_DIR)/test_ring.py                   \
    $(TEST_DIR)/test_ring_cyl.py               \
    $(TEST_DIR)/test_simulation.py             \
    $(TEST_DIR)/test_special_kz.py             \
    $(TEST_DIR)/test_source.py                 \
    $(TEST_DIR)/test_stop_when_flux_decayed.py \
    $(TEST_DIR)/test_subpixel_3d.py            \
    $(TEST_DIR)/test_timing_measurements.py    \
    $(TEST_DIR)/test_user_defined_material.py  \
    $(TEST_DIR)/test_verbosity_mgr.py          \
    $(TEST_DIR)/test_visualization.py          \
    $(MPBPYTEST)

# data files read by the tests
TEST_DATA_FILES =                                           \
    $(TEST_DIR)/data/amp_func_file.h5                       \
    $(TEST_DIR)/data/bragg-h.k01.b08.y.tm.h5                \
    $(TEST_DIR)/data/cavity-farfield-4-freqs.h5             \
    $(TEST_DIR)/data/cavity-farfield.h5                     \
    $(TEST_DIR)/data/cavity_arrayslice_1d.npy               \
    $(TEST_DIR)/data/cavity_arrayslice_2d.npy               \
    $(TEST_DIR)/data/converted-diamond-dpwr.k06.b05.h5      \
    $(TEST_DIR)/data/converted-tri-rods-e.k11.b08.tm.h5     \
    $(TEST_DIR)/data/convex_prism_vertices.npz              \
    $(TEST_DIR)/data/diamond-dpwr.k06.b05.h5                \
    $(TEST_DIR)/data/eps_input_file_test.h5                 \
    $(TEST_DIR)/data/hole-slab-h.k01.b09.z.zeven.h5         \
    $(TEST_DIR)/data/line-defect-e.k04.b12.z.tm.h5          \
    $(TEST_DIR)/data/mpa_unfilter_design.npy                \
    $(TEST_DIR)/data/mpa_unfilter_design_target.npy         \
    $(TEST_DIR)/data/nonconvex_prism_vertices1.npz          \
    $(TEST_DIR)/data/nonconvex_prism_vertices2.npz          \
    $(TEST_DIR)/data/nonconvex_prism_vertices3.npz          \
    $(TEST_DIR)/data/strip-flux.v.k01.b04.x.h5              \
    $(TEST_DIR)/data/subpixel_avg-epsilon.h5                \
    $(TEST_DIR)/data/tri-rods-e.k11.b08.z.tm-r-m3-n32.h5    \
    $(TEST_DIR)/data/tri-rods-e.k11.b08.z.tm.h5             \
    $(TEST_DIR)/data/tri-rods-epsilon-r-m3-n32.h5           \
    $(TEST_DIR)/data/tutorial-b.k16.b08.te.h5               \
    $(TEST_DIR)/data/tutorial-C.k16.b08.te.h5               \
    $(TEST_DIR)/data/tutorial-d.k01.b02.te.h5               \
    $(TEST_DIR)/data/tutorial-d.k16.b02.te.h5               \
    $(TEST_DIR)/data/tutorial-d.k16.b08.te.h5               \
    $(TEST_DIR)/data/tutorial-e.k16.b08.te.h5               \
    $(TEST_DIR)/data/tutorial-e.k16.b08.z.tm.h5             \
    $(TEST_DIR)/data/tutorial-epsilon.h5                    \
    $(TEST_DIR)/data/tutorial-flux.v.k11.b08.te.h5          \
    $(TEST_DIR)/data/tutorial-h.k16.b08.te.h5               \
    $(TEST_DIR)/data/tutorial-mu.h5                         \
    $(TEST_DIR)/data/tutorial-te-eigenvectors-3-3.h5        \
    $(TEST_DIR)/data/tutorial-te-eigenvectors-8-1.h5        \
    $(TEST_DIR)/data/tutorial-te-eigenvectors.h5            \
    $(TEST_DIR)/data/tutorial-tot.rpwr.k16.b08.te.h5

# only the *.py examples are distributed: the notebooks (*.ipynb) and
# media/ are omitted to keep the tarball small
EXAMPLES =                                                            \
    $(EXAMPLE_DIR)/3rd-harm-1d.py                                     \
    $(EXAMPLE_DIR)/absorbed_power_density.py                          \
    $(EXAMPLE_DIR)/absorber-1d.py                                     \
    $(EXAMPLE_DIR)/adjoint_optimization/binary_grating_levelset.py    \
    $(EXAMPLE_DIR)/adjoint_optimization/mode_converter.py             \
    $(EXAMPLE_DIR)/adjoint_optimization/multilayer_opt.py             \
    $(EXAMPLE_DIR)/antenna-radiation.py                               \
    $(EXAMPLE_DIR)/antenna_pec_ground_plane.py                        \
    $(EXAMPLE_DIR)/antenna_pec_ground_plane_1D.py                     \
    $(EXAMPLE_DIR)/bend-flux.py                                       \
    $(EXAMPLE_DIR)/bent-waveguide.py                                  \
    $(EXAMPLE_DIR)/binary_grating.py                                  \
    $(EXAMPLE_DIR)/binary_grating_n2f.py                              \
    $(EXAMPLE_DIR)/binary_grating_oblique.py                          \
    $(EXAMPLE_DIR)/binary_grating_phasemap.py                         \
    $(EXAMPLE_DIR)/cavity-farfield.py                                 \
    $(EXAMPLE_DIR)/cavity_arrayslice.py                               \
    $(EXAMPLE_DIR)/cherenkov-radiation.py                             \
    $(EXAMPLE_DIR)/chirped_pulse.py                                   \
    $(EXAMPLE_DIR)/coupler.py                                         \
    $(EXAMPLE_DIR)/cyl-ellipsoid.py                                   \
    $(EXAMPLE_DIR)/cylinder_cross_section.py                          \
    $(EXAMPLE_DIR)/differential_cross_section.py                      \
    $(EXAMPLE_DIR)/diffracted_planewave.py                            \
    $(EXAMPLE_DIR)/dipole_in_vacuum_1D.py                             \
    $(EXAMPLE_DIR)/dipole_in_vacuum_cyl_off_axis.py                   \
    $(EXAMPLE_DIR)/dipole_in_vacuum_cyl_on_axis.py                    \
    $(EXAMPLE_DIR)/disc_extraction_efficiency.py                      \
    $(EXAMPLE_DIR)/disc_radiation_pattern.py                          \
    $(EXAMPLE_DIR)/edge_emitter_2D.py                                 \
    $(EXAMPLE_DIR)/edge_emitter_3D.py                                 \
    $(EXAMPLE_DIR)/eps_fit_lorentzian.py                              \
    $(EXAMPLE_DIR)/extraction_eff_ldos.py                             \
    $(EXAMPLE_DIR)/faraday-rotation.py                                \
    $(EXAMPLE_DIR)/finite_grating.py                                  \
    $(EXAMPLE_DIR)/gaussian-beam.py                                   \
    $(EXAMPLE_DIR)/grating2d_triangular_lattice.py                    \
    $(EXAMPLE_DIR)/holey-wvg-bands.py                                 \
    $(EXAMPLE_DIR)/holey-wvg-cavity.py                                \
    $(EXAMPLE_DIR)/material-dispersion.py                             \
    $(EXAMPLE_DIR)/metal-cavity-ldos.py                               \
    $(EXAMPLE_DIR)/metasurface_lens.py                                \
    $(EXAMPLE_DIR)/mie_scattering.py                                  \
    $(EXAMPLE_DIR)/mode-decomposition.py                              \
    $(EXAMPLE_DIR)/mode_coeff_phase.py                                \
    $(EXAMPLE_DIR)/mpb_bragg.py                                       \
    $(EXAMPLE_DIR)/mpb_bragg_sine.py                                  \
    $(EXAMPLE_DIR)/mpb_data_analysis.py                               \
    $(EXAMPLE_DIR)/mpb_diamond.py                                     \
    $(EXAMPLE_DIR)/mpb_hole_slab.py                                   \
    $(EXAMPLE_DIR)/mpb_honey_rods.py                                  \
    $(EXAMPLE_DIR)/mpb_line_defect.py                                 \
    $(EXAMPLE_DIR)/mpb_sq_rods.py                                     \
    $(EXAMPLE_DIR)/mpb_strip.py                                       \
    $(EXAMPLE_DIR)/mpb_tri_holes.py                                   \
    $(EXAMPLE_DIR)/mpb_tri_rods.py                                    \
    $(EXAMPLE_DIR)/mpb_tutorial.py                                    \
    $(EXAMPLE_DIR)/multilevel-atom.py                                 \
    $(EXAMPLE_DIR)/oblique-planewave.py                               \
    $(EXAMPLE_DIR)/oblique-source.py                                  \
    $(EXAMPLE_DIR)/parallel-wvgs-force.py                             \
    $(EXAMPLE_DIR)/parallel-wvgs-mpb.py                               \
    $(EXAMPLE_DIR)/perturbation_theory.py                             \
    $(EXAMPLE_DIR)/perturbation_theory_2d.py                          \
    $(EXAMPLE_DIR)/phase_in_material.py                               \
    $(EXAMPLE_DIR)/planar_cavity_ldos.py                              \
    $(EXAMPLE_DIR)/plot_radiation_pattern_dipole.py                   \
    $(EXAMPLE_DIR)/point_dipole_cyl.py                                \
    $(EXAMPLE_DIR)/polarization_grating.py                            \
    $(EXAMPLE_DIR)/pw-source.py                                       \
    $(EXAMPLE_DIR)/refl-angular-kz2d.py                               \
    $(EXAMPLE_DIR)/refl-angular.py                                    \
    $(EXAMPLE_DIR)/refl-quartz.py                                     \
    $(EXAMPLE_DIR)/ring-cyl.py                                        \
    $(EXAMPLE_DIR)/ring-mode-overlap.py                               \
    $(EXAMPLE_DIR)/ring.py                                            \
    $(EXAMPLE_DIR)/ring_gds.py                                        \
    $(EXAMPLE_DIR)/solve-cw.py                                        \
    $(EXAMPLE_DIR)/stochastic_emitter.py                              \
    $(EXAMPLE_DIR)/stochastic_emitter_line.py                         \
    $(EXAMPLE_DIR)/stochastic_emitter_reciprocity.py                  \
    $(EXAMPLE_DIR)/straight-waveguide.py                              \
    $(EXAMPLE_DIR)/waveguide_crossing.py                              \
    $(EXAMPLE_DIR)/wvg-src.py                                         \
    $(EXAMPLE_DIR)/zone_plate.py

if WITH_COVERAGE
  PY_LOG_COMPILER = coverage run -p --source=$(top_srcdir)/python/meep --omit=$(top_srcdir)/python/tests/*,$(top_srcdir)/python/examples/*,$(top_srcdir)/python/meep/materials.py
else
  PY_LOG_COMPILER = $(RUNCODE) $(PYTHON)
endif

TEST_EXTENSIONS = .py
TESTS_ENVIRONMENT = export PYTHONPATH="`$(PYPATH_CONV) $(abs_top_builddir)/python`$(PYPATH_SEP)$$PYTHONPATH";

######################################################################
# instructions for building the _meep binary library that implements
# the low-level interface between python codes and the libmeep C++ library
######################################################################
BUILT_SOURCES = meep-python.cxx __init__.py meep.py
# files are listed individually (rather than whole directories) so that
# untracked files in tests/ and examples/ are not included in the tarball.
# The MPI- and MPB-dependent tests are listed separately because $(TESTS)
# depends on the configure options.  coupler.gds is read by examples/coupler.py.
EXTRA_DIST = $(BUILT_SOURCES) meep.i meep-python.hpp typemap_utils.cpp materials.py \
             $(TESTS) $(ADJOINT_TESTS) $(MATERIAL_DISPERSION_TEST) $(MPB_TESTS)     \
             $(TEST_DIR)/utils.py $(TEST_DATA_FILES) $(EXAMPLES) $(EXAMPLE_DIR)/coupler.gds
CLEANFILES = .coverage
MAINTAINERCLEANFILES = $(BUILT_SOURCES)

AM_CPPFLAGS = -I$(top_srcdir)/src               \
              -I$(top_srcdir)/libpympb		\
              -I$(top_builddir) # for config.h

LIBMEEP = $(top_builddir)/src/libmeep.la

SWIG_WARN_FLAGS = -Wno-unused-function -Wno-unused-parameter -Wno-missing-field-initializers

_meep_la_SOURCES = meep-python.cxx
_meep_la_LIBADD = $(LIBMEEP) $(PYTHON_LIBS)
_meep_la_LDFLAGS = -module $(PYMODULE_LDFLAGS) -avoid-version
_meep_la_CPPFLAGS = $(PYTHON_INCLUDES) $(AM_CPPFLAGS) $(SWIG_WARN_FLAGS)
_meep_la_CXXFLAGS = -O1

pkgpyexec_LTLIBRARIES = _meep.la

######################################################################
# instructions for building the _mpb binary library that implements
# the low-level interface between python codes and the libmpb C++ library
######################################################################
if WITH_MPB
  BUILT_SOURCES += mpb-python.cxx mpb.py
  EXTRA_DIST += solver.py mpb_data.py
  MPB_LA = _mpb.la
  pympbdir = $(pkgpythondir)/mpb
  pympb_PYTHON = mpb.py solver.py mpb_data.py
  pympb_LTLIBRARIES = _mpb.la
  _mpb_la_SOURCES = mpb-python.cxx
  _mpb_la_LIBADD = $(PYTHON_LIBS) $(top_builddir)/libpympb/libpympb.la $(LIBMEEP)
  _mpb_la_LDFLAGS = -module $(PYMODULE_LDFLAGS) -avoid-version
  _mpb_la_CPPFLAGS = $(PYTHON_INCLUDES) $(AM_CPPFLAGS) $(SWIG_WARN_FLAGS)
  _mpb_la_CXXFLAGS = -O1

install-data-hook:
	mv $(DESTDIR)$(pkgpythondir)/mpb/mpb.py $(DESTDIR)$(pkgpythondir)/mpb/__init__.py

uninstall-hook:
	rm -f $(DESTDIR)$(pkgpythondir)/mpb/__init__.py

endif

##################################################
# instructions for running SWIG to generate the
# source files for the previous step
##################################################
if MAINTAINER_MODE

SWIG_VERSION = $(shell $(SWIG) -version | grep Version | awk '{print $$3}')
MEEP_SWIG_SRC = meep.i numpy.i vec.i

HPPFILES=                                    \
 $(top_srcdir)/src/meep.hpp                  \
 $(top_srcdir)/src/meep_internals.hpp        \
 $(top_srcdir)/src/bicgstab.hpp              \
 $(top_srcdir)/src/meep/vec.hpp              \
 $(top_srcdir)/src/meep/mympi.hpp            \
 $(top_srcdir)/src/meepgeom.hpp              \
 $(top_srcdir)/src/material_data.hpp         \
 $(top_srcdir)/src/adjust_verbosity.hpp      \
 meep-python.hpp

if MEEP_SWIG_PYTHON_THREADS
SWIG_MEEP_FLAGS = -threads
endif # MEEP_SWIG_PYTHON_THREADS

meep-python.cxx: $(MEEP_SWIG_SRC) $(HPPFILES)
	$(SWIG) -Wextra $(SWIG_MEEP_FLAGS) $(AM_CPPFLAGS) -outdir $(builddir) -c++ -nofastunpack -python -o $@ $(srcdir)/meep.i

if WITH_MPB
MPB_SWIG_SRC = mpb.i
mpb-python.cxx: $(MPB_SWIG_SRC) $(top_srcdir)/libpympb/pympb.hpp $(top_srcdir)/src/meepgeom.hpp $(top_srcdir)/src/material_data.hpp $(top_srcdir)/src/material_data.cpp
	$(SWIG) -Wextra $(AM_CPPFLAGS) $(PYMPBINCLUDE) -outdir $(builddir) -c++ -nofastunpack -python -o $@ $(srcdir)/mpb.i
mpb.py: mpb-python.cxx
MPB_PY = mpb.py
endif # WITH_MPB

meep.py: meep-python.cxx

__init__.py: meep.py $(MPB_PY)
	cp $< $@
	echo "__version__ = '@VERSION@'" >> $@
	if [[ "${SWIG_VERSION}" = 3.0.12 ]]; then \
		sed -i.bak '/^if _swig_python_version_info >= (2, 7, 0):/,/^else:/d' $@; \
		sed -i.bak 's/    import _meep/from . import _meep/' $@; \
	fi
if WITH_MPB
	if [[ "${SWIG_VERSION}" == 3.0.12 ]]; then \
		sed -i.bak '/^if _swig_python_version_info >= (2, 7, 0):/,/^else:/d' mpb.py; \
		sed -i.bak 's/    import _mpb/from . import _mpb/' mpb.py; \
	fi
endif

INIT_PY = __init__.py
MPB_PY_PATH = $(builddir)/mpb.py

else # not maintainer mode

INIT_PY = $(srcdir)/__init__.py
MPB_PY_PATH = $(srcdir)/mpb.py

endif # MAINTAINER_MODE

######################################################################
# specification of python source files to be byte-compiled at installation
######################################################################
HL_IFACE =                  \
    $(srcdir)/binary_partition_utils.py       \
    $(srcdir)/chunk_balancer.py       \
    $(srcdir)/geom.py       \
    $(srcdir)/simulation.py \
    $(srcdir)/source.py     \
    $(srcdir)/timing_measurements.py       \
    $(srcdir)/visualization.py  \
    $(srcdir)/materials.py  \
    $(srcdir)/verbosity_mgr.py

pkgpython_PYTHON = __init__.py $(HL_IFACE)

adjointdir = $(pkgpythondir)/adjoint
adjoint_PYTHON  = $(srcdir)/adjoint/__init__.py              \
                  $(srcdir)/adjoint/basis.py                 \
                  $(srcdir)/adjoint/objective.py             \
                  $(srcdir)/adjoint/optimization_problem.py  \
                  $(srcdir)/adjoint/filters.py               \
                  $(srcdir)/adjoint/filter_source.py         \
                  $(srcdir)/adjoint/connectivity.py         \
                  $(srcdir)/adjoint/unfilter_design.py        \
                  $(srcdir)/adjoint/wrapper.py               \
                  $(srcdir)/adjoint/utils.py

PY_PKG_FILES = $(INIT_PY) $(HL_IFACE)

meep: _meep.la $(MPB_LA) __init__.py $(HL_IFACE) $(adjoint_PYTHON)
	mkdir -p meep/adjoint
	cp $(PY_PKG_FILES) meep
	cp .libs/_meep$(PYMODULE_BUILT_EXT) meep/_meep$(PYMODULE_EXT)
	cp $(adjoint_PYTHON) meep/adjoint
if WITH_MPB
	mkdir -p meep/mpb
	cp .libs/_mpb$(PYMODULE_BUILT_EXT) meep/mpb/_mpb$(PYMODULE_EXT)
	cp $(MPB_PY_PATH) meep/mpb/__init__.py
	cp $(srcdir)/solver.py $(srcdir)/mpb_data.py meep/mpb
endif

all-local: meep

check-adjoint:
	$(MAKE) check TESTS="$(ADJOINT_TESTS)"

clean-local:
	rm -rf meep __init__.py.bak

distclean-local:
	rm -f *.h5
	rm -f *.mp4
	rm -f *.gif
