# Meep Release Notes

## Meep 1.26.0

3/9/2023

* Improved Gaussian beam source in 2D ([#2333]).

* Support for returning the number of timesteps elapsed in simulation ([#2337]).

* Bug fix for fields update in cylindrical coordinates ([#2382]).

* Bug fix for PMLs in cylindrical coordinates ([#2383]).

* Bug fix in amplitude function of eigenmode source ([#2394]).

* The `doc/` directory (the manual) is no longer included in the release tarball to save space, since most people will view this online.  (It is still in the git repository.)

* Various improvements and minor bug fixes ([#2321], [#2349], [#2371], [#2380], [#2390], [#2413]), and additional unit tests and documentation ([#2314], [#2360], [#2364], [#2365], [#2387], [#2395], [#2402]).

## Meep 1.25.0

11/17/2022

* Support for connectivity constraints in adjoint solver ([#2207]).

* Support for animation in topology optimization ([#2186]).

* Support for `load_minus_flux` in adjoint solver ([#2271]).

* Support colorbars in `plot2D` ([#2289]).

* Support for `plot3D` ([#2305]).

* Various improvements and bug fixes ([#2176], [#2179], [#2190], [#2194], [#2202], [#2203], [#2208], [#2251], [#2253], [#2264], [#2290]), and additional unit tests and documentation.

## Meep 1.24.0

7/21/2022

* Support for adjoint gradients of local density of states (LDOS) ([#2077]).

* Improvements to memory usage of adjoint solver ([#1855]).

* Various bugfixes ([#1959], [#2044], [#2066], [#2073], [#2079], [#2091], [#2095], [#2114]) and additional unit tests ([#2032], [#2049], [#2053], [#2076], [#2082]).

## Meep 1.23.0

4/6/2022

* Support for termination condition function based on the field energy in the entire cell ([#2021]).

* Support for mode decomposition for 2d cell with out of plane wavevector ([#1968]).

* Type annotations for the `Simulation` class ([#1919]).

* Various improvements ([#1821], [#1895], [#2027]), bug fixes ([#1955], [#2016]), and additional documentation ([#2005]).

## Meep 1.22.0

1/11/2022

* Support for anisotropic materials and subpixel smoothing in adjoint gradients ([#1780], [#1801]).

* Support for damping property in `MaterialGrid` ([#1804]).

* Support for cylindrical coordinates in `plot2D` ([#1873]).

* Various bugfixes ([#1796], [#1830], [#1849], [#1860]), performance enhancements ([#1826], [#1839], [#1871], [#1872], [#1877]), and additional unit tests and documentation.

## Meep 1.21.0

10/12/2021

* Support for checkpointing structure and fields ([#1715], [#1738]).

* Support for multithreading via OpenMP ([#1628]).

* Support for cache-oblivious loop tiling of the E-from-D and H-from-B field updates ([#1655]).

* Support for load balancing parallel simulations using timing measurements ([#1775]).

* Support for automatic DFT decimation ([#1732]).

* Support for adjoint gradients in cylindrical coordinates ([#1749]).

* Revamped convergence criteria for DFT fields ([#1740]).

* Various bugfixes ([#1593], [#1769]), performance enhancements ([#1730], [#1763]), minor improvements, and additional documentation.

## Meep 1.20.0

8/11/2021

* Support for decimation of the DFT time-series updates ([#1684], [#1720], [#1722]).

* Support for optional single-precision floating point for the DFT fields arrays ([#1675]).

* Support for cache-oblivious loop tiling of the step-curl field updates ([#1655]).

* Performance improvements in chunk-to-chunk communication ([#1656], [#1721]).

* Code coverage for Python API via GitHub Actions ([#1651]).

* Various bugfixes ([#1692], [#1704]), minor improvements, and additional documentation.

## Meep 1.19.0

7/6/2021

* Support for optional single-precision floating point for fields arrays ([#1544]).

* Optional subpixel smoothing for `MaterialGrid` using analytical rather than numerical-quadrature framework ([#1539], [#1568]).

* Improvements to user-specified chunk layout as `BinaryPartition` object ([#1577]).

* JAX wrapper for `Simulation` object and adjoint solver ([#1569]).

* Continuous integration (CI) via GitHub Actions ([#1599], [#1608], [#1618], [#1623]).

* Remove MPI synchronization barrier at each timestep from connection of chunks ([#1635]).

* Various bugfixes ([#1546], [#1574], [#1588], [#1591], [#1592], [#1598], [#1599], [#1601], [#1603], [#1606], [#1634], [#1652]), additional documentation, and tests.

## Meep 1.18.0

3/26/2021

* New `get_epsilon_grid` function for evaluating ε on user-specified grid with arbitrary resolution ([#1522]).

* Support for user-specified chunk layouts for manual control over load-balancing ([#1528]).

* `MaterialGrid` `design_parameters` is renamed to `weights` and `U_SUM` is renamed to `U_MEAN` ([#1512]).

* Performance improvement in chunk division ([#1499]).

* Various bugfixes ([#1487], [#1515], [#1519], [#1521], [#1527]), additional documentation, and tests.

## Meep 1.17.1

1/8/2021

* Fix accidental breakage of the adjoint routines ([#1464]).

## Meep 1.17

1/4/2021

* `get_array_slice` now does continuous interpolation as the slice position is moved ([#1456]).

* New `contour` option for contour-plotting in `plot2D` ([#1437]).

* Adjoint optimization of `near2far` transformations ([#1417]).

* `get_array_metadata` is now consistent between array slices and DFT slices ([#1456]), no longer leaks memory ([#1447]), and returns numpy arrays rather than tuples ([#1458]).

* Bugfixes in adjoint-optimization filters ([#1427]).

## Meep 1.16.1

10/20/2020

* Bugfix in adjoint code ([#1403]).

## Meep 1.16.0

10/6/2020

* Gaussian beam source feature ([#1303] and [#1310]).

* New API for specifying planewave diffraction orders for eigenmode sources
  and coefficients ([#1316]).

* More accurate gradients in adjoint code ([#1285]).

* Simpler Python API for outputting ε or μ at a given frequency ([#1374]).

* `--with-libctl-dir` option of `configure` now accepts simply the installation `prefix` in addition to `prefix/share/libctl` ([#1286]).

* Less verbose mode-solver output from MPB ([#1302], [#1388]), and new
  `meep.verbosity` option in Python ([#1349]).

* Bug fix for single-point DFT monitor ([#1333]).

## Meep 1.15.0

7/8/2020

* Minimum-lengthscale filters for adjoint optimization ([#1205]).

* Python API documentation in docstrings ([#1240]).

* `MaterialGrid` material type in Python to interpolate an array of material
  values as the "material" of an object, especially for topology optimization ([#1242]).

* `merge_subgroup_data` Python function for coordinating parallel
  computations ([#1192]).

* Eigenmode sources now ensure that the source has the same
  frequency as the mode ([#1218]).

* Performance improvements to eigenmode sources ([#1233], [#1244], [#1257]).

## Meep 1.14.0

4/17/2020

* New adjoint solver for density-based topology optimization, including
  filtering, automatic differentiation, and other frequencies ([#1167]).

* DFT functions now allow you to pass an arbitrary array of frequencies, instead
  of being limited to equally spaced frequencies ([#1154] and [#1156]).

* Experimental shift-and-invert frequency-domain eigensolver ([#1158]).

* Renamed `omega` parameter to `frequency` at some places in the Python API,
  for consistency ([#1171]), and `dft_fields` object now takes `fcen` and `df` instead
  of `freq_min` and `freq_max` in Python.

* Support for SWIG 4.0 ([#1159]), and various other minor fixes.

## Meep 1.13.1

2/25/2020

* Avoid writing to source directory in remaining tests ([#1132]).

## Meep 1.13.0

2/19/2020

* Optional parameter `omega` for `output-epsilon` and similar functions,
  allowing the complex ε and μ at a given frequency to be outputted ([#1112], following [#919]).

* `near2far` computation now supports cylindrical coordinates ([#1090]).

* Experimental support for slanted prisms (requires libctl 4.5)
  via `sidewall_angle` parameter to prism objects ([#1129]).

* New `yee_grid=False` optional argument to `add_dft_fields`; by passing `True`
  one can compute the fields on the original Yee grid ([#1095]).

* New function `meep::make_output_directory()` to make a temporary
  directory (in `TMPDIR` or similar) and `meep::delete_directory(path)`
  to perform recursive deletion (like `rm -rf`).  These are now
  used in tests to avoid writing to the source directory ([#1121], [#1122] and [#1126]).

* Jupyter notebooks now show a graphical progress bar during simulations ([#1078]).

* `kz-2d` option in Scheme, mirroring Python `kz_2d` ([#1062]).

* Various bugfixes, documentation supplements, and other minor improvements.

## Meep 1.12.0

11/12/19

  * Faster 2d simulations with nonzero `kz` via the `kz_2d` option ([#1047]).

  * New Meep `verbosity` option superseding `quiet` and `verbose` flags ([#994]).

  * Output now only shows ≤ 10 geometric objects by default ([#1002]).

  * Performance improvements for `split_chunks_evenly=False`.

  * Fixed memory leaks ([#1041], [#1042]).

## Meep 1.11.0

7/29/19

  * Experimental support for gyrotropic media including magneto-optical effects ([#863]).

  * Mode decomposition for oblique waveguides ([#940], [#945]) and dispersive materials ([#919]).

  * Accept tuples in place of Vector3 arguments ([#960]).

  * Capture C++ error messages in Python notebooks ([#953]).

  * Automatically abort simulation if the fields blow up ([#922]).

  * Print additional timing statistics ([#927], [#952]).

  * Various small bugfixes and documentation improvements.

## Meep 1.10.0

6/5/19

  * New Python functions for simple visualization of the simulation domain ([#872]).

  * Capture Meep and MPB output in Python notebooks ([#891], [#894])

  * Add optional `meep.quiet()` parameter to the Python interface ([#876]).

  * Python evaluation of materials ε(ω) and μ(ω) ([#862]).

  * Experimental multithreading support for near2far calculation ([#868]) and other speedups ([#869]).

  * Add `stop_after_walltime` and `stop_on_interrupt` in Python ([#860]).

  * GDSII file introspection ([#817]).

  * Various small bugfixes and documentation improvements.


## Meep 1.9.0

4/17/19

  * Adjoint solver to compute sensitivity of solution to material perturbations ([#795]).

  * Experimental `do_averaging` feature for user-defined material functions ([#771], [#791]).

  * Periodic boundaries support in `near2far` via `nperiods` option ([#769], [#789]).

  * Capture more output in Python notebooks ([#785], [#807]).

  * `dft-energy` feature ([#744], [#747]).

  * Eigenmode sources are normalized to unit power ([#728]).

  * Fix interpolation of DFT slice output ([#787]).

  * Bug fix in `run-k-points` ([#779]).

  * Eigenmode sources for negative angles ([#752]).

  * Various other minor bugfixes, build fixes, documentation improvements, tutorials, etcetera.

## Meep 1.8.0

2/13/19

  * libctl 4.2 is required

  * Add `--without-scheme` flag to `./configure` ([#705])

  * Improve error messages in Python interface ([#699])

  * Allow `kguess` to specify MPB lattice vector for launching oblique waveguide modes ([#675])

  * Allow user materials when checking for conductivity ([#689])

  * Add `split_chunks_evenly` flag to `Simulation` constructor. Setting to `False` will improve parallel simulation performance by dividing chunks based on work instead of size ([#681])

  * Added `Simulation.visualize_chunks()` to visualize the chunk layout ([#671])

  * Improved stability of lorentzian susceptibility ([#666])

  * Get array metadata for `get_array` and `get_dft_array` ([#655])

  * Add ability to get a source slice as a numpy array ([#652])

  * Fixed performance issues in ModeSolver.find_k ([#644])

  * Add `force_all_components` flag to `Simulation` constructor ([#631])

  * libmeepgeom was merged into libmeep ([#630])

  * Expose `run_k_point` to access more Harminv data ([#626])

  * Various other bug fixes, documentation improvements, etc.

## Meep 1.7

11/16/18

 * Add `transform` method to `meep.Medium` ([#603]).

 * Read epsilon input from a numpy array when passed to a `Simulation` as `default_material`([#593]).

 * Support `geometry_center` in Python ([#599]).

 * Add Python `Ldos` class ([#581]).

 * Compute Fourier-transformed fields (e.g. fluxes) in `solve_cw` ([#570]).

 * Enable builds without MPB ([#558]).

 * Add birefringent materials to materials library ([#559]).

 * Print dominant planewave in `get_eigenmode` ([#531]).

 * Python API for GDSII regions ([#518])

 * Multilevel atom susceptibilities for Python and Scheme ([#500]).

 * Fix bug in `get_eigenmode_coefficients` for 2d cell with non-zero kz ([#602]).

 * Fix sync of eigenmode calculation when no mode is found ([#596]).

 * Fix memory leak in `get_dft_array` ([#577]).

 * Use same MPB phase on all processes, fixing bug with eigenmodes and multiprocessing ([#578]).

 * Fix memory leaks in `get_eigenmode` ([#558]).

 * Various other bug fixes, documentation improvements, etc.

## Meep 1.6

9/7/2018

 * Python interface to import GDSII files ([#392]).

 * New binary grating tutorial ([#376]).

 * Get source amplitude from HDF5 file ([#388]).

 * get_eigenmode_coefficients now returns group velocity and kpoints ([#396]).

 * New tutorial for visualizing 3d structures ([#416]).

 * Mode decomposition feature supports symmetries ([#417]).

 * Support for Guile >= 2.0.12 ([#419]). Merged upstream to SWIG repo ([#1288]).

 * Python get_eigenmode function and EigenmodeData class ([#422])

 * Symmetry support in dft arrays ([#427]).

 * Python 3.7 support ([#456]).

 * get-eigenmode-coefficients added to Scheme API ([#477]).

 * materials_library.py now part of Python package (e.g., from meep.materials import Al) ([#479]).

 * materials-library.scm automatically available in Meep scripts ([#483]).

 * Structure dump/load feature now supports dispersive materials ([#454]).

 * Various bug fixes, documentation improvements, etc.

## Meep 1.5

6/7/2018

 * Python interface to MPB ([#191] etc.).

 * Mode decomposition: given a DFT flux plane, decompose the fields
   at each frequency into the power in each mode of a waveguide
   or similar ([#192], [#248], etc.).

 * DFT slices: output Fourier-transformed fields in any given
   region of space ([#259]).

 * New `prism` geometric-object type for polygonal prisms ([#341], [#345])
   for upcoming GDSII import ([#357]).  Libctl 4.1.0 is required.

 * Structure dump/load feature to rapidly load in a geometry
   from a previous calculation ([#261], [#266]).

 * Susceptibilities are now supported in user-defined materials
   in Python ([#203], [#305]).

 * 64-bit support for extremely large computations ([#193]).

 * Various bug fixes, documentation improvements, etc.

## Meep 1.4.3

2/1/2018

 * Allow `meep` Python module to be imported without setting `PYTHONPATH` ([#189]).

## Meep 1.4.2

1/26/2018

  * Build fix for Python due to missing file ([#184]).

## Meep 1.4.1

1/19/2018

  * Minor packaging fixes.

## Meep 1.4

1/18/2018

  * Full-featured Python interface.

  * Migrated documentation to github/markdown/readthedocs ([#55]).

  * New feature to get slice as array in C++ and Python APIs ([#96], [#105]).

  * `libmeepgeom` library to allow C++ users to access geometric-object
    API ([#56]).

  * Removed overly conservative stability check for Lorentzian
    susceptibilities ([#150]).

  * Corrected small error in frequency interval for `dft-ldos` ([#40]).

  * Bug fixes in near-to-farfield spectra ([#21]), eigenmode source ([#20]),
    and LDOS ([#40]).

## Meep 1.3

31 March 2015.

  * New near-to-far-field functionality: given a bounding surface,
    automatically computes the Fourier-transformed field in any
    desired grid of "far-field" points arbitrarily far away.

  * Compatibility with Harminv 1.4 (fixes issue [#13]: ppc64 portability).

  * Fix compilation with latest C++ standard (e.g. on OS X 10.9).

  * Bug fix in CW solver convergence test; thanks to
    Wu Chuanren and @FilipDominec for the bug report.

  * Build fix for Fedora 21 (thanks to Dean Brettle) (issue [#14]).

## Meep 1.2.1

2 April 2014.

  * Added new absorber type, as an alternative to PML, which simply
    provides a scalar conductivity gradient for cases where PML fails.

  * Fixed bug which sometimes prevented dispersive materials from being
    used in PML regions.

  * Some fixes to BLAS/LAPACK linking.

  * Bug fixes in LDOS computation.

  * Work around gcc bug #54498, which caused a spurious PML test
    failure with gcc 4.7 and 4.7.1; thanks to Brahmanand Jogai and
    Thorsten Alteholz for the bug reports.

## Meep 1.2

20 July 2012.

  * Fixed to work with Guile version 2.x (older versions still work);
    requires libctl 3.2 or later.

  * Added `epsilon-input-file` feature to read a scalar dielectric function
    from an HDF5 file (similar to MPB).

  * Support for anisotropic dispersive materials (tensor sigma parameter).

  * Support for Drude dispersion model.  New syntax is
    `make drude-susceptibility`, `make lorentzian-susceptibility`, etc.
    (old `make polarizability` is still supported for backwards compatibility).

  * Support for "thermal" dispersive materials which include noise
    term in the polarization.

  * Added `dft-ldos` feature for efficient LDOS-spectrum computation.

  * Documented stress-tensor (force) spectrum computation feature.

  * Added `mean-stretch` property of PML (defaults to 1), to support
    real coordinate stretching for damping evanescent modes.

  * Support for eigenmode-source feature using upcoming MPB release.

  * Various small bugfixes.

## Meep 1.1.2

31 August 2009.

  * Added `make check` test (in 2D_convergence) for new `special-kz?`
    feature (for computing out-of-plane modes in 2d more efficiently).

  * Fix typo preventing Casimir calculations from running for periodic
    problems.

## Meep 1.1.1

24 August 2009.

  * Fixed release bug preventing Casimir calculation from running.

## Meep 1.1

20 August 2009.

  * Meep's PML is now a true PML for arbitrary anisotropic, dispersive,
    and conducting media.  (Now uses a slightly unconventional
    reformulation of PML described at ab-initio.mit.edu/meep/pml-meep.pdf)

  * Fixed bug which caused anisotropic non-diagonal mu to be unstable.

  * Fix compilation failure with gcc 4.4 due to missing cstdio
    header (thanks to Linran Fan and Bin Shao for the bug reports).

  * C++ interface: volume was renamed to grid_volume and geometric_volume
    was renamed to volume, to better reflect their respective roles.

  * Added `accurate-fields-near-cylorigin?` option to have more
    accurate fields near the r=0 origin for large m in cylindrical
    coordinates, at the expense of requiring a smaller Courant factor.
    (Default is `false`, corresponding to behavior in older Meep versions.)

  * In 2d computational cells, added much more efficient support for
    exp(ikz) z-dependence, enabled by new `special-kz?` input variable
    (default is `false` since it only works in 2d and is a little subtle
     for real fields).

  * Includes preliminary new features to aid in computation of
    optical forces (both classical and quantum Casimir forces);
    further documentation pending more testing.

  * Removed obsolete `doc` directory (all documentation is on the website
    these days).

  * Small performance improvements in Lorentzian dispersion handling.

  * Fix configure script failure when cross-compiling.

  * Fix compilation failure with MPICH.

## Meep 1.0.3

5 June 2009.

  * Allow `GUILE_CONFIG` environment variable to override location of
    `guile-config` program in `configure` script; this is useful when
    cross-compiling.

## Meep 1.0.2

2 June 2009.

  * Correct superficial `make check` failure on 32-bit x86 machines
    with gcc 4.3.x, due to slight impact on floating-point rounding
    by automatic SSE/SSE2 vectorization; thanks to Silviu Popescu for
    the bug report.

  * Correct superficial `make check` failure when compiling under icc.

Meep 1.0.1

28 May 2009.

  * Enable correct operation and passed test suite when `MEEP_SINGLE`
    (single-precision) mode is enabled in meep.hpp; thanks to
    Seyoon Kim for the bug reports.

  * Use new automake features to have less-verbose build output by
    default (you can build in verbose mode by `make V=1`), and
    running all test programs then reporting which ones failed
    instead of stopping at the first failure.

  * Fix superficial failure in 2D_convergence test under gcc 3.4.6;
    thanks to Alex Prengel for the bug report.

  * Fix failure in flux test under gcc 4.3.1 in some cases; thanks
    to Alex Prengel for the bug report.

  * Fix compilation problem with gcc 4.4, correcting Debian bug #505002.

## Meep 1.0

28 April 2009.

  * New timestepping scheme for off-diagonal anisotropic epsilon and
    mu, based on technique by Werner and Cary [ J. Comp. Phys. 226,
    1085 (2007) ], that improves FDTD stability when anisotropy is
    present (such as when subpixel averaging is used on isotropic media).

  * Scheme user interface now supports user-specified anisotropic
    (real-symmetric) epsilon and mu (via epsilon-diag, epsilon-offdiag,
    mu-diag, and mu-offdiag parameters, similar to MPB).  Accurate
    subpixel averaging of anisotropic media based on the method by
    Kottke, Farjadpour, & Johnson [ Phys. Rev. E. 77, 036611 (2008) ].

  * Anisotropic dispersive materials are now supported, although
    currently the dispersive part of the epsilon/mu tensor must be
    diagonal, via the new sigma-diag parameter of polarizability.
    (The corresponding C++ interface has also removed delta_epsilon.)

  * The delta-epsilon parameter of polarizability has been removed;
    you should use sigma instead.

  * New `fields::integrate2` function (and corresponding Scheme function
    `integrate2-field-function`) to perform integrations involving two
    simulations with the same computational cell (e.g. field-overlap
    calculations for coupled-mode theory).

  * In the Scheme interface, subpixel averaging is not used for
    user-specified material-function types; you only get subpixel
    averaging for the standard shapes (blocks, cylinders, etcetera).

  * Haskell code-generation is no longer used, and hsrc directory is
    removed.  Bitrotted and undocumented (hence unused)
    saturable-absorber feature has been removed, along with
    energy-saturation parameter of polarizability.

  * Some bug-fixes to test programs that made them overly sensitive
    to roundoff errors and possibly fail depending on the compiler.
    (New `fields::round_time` and meep-round-time functions to round
     times to single-precision, useful for robust time comparisons.)

## Meep 0.20.4

17 March 2009.

  * Bug fix in cylindrical code, which caused it to blow up
    in some circumstances for nonzero m.

  * Bug fix: non-integrated sources with conductivity are now
    second-order accurate, thanks to Alejandro Rodriguez.

  * Bug fix in writing strings with parallel HDF5, thanks to Zheng Li
    for the bug report.

  * Check that PML parameters are sensible (e.g. that total PML thickness
    is no greater than cell thickness) to avoid common mistakes.

  * New extra-materials input variable, so that you no longer
    have to use "dummy objects" to specify the existence of
    some materials when using material-function types.

## Meep 0.20.3

24 July 2008.

  * Fixed circular dependency in Makefile, which caused problems with some
    versions of make; thanks to Kaoru Narita for the bug report.

## Meep 0.20.2

21 July 2008.

  * Fixed incompatibility with Guile 1.6.x or earlier; thanks to the bug
    report by Andreas Unger.

## Meep 0.20.1

20 July 2008.

  * Improved handling of nested synchronized-magnetic calls.

  * Bug fix: parallel builds (`make -j`) should now work.

  * Bug fix: pkg-config file was incorrectly installed for MPI version;
    thanks to Majid Sodagar for the bug report.

## Meep 0.20

19 July 2008.

  * Support for user-specified permeability (mu).  Renamed `dielectric`
    to `medium` in libctl interface, new `mu` property and new
    output-bfield and output-mu functions, and new `Permeability` and `Bx`
    etc. field types.

  * Support for user-specified electric and/or magnetic conductivities.
    These are especially useful to add a desired dissipation loss
    (an imaginary part of epsilon/mu) in a narrow bandwidth, without
    messing around with Lorentzian dispersive materials.

  * Add predefined perfect-magnetic-conductor (mu = -infinity) material,
    along with perfect-electric-conductor (eps = -infinity).

  * Added synchronized-magnetic step function to allow step functions
    to run with the electric and magnetic fields synchronized in time
    to second-order accuracy.

  * New PML implementation (UPML instead of split-field), should
    have lower reflection in many cases.

  * User-specified PML profile and asymptotic reflection.

  * Internally, all timestepping code is now handwritten (and much shorter)
    rather than old verbose Haskell-generated code; this should make
    it easier to add new features.

  * Add support for non-integrated current sources, if the is-integrated?
    property of the current is set to false; this is now the default,
    to make handling of E and H sources more similar and intuitive.

  * Work with HDF5 1.8 (which previously would not compile unless you
    manually set a preprocessor flag, due to API changes).

  * Check for ctl.h in /usr/include/ctl/ctl.h (default in Fedora),
    and check for libctl in /usr/share/libctl3 (default in Debian & Ubuntu).

  * Bug fix: fixed relative phase of E and H sources (which were off
    from one another by half a timestep); thanks to M. Megens for bug report.

  * Bug fix: make sure h5 filenames have unique timestep for cases where
    dt is very small or very large.

## Meep 0.10.1

13 Nov. 2007.

  * Bug fix in flux_in_box, which accidentally returned the flux
    multiplied by the number of processors, instead of the flux.

  * Bug fix in epsilon averaging for structures including metals
    (`epsilon < 0`), fixing an instability.

  * Bug fix in output-png when running in parallel (removing race condition).

  * Fixed bug that disabled subpixel averaging for dimensions=1 (thanks
    to Mischa Megens for the bug report).

  * Fixed bug that caused output-tot-pwr to stop Meep with an error message;
    thanks to Vyacheslav Sokolov for the bug report.

  * Make `at-every` step functions less susceptible to rounding errors;
    thanks to L. Le Guyader for the bug report.

  * Fixed bug in dispersive media that wasted memory on parallel machines
    (the polarization memory was not parallelized); thanks to J. L. Silva
    for the bug report.

  * Bug fix in output-png+h5, thanks to a report by Chad Husko.

  * Fixed several deadlocks that could occur when the parallel Meep
    is used with a serial HDF5 library (we continue to recommend
    using the parallel HDF5 library with parallel Meep, however).
    Thanks in part to Lingling Tang for his bug report.

  * For maintainer-mode, improved detection of Haskell package names;
    thanks to Liang Huo for the bug report.

## Meep 0.10

21 Aug. 2006.

  * `eps-averaging?` is now turned on by default (in libctl interface),
    using much-improved algorithm by Ardavan Farjadpour.  This greatly
    improves accuracy, and also allows continuous tuning of geometric
    parameters.  (See our upcoming paper in Optics Lett., with a preprint
    linked on the web site.)  New input variables subpixel-tol and
    subpixel-maxeval to control the accuracy of the subpixel averaging.

  * Support for chi2 (Pockels) as well as chi3 (Kerr) nonlinearities.

  * Symmetries no longer require the cell size to be an even number of
    pixels.  Previously, Meep exited with an error in this case, whereas
    now it simply adds an extra pixel to the cell size as needed.

  * New with-prefix step function to allow you to use a different
    filename-prefix for selected outputs.

  * New feature for output-png: built-in shell variable $EPS that refers
    to the last-output epsilon .h5 file, which you can use to easily
    add dielectric contours/overlays to the field output image.

  * Added output-png+h5 function that outputs both .png and .h5 files.

  * New functions flux-in-box, electric-energy-in-box, magnetic-energy-in-box,
    and field-energy-in-box (convenience wrappers around C++ functions).

  * Bug fix in Kerr nonlinearity - chi3 was accidentally scaled by epsilon^4
    factor.

  * Bug fix: if you specified three or more symmetries, at most two
    symmetries were used (ignoring the rest).

  * Bug fix in rotate2 symmetry, which wasn't working correctly.

  * Bug fix in add-flux for multiple flux regions, thanks to K. Choi.

  * Bug fix in harminv where it wouldn't allow you to call harminv more
    than once for the same run loop; thanks to Aristos Karalis.

  * Bug fix in save-flux/load-flux that prevented it from working properly
    without output directories, thanks to Karl Koch.

  * Fixed abort that sometimes occurred due to rounding when the source
    was the same width as the cell (thanks to G. J. Parker).

  * Fixed minor build problems on Cygwin, SGI, and other systems,
    thanks to Christopher Kang, Robyn Landers, Florencio Garcia, and others.

## Meep 0.9

1 Apr. 2006.

  * Initial public release.

<!--- generated links: -->
[#13]: https://github.com/NanoComp/meep/issues/13
[#14]: https://github.com/NanoComp/meep/issues/14
[#20]: https://github.com/NanoComp/meep/issues/20
[#21]: https://github.com/NanoComp/meep/issues/21
[#40]: https://github.com/NanoComp/meep/issues/40
[#55]: https://github.com/NanoComp/meep/issues/55
[#56]: https://github.com/NanoComp/meep/issues/56
[#96]: https://github.com/NanoComp/meep/issues/96
[#105]: https://github.com/NanoComp/meep/issues/105
[#150]: https://github.com/NanoComp/meep/issues/150
[#184]: https://github.com/NanoComp/meep/issues/184
[#189]: https://github.com/NanoComp/meep/issues/189
[#191]: https://github.com/NanoComp/meep/issues/191
[#192]: https://github.com/NanoComp/meep/issues/192
[#193]: https://github.com/NanoComp/meep/issues/193
[#203]: https://github.com/NanoComp/meep/issues/203
[#248]: https://github.com/NanoComp/meep/issues/248
[#259]: https://github.com/NanoComp/meep/issues/259
[#261]: https://github.com/NanoComp/meep/issues/261
[#266]: https://github.com/NanoComp/meep/issues/266
[#305]: https://github.com/NanoComp/meep/issues/305
[#341]: https://github.com/NanoComp/meep/issues/341
[#345]: https://github.com/NanoComp/meep/issues/345
[#357]: https://github.com/NanoComp/meep/issues/357
[#376]: https://github.com/NanoComp/meep/issues/376
[#388]: https://github.com/NanoComp/meep/issues/388
[#392]: https://github.com/NanoComp/meep/issues/392
[#396]: https://github.com/NanoComp/meep/issues/396
[#416]: https://github.com/NanoComp/meep/issues/416
[#417]: https://github.com/NanoComp/meep/issues/417
[#419]: https://github.com/NanoComp/meep/issues/419
[#422]: https://github.com/NanoComp/meep/issues/422
[#427]: https://github.com/NanoComp/meep/issues/427
[#454]: https://github.com/NanoComp/meep/issues/454
[#456]: https://github.com/NanoComp/meep/issues/456
[#477]: https://github.com/NanoComp/meep/issues/477
[#479]: https://github.com/NanoComp/meep/issues/479
[#483]: https://github.com/NanoComp/meep/issues/483
[#500]: https://github.com/NanoComp/meep/issues/500
[#518]: https://github.com/NanoComp/meep/issues/518
[#531]: https://github.com/NanoComp/meep/issues/531
[#558]: https://github.com/NanoComp/meep/issues/558
[#559]: https://github.com/NanoComp/meep/issues/559
[#570]: https://github.com/NanoComp/meep/issues/570
[#577]: https://github.com/NanoComp/meep/issues/577
[#578]: https://github.com/NanoComp/meep/issues/578
[#581]: https://github.com/NanoComp/meep/issues/581
[#593]: https://github.com/NanoComp/meep/issues/593
[#596]: https://github.com/NanoComp/meep/issues/596
[#599]: https://github.com/NanoComp/meep/issues/599
[#602]: https://github.com/NanoComp/meep/issues/602
[#603]: https://github.com/NanoComp/meep/issues/603
[#626]: https://github.com/NanoComp/meep/issues/626
[#630]: https://github.com/NanoComp/meep/issues/630
[#631]: https://github.com/NanoComp/meep/issues/631
[#644]: https://github.com/NanoComp/meep/issues/644
[#652]: https://github.com/NanoComp/meep/issues/652
[#655]: https://github.com/NanoComp/meep/issues/655
[#666]: https://github.com/NanoComp/meep/issues/666
[#671]: https://github.com/NanoComp/meep/issues/671
[#675]: https://github.com/NanoComp/meep/issues/675
[#681]: https://github.com/NanoComp/meep/issues/681
[#689]: https://github.com/NanoComp/meep/issues/689
[#699]: https://github.com/NanoComp/meep/issues/699
[#705]: https://github.com/NanoComp/meep/issues/705
[#728]: https://github.com/NanoComp/meep/issues/728
[#744]: https://github.com/NanoComp/meep/issues/744
[#747]: https://github.com/NanoComp/meep/issues/747
[#752]: https://github.com/NanoComp/meep/issues/752
[#769]: https://github.com/NanoComp/meep/issues/769
[#771]: https://github.com/NanoComp/meep/issues/771
[#779]: https://github.com/NanoComp/meep/issues/779
[#785]: https://github.com/NanoComp/meep/issues/785
[#787]: https://github.com/NanoComp/meep/issues/787
[#789]: https://github.com/NanoComp/meep/issues/789
[#791]: https://github.com/NanoComp/meep/issues/791
[#795]: https://github.com/NanoComp/meep/issues/795
[#807]: https://github.com/NanoComp/meep/issues/807
[#817]: https://github.com/NanoComp/meep/issues/817
[#860]: https://github.com/NanoComp/meep/issues/860
[#862]: https://github.com/NanoComp/meep/issues/862
[#863]: https://github.com/NanoComp/meep/issues/863
[#868]: https://github.com/NanoComp/meep/issues/868
[#869]: https://github.com/NanoComp/meep/issues/869
[#872]: https://github.com/NanoComp/meep/issues/872
[#876]: https://github.com/NanoComp/meep/issues/876
[#891]: https://github.com/NanoComp/meep/issues/891
[#894]: https://github.com/NanoComp/meep/issues/894
[#919]: https://github.com/NanoComp/meep/issues/919
[#922]: https://github.com/NanoComp/meep/issues/922
[#927]: https://github.com/NanoComp/meep/issues/927
[#940]: https://github.com/NanoComp/meep/issues/940
[#945]: https://github.com/NanoComp/meep/issues/945
[#952]: https://github.com/NanoComp/meep/issues/952
[#953]: https://github.com/NanoComp/meep/issues/953
[#960]: https://github.com/NanoComp/meep/issues/960
[#994]: https://github.com/NanoComp/meep/issues/994
[#1002]: https://github.com/NanoComp/meep/issues/1002
[#1041]: https://github.com/NanoComp/meep/issues/1041
[#1042]: https://github.com/NanoComp/meep/issues/1042
[#1047]: https://github.com/NanoComp/meep/issues/1047
[#1062]: https://github.com/NanoComp/meep/issues/1062
[#1078]: https://github.com/NanoComp/meep/issues/1078
[#1090]: https://github.com/NanoComp/meep/issues/1090
[#1095]: https://github.com/NanoComp/meep/issues/1095
[#1112]: https://github.com/NanoComp/meep/issues/1112
[#1121]: https://github.com/NanoComp/meep/issues/1121
[#1122]: https://github.com/NanoComp/meep/issues/1122
[#1126]: https://github.com/NanoComp/meep/issues/1126
[#1129]: https://github.com/NanoComp/meep/issues/1129
[#1132]: https://github.com/NanoComp/meep/issues/1132
[#1154]: https://github.com/NanoComp/meep/issues/1154
[#1156]: https://github.com/NanoComp/meep/issues/1156
[#1158]: https://github.com/NanoComp/meep/issues/1158
[#1159]: https://github.com/NanoComp/meep/issues/1159
[#1167]: https://github.com/NanoComp/meep/issues/1167
[#1171]: https://github.com/NanoComp/meep/issues/1171
[#1192]: https://github.com/NanoComp/meep/issues/1192
[#1205]: https://github.com/NanoComp/meep/issues/1205
[#1218]: https://github.com/NanoComp/meep/issues/1218
[#1233]: https://github.com/NanoComp/meep/issues/1233
[#1240]: https://github.com/NanoComp/meep/issues/1240
[#1242]: https://github.com/NanoComp/meep/issues/1242
[#1244]: https://github.com/NanoComp/meep/issues/1244
[#1257]: https://github.com/NanoComp/meep/issues/1257
[#1285]: https://github.com/NanoComp/meep/issues/1285
[#1286]: https://github.com/NanoComp/meep/issues/1286
[#1288]: https://github.com/NanoComp/meep/issues/1288
[#1302]: https://github.com/NanoComp/meep/issues/1302
[#1303]: https://github.com/NanoComp/meep/issues/1303
[#1310]: https://github.com/NanoComp/meep/issues/1310
[#1316]: https://github.com/NanoComp/meep/issues/1316
[#1333]: https://github.com/NanoComp/meep/issues/1333
[#1349]: https://github.com/NanoComp/meep/issues/1349
[#1374]: https://github.com/NanoComp/meep/issues/1374
[#1388]: https://github.com/NanoComp/meep/issues/1388
[#1403]: https://github.com/NanoComp/meep/issues/1403
[#1417]: https://github.com/NanoComp/meep/issues/1417
[#1427]: https://github.com/NanoComp/meep/issues/1427
[#1437]: https://github.com/NanoComp/meep/issues/1437
[#1447]: https://github.com/NanoComp/meep/issues/1447
[#1456]: https://github.com/NanoComp/meep/issues/1456
[#1458]: https://github.com/NanoComp/meep/issues/1458
[#1464]: https://github.com/NanoComp/meep/issues/1464
[#1487]: https://github.com/NanoComp/meep/issues/1487
[#1499]: https://github.com/NanoComp/meep/issues/1499
[#1512]: https://github.com/NanoComp/meep/issues/1512
[#1515]: https://github.com/NanoComp/meep/issues/1515
[#1519]: https://github.com/NanoComp/meep/issues/1519
[#1521]: https://github.com/NanoComp/meep/issues/1521
[#1522]: https://github.com/NanoComp/meep/issues/1522
[#1527]: https://github.com/NanoComp/meep/issues/1527
[#1528]: https://github.com/NanoComp/meep/issues/1528
[#1539]: https://github.com/NanoComp/meep/issues/1539
[#1544]: https://github.com/NanoComp/meep/issues/1544
[#1546]: https://github.com/NanoComp/meep/issues/1546
[#1568]: https://github.com/NanoComp/meep/issues/1568
[#1569]: https://github.com/NanoComp/meep/issues/1569
[#1574]: https://github.com/NanoComp/meep/issues/1574
[#1577]: https://github.com/NanoComp/meep/issues/1577
[#1588]: https://github.com/NanoComp/meep/issues/1588
[#1591]: https://github.com/NanoComp/meep/issues/1591
[#1592]: https://github.com/NanoComp/meep/issues/1592
[#1593]: https://github.com/NanoComp/meep/issues/1593
[#1598]: https://github.com/NanoComp/meep/issues/1598
[#1599]: https://github.com/NanoComp/meep/issues/1599
[#1601]: https://github.com/NanoComp/meep/issues/1601
[#1603]: https://github.com/NanoComp/meep/issues/1603
[#1606]: https://github.com/NanoComp/meep/issues/1606
[#1608]: https://github.com/NanoComp/meep/issues/1608
[#1618]: https://github.com/NanoComp/meep/issues/1618
[#1623]: https://github.com/NanoComp/meep/issues/1623
[#1628]: https://github.com/NanoComp/meep/issues/1628
[#1634]: https://github.com/NanoComp/meep/issues/1634
[#1635]: https://github.com/NanoComp/meep/issues/1635
[#1651]: https://github.com/NanoComp/meep/issues/1651
[#1652]: https://github.com/NanoComp/meep/issues/1652
[#1655]: https://github.com/NanoComp/meep/issues/1655
[#1656]: https://github.com/NanoComp/meep/issues/1656
[#1675]: https://github.com/NanoComp/meep/issues/1675
[#1684]: https://github.com/NanoComp/meep/issues/1684
[#1692]: https://github.com/NanoComp/meep/issues/1692
[#1704]: https://github.com/NanoComp/meep/issues/1704
[#1715]: https://github.com/NanoComp/meep/issues/1715
[#1720]: https://github.com/NanoComp/meep/issues/1720
[#1721]: https://github.com/NanoComp/meep/issues/1721
[#1722]: https://github.com/NanoComp/meep/issues/1722
[#1730]: https://github.com/NanoComp/meep/issues/1730
[#1732]: https://github.com/NanoComp/meep/issues/1732
[#1738]: https://github.com/NanoComp/meep/issues/1738
[#1740]: https://github.com/NanoComp/meep/issues/1740
[#1749]: https://github.com/NanoComp/meep/issues/1749
[#1763]: https://github.com/NanoComp/meep/issues/1763
[#1769]: https://github.com/NanoComp/meep/issues/1769
[#1775]: https://github.com/NanoComp/meep/issues/1775
[#1780]: https://github.com/NanoComp/meep/issues/1780
[#1796]: https://github.com/NanoComp/meep/issues/1796
[#1801]: https://github.com/NanoComp/meep/issues/1801
[#1804]: https://github.com/NanoComp/meep/issues/1804
[#1821]: https://github.com/NanoComp/meep/issues/1821
[#1826]: https://github.com/NanoComp/meep/issues/1826
[#1830]: https://github.com/NanoComp/meep/issues/1830
[#1839]: https://github.com/NanoComp/meep/issues/1839
[#1849]: https://github.com/NanoComp/meep/issues/1849
[#1855]: https://github.com/NanoComp/meep/issues/1855
[#1860]: https://github.com/NanoComp/meep/issues/1860
[#1871]: https://github.com/NanoComp/meep/issues/1871
[#1872]: https://github.com/NanoComp/meep/issues/1872
[#1873]: https://github.com/NanoComp/meep/issues/1873
[#1877]: https://github.com/NanoComp/meep/issues/1877
[#1895]: https://github.com/NanoComp/meep/issues/1895
[#1919]: https://github.com/NanoComp/meep/issues/1919
[#1955]: https://github.com/NanoComp/meep/issues/1955
[#1959]: https://github.com/NanoComp/meep/issues/1959
[#1968]: https://github.com/NanoComp/meep/issues/1968
[#2005]: https://github.com/NanoComp/meep/issues/2005
[#2016]: https://github.com/NanoComp/meep/issues/2016
[#2021]: https://github.com/NanoComp/meep/issues/2021
[#2027]: https://github.com/NanoComp/meep/issues/2027
[#2032]: https://github.com/NanoComp/meep/issues/2032
[#2044]: https://github.com/NanoComp/meep/issues/2044
[#2049]: https://github.com/NanoComp/meep/issues/2049
[#2053]: https://github.com/NanoComp/meep/issues/2053
[#2066]: https://github.com/NanoComp/meep/issues/2066
[#2073]: https://github.com/NanoComp/meep/issues/2073
[#2076]: https://github.com/NanoComp/meep/issues/2076
[#2077]: https://github.com/NanoComp/meep/issues/2077
[#2079]: https://github.com/NanoComp/meep/issues/2079
[#2082]: https://github.com/NanoComp/meep/issues/2082
[#2091]: https://github.com/NanoComp/meep/issues/2091
[#2095]: https://github.com/NanoComp/meep/issues/2095
[#2114]: https://github.com/NanoComp/meep/issues/2114
[#2176]: https://github.com/NanoComp/meep/issues/2176
[#2179]: https://github.com/NanoComp/meep/issues/2179
[#2186]: https://github.com/NanoComp/meep/issues/2186
[#2190]: https://github.com/NanoComp/meep/issues/2190
[#2194]: https://github.com/NanoComp/meep/issues/2194
[#2202]: https://github.com/NanoComp/meep/issues/2202
[#2203]: https://github.com/NanoComp/meep/issues/2203
[#2207]: https://github.com/NanoComp/meep/issues/2207
[#2208]: https://github.com/NanoComp/meep/issues/2208
[#2251]: https://github.com/NanoComp/meep/issues/2251
[#2253]: https://github.com/NanoComp/meep/issues/2253
[#2264]: https://github.com/NanoComp/meep/issues/2264
[#2271]: https://github.com/NanoComp/meep/issues/2271
[#2289]: https://github.com/NanoComp/meep/issues/2289
[#2290]: https://github.com/NanoComp/meep/issues/2290
[#2305]: https://github.com/NanoComp/meep/issues/2305
[#2314]: https://github.com/NanoComp/meep/issues/2314
[#2321]: https://github.com/NanoComp/meep/issues/2321
[#2333]: https://github.com/NanoComp/meep/issues/2333
[#2337]: https://github.com/NanoComp/meep/issues/2337
[#2349]: https://github.com/NanoComp/meep/issues/2349
[#2360]: https://github.com/NanoComp/meep/issues/2360
[#2364]: https://github.com/NanoComp/meep/issues/2364
[#2365]: https://github.com/NanoComp/meep/issues/2365
[#2371]: https://github.com/NanoComp/meep/issues/2371
[#2380]: https://github.com/NanoComp/meep/issues/2380
[#2382]: https://github.com/NanoComp/meep/issues/2382
[#2383]: https://github.com/NanoComp/meep/issues/2383
[#2387]: https://github.com/NanoComp/meep/issues/2387
[#2390]: https://github.com/NanoComp/meep/issues/2390
[#2394]: https://github.com/NanoComp/meep/issues/2394
[#2395]: https://github.com/NanoComp/meep/issues/2395
[#2402]: https://github.com/NanoComp/meep/issues/2402
[#2413]: https://github.com/NanoComp/meep/issues/2413
