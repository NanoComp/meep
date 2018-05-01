# Meep Release Notes

## Meep 1.5 alpha

 * Python interface to MPB (#191 etc.).

 * Mode decomposition: given a DFT flux plane, decompose the fields
   at each frequency into the power in each mode of a waveguide
   or similar (#192, #248, etc.).

 * DFT slices: output Fourier-transformed fields in any given
   region of space (#259).

 * Structure dump/load feature to rapidly load in a geometry
   from a previous calculation (#261, #266).

 * Susceptibilities are now supported in user-defined materials
   in Python (#203, #305).

 * 64-bit support for extremely large computations (#193).

 * Various bug fixes, documentation improvements, etc.

## Meep 1.4.3

2/1/2018

 * Allow `meep` Python module to be imported without setting `PYTHONPATH` (#189).

## Meep 1.4.2

1/26/2018

  * Build fix for Python due to missing file (#184).

## Meep 1.4.1

1/19/2018

  * Minor packaging fixes.

## Meep 1.4

1/18/2018

  * Full-featured Python interface.

  * Migrated documentation to github/markdown/readthedocs (#55).

  * New feature to get slice as array in C++ and Python APIs (#96, #105).

  * `libmeepgeom` library to allow C++ users to access geometric-object
    API (#56).

  * Removed overly conservative stability check for Lorentzian
    susceptibilities (#150).

  * Corrected small error in frequency interval for `dft-ldos` (#40).

  * Bug fixes in near-to-farfield spectra (#21), eigenmode source (#20),
    and LDOS (#40).

## Meep 1.3

31 March 2015.

  * New near-to-far-field functionality: given a bounding surface,
    automatically computes the Fourier-transformed field in any 
    desired grid of "far-field" points arbitrarily far away.

  * Compatibility with Harminv 1.4 (fixes issue #13: ppc64 portability).

  * Fix compilation with latest C++ standard (e.g. on OS X 10.9).

  * Bug fix in CW solver convergence test; thanks to 
    Wu Chuanren and @FilipDominec for the bug report.

  * Build fix for Fedora 21 (thanks to Dean Brettle) (issue #14).

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
