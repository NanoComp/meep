![](doc/docs/images/Meep-banner.png)

[![Latest Docs](https://readthedocs.org/projects/meep/badge/?version=latest)](http://meep.readthedocs.io/en/latest/)
[![Build Status](https://travis-ci.org/NanoComp/meep.svg?branch=master)](https://travis-ci.org/NanoComp/meep)
[![Coverage Status](https://coveralls.io/repos/github/stevengj/meep/badge.svg?branch=master)](https://coveralls.io/github/stevengj/meep?branch=master)
![Python versions 2.7–3.6](https://img.shields.io/badge/python-2.7%2C%203.4%2C%203.5%2C%203.6-brightgreen.svg)

**Meep** is a free and open-source software package for [electromagnetics](https://en.wikipedia.org/wiki/Electromagnetism) simulation via the [finite-difference time-domain](https://en.wikipedia.org/wiki/Finite-difference_time-domain_method) (FDTD) method spanning a broad range of applications.

## Key Features

-   **Free and open-source software** under the [GNU GPL](https://en.wikipedia.org/wiki/GNU_General_Public_License).
-   Complete **scriptability** via [Python](https://meep.readthedocs.io/en/latest/Python_Tutorials/Basics/), [Scheme](https://meep.readthedocs.io/en/latest/Scheme_Tutorials/Basics), or [C++](https://meep.readthedocs.io/en/master/C++_Tutorial/) APIs.
-   Simulation in **1d, 2d, 3d**, and **cylindrical** coordinates.
-   Distributed memory **parallelism** on any system supporting [MPI](https://en.wikipedia.org/wiki/MPI).
-   Portable to any Unix-like operating system such as [Linux](https://en.wikipedia.org/wiki/Linux), [macOS](https://en.wikipedia.org/wiki/macOS), and [FreeBSD](https://en.wikipedia.org/wiki/FreeBSD).
-   **Precompiled binary packages** of official releases and nightly builds via [Conda](https://meep.readthedocs.io/en/latest/Installation/#conda-packages).
-   Arbitrary **anisotropic** electric permittivity ε and magnetic permeability μ, along with **dispersive** ε(ω) and μ(ω) including loss/gain, **nonlinear** (Kerr & Pockels) dielectric and magnetic materials, electric/magnetic **conductivities** σ, **saturable** gain/absorption, and **gyrotropic** media (magneto-optical effects).
-   **Perfectly-matched layer** (**PML**) absorbing boundaries as well as **Bloch-periodic** and perfect-conductor boundary conditions.
-   Exploitation of **symmetries** to reduce the computation size, including even/odd mirror planes and 90°/180° rotations.
-   [Subpixel smoothing](https://meep.readthedocs.io/en/latest/Subpixel_Smoothing/) for improving accuracy and shape optimization.
-   Arbitrary current sources including a [mode launcher](https://meep.readthedocs.io/en/latest/Python_Tutorials/Eigenmode_Source/).
-   Frequency-domain solver for finding the response to a [continuous-wave](https://en.wikipedia.org/wiki/Continuous_wave) (CW) source.
-   ε/μ and field import/export in the [HDF5](https://en.wikipedia.org/wiki/HDF5) data format.
-   **GDSII** file import for planar geometries.
-   **Materials library** containing predefined broadband, complex refractive indices.
-   Field analyses including [Poynting flux](https://meep.readthedocs.io/en/latest/Python_Tutorials/Basics/#transmittance-spectrum-of-a-waveguide-bend), [mode decomposition](https://meep.readthedocs.io/en/latest/Python_Tutorials/Mode_Decomposition/) (for [S-parameters](https://meep.readthedocs.io/en/latest/Python_Tutorials/GDSII_Import/)), [energy density](https://meep.readthedocs.io/en/latest/Python_User_Interface/#energy-density-spectra), [near to far transformations](https://meep.readthedocs.io/en/latest/Python_Tutorials/Near_to_Far_Field_Spectra/), [frequency extraction](https://meep.readthedocs.io/en/latest/Python_Tutorials/Basics/#modes-of-a-ring-resonator), [local density of states](https://meep.readthedocs.io/en/latest/Python_Tutorials/Local_Density_of_States/) (LDOS), [modal volume](https://meep.readthedocs.io/en/latest/Python_User_Interface/#field-computations), [Maxwell stress tensor](https://meep.readthedocs.io/en/latest/Python_Tutorials/Optical_Forces/), [arbitrary functions](https://meep.readthedocs.io/en/latest/Field_Functions/); completely programmable.
-   [Adjoint solver](https://meep.readthedocs.io/en/latest/Python_Tutorials/AdjointSolver) for **sensitivity analysis** and **automated design optimization**.
-   [Visualization routines](https://meep.readthedocs.io/en/latest/Python_User_Interface/#data-visualization) for the simulation domain involving geometries, fields, boundary layers, sources, and monitors.

## Citing Meep

We kindly request that you cite the following paper in any published work for which you used Meep:

- A. Oskooi, D. Roundy, M. Ibanescu, P. Bermel, J.D. Joannopoulos, and S.G. Johnson, [MEEP: A flexible free-software package for electromagnetic simulations by the FDTD method](http://dx.doi.org/doi:10.1016/j.cpc.2009.11.008), Computer Physics Communications, Vol. 181, pp. 687-702, 2010 ([pdf](http://ab-initio.mit.edu/~oskooi/papers/Oskooi10.pdf)).


## Documentation

See the [manual on readthedocs](https://meep.readthedocs.io/en/latest) for the latest documentation.
