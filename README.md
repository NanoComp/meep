![](doc/docs/images/Meep-banner.png)

[![Latest Docs](https://readthedocs.org/projects/meep/badge/?version=latest)](http://meep.readthedocs.io/en/latest/)
[![Build Status](https://travis-ci.org/stevengj/meep.svg?branch=master)](https://travis-ci.org/stevengj/meep)
![Python versions 2.7–3.6](https://img.shields.io/badge/python-2.7%2C%203.4%2C%203.5%2C%203.6-brightgreen.svg)

**Meep** is a free and open-source software package for simulating electromagnetic systems via the [finite-difference time-domain](https://en.wikipedia.org/wiki/Finite-difference_time-domain_method) (FDTD) method. Meep is an acronym for *MIT Electromagnetic Equation Propagation*.

**Features**

-   **Free and open-source software** under the [GNU GPL](https://en.wikipedia.org/wiki/GNU_General_Public_License).
-   Complete **scriptability** via [Python](https://meep.readthedocs.io/en/latest/Python_Tutorials/Basics/), [Scheme](https://meep.readthedocs.io/en/latest/Scheme_Tutorials/Basics), or [C++](https://meep.readthedocs.io/en/master/C++_Tutorial/) APIs.
-   Simulation in **1d, 2d, 3d**, and **cylindrical** coordinates.
-   Distributed memory **parallelism** on any system supporting the [MPI](https://en.wikipedia.org/wiki/MPI) standard.
-   Portable to any Unix-like operating system such as [Linux](https://en.wikipedia.org/wiki/Linux), [macOS](https://en.wikipedia.org/wiki/macOS), and [FreeBSD](https://en.wikipedia.org/wiki/FreeBSD).
-   **Precompiled binary packages** of official releases and nightly builds via [Conda](https://meep.readthedocs.io/en/latest/Installation/#conda-packages).
-   Arbitrary **anisotropic** electric permittivity ε and magnetic permeability μ, along with **dispersive** ε(ω) and μ(ω) including loss/gain, **nonlinear** (Kerr & Pockels) dielectric and magnetic materials, and electric/magnetic **conductivities** σ.
-   **PML** absorbing boundaries as well as **Bloch-periodic** and perfect-conductor boundary conditions.
-   Exploitation of **symmetries** to reduce the computation size, including even/odd mirror planes and 90°/180° rotations.
-   Arbitrary current sources including a guided-mode launcher.
-   Frequency-domain solver for finding the response to a [continuous-wave (CW)](https://en.wikipedia.org/wiki/Continuous_wave) source.
-   ε/μ and field import/export in the [HDF5](https://en.wikipedia.org/wiki/HDF5) data format.
-   **GDSII** file import for planar geometries.
-   **Materials library** containing predefined broadband, complex refractive indices.
-   Field analyses including Poynting flux, mode decomposition, near to far transformations, frequency extraction, local density of states (LDOS), modal volume, Maxwell stress tensor, arbitrary functions; completely programmable.

# Documentation

See the [manual on readthedocs](https://meep.readthedocs.io/en/latest) for the latest documentation.


