![](doc/docs/images/Meep-banner.png)

[![Latest Docs](https://readthedocs.org/projects/meep/badge/?version=latest)](http://meep.readthedocs.io/en/latest/)
[![Build Status](https://travis-ci.org/stevengj/meep.svg?branch=master)](https://travis-ci.org/stevengj/meep)

# Meep

**Meep** is a free and open-source software package for simulating electromagnetic systems via the [finite-difference time-domain](https://en.wikipedia.org/wiki/Finite-difference_time-domain_method) (FDTD) method. Meep is an acronym which officially stands for *MIT Electromagnetic Equation Propagation*. Features include:

-   **Free and open-source software** under the [GNU GPL](https://en.wikipedia.org/wiki/GNU_General_Public_License).
-   Complete **scriptability** via [Python](http://meep.readthedocs.io/en/latest/Python_Tutorials/Basics/), [Scheme](http://meep.readthedocs.io/en/latest/Scheme_Tutorials/Basics), or [C++](C++) APIs.
-   Simulation in **1d, 2d, 3d**, and **cylindrical** coordinates.
-   Distributed memory **parallelism** on any system supporting the [MPI](https://en.wikipedia.org/wiki/MPI) standard.
-   Portable to any Unix-like operating system such as [Linux](https://en.wikipedia.org/wiki/Linux), [macOS](https://en.wikipedia.org/wiki/macOS), and and [FreeBSD](https://en.wikipedia.org/wiki/FreeBSD).
-   Arbitrary **anisotropic** electric permittivity ε and magnetic permeability μ, along with **dispersive** ε(ω) and μ(ω) including loss/gain, **nonlinear** (Kerr & Pockels) dielectric and magnetic materials, and electric/magnetic **conductivities** σ.
-   **PML** absorbing boundaries and/or perfect conductor and/or **Bloch-periodic** boundary conditions.
-   Exploitation of **symmetries** to reduce the computation size — even/odd mirror symmetries and 90°/180° rotations.
-   Field output in the [HDF5](https://en.wikipedia.org/wiki/HDF5) standard scientific data format.
-   Arbitrary current sources including a guided-mode launcher.
-   Frequency-domain solver for finding the response to a continuous-wave source.
-   Field analyses including flux spectra, near to far field transformations, modal decomposition, frequency extraction, local density of states, modal volume, Maxwell stress tensor, arbitrary functions; completely programmable.

# Documentation

See the [manual on readthedocs](http://meep.readthedocs.io/en/latest) for the latest documentation.


