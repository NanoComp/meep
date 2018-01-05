![](doc/docs/images/Meep-banner.png)

[![Latest Docs](https://readthedocs.org/projects/meep/badge/?version=latest)](http://meep.readthedocs.io/en/latest/)
[![Build Status](https://travis-ci.org/stevengj/meep.svg?branch=master)](https://travis-ci.org/stevengj/meep)

# Meep

**Meep** is a free finite-difference time-domain (FDTD) simulation software package to model electromagnetic systems. Meep is an acronym which officially stands for *MIT Electromagnetic Equation Propagation*. Its features include:

-   **Free and open-source software** under the [GNU GPL](https://en.wikipedia.org/wiki/GNU_General_Public_License).
-   Complete **scriptability** via [Python](http://meep.readthedocs.io/en/latest/Python_Tutorials/Basics/), [Scheme](http://meep.readthedocs.io/en/latest/Scheme_Tutorials/Basics), or [C++](C++).
-   Simulation in **1d, 2d, 3d**, and **cylindrical** coordinates.
-   Distributed memory **parallelism** on any system supporting the [MPI](https://en.wikipedia.org/wiki/MPI) standard. Portable to any Unix-like operating system such as [Linux](https://en.wikipedia.org/wiki/Linux) and [macOS](https://en.wikipedia.org/wiki/macOS).
-   Arbitrary **anisotropic** electric permittivity ε and magnetic permeability μ, along with **dispersive** ε(ω) and μ(ω) including loss/gain and **nonlinear** (Kerr & Pockels) dielectric and magnetic materials, and electric/magnetic **conductivities** σ.
-   **PML** absorbing boundaries and/or perfect conductor and/or **Bloch-periodic** boundary conditions.
-   Exploitation of **symmetries** to reduce the computation size — even/odd mirror symmetries and 90°/180° rotations.
-   Field output in the [HDF5](https://en.wikipedia.org/wiki/HDF5) standard scientific data format, supported by many visualization tools.
-   Arbitrary material and source distributions.
-   Field analyses including flux spectra, Maxwell stress tensor, frequency extraction, local density of states, arbitrary functions, near to far field transformations; completely programmable.

# Documentation

See the [Meep manual on readthedocs](http://meep.readthedocs.io/en/latest) for the latest documentation.

# Download

The latest official release can be found at the [Meep download page](http://meep.readthedocs.io/en/latest/Download/).

# Git Source

To compile directly from a clone of the [Meep git repository](https://github.com/stevengj/meep), you need to run
```
sh autogen.sh
make
```
in the cloned directory in order to generate the necessary files. You will need GNU autotools and SWIG installed, among other packages.
