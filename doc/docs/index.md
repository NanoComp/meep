---
# Meep
---

<center>
![](images/Meep-banner.png)
</center>

 **Meep** is a free finite-difference time-domain (FDTD) simulation software package to model electromagnetic systems. Meep is an acronym which officially stands for *MIT Electromagnetic Equation Propagation*. Its features include:

-   **Free software** under the [GNU GPL](https://en.wikipedia.org/wiki/GNU_General_Public_License).
-   Complete **scriptability** via [Python](Python_Tutorials/Basics), [Scheme](Scheme_Tutorials/Basics), or [C++](C++_Tutorial).
-   Simulation in **1d, 2d, 3d**, and **cylindrical** coordinates.
-   Distributed memory **parallelism** on any system supporting the [MPI](https://en.wikipedia.org/wiki/MPI) standard. Portable to any Unix-like operating system such as [Linux](https://en.wikipedia.org/wiki/Linux) and [MacOS](https://en.wikipedia.org/wiki/MacOS).
-   Arbitrary **anisotropic** electric permittivity $\varepsilon$ and magnetic permeability $\mu$, along with **dispersive** $\varepsilon(\omega)$ and $\mu(\omega)$ including loss/gain and **nonlinear** (Kerr & Pockels) dielectric and magnetic materials, and electric/magnetic **conductivities** $\sigma$.
-   **PML** absorbing boundaries and/or perfect conductor and/or **Bloch-periodic** boundary conditions.
-   Exploitation of **symmetries** to reduce the computation size &mdash; even/odd mirror symmetries and 90°/180° rotations.
-   Field output in the [HDF5](https://en.wikipedia.org/wiki/HDF5) standard scientific data format, supported by many visualization tools.
-   Arbitrary material and source distributions.
-   Field analyses including flux spectra, Maxwell stress tensor, frequency extraction, local density of states, arbitrary functions, near to far field transformations; completely programmable.

Time-Domain Simulation
----------------------

A time-domain electromagnetic simulation simply takes [Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell's_equations) and evolves them over time within some finite computational region, essentially performing a kind of **numerical experiment**. This can be used to calculate a wide variety of useful quantities, but major applications include:

-   **Transmission and Reflection Spectra** &mdash; by Fourier-transforming the response to a short pulse, a single simulation can yield the scattering amplitudes over a wide spectrum of frequencies.
-   **Resonant Modes and Frequencies** &mdash; by analyzing the response of the system to a short pulse, one can extract the frequencies, decay rates, and field patterns of the harmonic modes of lossy and lossless systems including waveguide and cavity modes.
-   **Field Patterns** (e.g. Green's functions) &mdash; in response to an arbitrary source, archetypically a [CW](https://en.wikipedia.org/wiki/Continuous_wave) (fixed-$\omega$) input.

Meep's scriptable interface makes it possible to combine many sorts of computations along with multi-parameter optimization etcetera in sequence or in parallel.

[Tutorial/Basics](Python_Tutorials/Basics.md) gives examples of all of these kinds of computations.

Download
--------

The latest development sources are available on [GitHub](https://github.com/stevengj/meep). The source tarballs are available on the [Download](Download.md) page. The release history is described in the [Release Notes](Release_Notes.md). The installation instructions can be found in the [Installation](Installation.md) page.

Documentation
-------------

See the navigation sidebar at left. In particular, the [Introduction](Introduction.md) and [Tutorial/Basics](Python_Tutorials/Basics.md) are the most important things to read. There is also an [FAQ](FAQ.md).

Please [cite Meep](Acknowledgements.md#referencing) in any publication for which you found it useful.

### Mailing Lists

Subscribe to the read-only [meep-announce mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/meep-announce) to receive notifications of updates and releases. Subscribe to the [meep-discuss mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/meep-discuss) for discussions about using Meep. The [meep-discuss archives](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/) includes all postings since 2006 spanning a large number and variety of discussion topics related to installation, setting up simulations, post-processing output, etc.

### Bug Reports and Feature Requests

For bug reports and feature requests, please [file a Meep GitHub issue](https://github.com/stevengj/meep/issues).

Acknowledgements
----------------

The Meep project is maintained by [Simpetus](http://www.simpetus.com) and the open-source community on [GitHub](https://github.com/stevengj/meep). Please see the [Acknowledgements](Acknowledgements.md) for a more complete listing of the project contributors.

Contacts and Feedback
---------------------

If you have questions or problems regarding Meep, you are encouraged to query the [mailing list](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/).

Professional consulting services for photonic design and modeling including development of turn-key simulation modules, training and technical support for getting up and running with Meep as well as free access to Meep in the public cloud via Amazon Web Services (AWS) are provided by [Simpetus](http://www.simpetus.com).
