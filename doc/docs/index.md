<center>
![](images/Meep-banner.png)
</center>

 **Meep** is a free and open-source software package for simulating electromagnetic systems via the [finite-difference time-domain](https://en.wikipedia.org/wiki/Finite-difference_time-domain_method) (FDTD) method. Meep is an acronym for *MIT Electromagnetic Equation Propagation*.

**Features**

-   **Free and open-source software** under the [GNU GPL](https://en.wikipedia.org/wiki/GNU_General_Public_License).
-   Complete **scriptability** via [Python](Python_Tutorials/Basics), [Scheme](Scheme_Tutorials/Basics), or [C++](C++_Tutorial) APIs.
-   Simulation in **1d, 2d, 3d**, and **cylindrical** coordinates.
-   Distributed memory **parallelism** on any system supporting [MPI](https://en.wikipedia.org/wiki/MPI).
-   Portable to any Unix-like operating system such as [Linux](https://en.wikipedia.org/wiki/Linux), [macOS](https://en.wikipedia.org/wiki/macOS), and [FreeBSD](https://en.wikipedia.org/wiki/FreeBSD).
-   Arbitrary **anisotropic**, electric permittivity ε and magnetic permeability μ, along with **dispersive** ε(ω) and μ(ω) including loss/gain, **nonlinear** (Kerr & Pockels) dielectric and magnetic materials, and electric/magnetic **conductivities** σ.
-   **Perfectly-matched layer** (**PML**) absorbing boundaries as well as perfect conductor and **Bloch-periodic** boundary conditions.
-   Exploitation of **symmetries** to reduce the computation size: even/odd mirror planes and 90°/180° rotations.
-   Field output in the [HDF5](https://en.wikipedia.org/wiki/HDF5) data format.
-   Arbitrary current sources including a guided-mode launcher.
-   Materials library containing list of predefined broadband, complex refractive indices.
-   Frequency-domain solver for finding the response to a continuous-wave source.
-   Field analyses including flux spectra, near to far transformations, modal decomposition, frequency extraction, local density of states, modal volume, Maxwell stress tensor, arbitrary functions; completely programmable.

Time-Domain Simulation
----------------------

A time-domain electromagnetic simulation simply evolves [Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell's_equations) over time within some finite computational volume, essentially performing a kind of **numerical experiment**. This can be used to calculate a wide variety of useful quantities. Major applications include:

-   **Transmittance and Reflectance Spectra** &mdash; by Fourier-transforming the response to a short pulse, a single simulation can yield the scattering amplitudes over a broadband spectrum.
-   **Resonant Modes and Frequencies** &mdash; by analyzing the response of the system to a short pulse, one can extract the frequencies, decay rates, and field patterns of the harmonic modes of lossy and lossless systems including waveguide and cavity modes.
-   **Field Patterns** (e.g. Green's functions) &mdash; in response to an arbitrary source via a [continuous-wave](https://en.wikipedia.org/wiki/Continuous_wave) (CW) input (fixed-ω).

Meep's scriptable interface makes it possible to combine many sorts of computations along with multi-parameter optimization etcetera in sequence or in parallel.

[Tutorial/Basics](Python_Tutorials/Basics.md) provides examples of all these kinds of computations.

Download
--------

The source repository is on [GitHub](https://github.com/stevengj/meep). Gzipped tarballs of stable versions are available in [Releases](https://github.com/stevengj/meep/releases). The release history is described in [NEWS](https://github.com/stevengj/meep/blob/master/NEWS.md). Installation instructions are in [Installation](Installation.md).

Documentation
-------------

See the left navigation sidebar. In particular, the [Introduction](Introduction.md), [Tutorial/Basics](Python_Tutorials/Basics.md), and [FAQ](FAQ.md) are the most important things to review.

### Mailing Lists

Subscribe to the read-only [meep-announce mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/meep-announce) to receive notifications of updates and releases. Subscribe to the [meep-discuss mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/meep-discuss) for discussions regarding using Meep. The [meep-discuss archives](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/) includes all postings since 2006 spanning a large number and variety of discussion topics related to installation, setting up simulations, post-processing output, etc.

### Bug Reports and Feature Requests

For bug reports and feature requests, please file a [GitHub issue](https://github.com/stevengj/meep/issues).

Acknowledgements
----------------

The Meep project is maintained by [Simpetus](http://www.simpetus.com) and the open-source developer community on [GitHub](https://github.com/stevengj/meep). [Acknowledgements](Acknowledgements.md) provides a complete listing of the project contributors.

Support and Feedback
---------------------

If you have questions or problems regarding Meep, you are encouraged to query the [mailing list](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/).

Professional consulting services for photonic design and modeling including development of custom, turn-key simulation modules, training, technical support, and access to Meep in the public cloud via Amazon Web Services (AWS) are provided by [Simpetus](http://www.simpetus.com).
