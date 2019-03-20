<center>
![](images/Meep-banner.png)
</center>

 **Meep** is a free and open-source software package for [electromagnetics](https://en.wikipedia.org/wiki/Electromagnetism) simulation via the [finite-difference time-domain](https://en.wikipedia.org/wiki/Finite-difference_time-domain_method) (**FDTD**) method spanning a broad range of applications.

**Key Features**

-   **Free and open-source software** under the [GNU GPL](https://en.wikipedia.org/wiki/GNU_General_Public_License).
-   Complete **scriptability** via [Python](Python_Tutorials/Basics), [Scheme](Scheme_Tutorials/Basics), or [C++](C++_Tutorial) APIs.
-   Simulation in **1d, 2d, 3d**, and **cylindrical** coordinates.
-   Distributed memory **parallelism** on any system supporting [MPI](https://en.wikipedia.org/wiki/MPI).
-   Portable to any Unix-like operating system such as [Linux](https://en.wikipedia.org/wiki/Linux), [macOS](https://en.wikipedia.org/wiki/macOS), and [FreeBSD](https://en.wikipedia.org/wiki/FreeBSD).
-   **Precompiled binary packages** of official releases and nightly builds via [Conda](Installation.md#conda-packages).
-   Arbitrary **anisotropic** electric permittivity ε and magnetic permeability μ, along with **dispersive** ε(ω) and μ(ω) including loss/gain, **nonlinear** (Kerr & Pockels) dielectric and magnetic materials, electric/magnetic **conductivities** σ, and **saturable** gain/absorption.
-   **Perfectly-matched layer** (**PML**) absorbing boundaries as well as **Bloch-periodic** and perfect-conductor boundary conditions.
-   Exploitation of **symmetries** to reduce the computation size, including even/odd mirror planes and 90°/180° rotations.
-   Arbitrary current sources including a [unidirectional mode launcher](Python_Tutorials/Eigenmode_Source.md).
-   Frequency-domain solver for finding the response to a [continuous-wave](https://en.wikipedia.org/wiki/Continuous_wave) (CW) source.
-   ε/μ and field import/export in the [HDF5](https://en.wikipedia.org/wiki/HDF5) data format.
-   **GDSII** file import for planar geometries.
-   **Materials library** containing predefined broadband, complex refractive indices.
-   Field analyses including [Poynting flux](Python_Tutorials/Basics.md#transmittance-spectrum-of-a-waveguide-bend), [mode decomposition](Python_Tutorials/Mode_Decomposition.md) (for [S-parameters](Python_Tutorials/GDSII_Import.md)), energy density, [near to far transformations](Python_Tutorials/Near_to_Far_Field_Spectra.md), [frequency extraction](Python_Tutorials/Basics.md#modes-of-a-ring-resonator), [local density of states](Python_Tutorials/Local_Density_of_States.md) (LDOS), [modal volume](Python_User_Interface.md#field-computations), [Maxwell stress tensor](Python_Tutorials/Optical_Forces.md), [arbitrary functions](Field_Functions.md); completely programmable.

Time-Domain Simulation
----------------------

A time-domain electromagnetic simulation simply evolves [Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell's_equations) over time within some finite computational volume, essentially performing a kind of **numerical experiment**. This can be used to calculate a wide variety of useful quantities. Major applications include:

-   **Transmittance and Reflectance Spectra** &mdash; by Fourier-transforming the response to a short pulse, a single simulation can yield the scattering amplitudes over a broadband spectrum.
-   **Resonant Modes and Frequencies** &mdash; by analyzing the response of the system to a short pulse, one can extract the frequencies, decay rates, and field patterns of the harmonic modes of lossy and lossless systems including waveguide and cavity modes.
-   **Field Patterns** (e.g. Green's functions) &mdash; in response to an arbitrary source via a [continuous-wave](https://en.wikipedia.org/wiki/Continuous_wave) (CW) input (fixed-ω).

Meep's scriptable interface makes it possible to combine many sorts of computations along with multi-parameter optimization in sequence or in parallel.

[Tutorial/Basics](Python_Tutorials/Basics.md) provides examples of the various kinds of computations.

Download
--------

The source repository is on [GitHub](https://github.com/NanoComp/meep). Gzipped tarballs of stable versions are in [Releases](https://github.com/NanoComp/meep/releases). The release history is in [NEWS](https://github.com/NanoComp/meep/blob/master/NEWS.md). Installation instructions are in [Installation](Installation.md).

Documentation
-------------

For a list of topics, see the left navigation sidebar. For first-time users, the [Introduction](Introduction.md), [Tutorial/Basics](Python_Tutorials/Basics.md), and [FAQ](FAQ.md) are the most important items to review.

This documentation is for the master branch of the [source repository](Download.md#github-source-repository).

### Mailing Lists

Subscribe to the read-only [meep-announce mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/meep-announce) to receive notifications of updates and releases. Subscribe to the [meep-discuss mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/meep-discuss) for discussions regarding using Meep. The [meep-discuss archives](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/) includes all postings since 2006 spanning a large number and variety of discussion topics related to installation, setting up simulations, post-processing output, etc. This list can also be accessed using a [newsgroup reader](https://en.wikipedia.org/wiki/List_of_Usenet_newsreaders) via the NNTP interface address: `news.gmane.org/gmane.comp.science.electromagnetism.meep.general`.

### Bug Reports and Feature Requests

For bug reports and feature requests, please file a [GitHub issue](https://github.com/NanoComp/meep/issues).

Acknowledgements
----------------

The Meep project is maintained by [Simpetus](http://www.simpetus.com) and the developer community on [GitHub](https://github.com/NanoComp/meep). [Acknowledgements](Acknowledgements.md) provides a complete listing of the project contributors.

Support and Feedback
---------------------

If you have questions or problems regarding Meep, you are encouraged to query the [mailing list](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/).

Professional consulting services for photonic design and modeling including development of custom, turn-key simulation modules, training, technical support, and access to Meep in the public cloud via Amazon Web Services (AWS) are provided by [Simpetus](http://www.simpetus.com).
