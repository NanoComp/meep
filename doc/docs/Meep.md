---
title: Meep
permalink: /Meep/
---


![440px|center|Meep logo banner](images/Meep-banner.png)

 __NOTOC__ **Meep** (or [MEEP](Meep_acronym_expansions.md)) is a free finite-difference time-domain (FDTD) simulation software package developed at MIT to model electromagnetic systems, along with our [MPB](http://ab-initio.mit.edu/wiki/index.php/MPB) eigenmode package. Its features include:

-   **Free software** under the [GNU GPL](https://en.wikipedia.org/wiki/GNU_General_Public_License).
-   Simulation in **1d, 2d, 3d**, and **cylindrical** coordinates.
-   Distributed memory **parallelism** on any system supporting the [MPI](https://en.wikipedia.org/wiki/MPI) standard. Portable to any Unix-like system ([GNU/Linux](https://en.wikipedia.org/wiki/Linux) is fine).
-   Arbitrary **anisotropic** electric permittivity ε and magnetic permeability μ, along with **dispersive** ε(ω) and μ(ω) (including loss/gain) and **nonlinear** (Kerr & Pockels) dielectric and magnetic materials, and electric/magnetic **conductivities** σ.
-   **PML** absorbing boundaries and/or perfect conductor and/or **Bloch-periodic** boundary conditions.
-   Exploitation of **symmetries** to reduce the computation size — even/odd mirror symmetries and 90°/180° rotations.
-   Complete **scriptability** — either via a [Scheme](https://en.wikipedia.org/wiki/Scheme_programming_language) scripting front-end (as in [libctl](http://ab-initio.mit.edu/wiki/index.php/Libctl) and [MPB](http://ab-initio.mit.edu/wiki/index.php/MPB)), or callable as a [C++](https://en.wikipedia.org/wiki/C_plus_plus) library; a [Python](https://en.wikipedia.org/wiki/Python_programming_language) interface is also available.
-   Field output in the [HDF5](https://en.wikipedia.org/wiki/HDF5) standard scientific data format, supported by many visualization tools.
-   Arbitrary material and source distributions.
-   Field analyses including flux spectra, Maxwell stress tensor, frequency extraction, local density of states and energy integrals, near to far field transformations; completely programmable.
-   Multi-parameter optimization, root-finding, integration, etcetera (via [libctl](http://ab-initio.mit.edu/wiki/index.php/Libctl)).

*Meep* officially stands for *MIT Electromagnetic Equation Propagation*, but we also have [several unofficial meanings](Meep_acronym_expansions.md) of the acronym.

Time-domain simulation
----------------------

A time-domain electromagnetic simulation simply takes [Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell's_equations) and evolves them over time within some finite computational region, essentially performing a kind of **numerical experiment**. This can be used to calculate a wide variety of useful quantities, but major applications include:

-   **Transmission and reflection spectra** — by Fourier-transforming the response to a short pulse, a single simulation can yield the scattering amplitudes over a wide spectrum of frequencies.
-   **Resonant modes and frequencies** — by analyzing the response of the system to a short pulse, one can extract the frequencies, decay rates, and field patterns of the harmonic modes of a system (including waveguide and cavity modes, and including losses).
-   **Field patterns** (e.g. Green's functions) in response to an arbitrary source, archetypically a [CW](https://en.wikipedia.org/wiki/Continuous_wave) (fixed-ω) input.

Using these results, one can then compute many other things, such as the local density of states (from the trace of the Green's function). Meep's scriptable interface makes it possible to combine many sorts of computations (along with multi-parameter optimization etcetera) in sequence or in parallel.

The [Meep manual](Meep_manual.md) gives examples of all of these kinds of computations.

Download
--------

Please see the [Meep Download](Meep_Download.md) page to get the latest version of Meep; the differences between versions are described in the [Meep release notes](Meep_release_notes.md). The installation instructions can be found in the [installation section](Meep_installation.md) of the [Meep manual](Meep_manual.md).

The latest development sources are available on [Github](https://github.com/stevengj/meep).

Documentation
-------------

See the [Meep manual](Meep_manual.md), and also the navigation sidebar at right. In particular, the [Meep Introduction](Meep_Introduction.md) and [Meep Tutorial](Meep_Tutorial.md) are the most important things to read. We also have a [Meep FAQ](Meep_FAQ.md).

Please [cite Meep](Citing_Meep.md) in any publication for which you found it useful.

### Mailing Lists

The Meep mailing lists and their archives are another source of information about Meep.

Subscribe to the read-only [meep-announce mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/meep-announce) to receive notifications of updates and releases. Subscribe to the unmoderated [meep-discuss mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/meep-discuss) for discussions about using Meep. Archives are available [here](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/). You can also read and post to the list via the [gmane.comp.science.electromagnetism.meep.general](news://news.gmane.org/gmane.comp.science.electromagnetism.meep.general) newsgroup from [Gmane](http://www.gmane.org/).

### Bug reports and feature requests

For bug reports and feature requests, please [file a Meep Github issue](https://github.com/stevengj/meep/issues).

Acknowledgements
----------------

Meep's active developers are Ardavan Oskooi and Steven G. Johnson. Please see the [Meep Acknowledgements](Meep_Acknowledgements.md) for a more complete listing of those to whom we are grateful.

Contacts and Feedback
---------------------

If you have questions or problems regarding Meep, you are encouraged to query the [mailing list](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/).

For professional consulting as well as free access to Meep in the public cloud via Amazon Web Services (AWS), see [Simpetus](http://www.simpetuscloud.com).

[Category:Meep](Meep.md)
