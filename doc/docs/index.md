---
# Meep
---

<center>
![](images/Meep-banner.png)
</center>

 **Meep** is a free finite-difference time-domain (FDTD) simulation software package to model electromagnetic systems. Meep is an acronym which officially stands for *MIT Electromagnetic Equation Propagation*. Its features include:

-   **Free software** under the [GNU GPL](https://en.wikipedia.org/wiki/GNU_General_Public_License).
-   Simulation in **1d, 2d, 3d**, and **cylindrical** coordinates.
-   Distributed memory **parallelism** on any system supporting the [MPI](https://en.wikipedia.org/wiki/MPI) standard. Portable to any Unix-like operating system such as [Linux](https://en.wikipedia.org/wiki/Linux) and [MacOS](https://en.wikipedia.org/wiki/MacOS).
-   Arbitrary **anisotropic** electric permittivity $\varepsilon$ and magnetic permeability $\mu$, along with **dispersive** $\varepsilon(\omega)$ and $\mu(\omega)$ including loss/gain and **nonlinear** (Kerr & Pockels) dielectric and magnetic materials, and electric/magnetic **conductivities** $\sigma$.
-   **PML** absorbing boundaries and/or perfect conductor and/or **Bloch-periodic** boundary conditions.
-   Exploitation of **symmetries** to reduce the computation size &mdash; even/odd mirror symmetries and 90°/180° rotations.
-   Complete **scriptability** &mdash; either via a [Scheme](Scheme_Tutorials/Basics) scripting front-end via [libctl](https://libctl.readthedocs.io), or callable as a [C++](C++_Tutorial) library. An official [Python](Python_Tutorials/Basics) interface is under development.
-   Field output in the [HDF5](https://en.wikipedia.org/wiki/HDF5) standard scientific data format, supported by many visualization tools.
-   Arbitrary material and source distributions.
-   Field analyses including flux spectra, Maxwell stress tensor, frequency extraction, local density of states and energy integrals, near to far field transformations; completely programmable.
-   Multi-parameter optimization, root-finding, integration, etcetera via [libctl](https://libctl.readthedocs.io) and/or [NLopt](https://nlopt.readthedocs.io).

Time-Domain Simulation
----------------------

A time-domain electromagnetic simulation simply takes [Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell's_equations) and evolves them over time within some finite computational region, essentially performing a kind of **numerical experiment**. This can be used to calculate a wide variety of useful quantities, but major applications include:

-   **Transmission and Reflection Spectra** — by Fourier-transforming the response to a short pulse, a single simulation can yield the scattering amplitudes over a wide spectrum of frequencies.
-   **Resonant Modes and Frequencies** — by analyzing the response of the system to a short pulse, one can extract the frequencies, decay rates, and field patterns of the harmonic modes of lossy and lossless systems including waveguide and cavity modes.
-   **Field Patterns** (e.g. Green's functions) in response to an arbitrary source, archetypically a [CW](https://en.wikipedia.org/wiki/Continuous_wave) (fixed-$\omega$) input.

Using these results, one can then compute many other things, such as the local density of states from the trace of the Green's function. Meep's scriptable interface makes it possible to combine many sorts of computations along with multi-parameter optimization etcetera in sequence or in parallel.

[Tutorial/Basics](Scheme_Tutorials/Basics.md) gives examples of all of these kinds of computations.

Download
--------

The latest development sources are available on [GitHub](https://github.com/stevengj/meep). The source tarballs are available on the [Download](Download.md) page. The release history is described in the [Release Notes](Release_Notes.md). The installation instructions can be found in the [Installation](Installation.md) page.

Documentation
-------------

See the navigation sidebar at left. In particular, the [Introduction](Introduction.md) and [Tutorial](Scheme_Tutorials/Basics.md) are the most important things to read. We also have an [FAQ](FAQ.md).

Please [cite Meep](Acknowledgements.md#referencing) in any publication for which you found it useful.

### Mailing Lists

Subscribe to the read-only [meep-announce mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/meep-announce) to receive notifications of updates and releases. Subscribe to the [meep-discuss mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/meep-discuss) for discussions about using Meep. Archives are available [here](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/). You can also read and post to the list via the [gmane.comp.science.electromagnetism.meep.general](news://news.gmane.org/gmane.comp.science.electromagnetism.meep.general) newsgroup from [Gmane](http://www.gmane.org/).

### Bug Reports and Feature Requests

For bug reports and feature requests, please [file a Meep GitHub issue](https://github.com/stevengj/meep/issues).

Acknowledgements
----------------

The Meep project is maintained by [Simpetus](http://www.simpetus.com) and the open-source community on [GitHub](https://github.com/stevengj/meep). Please see the [Acknowledgements](Acknowledgements.md) for a more complete listing of those to whom we are grateful.

Contacts and Feedback
---------------------

If you have questions or problems regarding Meep, you are encouraged to query the [mailing list](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/).

Professional consulting services as well as free access to Meep in the public cloud via Amazon Web Services (AWS) are provided by [Simpetus](http://www.simpetus.com).
