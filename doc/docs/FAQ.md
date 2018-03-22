---
# FAQ
---

The following are frequently asked questions.

[TOC]

General
-------

### What is Meep?

Meep is a [free and open-source](https://en.wikipedia.org/wiki/Free_and_open-source_software) software package for simulating electromagnetic systems via the [finite-difference time-domain](https://en.wikipedia.org/wiki/Finite-difference_time-domain_method) (FDTD) method. Meep is an acronym for *MIT Electromagnetic Equation Propagation*.

### Who are the developers of Meep?

Meep was originally developed as part of graduate research at MIT. The project is now being maintained by [Simpetus](http://www.simpetus.com) and the open-source developer community on [GitHub](https://github.com/stevengj/meep).

### Where can I ask questions regarding Meep?

There is a public [mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/meep-discuss) for users to discuss issues pertaining to setting up simulations, post-processing output, installation, etc. A good place to start is the [list archives](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/) which includes all postings (6000+) since 2006 spanning a variety of discussion topics. Bug reports and new feature requests should be filed as a [GitHub issue](https://github.com/stevengj/meep/issues).

[Simpetus](http://www.simpetus.com), a company started by Meep's developers, provides professional consulting services for photonic design and modeling including development of turn-key simulation modules as well as training and technical support for getting up and running with Meep.

### How can I contribute to the Meep project?

[Pull requests](https://github.com/stevengj/meep/pulls) involving bug fixes, new features, and general improvements are welcome and can be made to the master branch on GitHub. This includes tweaks, revisions, and updates to this documentation which is also part of the [source repository](https://github.com/stevengj/meep/tree/master/doc).

### Is there a technical reference on Meep?

The technical details of Meep's inner workings are described in the peer-reviewed publication [MEEP: A flexible free-software package for electromagnetic simulations by the FDTD method](http://dx.doi.org/doi:10.1016/j.cpc.2009.11.008), Computer Physics Communications, Vol. 181, pp. 687-702, 2010 ([pdf](http://ab-initio.mit.edu/~oskooi/papers/Oskooi10.pdf)). Additional information is provided in the book [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707) in Chapters 4 ("Electromagnetic Wave Source Conditions"), 5 ("Rigorous PML Validation and a Corrected Unsplit PML for Anisotropic Dispersive Media"), 6 ("Accurate FDTD Simulation of Discontinuous Materials by Subpixel Smoothing"), and 20 ("MEEP: A Flexible Free FDTD Software Package"). [Lecture presentation](https://www.youtube.com/watch?v=9CA949csYvM) and [slides](http://ab-initio.mit.edu/~ardavan/stuff/IEEE_Photonics_Society_SCV3.pdf) are also available.

Installation
------------

### Where can I install Meep?

Meep runs on any Unix-like operating system, such as Linux and macOS, from notebooks to desktops to supercomputers. [Conda packages](Installation/#conda-packages) are available for Linux and macOS. Meep can also be installed on Windows using the open-source [Cygwin](https://en.wikipedia.org/wiki/Cygwin) Unix-compatibility environment. See [Installation](Installation) for details.

Installing Meep from source code requires some understanding of Unix, especially to install the various prerequisites. Installation shell scripts are available for [Ubuntu 16.04](http://ab-initio.mit.edu/~oskooi/meep_discuss/build_meep_python_mpi.sh) and [macOS Sierra](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/msg05811.html).

Meep is also available preinstalled on Ubuntu on Amazon Web Services (AWS) Elastic Compute Cloud (EC2) as a free [Amazon Machine Image (AMI)](https://aws.amazon.com/marketplace/pp/B01KHWH0AS). To access this AMI, follow these [instructions](http://www.simpetus.com/launchsims.html).

### Guile is installed, but configure complains that it can't find `guile`

With most Linux distributions as well as Cygwin, packages like [Guile](http://www.gnu.org/software/guile) are split into two parts: a `guile` package that just contains the libraries and executables, and a `guile-dev` or `guile-devel` package that contains the header files and other things needed to compile programs using Guile. Usually, the former is installed by default but the latter is not. You need to install both, which means that you probably need to install `guile-dev`. Similarly for any other library packages needed by Meep.

Physics
-------

### How does the current amplitude relate to the resulting field amplitude?

There is no simple formula relating the input current amplitude (**J** in Maxwell's equations) to the resulting fields (**E**) etcetera, even at the same point as the current. The exact same current will produce a different field and radiate a different total power depending upon the surrounding materials/geometry, and depending on the frequency. This is a physical consequence of the geometry's effect on the local density of states; it can also be thought of as feedback from reflections on the source. As a simple example, if you put a current source inside a perfect electric conductor, the resulting field will be zero. As another example, the frequency-dependence of the radiated power in vacuum is part of the reason why the sky is blue.

See also Section 4.4 ("Currents and Fields: The Local Density of States") in [Chapter 4](http://arxiv.org/abs/arXiv:1301.5366) ("Electromagnetic Wave Source Conditions") of the book [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707).

If you are worried about this, then you are probably setting up your calculation in the wrong way. Especially in linear materials, the absolute magnitude of the field is useless; the only meaningful quantities are dimensionless ratios like the fractional transmission: the transmitted power relative to the transmitted power in some reference calculation. Almost always, you want to perform two calculations, one of which is a reference, and compute the ratio of a result in one calculation to the result in the reference. For nonlinear calculations, see [Units and Nonlinearity](Units_and_Nonlinearity.md).

### How do I set the imaginary part of ε?

If you only care about the imaginary part of $\varepsilon$ in a narrow bandwidth around some frequency $\omega$, you should set it by using the electric [conductivity](Materials/#conductivity-and-complex). If you care about the imaginary part of $\varepsilon$ over a broad bandwidth, then for any physical material the imaginary part will be frequency-dependent and you will have to fit the data to a [Drude-Lorentzian susceptibility profile](Materials#material-dispersion).

Meep doesn't implement a frequency-independent complex $\varepsilon$. Not only is this not physical, but it also leads to both exponentially decaying and exponentially growing solutions in Maxwell's equations from positive- and negative-frequency Fourier components, respectively. Thus, it cannot be simulated in the time domain.

### Why does my simulation diverge if ε &lt; 0?

Maxwell's equations have exponentially growing solutions for a frequency-independent negative $\varepsilon$. For any physical medium with negative $\varepsilon$, there must be dispersion, and you must likewise use dispersive materials in Meep to have negative $\varepsilon$ at some desired frequency. The requirement of dispersion to obtain negative $\varepsilon$ follows from the Kramers–Kronig relations, and also follows from thermodynamic considerations that the energy in the electric field must be positive; see, for example, the book [Electrodynamics of Continuous Media](https://www.amazon.com/Electrodynamics-Continuous-Media-Second-Theoretical/dp/0750626348) by Landau, Pitaevskii, and Lifshitz. At an even more fundamental level, it can be derived from [passivity constraints](http://arxiv.org/abs/arXiv:1405.0238).

If you solve Maxwell's equations in a homogeneous-epsilon material at some real wavevector **k**, you get a dispersion relation $\omega^2 = c^2 |\mathbf{k}|^2 / \varepsilon$. If $\varepsilon$ is positive, there are two real solutions $\omega^2 = \pm c |\mathbf{k}| / \sqrt{\varepsilon}$, giving oscillating solutions. If $\varepsilon$ is negative, there are two *imaginary* solutions $\mu$, corresponding to exponentially decaying and *exponentially growing solutions* from any current source. These solutions can always be spatially decomposed into a superposition of real-**k** values via a spatial Fourier transform.

If you do a simulation of any kind in the time domain (not just FDTD), you pretty much can't avoid exciting both the decaying and the growing solutions. This is *not* a numerical instability, it is a real solution of the underlying equations for an unphysical material.

See [Materials](Materials/#material-dispersion) for how to include dispersive materials which can have negative $\varepsilon$ and loss.

If you have negative $\varepsilon$ *and* negative $\mu$ *everywhere*, the case of a negative-index material, then the simulation is fine. However at the boundary between negative- and positive-index materials, you will encounter instabilities: because of the way Maxwell's equations are discretized in FDTD, the $\varepsilon$ and $\mu$ are discretized on different spatial grids, so you will get a half-pixel or so of $\varepsilon\mu$ &lt; 0 at the boundary between negative and positive indices, which will cause the simulation to diverge. But of course, any physical negative-index metamaterial also involves dispersion.

Note also that, as a consequence of the above analysis, $\varepsilon$ must go to a positive value in the $\omega\to\pm\infty$ limit to get non-diverging solutions of Maxwell's equations. So the $\varepsilon_\infty$ in your [dispersion model](Materials.md) must be positive.

### Why are there strange peaks in my reflection/transmission spectrum when modeling planar or periodic structures?

Modeling flat/planar structures typically requires a 1d computational cell and periodic structures a single unit cell. You may be using a higher-dimensional cell with multiple periods (a supercell) which introduces unwanted additional modes due to band folding. For more details, see Section 4.6 ("Sources in Supercells") in [Chapter 4](http://arxiv.org/abs/arXiv:1301.5366) ("Electromagnetic Wave Source Conditions") of the book [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707). Note that a 1d cell must be along the $z$ direction with only the $E_x$ and $H_y$ field components permitted.

Usage
-----

### Is there a Python interface?

An official Python interface was released in January 2018 with version 1.4. An unofficial [Python interface](https://github.com/FilipDominec/python-meep-utils) has been developed independently by researchers at the Institute of Physics at the Czech Academy of Sciences and Ghent University. Unfortunately, this interface has a number of shortcomings including missing support for geometric objects, lack of high-level abstractions for low-level functionality, and limited documentation. The official interface addresses all these issues.

### What are the different ways to define the material geometry?

There are currently three ways to define the material geometry via: (1) the [`GeometricObject`](Python_User_Interface/#geometricobject) (Python) or [`geometric-object`](Scheme_User_Interface/#geometric-object) (Scheme) class used to specify a collection of shapes including spheres, cylinders, cones, blocks, and ellipsoids, (2) `material_function` (Python) or `material-function` (Scheme) used to define an arbitrary function, or (3) importing the scalar, real-valued, frequency-independent permittivity from an HDF5 file via the `epsilon_input_file` (Python) or `epsilon-input-file` (Scheme) input parameter. Combinations of (1) and (2) are allowed but not (3).

### Does Meep support importing GDSII files?

Not currently, but work is underway to add support for this feature with expected release in mid 2018. Importing [GDSII](https://en.wikipedia.org/wiki/GDSII) files will facilitate the simulation of 2d/planar structures which are fabricated using semiconductor foundries. Also, this feature will enable Meep's plug-and-play capability with [electronic design automation](https://en.wikipedia.org/wiki/Electronic_design_automation) (EDA) circuit-layout editors (e.g., Cadence Virtuoso Layout, Silvaco Expert, KLayout, etc.). EDA is used for the synthesis and verification of large and complex integrated circuits.

### Why doesn't turning off subpixel averaging work?

By default, when Meep assigns a dielectric constant $\varepsilon$ or $\mu$ to each pixel, it uses a carefully designed average of the $\varepsilon$ values within that pixel. This subpixel averaging generally improves the accuracy of the simulation &mdash; perhaps counter-intuitively, for geometries with discontinous $\varepsilon$ it is *more* accurate (i.e. closer to the exact Maxwell result for the *discontinuous* case) to do the simulation with the subpixel-averaged (*smoothed*) $\varepsilon$, as long as the averaging is done properly. For details, see Section 3 ("Interpolation and the illusion of continuity") of [Computer Physics Communications, Vol. 181, pp. 687-702, 2010](http://ab-initio.mit.edu/~oskooi/papers/Oskooi10.pdf).

Still, there are times when, for whatever reason, you might not want this feature. For example, if your accuracy is limited by other issues, or if you want to skip the wait at the beginning of the simulation for it do to the averaging. In this case, you can disable the subpixel averaging by setting `Simulation.eps_averaging = False` (Python) or `(set! eps-averaging? false)` (Scheme). See the [User Interface](Python_User_Interface.md).

Even if you disable the subpixel averaging, however, when you output the dielectric function to a file and plot it, you may notice that there are some pixels with intermediate $\varepsilon$ values, right at the boundary between two materials. This has a completely different source. Internally, Meep's simulation is performed on a [Yee grid](Yee_Lattice.md), in which every field component is stored on a slightly different grid which are offset from one another by half-pixels, and the $\varepsilon$ values are also stored on this Yee grid. For output purposes, however, it is more user-friendly to output all fields etcetera on the same grid at the center of each pixel, so all quantities are interpolated onto this grid for output. Therefore, even though the internal $\varepsilon$ values are indeed discontinuous when you disable subpixel averaging, the *output* file will still contain some "averaged" values at interfaces due to the interpolation from the Yee grid to the center-pixel grid.

### How to set up an oblique planewave source?

A planewave incident at any angle can be generated by setting the amplitude function of a 1d/line source (for a 2d computational cell) or 2d/planar source (for a 3d cell). This is discussed on the [mailing list](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/msg00692.html). Examples are provided in [Python](https://github.com/stevengj/meep/blob/master/python/examples/pw-source.py) and [Scheme](https://github.com/stevengj/meep/blob/master/examples/pw-source.ctl). Note: the oblique planewave is incident at the given angle for only a single frequency component. For accuracy involving broadband calculations, this will typically require splitting up the spectrum into subintervals (requiring multilple simulations) and recombining the results in post processing. For more details, refer to Section 4.5 ("Efficiency Frequency-Angle Coverage") in [Chapter 4](https://arxiv.org/abs/1301.5366) ("Electromagnetic Wave Source Conditions") of the book [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707).

### What is Meep's frequency-domain solver and how does it work? 

Meep contains a [frequency-domain solver](Python_User_Interface/#frequency-domain-solver) that directly computes the fields produced in a geometry in response to a [continuous-wave (CW) source](https://en.wikipedia.org/wiki/Continuous_wave), using an [iterative linear solver](https://en.wikipedia.org/wiki/Iterative_method) instead of time-stepping. This is possible because the FDTD timestep can be used to directly plug a frequency-domain problem into an iterative linear solver. The frequency-domain response can often be determined using many fewer timesteps while exploiting the FDTD code almost without modification. For details, see Section 5.3 ("Frequency-domain solver") of [Computer Physics Communications, Vol. 181, pp. 687-702, 2010](http://ab-initio.mit.edu/~oskooi/papers/Oskooi10.pdf).

This means that all of the features from the time-domain solver (e.g., arbitrary materials, symmetries, subpixel averaging, parallelization, etc.) are also available as a frequency-domain solver. For certain problems, such as cavities (e.g., ring resonators) with long-lived resonant modes, the frequency-domain solver converges much faster than the straightforward approach of simply running a long simulation until transients have disappeared. Another benefit is that an arbitrary complex refractive index can be specified directly using the [electric conductivity](Materials/#conductivity-and-complex) without having to fit the data to a sum of [Lorentzian-Drude susceptibility terms](Materials/#material-dispersion).

Examples are provided in [Tutorials/Frequency-Domain Solver](Python_Tutorials/Frequency_Domain_Solver/).

### Is there a materials library?

A materials library is available containing 11 commonly used metals in optoelectronic devices: Ag, Au, Cu, Al, Be, Cr, Ni, Pd, Pt, Ti, W. Additional information is provided in [Materials](Materials/#materials-library).

### Should I expect linear [speedup](https://en.wikipedia.org/wiki/Speedup) from the parallel Meep?

For a given computational grid, Meep divides the grid points roughly equally among the processors,
and each process is responsible for all computations involving its "own" grid points (computing
ε from the materials, timestepping the fields, accumulating Fourier transforms, computing far fields, etcetera).
How much speedup this parallelization translates into depends on a number of factors, especially:

* The ratio of communications to computation, and the speed of your network.  During timestepping, each processor needs to communicate neighboring grid points with other processors, and if you have too few grid points per processor (or your network is too slow) then the cost of this communication could overwhelm the computational gains.
* [Load balancing](https://en.wikipedia.org/wiki/Load_balancing_(computing)): different portions of the grid may be more expensive than other portions, causing processors in the latter portions to sit idle while a few processors work on the expensive regions.   (For example, setting up the materials at the beginning is more expensive in regions with lots of objects or interfaces.  Timestepping is more expensive in regions with Fourier-transformed flux planes.   Computing far fields only uses the processors where the corresponding near fields are located.)
* If you write lots of fields to files, the parallel I/O speed (which depends on your network, filesystem, etc) may dominate.

In general, you will need large simulations to benefit from lots of processors.   A rule of thumb is to keep doubling the number of processors until you no longer see much speedup.
