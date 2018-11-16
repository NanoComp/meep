---
# FAQ
---

The following are frequently asked questions.

[TOC]

General
-------

### What is Meep?

Meep is a [free and open-source](https://en.wikipedia.org/wiki/Free_and_open-source_software) software package for [electromagnetics](https://en.wikipedia.org/wiki/Electromagnetism) simulation via the [finite-difference time-domain](https://en.wikipedia.org/wiki/Finite-difference_time-domain_method) (FDTD) method. The name Meep is an acronym for *MIT Electromagnetic Equation Propagation*.

### Who are the developers of Meep?

Meep was originally developed as part of graduate research at MIT. The project is now being maintained by [Simpetus](http://www.simpetus.com) and the open-source developer community on [GitHub](https://github.com/stevengj/meep).

### Where can I ask questions regarding Meep?

There is a public [mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/meep-discuss) for users to discuss issues pertaining to setting up simulations, post-processing output, installation, etc. A good place to start is the [list archives](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/) which includes all postings (6000+) since 2006 spanning a variety of topics. Bug reports and new feature requests should be filed as a [GitHub issue](https://github.com/stevengj/meep/issues).

### Are professional consulting services available?

Yes. [Simpetus](http://www.simpetus.com), a company started by Meep's developers and maintainers, provides professional consulting services for photonic design and modeling including development of turn-key simulation modules as well as training and technical support for getting up and running with Meep.

### How can I contribute to the Meep project?

[Pull requests](https://github.com/stevengj/meep/pulls) involving bug fixes, new features, and general improvements are welcome and can be made to the master branch on GitHub. This includes tweaks, revisions, and updates to this documentation which is also part of the [source repository](https://github.com/stevengj/meep/tree/master/doc).

### Is there a technical reference for Meep?

Yes. The technical details of Meep's inner workings are described in the peer-reviewed publication [MEEP: A flexible free-software package for electromagnetic simulations by the FDTD method](http://dx.doi.org/doi:10.1016/j.cpc.2009.11.008), Computer Physics Communications, Vol. 181, pp. 687-702, 2010 ([pdf](http://ab-initio.mit.edu/~oskooi/papers/Oskooi10.pdf)). Additional information is provided in the book [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707) in Chapters 4 ("Electromagnetic Wave Source Conditions"), 5 ("Rigorous PML Validation and a Corrected Unsplit PML for Anisotropic Dispersive Media"), 6 ("Accurate FDTD Simulation of Discontinuous Materials by Subpixel Smoothing"), and 20 ("MEEP: A Flexible Free FDTD Software Package"). A [video presentation](https://www.youtube.com/watch?v=9CA949csYvM) and [slides](http://ab-initio.mit.edu/~ardavan/stuff/IEEE_Photonics_Society_SCV3.pdf) as well as a [podcast](http://www.rce-cast.com/Podcast/rce-118-meep.html) are also available.

### Where can I find a list of projects which have used Meep?

For a list of more than 2500 published works which have used Meep, see the [Google Scholar citation page](https://scholar.google.com/scholar?hl=en&q=meep+software) as well as that for the [technical reference](https://scholar.google.com/scholar?cites=17712807607104508775) and also the [subpixel smoothing reference](https://scholar.google.com/scholar?cites=410731148689673259).

### Can I access Meep in the public cloud?

Yes. Meep is available preinstalled on Ubuntu on Amazon Web Services (AWS) Elastic Compute Cloud (EC2) as a free [Amazon Machine Image (AMI)](https://aws.amazon.com/marketplace/pp/B01KHWH0AS). To access this AMI, follow these [instructions](http://www.simpetus.com/launchsims.html).

Installation
------------

### Where can I install Meep?

Meep runs on any Unix-like operating system, such as Linux, macOS, and FreeBSD, from notebooks to desktops to supercomputers. [Conda packages](Installation.md#conda-packages) of the latest released version are available for Linux and macOS. There are also Conda packages of [nightly development builds](Installation.md#nightly-builds) which can be used to experiment with new features. Installing Meep from the source code requires some understanding of Unix, especially to install the various dependencies. Installation shell scripts are available for [Ubuntu 16.04 and 18.04](Build_From_Source.md#building-from-source) and [macOS Sierra](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/msg05811.html).

### Can I install Meep on Windows machines?

Yes. For Windows 10, you can install the [Ubuntu terminal](https://www.microsoft.com/en-us/p/ubuntu/9nblggh4msv6) as an app and then follow the instructions for [obtaining the Conda packages](Installation.md#conda-packages) or [building from source](Build_From_Source.md#building-from-source). For Windows 8 and older versions, you can use the free Unix-compatibility environment [Cygwin](http://www.cygwin.org/) following these [instructions](http://novelresearch.weebly.com/installing-meep-in-windows-8-via-cygwin.html).

### Are there precompiled binary packages for Ubuntu?

Yes. These packages can be obtained via the package manager [APT](https://en.wikipedia.org/wiki/APT_(Debian)) as described in [Download](Download.md#precompiled-packages-for-ubuntu). However, the current packages are for version 1.3 which is out of date. Up to date packages are being prepared. In the meantime, you can use the [Conda packages](Installation.md#conda-packages) which contain the official releases as well as nightly builds of the source repository.

### Guile is installed, but configure complains that it can't find `guile`

With most Linux distributions as well as Cygwin, packages like [Guile](http://www.gnu.org/software/guile) are split into two parts: a `guile` package that just contains the libraries and executables, and a `guile-dev` or `guile-devel` package that contains the header files and other things needed to compile programs using Guile. Usually, the former is installed by default but the latter is not. You need to install both, which means that you probably need to install `guile-dev`. Similarly for any other library packages needed by Meep.

Physics
-------

### How does the current amplitude relate to the resulting field amplitude?

There is no simple formula relating the input current amplitude (**J** in Maxwell's equations) to the resulting fields (**E**) etcetera, even at the same point as the current. The exact same current will produce a different field and radiate a different total power depending upon the surrounding materials/geometry, and depending on the frequency. This is a physical consequence of the geometry's effect on the local density of states; it can also be thought of as feedback from reflections on the source. As a simple example, if you put a current source inside a perfect electric conductor, the resulting field will be zero. As another example, the frequency-dependence of the radiated power in vacuum is part of the reason why the sky is blue.

See also Section 4.4 ("Currents and Fields: The Local Density of States") in [Chapter 4](http://arxiv.org/abs/arXiv:1301.5366) ("Electromagnetic Wave Source Conditions") of the book [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707).

If you are worried about this, then you are probably setting up your calculation in the wrong way. Especially in linear materials, the absolute magnitude of the field is useless; the only meaningful quantities are dimensionless ratios like the fractional transmittance: the transmitted power relative to the transmitted power in some reference calculation. Almost always, you want to perform two calculations, one of which is a reference, and compute the ratio of a result in one calculation to the result in the reference. For nonlinear calculations, see [Units and Nonlinearity](Units_and_Nonlinearity.md).

### How do I set the imaginary part of ε?

If you only care about the imaginary part of ε in a narrow bandwidth around some frequency ω, you should set it by using the electric [conductivity](Materials.md#conductivity-and-complex). If you care about the imaginary part of ε over a broad bandwidth, then for any physical material the imaginary part will be frequency-dependent and you will have to fit the data to a [Drude-Lorentz susceptibility model](Materials.md#material-dispersion).

Meep doesn't implement a frequency-independent complex ε. Not only is this not physical, but it also leads to both exponentially decaying and exponentially growing solutions in Maxwell's equations from positive- and negative-frequency Fourier components, respectively. Thus, it cannot be simulated in the time domain.

### Why does my simulation diverge if ε &lt; 0?

Maxwell's equations have exponentially growing solutions for a frequency-independent negative ε. For any physical medium with negative ε, there must be dispersion, and you must likewise use dispersive materials in Meep to obtain negative ε at some desired frequency. The requirement of dispersion to obtain negative ε follows from the [Kramers–Kronig relations](https://en.wikipedia.org/wiki/Kramers%E2%80%93Kronig_relations), and also follows from thermodynamic considerations that the energy in the electric field must be positive. For example, see [Electrodynamics of Continuous Media](https://www.amazon.com/Electrodynamics-Continuous-Media-Second-Theoretical/dp/0750626348) by Landau, Pitaevskii, and Lifshitz. At an even more fundamental level, it can be derived from passivity constraints as shown in [Physical Review A, Vol. 90, 023847, 2014](http://arxiv.org/abs/arXiv:1405.0238).

If you solve Maxwell's equations in a homogeneous-epsilon material at some real wavevector **k**, you get a dispersion relation $\omega^2 = c^2 |\mathbf{k}|^2 / \varepsilon$. If ε is positive, there are two real solutions $\omega = \pm c |\mathbf{k}| / \sqrt{\varepsilon}$, giving oscillating solutions. If ε is negative, there are two imaginary solutions corresponding to exponentially decaying and exponentially growing solutions from any current source. These solutions can always be spatially decomposed into a superposition of real-**k** values via a spatial Fourier transform.

If you do a simulation of any kind in the time domain (not just FDTD), you pretty much can't avoid exciting both the decaying and the growing solutions. This is *not* a numerical instability, it is a real solution of the underlying equations for an unphysical material.

See [Materials](Materials.md#material-dispersion) for how to include dispersive materials which can have negative ε and loss.

If you have negative ε *and* negative μ *everywhere*, the case of a negative-index material, then the simulation is fine. However at the boundary between negative- and positive-index materials, you will encounter instabilities: because of the way Maxwell's equations are discretized in FDTD, the ε and μ are discretized on different spatial grids, so you will get a half-pixel or so of εμ &lt; 0 at the boundary between negative and positive indices, which will cause the simulation to diverge. But of course, any physical negative-index metamaterial also involves dispersion.

Note also that, as a consequence of the above analysis, ε must go to a positive value in the ω $\to\pm\infty$ limit to get non-diverging solutions of Maxwell's equations. So the ε$_\infty$ in your [dispersion model](Materials.md#material-dispersion) must be positive.

### Why are there strange peaks in my reflectance/transmittance spectrum when modeling planar or periodic structures?

Firstly, the simulation run time may be too short and your results may have yet to be sufficiently [converged](#checking-convergence). Secondly, modeling flat/planar structures typically requires a 1d cell and periodic structures a single unit cell in 2d/3d. You may be using a higher-dimensional cell with multiple periods (a supercell) which introduces unwanted additional modes due to band folding. For more details, see Section 4.6 ("Sources in Supercells") in [Chapter 4](http://arxiv.org/abs/arXiv:1301.5366) ("Electromagnetic Wave Source Conditions") of the book [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707). Note that a 1d cell must be along the $z$ direction with only the E<sub>x</sub> and H<sub>y</sub> field components permitted.

Usage
-----

### Is there a Python interface?

Yes. An official Python interface was released in January 2018 with version 1.4. An unofficial [Python interface](https://www.fzu.cz/~dominecf/meep/), which is **not** compatible with the official version, has been developed independently by researchers at the Institute of Physics at the Czech Academy of Sciences and Ghent University, and maintained by [Filip Dominec](https://github.com/FilipDominec/python-meep-utils). Unfortunately, this interface has several shortcomings including missing support for geometric objects, lack of high-level abstractions for low-level functionality, and limited documentation. The official interface addresses all these issues.

### What are the different ways to define the material geometry?

There are four ways to define the material geometry: (1) the [`GeometricObject`](Python_User_Interface.md#geometricobject) (Python) or [`geometric-object`](Scheme_User_Interface.md#geometric-object) (Scheme) class used to specify a collection of predefined shapes including prisms, spheres, cylinders, cones, blocks, and ellipsoids, (2) `material_function` (Python) or `material-function` (Scheme) used to define any arbitrary function (i.e., for a given position in the computational cell, return the ε/μ at that point), (3) importing the scalar, real-valued, frequency-independent permittivity from an HDF5 file via the `epsilon_input_file` (Python) or `epsilon-input-file` (Scheme) input parameter, or (4) importing planar geometries from a [GDSII file](Python_User_Interface.md#gdsii-support). Combinations of (1) and (2) are allowed but not (3).

### Is there a materials library?

Yes. A materials library is available containing [crystalline silicon](https://en.wikipedia.org/wiki/Crystalline_silicon) (c-Si), [amorphous silicon](https://en.wikipedia.org/wiki/Amorphous_silicon) (a-Si) including also hydrogenated form, [silicon dioxide](https://en.wikipedia.org/wiki/Silicon_dioxide) (SiO<sub>2</sub>), [indium tin oxide](https://en.wikipedia.org/wiki/Indium_tin_oxide) (ITO), [alumina](https://en.wikipedia.org/wiki/Aluminium_oxide) (Al<sub>2</sub>O<sub>3</sub>), [gallium arsenide](https://en.wikipedia.org/wiki/Gallium_arsenide) (GaAs), [gallium nitride](https://en.wikipedia.org/wiki/Gallium_nitride) (GaN), [aluminum arsenide](https://en.wikipedia.org/wiki/Aluminium_arsenide) (AlAs), [aluminum nitride](https://en.wikipedia.org/wiki/Aluminium_nitride) (AlN), [borosilicate glass](https://en.wikipedia.org/wiki/Borosilicate_glass) (BK7), [fused quartz](https://en.wikipedia.org/wiki/Fused_quartz), [silicon nitride](https://en.wikipedia.org/wiki/Silicon_nitride) (Si<sub>3</sub>N<sub>4</sub>), [germanium](https://en.wikipedia.org/wiki/Germanium) (Ge), [indium phosphide](https://en.wikipedia.org/wiki/Indium_phosphide) (InP), [lithium niobate](https://en.wikipedia.org/wiki/Lithium_niobate) (LiNbO<sub>3</sub>), as well as 11 elemental metals: [silver](https://en.wikipedia.org/wiki/Silver) (Ag), [gold](https://en.wikipedia.org/wiki/Gold) (Au), [copper](https://en.wikipedia.org/wiki/Copper) (Cu), [aluminum](https://en.wikipedia.org/wiki/Aluminium) (Al), [berylium](https://en.wikipedia.org/wiki/Beryllium) (Be), [chromium](https://en.wikipedia.org/wiki/Chromium) (Cr), [nickel](https://en.wikipedia.org/wiki/Nickel) (Ni), [palladium](https://en.wikipedia.org/wiki/Palladium) (Pd), [platinum](https://en.wikipedia.org/wiki/Platinum) (Pt), [titanium](https://en.wikipedia.org/wiki/Titanium) (Ti), and [tungsten](https://en.wikipedia.org/wiki/Tungsten) (W). Additional information is provided in [Materials](Materials.md#materials-library).

### How do I import n and k values into Meep?

You can import any arbitrary complex permittivity profile via n and k values into Meep by fitting the wavelength- or frequency-dependent data to a sum of Drude-Lorentz polarizability terms as described in [Materials](Materials.md#material-dispersion). In general, you have to use nonlinear optimization to do the fit (e.g., to minimize the sum-of-squares errors or whatever error criterion you prefer). Enough Lorentzians should form a complete basis, so you should be able to fit any function given enough Lorentzians. A wavelength-dependent, purely-real permittivity (i.e., with no loss) which can be represented using the [Sellmeier equation](https://en.wikipedia.org/wiki/Sellmeier_equation) can be directly transferred to the Lorentz model using a simple substitution of variables. This is also described in [Materials](Materials.md#sellmeier-coefficients). Note that Meep only does subpixel averaging of the nondispersive part of ε (and μ).

### Why are the fields blowing up in my simulation?

Instability in the fields is likely due to one of three causes: (1) [PML](Python_User_Interface.md#pml) overlapping dispersive materials based on a [Drude-Lorentzian susceptibility](Python_User_Interface.md#lorentziansusceptibility) in the presence of [backward-wave modes](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.79.065601) (fix: replace the PML with an [Absorber](Python_User_Interface.md#absorber)), (2) the frequency of a Lorentzian susceptibility term is too high relative to the grid discretization (fix: increase the resolution and/or turn off subpixel smothing and/or reduce the Courant factor), or (3) a material with a [wavelength-independent negative real permittivity](#why-does-my-simulation-diverge-if-0) (fix: [fit the permittivity to a broadband Drude-Lorentzian susceptibility](#how-do-i-import-n-and-k-values-into-meep)).

### Does Meep support importing GDSII files?

Yes. The [`get_GDSII_prisms`](Python_User_Interface.md#gdsii-support) routine is used to import [GDSII](https://en.wikipedia.org/wiki/GDSII) files. See [Tutorial/GDSII Import](Python_Tutorials/GDSII_Import.md) for an example. This feature facilitates the simulation of 2d/planar structures which are fabricated using semiconductor foundries. Also, it enables Meep's plug-and-play capability with [electronic design automation](https://en.wikipedia.org/wiki/Electronic_design_automation) (EDA) circuit-layout editors (e.g., Cadence Virtuoso Layout, Silvaco Expert, KLayout, etc.). EDA is used for the synthesis and verification of large and complex integrated circuits.

### Checking convergence

In any computer simulation like Meep, you should check that your results are *converged* with respect to any approximation that was made.   There is no simple formula that will tell you in advance exactly how much resolution (etc.) is required for a given level of accuracy; the most reliable procedure is to simply double the resolution and verify that the answers you care about don't change to your desired tolerance.   Useful things to check (ideally by doubling) in this way are: **resolution**, **run time** (for Fourier spectra), **PML thickness**.   

Meep's [subpixel smoothing](Introduction.md#the-illusion-of-continuity) often improves the rate of convergence and makes convergence a smoother function of resolution.  However, subpixel smoothing does not occur for dispersive materials or user-defined material functions ε(x) (instead of the built-in geometric objects).

For flux calculations involving pulsed (i.e., Gaussian) sources, it is important to run the simulation long enough to ensure that all the fields have sufficiently decayed away (i.e., due to absorption by the PMLs, etc). Terminating the simulation prematurely will result in the Fourier-transformed fields, which are being accumulated during the time stepping (as explained in [Introduction](Introduction.md#transmittancereflectance-spectra)), to not be fully converged. Convergence of the fields is typically achieved by lowering the `decay_by` parameter in the `stop_when_fields_decayed` [run function](Python_User_Interface.md#run-functions).  Alternatively, you can explicitly set the run time to some numeric value that you repeatedly double, instead of using the field decay.  Sometimes it is also informative to double the `cutoff` parameter of sources to increase their smoothness (reducing the amplitude of long-lived high-frequency modes).

### Why are my reflectance/transmittance values less than zero and/or greater than one?

There are three possible explanations: (1) the normalization and the scattering runs are not comparable because e.g., the sources or flux monitors are not in the same position within the structure, (2) the [run time is not long enough](#checking-convergence) and hence all of the flux is not being collected in either or both runs, or (3) the flux is being computed at a frequency which is too far away from the center of the source bandwidth; in such cases the flux values are too small and may be dominated by rounding errors.

### Why doesn't turning off subpixel averaging work?

By default, when Meep assigns a dielectric constant ε or μ to each pixel, it uses a carefully designed average of the ε values within that pixel. This subpixel averaging generally improves the accuracy of the simulation &mdash; perhaps counter-intuitively, for geometries with discontinuous ε it is *more* accurate (i.e. closer to the exact Maxwell result for the *discontinuous* case) to do the simulation with the subpixel-averaged (*smoothed*) ε, as long as the averaging is done properly. For details, see Section 3 ("Interpolation and the illusion of continuity") of [Computer Physics Communications, Vol. 181, pp. 687-702, 2010](http://ab-initio.mit.edu/~oskooi/papers/Oskooi10.pdf).

Still, there are times when, for whatever reason, you might not want this feature. For example, if your accuracy is limited by other issues, or if you want to skip the wait at the beginning of the simulation for it do to the averaging. In this case, you can disable the subpixel averaging by setting `Simulation.eps_averaging = False` (Python) or `(set! eps-averaging? false)` (Scheme). For more details, see [Python User Interface](Python_User_Interface.md).

Even if you disable the subpixel averaging, however, when you output the dielectric function to a file and plot it, you may notice that there are some pixels with intermediate ε values, right at the boundary between two materials. This has a completely different source. Internally, Meep's simulation is performed on a [Yee grid](Yee_Lattice.md), in which every field component is stored on a slightly different grid which are offset from one another by half-pixels, and the ε values are also stored on this Yee grid. For output purposes, however, it is more user-friendly to output all fields etcetera on the same grid at the center of each pixel, so all quantities are interpolated onto this grid for output. Therefore, even though the internal ε values are indeed discontinuous when you disable subpixel averaging, the output file will still contain some "averaged" values at interfaces due to the interpolation from the Yee grid to the center-pixel grid.

### Can subpixel averaging be applied to dispersive materials?

Meep only does subpixel averaging of the *nondispersive* part of ε and μ. The dispersive part is not averaged at all.  This means that any sharp interfaces between dispersive materials will dominate the error, and you will probably get only first-order convergence, the same as if you do no subpixel averaging at all. It is possible that the subpixel averaging may still improve the constant factor in the convergence if not the asymptotic convergence rate, if you also have a lot of interfaces between nondispersive materials or if the dispersion is small (i.e., if ε is close to ε<sub>&#8734;</sub> over your bandwidth). On the other hand, if the dispersion is large and most of your interfaces are between large-dispersion materials, then subpixel averaging may not help at all and you might as well turn it off (which may improve [stability](#why-are-the-fields-blowing-up-in-my-simulation)). Generally, the subpixel averaging will not degrade accuracy though it will affect performance.

### How do I set up an oblique planewave source?

A planewave incident at any angle can be generated by typically setting the amplitude function of a 1d/line source (for a 2d computational cell) or 2d/planar source (for a 3d cell) and the Bloch-periodic boundary conditions. Tutorial examples are provided for [Python](Python_Tutorials/Basics.md#angular-reflectance-spectrum-of-a-planar-interface) and [Scheme](Scheme_Tutorials/Basics.md#angular-reflectance-spectrum-of-a-planar-interface). Additional examples are available for [Python](https://github.com/stevengj/meep/blob/master/python/examples/pw-source.py) and [Scheme](https://github.com/stevengj/meep/blob/master/scheme/examples/pw-source.ctl). See also [Tutorial/Mode Decomposition](Python_Tutorials/Mode_Decomposition.md#reflectance-and-transmittance-spectra-for-planewave-at-oblique-incidence). This topic is discussed on the [mailing list](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/msg00692.html). Note: the oblique planewave is incident at the given angle for only a *single* frequency component. For more details, refer to Section 4.5 ("Efficient Frequency-Angle Coverage") in [Chapter 4](https://arxiv.org/abs/1301.5366) ("Electromagnetic Wave Source Conditions") of the book [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707).

### What is Meep's frequency-domain solver and how does it work? 

Meep contains a [frequency-domain solver](Python_User_Interface.md#frequency-domain-solver) that directly computes the fields produced in a geometry in response to a [continuous-wave (CW) source](https://en.wikipedia.org/wiki/Continuous_wave), using an [iterative linear solver](https://en.wikipedia.org/wiki/Iterative_method) instead of time-stepping. This is possible because the FDTD timestep can be used to formulate a frequency-domain problem via an iterative linear solver. The frequency-domain response can often be determined using many fewer timesteps while exploiting the FDTD code almost without modification. For details, see Section 5.3 ("Frequency-domain solver") of [Computer Physics Communications, Vol. 181, pp. 687-702, 2010](http://ab-initio.mit.edu/~oskooi/papers/Oskooi10.pdf).

This means that all of the features from the time-domain solver (e.g., arbitrary materials, symmetries, subpixel averaging, parallelization, etc.) are also available as a frequency-domain solver. For certain problems, such as cavities (e.g., ring resonators) with long-lived resonant modes, the frequency-domain solver converges much faster than the straightforward approach of simply running a long simulation until transients have disappeared. Another benefit is that an arbitrary, complex, refractive index can be specified directly using the [electric conductivity](Materials.md#conductivity-and-complex) without having to fit the data to a sum of [Drude-Lorentz susceptibility terms](Materials.md#material-dispersion).

For examples, see [Tutorial/Frequency-Domain Solver](Python_Tutorials/Frequency_Domain_Solver.md) and [Tutorial/Mode Decomposition](Python_Tutorials/Mode_Decomposition.md#reflectance-and-transmittance-spectra-for-planewave-at-oblique-incidence).

### Should I expect linear [speedup](https://en.wikipedia.org/wiki/Speedup) from the parallel Meep?

For a given computational grid, Meep divides the grid points roughly equally among the processors,
and each process is responsible for all computations involving its "own" grid points (computing
ε from the materials, timestepping the fields, accumulating Fourier transforms, computing far fields, etcetera).
How much speedup this parallelization translates into depends on a number of factors, especially:

* The ratio of communications to computation, and the speed of your network. During timestepping, each processor needs to communicate neighboring grid points with other processors, and if you have too few grid points per processor (or your network is too slow) then the cost of this communication could overwhelm the computational gains.
* [Load balancing](https://en.wikipedia.org/wiki/Load_balancing_(computing)): different portions of the grid may be more expensive than other portions, causing processors in the latter portions to sit idle while a few processors work on the expensive regions. For example, setting up the materials at the beginning is more expensive in regions with lots of objects or interfaces. Timestepping is more expensive in regions with Fourier-transformed flux planes. Computing far fields only uses the processors where the corresponding near fields are located.
* If you write lots of fields to files, the parallel I/O speed (which depends on your network, filesystem, etc) may dominate.

In general, you will need large simulations to benefit from lots of processors. A rule of thumb is to keep doubling the number of processors until you no longer see much speedup.

### Why does the amplitude of my point dipole source increase with resolution?

The field from a point source is singular &mdash; it blows up as you approach the source. At any finite resolution, this singularity is truncated to a finite value by the discretization but the peak field at the source location increases as you increase the resolution.

### How does Meep deal with numerical dispersion?

Numerical dispersion can be analyzed and quantified analytically for a homogeneous medium. For details, see e.g., Chapter 4 ("Numerical Dispersion and Stability") of [Computational Electrodynamics: The Finite Difference Time-Domain Method (3rd edition)](https://www.amazon.com/Computational-Electrodynamics-Finite-Difference-Time-Domain-Method/dp/1580538320). However, in practice numerical dispersion is rarely the dominant source of error in FDTD calculations which almost always involve material inhomogeneities that give rise to much larger errors. Similar to other errors associated with the finite resolution, numerical dispersion decreases with resolution, so you can deal with it by increasing the resolution until convergence is obtained to the desired accuracy. In particular, the errors from numerical dispersion vary *quadratically* with resolution (in the ordinary center-difference FDTD scheme). On the other hand, the errors introduced by discretization of material interfaces go *linearly* with the resolution, so they are almost always dominant. Meep can partially correct for these errors using [subpixel averaging](Introduction.md#the-illusion-of-continuity).

### How do I compute S-parameters?

Meep contains a [mode-decomposition feature](Mode_Decomposition) which can be used to compute complex [S-parameters](https://en.wikipedia.org/wiki/Scattering_parameters). An example is provided for a two-port network in [Tutorial/Mode Decomposition](Python_Tutorials/Mode_Decomposition.md#reflectance-of-a-waveguide-taper).

### `Harminv` is unable to find the resonant modes of my structure

There are four possible explanations for why [`Harminv`](Python_User_Interface.md#harminv) could not find the resonant modes: (1) the run time was not long enough and the decay rate of the mode is so small that the `Harminv` data was mainly noise, (2) the `Harminv` call was not wrapped in [`after_sources`](Python_User_Interface.md#controlling-when-a-step-function-executes); if `Harminv` overlaps sources it will get confused because the sources are not exponentially decaying fields, (3) the `Harminv` monitor point is near where the mode has a node (e.g., in a symmetry plane), or (4) there are field instabilities where the fields are actually [blowing up](#why-are-the-fields-blowing-up-in-my-simulation) slowly; this may result in `Harminv` returning a *negative* quality factor. For calculations involving `Harminv`, it is always preferrable to run with a narrower bandwidth source around the frequency of interest and/or for a longer time.

Note: any real-valued signal consists of both positive and negative frequency components (with complex-conjugate amplitudes) in a Fourier domain decomposition into complex exponentials. `Harminv` usually is set up to find just one sign of the frequency, but occasionally converges to a negative-frequency component as well; these are just as meaningful as the positive frequencies.

### When outputting the dielectric function to a file, I don't see any dispersive materials

Only the real, frequency-independent part of ε/μ is written to an HDF5 file. As an example, many of the dispersive materials in the [materials library](Materials.md#materials-library) which have a broadband, complex, refractive index will appear as ε=1 in the output file. Thus, in order to verify the material geometry during debugging using visualization tools, etc., you may have to artificially adjust the `epsilon` value.

### Does Meep support a non-uniform grid?

No. Meep does not support grids with spatially varying resolution. One possible approach, which does not require changes to the underlying code and is not yet implemented, is to use a coordinate transformation to selectively increase the resolution in a given region of the computational cell. This is possible using transformation optics which involves a change of materials. For more details, see the notes [Coordinate Transformation and Invariance in Electromagnetism](http://math.mit.edu/~stevenj/18.369/coordinate-transform.pdf) and the [notes by Felix Schwarz on variable resolution in Meep](https://github.com/fesc3555/meep_variable_resolution) using this technique.

### How do I visualize the structure and fields in 3d?

You can use [Mayavi](http://docs.enthought.com/mayavi/mayavi/index.html). For an example, see [Tutorial/Basics](Python_Tutorials/Basics.md#visualizing-3d-structures).

### How do I compute the Poynting flux in an arbitrary direction?

The flux plane defined using [`FluxRegion`](Python_User_Interface.md#fluxregion) is always perpendicular to one of the coordinate axes. You can also integrate the Poynting vector over a volume instead of a surface. However, for a given `FluxRegion`, you can use its `direction` property to integrate a different component of the Poynting vector. Thus, in order to compute the Poynting flux in an arbitrary direction, you could specify three *overlapping* flux planes, each of which computes the integral of the Poynting vector in the *x*, *y*, and *z* directions.  Then a linear combination of these will give the Poynting flux in any desired direction.

### Can Meep be used to investigate lasing phenomena?

Yes. More specifically, Meep can be used to model saturable gain and absorption via multilevel atomic susceptibility. This feature may be used to investigate optically-pumped lasing phenomena such as [Raman lasers](https://en.wikipedia.org/wiki/Raman_laser). For details, see [Materials/Saturable Gain and Absorption](Materials.md#saturable-gain-and-absorption).

### Does Meep support adjoint-based optimization?

Not currently but work is underway to add support for this feature with expected release in early 2019.