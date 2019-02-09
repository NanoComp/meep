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

Meep was originally developed as part of graduate research at MIT. The project is now being maintained by [Simpetus](http://www.simpetus.com) and the developer community on [GitHub](https://github.com/NanoComp/meep).

### Where can I ask questions regarding Meep?

There is a public [mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/meep-discuss) for users to discuss issues pertaining to setting up simulations, post-processing output, installation, etc. A useful place to start is the [list archives](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/) which includes all postings (6000+) since 2006 spanning a variety of topics. Bug reports and new feature requests should be filed as a [GitHub issue](https://github.com/NanoComp/meep/issues).

### Are professional consulting services available?

Yes. [Simpetus](http://www.simpetus.com), a company started by Meep's developers and maintainers, provides professional consulting services for photonic design and modeling including development of turn-key simulation modules as well as training and technical support for getting up and running with Meep.

### How can I contribute to the Meep project?

[Pull requests](https://github.com/NanoComp/meep/pulls) involving bug fixes, new features, and general improvements are welcome and can be made to the master branch on GitHub. This includes tweaks, revisions, and updates to this documentation, generated from [markdown](https://en.wikipedia.org/wiki/Markdown), which is also part of the [source repository](https://github.com/NanoComp/meep/tree/master/doc).

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

Yes. For Windows 10, you can install the [Ubuntu terminal](https://www.microsoft.com/en-us/p/ubuntu/9nblggh4msv6) as an app which is based on the [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about) framework and then follow the instructions for [obtaining the Conda packages](Installation.md#conda-packages) or [building from source](Build_From_Source.md#building-from-source). For Windows 8 and older versions, you can use the free Unix-compatibility environment [Cygwin](http://www.cygwin.org/) following these [instructions](http://novelresearch.weebly.com/installing-meep-in-windows-8-via-cygwin.html).

### Are there precompiled binary packages for Ubuntu?

Yes. Ubuntu and Debian packages can be obtained via the package manager [APT](https://en.wikipedia.org/wiki/APT_(Debian)) as described in [Download](Download.md#precompiled-packages-for-ubuntu). However, the Meep packages for Ubuntu 16.04 ([serial](https://packages.ubuntu.com/xenial/meep) and [parallel](https://packages.ubuntu.com/xenial/meep-openmpi)) and 18.04 ([serial](https://packages.ubuntu.com/bionic/meep) and [parallel](https://packages.ubuntu.com/bionic/meep-openmpi)) are for [version 1.3](https://github.com/NanoComp/meep/releases) (March 2015) which is out of date. The Meep package for Ubuntu is in the process of being updated and will likely appear in Ubuntu 19.10 as derived from the [unstable Debian package](https://packages.debian.org/unstable/meep). In the meantime, since the [Scheme interface](Scheme_User_Interface.md) is no longer being supported and has been replaced by the [Python interface](Python_User_Interface.md), you can use the [Conda packages](Installation.md#conda-packages) which contain the official releases as well as nightly builds of the master branch of the source repository.

### Guile is installed, but configure complains that it can't find `guile`

With most Linux distributions as well as Cygwin, packages like [Guile](http://www.gnu.org/software/guile) are split into two parts: a `guile` package that just contains the libraries and executables, and a `guile-dev` or `guile-devel` package that contains the header files and other things needed to compile programs using Guile. Usually, the former is installed by default but the latter is not. You need to install both, which means that you probably need to install `guile-dev`. Similarly for any other library packages needed by Meep.

Physics
-------

### How does the current amplitude relate to the resulting field amplitude?

There is no simple formula relating the input current amplitude (**J** in Maxwell's equations) to the resulting fields (**E**) etcetera, even at the same point as the current. The exact same current will produce a different field and radiate a different total power depending upon the surrounding materials/geometry, and depending on the frequency. This is a physical consequence of the geometry's effect on the local density of states; it can also be thought of as feedback from reflections on the source. A classic example is an antenna in front of a ground plane, which radiates very different amounts of power depending on the distance between the antenna and the plane (half wavelength vs. quarter wavelength, for example). Alternatively, if you put a current source inside a perfect electric conductor, the resulting field will be zero. Also, as the frequency of the current increases, the amplitude of the resulting field will also increase. This is related to why the sky is blue: scattered power increases with frequency (alternatively the density of states increases as the frequency to the d-1 power in d dimensions).

For a mathematical description, see Section 4.4 ("Currents and Fields: The Local Density of States") in [Chapter 4](http://arxiv.org/abs/arXiv:1301.5366) ("Electromagnetic Wave Source Conditions") of [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707).

If you are worried about this, then you are probably setting up your calculation in the wrong way. Especially in linear materials, the absolute magnitude of the field is useless; the only meaningful quantities are dimensionless ratios like the fractional transmittance: the transmitted power relative to the transmitted power in some reference calculation. Almost always, you want to perform two calculations, one of which is a reference, and compute the ratio of a result in one calculation to the result in the reference. For nonlinear calculations, see [Units and Nonlinearity](Units_and_Nonlinearity.md).

### How is the source current defined?

The source current in Meep is defined as a [free charge current **J** in Maxwell's equations](Introduction.md#maxwells-equations). Meep does not simulate the driving force behind this free charge current, nor does the current have to be placed in a conductor. Specifying a current means that somehow you are shaking a [charge](https://en.wikipedia.org/wiki/Electric_charge) at that point (by whatever means, Meep doesn't care) and you want to know the [resulting fields](#how-does-the-current-amplitude-relate-to-the-resulting-field-amplitude).

In the [interface](Python_User_Interface.md#source), the source currents are labeled E<sub>x</sub> or H<sub>y</sub> etc. according to what components of the electric/magnetic fields they correspond to.

### How do I compute the local density of states (LDOS) in a lossy material?

If you put a point source *inside* a lossy material (e.g., a [Lorentz-Drude metal](Materials.md#material-dispersion)), then the power expended by a dipole diverges as you increase the resolution. In the limit of infinite resolution, infinite power is absorbed. LDOS is not well defined for points inside of a lossy material.

However, LDOS is perfectly well defined for points *outside* of a lossy material. For example, you can choose a point outside of a lossy object and calculate the [LDOS](Python_User_Interface.md#ldos-spectra), and it will converge to a finite value as you increase the resolution. For an example, see [Tutorial/Local Density of States](../Python_Tutorials/Local_Density_of_States/).

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

There are two possible explanations: (1) the simulation run time may be too short and your results have not sufficiently [converged](#checking-convergence), or (2) you may be using a higher-dimensional cell with multiple periods (a supercell) which introduces unwanted additional modes due to band folding. Modeling flat/planar structures typically requires a 1d cell and periodic structures a single unit cell in 2d/3d. For more details, see Section 4.6 ("Sources in Supercells") in [Chapter 4](http://arxiv.org/abs/arXiv:1301.5366) ("Electromagnetic Wave Source Conditions") of [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707). Note that a 1d cell must be along the $z$ direction with only the E<sub>x</sub> and H<sub>y</sub> field components permitted.

### How do I model the solar radiation spectrum?

For simulations involving [solar radiation](https://en.wikipedia.org/wiki/Sunlight#Surface_illumination), including the [air mass](https://en.wikipedia.org/wiki/Air_mass_(solar_energy)), the [reflectance/transmittance spectra](Introduction.md#transmittancereflectance-spectra) is computed as normal. Since typical solar-cell problems are linear, the reflected or transmitted power can then be obtained by simply multiplying the reflectance or transmittance by the solar spectrum.

Usage
-----

### Is there a Python interface?

Yes. An official [Python interface](Python_User_Interface.md) was released in [version 1.4](https://github.com/NanoComp/meep/releases) and replaces the [Scheme interface](Scheme_User_Interface.md) which is no longer being supported. An unofficial [Python interface](https://www.fzu.cz/~dominecf/meep/), which predates and is **incompatible** with the official version, has been developed independently by researchers at the Institute of Physics at the Czech Academy of Sciences and Ghent University, and maintained by [Filip Dominec](https://github.com/FilipDominec/python-meep-utils). Unfortunately, this interface has several shortcomings including missing support for geometric objects, lack of high-level abstractions for low-level functionality, and limited documentation. The official interface addresses all these issues.

### What are the different ways to define a structure?

There are four ways to define a structure: (1) the [`GeometricObject`](Python_User_Interface.md#geometricobject) (Python) or [`geometric-object`](Scheme_User_Interface.md#geometric-object) (Scheme) class used to specify a collection of predefined shapes including `Prism`, `Sphere`, `Cylinder`, `Cone`, `Block`, and `Ellipsoid`, (2) `material_function` (Python) or `material-function` (Scheme) used to define an arbitrary function: for a given position in the cell, return the ε/μ at that point, (3) import the scalar, real-valued, frequency-independent permittivity from an HDF5 file via the `epsilon_input_file` (Python) or `epsilon-input-file` (Scheme) input parameter, or (4) import planar geometries from a [GDSII file](Python_User_Interface.md#gdsii-support). Combinations of (1), (2), and (4) are allowed but not (3).

### Is there a materials library?

Yes. A materials library is available containing [crystalline silicon](https://en.wikipedia.org/wiki/Crystalline_silicon) (c-Si), [amorphous silicon](https://en.wikipedia.org/wiki/Amorphous_silicon) (a-Si) including the hydrogenated form, [silicon dioxide](https://en.wikipedia.org/wiki/Silicon_dioxide) (SiO<sub>2</sub>), [indium tin oxide](https://en.wikipedia.org/wiki/Indium_tin_oxide) (ITO), [alumina](https://en.wikipedia.org/wiki/Aluminium_oxide) (Al<sub>2</sub>O<sub>3</sub>), [gallium arsenide](https://en.wikipedia.org/wiki/Gallium_arsenide) (GaAs), [gallium nitride](https://en.wikipedia.org/wiki/Gallium_nitride) (GaN), [aluminum arsenide](https://en.wikipedia.org/wiki/Aluminium_arsenide) (AlAs), [aluminum nitride](https://en.wikipedia.org/wiki/Aluminium_nitride) (AlN), [borosilicate glass](https://en.wikipedia.org/wiki/Borosilicate_glass) (BK7), [fused quartz](https://en.wikipedia.org/wiki/Fused_quartz), [silicon nitride](https://en.wikipedia.org/wiki/Silicon_nitride) (Si<sub>3</sub>N<sub>4</sub>), [germanium](https://en.wikipedia.org/wiki/Germanium) (Ge), [indium phosphide](https://en.wikipedia.org/wiki/Indium_phosphide) (InP), [lithium niobate](https://en.wikipedia.org/wiki/Lithium_niobate) (LiNbO<sub>3</sub>), as well as 11 elemental metals: [silver](https://en.wikipedia.org/wiki/Silver) (Ag), [gold](https://en.wikipedia.org/wiki/Gold) (Au), [copper](https://en.wikipedia.org/wiki/Copper) (Cu), [aluminum](https://en.wikipedia.org/wiki/Aluminium) (Al), [berylium](https://en.wikipedia.org/wiki/Beryllium) (Be), [chromium](https://en.wikipedia.org/wiki/Chromium) (Cr), [nickel](https://en.wikipedia.org/wiki/Nickel) (Ni), [palladium](https://en.wikipedia.org/wiki/Palladium) (Pd), [platinum](https://en.wikipedia.org/wiki/Platinum) (Pt), [titanium](https://en.wikipedia.org/wiki/Titanium) (Ti), and [tungsten](https://en.wikipedia.org/wiki/Tungsten) (W). Additional information is provided in [Materials](Materials.md#materials-library).

### How do I import n and k values into Meep?

You can import any arbitrary complex permittivity profile via n and k values into Meep by fitting the wavelength- or frequency-dependent data to a sum of Drude-Lorentz polarizability terms as described in [Materials](Materials.md#material-dispersion). In general, you have to use nonlinear optimization to do the fit (e.g., to minimize the sum-of-squares errors or whatever error criterion you prefer). Enough Lorentzians should form a complete basis, so you should be able to fit any function given enough Lorentzians. A wavelength-dependent, purely-real permittivity (i.e., with no loss) which can be represented using the [Sellmeier equation](https://en.wikipedia.org/wiki/Sellmeier_equation) can be directly [transferred to the Lorentz model using a simple substitution of variables](Materials.md#sellmeier-coefficients). Note that Meep only does [subpixel averaging of the nondispersive part of ε (and μ)](#can-subpixel-averaging-be-applied-to-dispersive-materials).

### Does Meep support gyrotropic materials?

No. Currently, Meep only supports anisotropic, real-symmetric, permittivity tensors. In the [magneto-optic effect](https://en.wikipedia.org/wiki/Magneto-optic_effect), an external magnetic field yields imaginary off-diagonal components of ε (with no absorption) which is not yet supported (issue [#60](https://github.com/NanoComp/meep/issues/60)).

### Why are the fields blowing up in my simulation?

Instability in the fields is likely due to one of four causes: (1) [PML](Python_User_Interface.md#pml) overlapping dispersive materials based on a [Drude-Lorentzian susceptibility](Python_User_Interface.md#lorentziansusceptibility) in the presence of [backward-wave modes](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.79.065601) (fix: replace the PML with an [Absorber](Python_User_Interface.md#absorber)), (2) the frequency of a Lorentzian susceptibility term is *too high* relative to the grid discretization (fix: increase the resolution and/or reduce the Courant factor), (3) a material with a [wavelength-independent negative real permittivity](#why-does-my-simulation-diverge-if-0) (fix: [fit the permittivity to a broadband Drude-Lorentzian susceptibility](#how-do-i-import-n-and-k-values-into-meep)), or (4) a grid voxel contains more than one dielectric interface (fix: turn off subpixel averaging).

### Does Meep support importing GDSII files?

Yes. The [`get_GDSII_prisms`](Python_User_Interface.md#gdsii-support) routine is used to import [GDSII](https://en.wikipedia.org/wiki/GDSII) files. See [Tutorial/GDSII Import](Python_Tutorials/GDSII_Import.md) for an example. This feature facilitates the simulation of 2d/planar structures which are fabricated using semiconductor foundries. Also, it enables Meep's plug-and-play capability with [electronic design automation](https://en.wikipedia.org/wiki/Electronic_design_automation) (EDA) circuit-layout editors (e.g., Cadence Virtuoso Layout, Silvaco Expert, KLayout, etc.). EDA is used for the synthesis and verification of large and complex integrated circuits.

### Checking convergence

In any computer simulation like Meep, you should check that your results are *converged* with respect to any approximation that was made. There is no simple formula that will tell you in advance exactly how much resolution (etc.) is required for a given level of accuracy; the most reliable procedure is to simply double the resolution and verify that the answers you care about don't change to your desired tolerance. Useful things to check (ideally by doubling) in this way are: **resolution**, **run time** (for Fourier spectra), **PML thickness**.

Meep's [subpixel smoothing](Introduction.md#the-illusion-of-continuity) often improves the rate of convergence and makes convergence a smoother function of resolution. However, subpixel smoothing does not occur for [dispersive materials](#can-subpixel-averaging-be-applied-to-dispersive-materials) or [user-defined material functions](#why-does-subpixel-averaging-take-so-long) ε(x) instead of the built-in geometric objects.

For flux calculations involving pulsed (i.e., Gaussian) sources, it is important to run the simulation long enough to ensure that all the transient fields have sufficiently decayed away (i.e., due to absorption by the PMLs, etc). Terminating the simulation prematurely will result in the Fourier-transformed fields, which are being accumulated during the time stepping (as explained in [Introduction](Introduction.md#transmittancereflectance-spectra)), to not be fully converged. Convergence of the fields is typically achieved by lowering the `decay_by` parameter in the `stop_when_fields_decayed` [run function](Python_User_Interface.md#run-functions).  Alternatively, you can explicitly set the run time to some numeric value that you repeatedly double, instead of using the field decay.  Sometimes it is also informative to double the `cutoff` parameter of sources to increase their smoothness (reducing the amplitude of long-lived high-frequency modes).

### Why are my reflectance/transmittance values less than zero and/or greater than one?

There are four possible explanations: (1) the normalization and the scattering runs are not comparable because e.g., the sources or flux monitors are not in the same position within the structure, (2) the [run time is not long enough](#checking-convergence) and hence all of the flux is not being collected in either or both runs, (3) the flux is being computed at a frequency which is too far away from the center of the source bandwidth; in such cases the flux values are too small and may be dominated by rounding errors, or (4) the source or flux monitor is positioned too close to the scatterer which [modifies the local density of states (LDOS)](#how-does-the-current-amplitude-relate-to-the-resulting-field-amplitude); for example, a source emits more power near a band edge and less power within a bandgap than the same source surrounded by many wavelengths of vacuum.

### How do I compute the modes of a non-orthogonal (i.e., triangular) lattice?

Meep does not support non-rectangular unit cells. To deal with a triangular lattice, you have to use a supercell. This will cause the band structure to be [folded](#why-are-there-strange-peaks-in-my-reflectancetransmittance-spectrum-when-modeling-planar-or-periodic-structures).  However, if you take your point source and replicate it according to the underlying triangular lattice vectors, with the right phase relationship according to the Bloch wavevector, then it should excite the folded bands only with very low amplitude as reported by [`Harminv`](Python_User_Interface.md#harminv). Also, for every `Harminv` point you put in, you should analyze the fields from the periodic copies of that point (with the periodicity of the underlying lattice). Then, reject any frequency that is not detected at *all* points, with an amplitude that is related by something close to the correct exp(ikx) phase.

In principle, the excitation of the folded bands would be exactly zero if you place your sources correctly in the supercell. However, this doesn't happen in FDTD because the finite grid spoils the symmetry slightly. It also means that the detection of folded bands will vary with resolution. 

For an example, see Section 4.6 ("Sources in Supercells") in [Chapter 4](http://arxiv.org/abs/arXiv:1301.5366) ("Electromagnetic Wave Source Conditions") of [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707).

### Why doesn't turning off subpixel averaging work?

By default, when Meep assigns a dielectric constant ε or μ to each pixel, it uses a carefully designed average of the ε values within that pixel. This subpixel averaging generally improves the accuracy of the simulation &mdash; perhaps counter-intuitively, for geometries with discontinuous ε it is *more* accurate (i.e. closer to the exact Maxwell result for the *discontinuous* case) to do the simulation with the subpixel-averaged (*smoothed*) ε, as long as the averaging is done properly. For details, see Section 3 ("Interpolation and the illusion of continuity") of [Computer Physics Communications, Vol. 181, pp. 687-702, 2010](http://ab-initio.mit.edu/~oskooi/papers/Oskooi10.pdf).

Still, there are times when, for whatever reason, you might not want this feature. For example, if your accuracy is limited by other issues, or if you want to skip the wait at the beginning of the simulation for it do to the averaging. In this case, you can disable the subpixel averaging by setting `Simulation.eps_averaging = False` (Python) or `(set! eps-averaging? false)` (Scheme). For more details, see [Python User Interface](Python_User_Interface.md).

Even if you disable the subpixel averaging, however, when you output the dielectric function to a file and visualize it, you may notice that there are some pixels with intermediate ε values, right at the boundary between two materials. This is due to a completely different reason. Internally, Meep's simulation is performed on a [Yee grid](Yee_Lattice.md), in which every field component is stored on a slightly different grid which are offset from one another by half-pixels, and the ε values are also stored on this Yee grid. For output purposes, however, it is more user-friendly to output all fields etcetera on the same grid at the center of each pixel, so all quantities are interpolated onto this grid for output. Therefore, even though the internal ε values are indeed discontinuous when you disable subpixel averaging, the output file will still contain some "averaged" values at interfaces due to the interpolation from the Yee grid to the center-pixel grid.   For the same reason, if `k_point` is set and the boundaries are Bloch-periodic, the permittivity function of the entire cell obtained via `get_epsilon` or `output_epsilon` will show that a little of the cell from one edge "leaks" over to the other edge. This is independent of PML and the way the structure is defined (i.e., using geometric objects or a material function, etc.). An example is shown in the figure below comparing `output_epsilon` for two cases involving with and without `k_point`. The discretization artifacts are highlighted.

<center>
![](images/output_epsilon_kpoint.png)
</center>

### Can subpixel averaging be applied to dispersive materials?

No. Meep only does subpixel averaging of the non-dispersive part of ε and μ. The dispersive part is not averaged at all.  This means that any sharp interfaces between dispersive materials will dominate the error, and you will probably get only first-order convergence, the same as if you do no subpixel averaging at all. It is possible that the subpixel averaging may still improve the constant factor in the convergence if not the asymptotic convergence rate, if you also have a lot of interfaces between non-dispersive materials or if the dispersion is small (i.e., if ε is close to ε<sub>&#8734;</sub> over your bandwidth). On the other hand, if the dispersion is large and most of your interfaces are between large-dispersion materials, then subpixel averaging may not help at all and you might as well turn it off (which may improve [stability](#why-are-the-fields-blowing-up-in-my-simulation)). Generally, the subpixel averaging will not degrade accuracy though it will affect performance.

### What is a good rule of thumb for the grid resolution?

At least 8 pixels per wavelength in the lossless dielectric material with the highest index. Resolving the [skin depth of metals](https://en.wikipedia.org/wiki/Skin_effect), which is typically tens of nanometers at optical frequencies, will require a pixel size of comparable dimensions since [subpixel averaging does not apply to dispersive materials](#can-subpixel-averaging-be-applied-to-dispersive-materials).

### Why does subpixel averaging take so long?

There are at least two possible reasons due to using: (1) a `material_function` to define a [`Medium`](Python_User_Interface.md#medium) object or (2) the [C++](C++_Tutorial) interface. Unlike either the [Python](Python_User_Interface/) or [Scheme](Scheme_User_Interface/) interfaces which are based on analytically computing the averaged permittivity for boundary voxels involving at most one [`GeometricObject`](Python_User_Interface.md#geometricobject) (e.g., `Sphere`, `Prism`, `Block`, etc.), the C++ interface computes these averages from the `material_function` using [numerical quadrature](https://en.wikipedia.org/wiki/Numerical_integration) if the parameter `use_anisotropic_averaging=true` is passed to the constructor of `set_materials`. This procedure involves calling the `material_function` many times for every voxel in the [structure object](C++_Developer_Information.md#data-structures-and-chunks) which can be slow due to the [SWIG](http://www.swig.org/) callbacks, particularly because the voxel density is repeatedly doubled until a given threshold tolerance (`subpixel_tol`) or maximum iteration number (`subpixel_maxeval`) is reached. Because of this discrepancy in the subpixel averaging, the results for the C++ and Python/Scheme interfaces may be slightly different at the same resolution. You can potentially speed up subpixel averaging by increasing `subpixel_tol` or decreasing `subpixel_maxeval`. Note that the slow callbacks may still be noticeable during the grid initialization even when subpixel averaging is turned off. Just remember that if you turn off subpixel averaging, it usually means that you may need to increase the grid resolution to obtain the same accuracy. You will have to determine how much accuracy you want to trade for time. Alternatively, in the C++ interface you can use the [`meepgeom.hpp`](https://github.com/NanoComp/meep/blob/master/src/meepgeom.hpp) routines to define your geometry in terms of blocks, cylinders, etcetera similar to Python and Scheme, with semi-analytical subpixel averaging.

### For calculations involving Fourier-transformed fields, why should the source be a pulse rather than a continuous wave?

A continuous-wave source ([ContinuousSource](Python_User_Interface.md#continuoussource)) produces fields which are not integrable: their Fourier transform will not converge as the run time of the simulation is increased because the source never terminates. The Fourier-transformed fields are therefore arbitrarily defined by the run time. This [windowing](https://en.wikipedia.org/wiki/Window_function) does different things to the normalization and scattering runs because the spectra are different in the two cases. In contrast, a pulsed source ([GaussianSource](Python_User_Interface.md#gaussiansource)) produces fields which are [L2](https://en.wikipedia.org/wiki/Norm_(mathematics)#Euclidean_norm)-integrable: their Fourier transform is well defined and convergent as long as the run time is sufficiently large and the [fields have decayed away](#checking-convergence).

### How do I compute the group velocity of a mode?

There are two possible approaches for computing the group velocity: (1) compute the [band diagram](https://en.wikipedia.org/wiki/Band_diagram) ω(**k**) using [`Harminv`](Python_User_Interface.md#harminv), fit it to a polynomial, and calculate its derivative using a [finite difference](https://en.wikipedia.org/wiki/Finite_difference), or (2) excite the mode using a narrowband pulse and compute the ratio of the flux to energy density.

### How do I set up an oblique planewave source?

A planewave incident at any angle can be generated by typically setting the amplitude function [`amp_func`](Python_User_Interface.md#source) of a 1d/line source for a 2d cell or 2d/planar source for a 3d cell, as well as the Bloch-periodic boundary conditions via `k_point`. For a 1d example, see [Tutorial/Basics](Python_Tutorials/Basics.md#angular-reflectance-spectrum-of-a-planar-interface). There is also a [Scheme version](Scheme_Tutorials/Basics.md#angular-reflectance-spectrum-of-a-planar-interface). For 2d, see [Tutorial/Mode Decomposition](Python_Tutorials/Mode_Decomposition.md#reflectance-and-transmittance-spectra-for-planewave-at-oblique-incidence). Additional examples are available for [Python](https://github.com/NanoComp/meep/blob/master/python/examples/pw-source.py) and [Scheme](https://github.com/NanoComp/meep/blob/master/scheme/examples/pw-source.ctl).

The oblique planewave is incident at a given angle for only a *single* frequency component of the source. Alternatively, each frequency of a pulsed source corresponds to a *different* angle. This is a fundamental feature of FDTD simulations and not of Meep per se. Thus, to simulate an incident planewave at multiple angles for a given frequency ω, you will need to do separate simulations involving different values of k (`k_point`) since each set of (k,ω) specifying the Bloch-periodic boundaries and the frequency of the source will produce a different angle of the planewave.

For more details, refer to Section 4.5 ("Efficient Frequency-Angle Coverage") in [Chapter 4](https://arxiv.org/abs/1301.5366) ("Electromagnetic Wave Source Conditions") of [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707).

### What is Meep's frequency-domain solver and how does it work? 

Meep contains a [frequency-domain solver](Python_User_Interface.md#frequency-domain-solver) that directly computes the steady-state fields produced in a geometry in response to a [continuous-wave (CW) source](https://en.wikipedia.org/wiki/Continuous_wave), using an [iterative linear solver](https://en.wikipedia.org/wiki/Iterative_method) instead of time-stepping. This is possible because the FDTD timestep can be used to formulate a frequency-domain problem via an iterative linear solver. The frequency-domain response can often be determined using many fewer timesteps while exploiting the FDTD code almost without modification. For details, see Section 5.3 ("Frequency-domain solver") of [Computer Physics Communications, Vol. 181, pp. 687-702, 2010](http://ab-initio.mit.edu/~oskooi/papers/Oskooi10.pdf).

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

Unless the computational parallelism outweighs the extra communications overhead, the parallel program will actually be *slower* than the serial one.  This means, for example, that even if you really have two or more physical processors you won't be able to benefit from parallelization until the problem is sufficiently large. In general, you will need large simulations to benefit from lots of processors. A rule of thumb is to keep doubling the number of processors until you no longer see much speedup.

### Why are simulations involving Fourier-transformed fields slow?

The [discrete time Fourier transform](https://en.wikipedia.org/wiki/Discrete-time_Fourier_transform) (DTFT) of the fields, which is necessary for computing the [Poynting flux](Python_User_Interface.md#flux-spectra), [local density of states](Python_User_Interface.md#ldos-spectra) (LDOS), [near to far field transformation](Python_User_Interface.md#near-to-far-field-spectra), etc., is accumulated at every time step for every point in the [flux region](Python_User_Interface.md#fluxregion). The DTFT computation is parallelized but only in the sense that each processor computes the DTFT fields at points in its own [chunk](Chunks_and_Symmetry.md) of the grid. If the division of the grid among processors into approximately equal-sized chunks allocates most of the points where the DTFT fields are computed to one processor, it is *not* going to parallelize.

Ideally, the parallelization should take the DTFT computation into account: the grid should be divided such that the flux region, if it is sufficiently expensive, is divided among the processors. The best approach would be some kind of adaptive [load-balancing](https://en.wikipedia.org/wiki/Load_balancing_(computing)). However, nothing like this is currently implemented. Also, it is not possible in Meep to customize how the cell is divided among the processes.

A simple approach to reduce the cost of the DTFT computation is to reduce the number of frequency points. If you need high frequency resolution in a certain bandwidth, consider adding a second flux region just for that bandwidth, with as many points as you need there, and use a smaller number of frequency points over a broad bandwidth.

Another approach might be to change your structure so that the flux region is closer to the center of the grid, which will increase the likelihood that the parallelization will divide the grid in such a way as to split the flux region between different processors.

### Does Meep support shared-memory parallelism?

Meep can run in parallel on a shared-memory machine using MPI. However, it doesn't yet take special advantage of shared memory using [multithreading](https://en.wikipedia.org/wiki/Thread_(computing)#Multithreading) (issue [#228](https://github.com/NanoComp/meep/issues/228)).

### Why does the amplitude of a dipole point source increase with resolution?

The field from a point source is singular &mdash; it blows up as you approach the source. At any finite resolution, this singularity is truncated to a finite value by the discretization but the peak field at the source location increases as you increase the resolution.

### What normalization convention does meep use for the fields of eigenmode sources?

An [eigenmode source](Python_User_Interface.md#eigenmodesource)
is a localized distribution of electric and magnetic currents 
$\mathbf{J}(\mathbf{x}), \mathbf{M}(\mathbf{x})$, 
confined to a cross-sectional line (2D) or plane (3D) through a
waveguide or similar geometry, with the property that the
spatial distribution of the
electric and magnetic fields radiated by the sources
exactly reproduces the spatial distribution of
one specific (quasi-)normal mode of the
geometry (and is thus orthogonal to all other modes).
This condition leaves unspecified an overall scale
factor---if a field distribution
$\{\mathbf{E}(\mathbf x), \mathbf{H}(\mathbf x)\}$
satisfies the single-mode condition,
then so does the scaled distribution
$\{\lambda \mathbf{E}(\mathbf x), \lambda \mathbf{H}(\mathbf x)\}$
for any arbitrary $\lambda.$ To pin down this
ambiguity, MEEP chooses the overall amplitude of eigenmode
sources to ensure that the total power flux
carried by the fields they produce---the integral
of the normal Poynting vector over the cross-sectional
line or plane---evaluates numerically to 1.

More specifically,
in a time-harmonic problem in which all sources and fields
have time dependence $e^{-i\omega_m t}$ (where $\omega$

In practice, this has the following ramifications for MEEP
calculations using eigenmode sources:

+ For [frequency-domain calculations](Python_User_Interface.md#frequency-domain-solver)
  involving an eigenmode source with a `ContinuousSrc` time envelope---corresponding
  to monochromatic fields and sources that oscillate forever at the eigenfrequency
  of the mode---the total flux 

+ On the other hand, for the typical case of *time-domain* calculations,
  in which the spatial current distributions
  $\mathbf{J}(\mathbf{x}), \mathbf{M}(\mathbf x)}$ of the eigenmode source
  are paired with a temporal envelope function $W(t)$---in practice
  either a 
  [gaussian](Python_User_Interface.md#gaussiansource) or a
  [user-specified custom envelope](Python_User_Interface.md#customsource)---the
  values reported by MEEP for all frequency-domain field amplitudes at
  frequency $\omega$ will include factors of $\widetwiddle{W}(\omega)$
  (the Fourier transform of the temporal envelope), while
  while field-bilinear quantities like Poynting vectors and power fluxes
  will include factors of $|\widetwiddle{W}(\omega)|^2.$

### How does Meep deal with numerical dispersion?

Numerical dispersion can be analyzed and quantified analytically for a homogeneous medium. For details, see e.g., Chapter 4 ("Numerical Dispersion and Stability") of [Computational Electrodynamics: The Finite Difference Time-Domain Method (3rd edition)](https://www.amazon.com/Computational-Electrodynamics-Finite-Difference-Time-Domain-Method/dp/1580538320). However, in practice numerical dispersion is rarely the dominant source of error in FDTD calculations which almost always involve material inhomogeneities that give rise to much larger errors. Similar to other errors associated with the finite grid resolution, numerical dispersion decreases with resolution, so you can deal with it by increasing the resolution until convergence is obtained to the desired accuracy. In particular, the errors from numerical dispersion vary *quadratically* with resolution (in the ordinary center-difference FDTD scheme). On the other hand, the errors introduced by discretization of material interfaces go *linearly* with the resolution, so they are almost always dominant. Meep can partially correct for these errors using [subpixel averaging](Introduction.md#the-illusion-of-continuity).

### How do I compute S-parameters?

Meep contains a [mode-decomposition feature](Mode_Decomposition.md) which can be used to compute complex-valued [S-parameters](https://en.wikipedia.org/wiki/Scattering_parameters). An example is provided for a [two-port network](https://en.wikipedia.org/wiki/Two-port_network#Scattering_parameters_(S-parameters)) based on a silicon directional coupler in [Tutorial/GDSII Import](/Python_Tutorials/GDSII_Import/). An additional example is available for a [waveguide mode converter](Python_Tutorials/Mode_Decomposition.md#reflectance-of-a-waveguide-taper).

### `Harminv` is unable to find the resonant modes of my structure

There are five possible explanations for why [`Harminv`](Python_User_Interface.md#harminv) could not find the resonant modes: (1) the run time was not long enough and the decay rate of the mode is so small that the `Harminv` data was mainly noise, (2) the `Harminv` call was not wrapped in [`after_sources`](Python_User_Interface.md#controlling-when-a-step-function-executes); if `Harminv` overlaps sources turning on and off it will get confused because the sources are not exponentially decaying fields, (3) the `Harminv` monitor is near the mode's nodal point (e.g., in a symmetry plane), (4) there are field instabilities where the fields are actually [blowing up](#why-are-the-fields-blowing-up-in-my-simulation); this may result in `Harminv` returning a negative [quality factor](https://en.wikipedia.org/wiki/Q_factor), or (5) the decay rate of the mode is *too high*; `Harminv` discards any modes which have a quality factor less than 50 (and thus `Harminv` won't find any modes inside [metals](Materials.md#material-dispersion)).

In order to resolve two closely-spaced modes, in general it is preferable to run with a narrow bandwidth source around the frequency of interest to excite/analyze as few modes as possible and/or increase the run time to improve the frequency resolution.

Note: any real-valued signal consists of both positive and negative frequency components (with complex-conjugate amplitudes) in a Fourier domain decomposition into complex exponentials. `Harminv` usually is set up to find just one sign of the frequency, but occasionally converges to a negative-frequency component as well; these are just as meaningful as the positive frequencies.

### How does `k_point` define the phase relation between adjacent unit cells?

If you set the `k_point` to any `meep.Vector3`, the structure will be periodic in **all** directions. (There is a lower-level `field::set_boundary` function that allows you to set individual boundary conditions independently, however.)

The value of the `k_point` determines the phase relation between the fields and sources in adjacent periodic unit cells. In general, if you have period (`Lx`,`Ly`) and you are looking at the (`n`,`m`) unit cell it has a phase of exp(2πi * (`kx` * `Lx` * `n` + `ky` * `Ly` * `m`)). For example, if you set the `k_point` to `meep.Vector3(0,0,0)`, that means the fields/sources are periodic: the phase is unity from one cell to the next. If you set the `k_point` to `meep.Vector3(1,0,0)` it means that there is a phase difference of exp(2πi * `Lx`) between adjacent cells in the *x* direction. This is known as a [Bloch wave](https://en.wikipedia.org/wiki/Bloch_wave).

Note: in any cell direction where there is a [PML](Perfectly_Matched_Layer.md), the boundary conditions are mostly irrelevant. For example, if there is a PML in front of a periodic boundary, the periodicity doesn't matter because the field will have decayed almost to zero by the time it "wraps around" to the other side of the cell.

### Should I include the 2π factor when defining the frequency or the wavevector?

No. Frequency inputs and outputs in Meep are the ordinary frequency `f`, not the angular frequency ω=2πf. Similarly, spatial wavevectors k (e.g. for Bloch-periodic boundary conditions) are specified without the 2π factor, so that the spatial dependence is exp(2πikx).
For example, if you specify a `frequency=0.3` in a source, then the time-dependence of the source is exp(-2πi0.3t), where time t is also in Meep units. Similarly, if you specify `k_point = meep.Vector3(0.4,0,0)` in the interface, then the phase factor between adjacent unit cells with period L in the x direction is exp(2πi0.4L).

### Can Meep simulate time-varying structures?

Yes. The most general method is to re-initialize the material at every timestep by calling `field::set_materials` or `set_materials_from_geometry` in C++, or `simulation.set_materials` in Python. However, this is potentially quite slow. One alternative is a function [`field::phase_in_material`](Python_User_Interface.md#field-computations) that allows you to linearly interpolate between two precomputed structures, gradually transitioning over a given time period; we hope to have a more general version of this functionality in the future (issue [#207](https://github.com/NanoComp/meep/issues/207)).

### How do I model a moving point charge?

You can use an instantaneous [`ContinuousSource`](Python_User_Interface.md#continuoussource) with large wavelength (or nearly-zero frequency). This is analogous to a [direct current](https://en.wikipedia.org/wiki/Direct_current). You will also need to create a [run function](Python_User_Interface.md#run-functions) which contains [`change_sources`](Python_User_Interface.md#reloading-parameters) and specify the `center` property of the point source to be time dependent. As an example, the following image demonstrates [Cherenkov radiation](https://en.wikipedia.org/wiki/Cherenkov_radiation) involving a moving point charge with [superluminal phase velocity](https://en.wikipedia.org/wiki/Faster-than-light#Phase_velocities_above_c) (see [examples/cherenkov-radiation.py](https://github.com/NanoComp/meep/blob/master/python/examples/cherenkov-radiation.py)).

<center>
![](images/cherenkov_radiation.png)
</center>

### When outputting the dielectric function to a file, I don't see any dispersive materials

Only the real, frequency-independent part of ε/μ is written to an HDF5 file. As an example, many of the dispersive materials in the [materials library](Materials.md#materials-library) which have a broadband, complex, refractive index will appear as ε=1 in the output file. Thus, in order to verify the material geometry during debugging using visualization tools, etc., you may have to artificially adjust the `epsilon` value.

### Does Meep support grids with non-uniform discretization?

No. Meep does not support non-orthogonal grids with spatially varying resolution. One possible approach, which does not require changes to the underlying code and is not yet implemented, is to use a coordinate transformation to selectively increase the resolution in a given region of the cell. This is possible using transformation optics which involves a change of materials: an arbitrary coordinate transformation can be mapped to Cartesian coordinates with transformed ε/μ. For more details, see the notes [Coordinate Transformation and Invariance in Electromagnetism](http://math.mit.edu/~stevenj/18.369/coordinate-transform.pdf) and [Variable Resolution in Meep](https://github.com/fesc3555/meep_variable_resolution) using this technique.

### Can Meep model electrostatic effects?

In principle, this corresponds to the limit as the frequency goes to zero or the wavelength goes to infinity.  However, a time-domain simulation is rather inefficient for such [electrostatic](https://en.wikipedia.org/wiki/Electrostatics) (or magnetostatic) calculation; this includes [lumped circuit models](https://en.wikipedia.org/wiki/Lumped_element_model) such as resistance, voltage, capacitance, etc. In this regime, you are usually much better off directly solving e.g. [Poisson's equation](https://en.wikipedia.org/wiki/Poisson%27s_equation#Electrostatics) to obtain the fields from a given charge distribution. There are many available Poisson solvers based on [finite](https://en.wikipedia.org/wiki/Finite_element_method) or [boundary](https://en.wikipedia.org/wiki/Boundary_element_method) element methods.  In Meep, probably the best you can do is to use a source with a very low frequency and a gradual turn-on specified by the `width` parameter of [`ContinuousSrc`](Python_User_Interface.md#continuoussource).

### How do I access or visualize the structure, fields, or sources in a subregion of the cell?

You can use the routines [`get_array`](Python_User_Interface.md#array-slices), `get_dft_array`, or [`get_source_slice`](Python_User_Interface.md#source-slices) to obtain the fields/sources and [`get_array_metadata`](Python_User_Interface.md#array-metadata) or `get_dft_array_metadata` to obtain information for the geometric slice. Visualization in 3d can be done with [Mayavi](http://docs.enthought.com/mayavi/mayavi/index.html). For an example, see [Tutorial/Basics](Python_Tutorials/Basics.md#visualizing-3d-structures).

To output the data to an HDF5 file, you can use the [`in_volume`](Python_User_Interface.md#modifying-hdf5-output) or `in_point` routines as part of your [run function](../Python_User_Interface/#run-functions). For example, to restrict the output to a line, you could use: `meep.in_volume(meep.Volume(center=meep.Vector3(0,0,0), size=meep.Vector3(10,0,0)), meep.output_dpwr)` which outputs ε|E|<sup>2</sup> along a line of length 10 in the x direction centered at (0,0,0). You can even wrap this statement in `to_appended("line.h5", ...)` to output the intensity along the line as a function of time to a 2d HDF5 dataset. This would enable you to plot intensity vs. time and space as a 2d color image.

### Can Meep be used to investigate lasing phenomena?

Yes. More specifically, Meep can be used to model saturable gain and absorption via multilevel atomic susceptibility. This feature may be used to investigate optically-pumped lasing phenomena such as [Raman lasers](https://en.wikipedia.org/wiki/Raman_laser). For details, see [Materials/Saturable Gain and Absorption](Materials.md#saturable-gain-and-absorption).

### Does Meep support adjoint-based optimization?

Not currently but work is underway to add support for this feature with expected release in early 2019 (issue [#600](https://github.com/NanoComp/meep/pull/600)).
