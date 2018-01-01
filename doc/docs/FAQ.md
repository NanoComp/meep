---
# FAQ
---

The following are frequently asked questions about Meep.

[TOC]

General
-------

### What is Meep?

Meep is a free/open-source software package for finite-difference time-domain (FDTD) simulation. Meep is an acronym for MIT Electromagnetic Equation Propagation.

### Who are the developers of Meep?

Meep was originally developed as part of graduate research at MIT. The project is now being maintained by [Simpetus](http://www.simpetus.com) and the open-source community on [GitHub](https://github.com/stevengj/meep).

### Where can I ask questions regarding Meep?

There is a public [mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/meep-discuss) for users to discuss issues pertaining to setting up simulations, post-processing output, installation, etc. A good place to start is the [list archives](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/) which includes all postings since 2006 spanning a large number and variety of discussion topics. Bug reports and new feature requests should be filed as a [GitHub issue](https://github.com/stevengj/meep/issues).

[Simpetus](http://www.simpetus.com), a company started by Meep's developers, provides professional consulting services for photonic design and modeling including development of turn-key simulation modules as well as training and technical support for getting up and running with Meep.

Installation
------------

### Where can I install Meep?

Meep runs on any Unix-like operating system, such as Linux and macOS, from notebooks to desktops to high-performance computing (HPC) clusters. Conda packages are available for Linux and macOS. Meep can also be installed on Windows using the open-source [Cygwin](https://en.wikipedia.org/wiki/Cygwin) Unix-compatibility environment. See [Installation](Installation) for details.

Installing Meep from source code requires some understanding of Unix, especially to install the various prerequisites. Installation shell scripts are available for [Ubuntu 16.04](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/msg05884.html) and [macOS Sierra](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/msg05811.html).

Meep is also available preinstalled on Ubuntu on Amazon Web Services (AWS) Elastic Compute Cloud (EC2) as a free [Amazon Machine Image (AMI)](https://aws.amazon.com/marketplace/pp/B01KHWH0AS). To access this AMI, follow these [instructions](http://www.simpetus.com/launchsims.html).

### Guile is installed, but configure complains that it can't find `guile`

With most Linux distributions as well as Cygwin, packages like [Guile](http://www.gnu.org/software/guile) are split into two parts: a `guile` package that just contains the libraries and executables, and a `guile-dev` or `guile-devel` package that contains the header files and other things needed to compile programs using Guile. Usually, the former is installed by default but the latter is not. You need to install both, which means that you probably need to install `guile-dev`. Similarly for any other library packages needed by Meep.

Physics
-------

### How does the current amplitude relate to the resulting field amplitude?

There is no simple formula relating the input current amplitude (**J** in Maxwell's equations) to the resulting fields (**E**) etcetera, even at the same point as the current. The exact same current will produce a different field and radiate a different total power depending upon the surrounding materials/geometry, and depending on the frequency. This is a physical consequence of the geometry's effect on the local density of states; it can also be thought of as feedback from reflections on the source. As a simple example, if you put a current source inside a perfect electric conductor, the resulting field will be zero. As another example, the frequency-dependence of the radiated power in vacuum is part of the reason why the sky is blue.

See also this book chapter on [Wave Source Conditions](http://arxiv.org/abs/arXiv:1301.5366).

If you are worried about this, then you are probably setting up your calculation in the wrong way. Especially in linear materials, the absolute magnitude of the field is useless; the only meaningful quantities are dimensionless ratios like the fractional transmission: the transmitted power relative to the transmitted power in some reference calculation. Almost always, you want to perform two calculations, one of which is a reference, and compute the ratio of a result in one calculation to the result in the reference. For nonlinear calculations, see [Units and Nonlinearity](Units_and_Nonlinearity.md).

### How do I set the imaginary part of ε?

If you only care about the imaginary part of $\varepsilon$ in a narrow bandwidth around some frequency $\omega$, you should set it by using the electric conductivity as described in [Materials](Materials/#conductivity-and-complex). If you care about the imaginary part of $\varepsilon$ over a broad bandwidth, then for any physical material the imaginary part will be frequency-dependent and you will have to fit it to Meep's dispersive-$\varepsilon$ support. See [Materials](Materials#material-dispersion).

Meep doesn't implement a frequency-independent complex $\varepsilon$. Not only is this not physical, but it also leads to both exponentially decaying and exponentially growing solutions in Maxwell's equations from positive- and negative-frequency Fourier components, respectively. Thus, it cannot be simulated in the time domain.

### Why does my simulation diverge if ε &lt; 0?

Maxwell's equations have exponentially growing solutions for a frequency-independent negative $\varepsilon$. For any physical medium with negative $\varepsilon$, there must be dispersion, and you must likewise use dispersive materials in Meep to have negative $\varepsilon$ at some desired frequency. The requirement of dispersion to obtain negative $\varepsilon$ follows from the Kramers–Kronig relations, and also follows from thermodynamic considerations that the energy in the electric field must be positive; see, for example, *Electrodynamics of Continuous Media* by Landau, Pitaevskii, and Lifshitz. At an even more fundamental level, it can be derived from [passivity constraints](http://arxiv.org/abs/arXiv:1405.0238).

If you solve Maxwell's equations in a homogeneous-epsilon material at some real wavevector **k**, you get a dispersion relation $\omega^2 = c^2 |\mathbf{k}|^2 / \varepsilon$. If $\varepsilon$ is positive, there are two real solutions $\omega^2 = \pm c |\mathbf{k}| / \sqrt{\varepsilon}$, giving oscillating solutions. If $\varepsilon$ is negative, there are two *imaginary* solutions $\mu$, corresponding to exponentially decaying and *exponentially growing solutions* from any current source. These solutions can always be spatially decomposed into a superposition of real-**k** values via a spatial Fourier transform.

If you do a simulation of any kind in the time domain (not just FDTD), you pretty much can't avoid exciting both the decaying and the growing solutions. This is *not* a numerical instability, it is a real solution of the underlying equations for an unphysical material.

See [Materials](Materials/#material-dispersion) for how to include dispersive materials which can have negative $\varepsilon$ and loss.

If you have negative $\varepsilon$ *and* negative $\mu$ *everywhere*, the case of a negative-index material, then the simulation is fine. However at the boundary between negative- and positive-index materials, you will encounter instabilities: because of the way Maxwell's equations are discretized in FDTD, the $\varepsilon$ and $\mu$ are discretized on different spatial grids, so you will get a half-pixel or so of $\varepsilon\mu$ &lt; 0 at the boundary between negative and positive indices, which will cause the simulation to diverge. But of course, any physical negative-index metamaterial also involves dispersion.

Note also that, as a consequence of the above analysis, $\varepsilon$ must go to a positive value in the $\omega\to\pm\infty$ limit to get non-diverging solutions of Maxwell's equations. So the $\varepsilon_\infty$ in your [dispersion model](Materials.md) must be positive.

Usage
-----

### Is there a Python interface?

An official Python interface is currently under development and will likely be released by the end of 2017. An unofficial [Python interface](https://github.com/FilipDominec/python-meep-utils) has been independently developed by researchers at the Institute of Physics at the Czech Academy of Sciences and Ghent University. Unfortunately, this interface has a number of shortcomings including missing support for geometric objects, lack of high-level abstractions for low-level functionality, and limited documentation. These will all be addressed in the official interface.

### Why doesn't turning off subpixel averaging work?

By default, when Meep assigns a dielectric constant $\varepsilon$ or $\mu$ to each pixel, it uses a carefully designed average of the $\varepsilon$ values within that pixel. This subpixel averaging generally improves the accuracy of the simulation &mdash; perhaps counter-intuitively, for geometries with discontinous $\varepsilon$ it is *more* accurate (i.e. closer to the exact Maxwell result for the *discontinuous* case) to do the simulation with the subpixel-averaged (*smoothed*) $\varepsilon$, as long as the averaging is done properly. For details, see the [reference publication](Acknowledgements/#referencing).

Still, there are times when, for whatever reason, you might not want this feature. For example, if your accuracy is limited by other issues, or if you want to skip the wait at the beginning of the simulation for it do to the averaging. In this case, you can disable the subpixel averaging by setting `eps_averaging = False` via the `Simulation` instance (Python) or `(set! eps-averaging? false)` (Scheme) in your script file. See the [User Interface](Python_User_Interface.md).

Even if you disable the subpixel averaging, however, when you output the dielectric function to a file and plot it, you may notice that there are some pixels with intermediate $\varepsilon$ values, right at the boundary between two materials. This has a completely different source. Internally, Meep's simulation is performed on a [Yee grid](Yee_Lattice.md), in which every field component is stored on a slightly different grid which are offset from one another by half-pixels, and the $\varepsilon$ values are also stored on this Yee grid. For output purposes, however, it is more user-friendly to output all fields etcetera on the same grid at the center of each pixel, so all quantities are interpolated onto this grid for output. Therefore, even though the internal $\varepsilon$ values are indeed discontinuous when you disable subpixel averaging, the *output* file will still contain some "averaged" values at interfaces due to the interpolation from the Yee grid to the center-pixel grid.

### Does Meep support importing GDSII files?

Not currently, but work is underway to add support for this feature (expected release: 2018). This will facilitate the fabrication of simulated, [planar 2d designs](https://en.wikipedia.org/wiki/GDSII) using semiconductor foundries. Also, this feature will enable plug-and-play capability with [electronic design automation](https://en.wikipedia.org/wiki/Electronic_design_automation) (EDA) circuit-layout editors (e.g., Cadence Virtuoso Layout, Silvaco Expert, KLayout, etc.). EDA is used for the synthesis and verification of large and complex integrated circuits.
