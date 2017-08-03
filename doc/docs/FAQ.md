---
# FAQ
---

The following are frequently asked questions about Meep.

General
-------

### What is Meep?

Meep is a free-software package for finite-difference time-domain (FDTD) simulation. Meep is an acronym for MIT Electromagnetic Equation Propagation.

### Who are the developers of Meep?

Meep was originally developed as part of graduate research at MIT. The project is now being maintained by [Simpetus](http://www.simpetuscloud.com) and the open-source community on [GitHub](https://github.com/stevengj/meep).

Installation
------------

### Where can I install Meep?

Meep runs on any Unix-like system, such as Linux and macOS, from notebooks to desktops to high-performance computing (HPC) clusters. Precompiled packages are available for Ubuntu. Meep can also be installed on Windows using the free [Cygwin](https://en.wikipedia.org/wiki/Cygwin) Unix-compatibility environment.

Installing Meep from source code requires some understanding of Unix, especially to install the various prerequisites. Step-by-step instructions are available for [Ubuntu 16.04](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/msg05819.html) and [macOS Sierra](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/msg05811.html).

Meep is also available preinstalled on Ubuntu on Amazon Web Services (AWS) Elastic Compute Cloud (EC2) as a free [Amazon Machine Image (AMI)](https://aws.amazon.com/marketplace/pp/B01KHWH0AS). To access this AMI, follow these [instructions](http://www.simpetuscloud.com/launchsims.html).

### Guile is installed, but configure complains that it can't find `guile`

With most Linux distributions as well as Cygwin, packages like [Guile](http://www.gnu.org/software/guile) are split into two parts: a `guile` package that just contains the libraries and executables, and a `guile-dev` or `guile-devel` package that contains the header files and other things needed to compile programs using Guile. Usually, the former is installed by default but the latter is not. You need to install both, which means that you probably need to install `guile-dev`. Similarly for any other library packages needed by Meep.

Physics
-------

### How does the current amplitude relate to the resulting field amplitude?

There is no simple formula relating the input current amplitude (**J** in Maxwell's equations) to the resulting fields (**E**) etcetera, even at the same point as the current. The exact same current will produce a different field and radiate a different total power depending upon the surrounding materials/geometry, and depending on the frequency. This is a physical consequence of the geometry's effect on the local density of states; it can also be thought of as feedback from reflections on the source. As a simple example, if you put a current source inside a perfect electric conductor, the resulting field will be zero. As another example, the frequency-dependence of the radiated power in vacuum is part of the reason why the sky is blue.

See also our online book chapter on [Wave Source Conditions](http://arxiv.org/abs/arXiv:1301.5366).

If you are worried about this, then you are probably setting up your calculation in the wrong way. Especially in linear materials, the absolute magnitude of the field is useless; the only meaningful quantities are dimensionless ratios like the fractional transmission: the transmitted power relative to the transmitted power in some reference calculation. Almost always, you want to perform two calculations, one of which is a reference, and compute the ratio of a result in one calculation to the result in the reference. For nonlinear calculations, see [Units and Nonlinearity](Units_and_Nonlinearity.md).

### How do I set the imaginary part of ε?

If you only care about the imaginary part of ε in a narrow bandwidth around some frequency ω, you should set it by using the electric conductivity as described in [Conductivity](Materialsmd#Conductivity). If you care about the imaginary part of ε over a broad bandwidth, then for any physical material the imaginary part will be frequency-dependent and you will have to fit it to Meep's dispersive-ε support. See [Material Dispersion](Materials.md#material-dispersion).

Meep doesn't implement a frequency-independent complex ε. Not only is this not physical, but it also leads to both exponentially decaying and exponentially growing solutions in Maxwell's equations from positive- and negative-frequency Fourier components, respectively. Thus, it cannot be simulated in the time domain.

### Why does my simulation diverge if ε &lt; 0?

Maxwell's equations have exponentially growing solutions for a frequency-independent negative ε. For any physical medium with negative ε, there must be dispersion, and you must likewise use dispersive materials in Meep to have negative ε at some desired frequency. The requirement of dispersion to obtain negative ε follows from the Kramers–Kronig relations, and also follows from thermodynamic considerations that the energy in the electric field must be positive; see, for example, *Electrodynamics of Continuous Media* by Landau, Pitaevskii, and Lifshitz. At an even more fundamental level, it can be derived from [passivity constraints](http://arxiv.org/abs/arXiv:1405.0238).

If you solve Maxwell's equations in a homogeneous-epsilon material at some real wavevector **k**, you get a dispersion relation $\omega^2 = c^2 |\mathbf{k}|^2 / \varepsilon$. If ε is positive, there are two real solutions $\omega^2 = \pm c |\mathbf{k}| / \sqrt{\varepsilon}$, giving oscillating solutions. If ε is negative, there are two *imaginary* solutions ω, corresponding to exponentially decaying and *exponentially growing solutions* from any current source. These solutions can always be spatially decomposed into a superposition of real-**k** values via a spatial Fourier transform.

If you do a simulation of any kind in the time domain (not just FDTD), you pretty much can't avoid exciting both the decaying and the growing solutions. This is *not* a numerical instability, it is a real solution of the underlying equations for an unphysical material.

See [Materials](Materials.md) for how to include dispersive materials which can have negative ε and loss.

If you have negative ε *and* negative μ *everywhere*, the case of a negative-index material, then the simulation is fine. However at the boundary between negative- and positive-index materials, you will encounter instabilities: because of the way Maxwell's equations are discretized in FDTD, the ε and μ are discretized on different spatial grids, so you will get a half-pixel or so of εμ &lt; 0 at the boundary between negative and positive indices, which will cause the simulation to diverge. But of course, any physical negative-index metamaterial also involves dispersion.

Note also that, as a consequence of the above analysis, ε must go to a positive value in the $\omega\to\pm\infty$ limit to get non-diverging solutions of Maxwell's equations. So the $\epsilon_\infty$ in your [dispersion model](Materials.md) must be positive.

Usage
-----

### Why doesn't turning off subpixel averaging work?

By default, when Meep assigns a dielectric constant ε or μ to each pixel, it uses a carefully designed average of the ε values within that pixel. This subpixel averaging generally improves the accuracy of the simulation—perhaps counter-intuitively, for geometries with discontinous ε it is *more* accurate (i.e. closer to the exact Maxwell result for the *discontinuous* case) to do the simulation with the subpixel-averaged (*smoothed*) ε, as long as the averaging is done properly. For details, see the [reference publication](License_and_Copyright.md#referencing).

Still, there are times when, for whatever reason, you might not want this feature. For example, if your accuracy is limited by other issues, or if you want to skip the wait at the beginning of the simulation for it do to the averaging. In this case, you can disable the subpixel averaging by doing `(set! eps-averaging? false)` in your control file. See the [User Interface](Scheme_User_Interface.md).

Even if you disable the subpixel averaging, however, when you output the dielectric function to a file and plot it, you may notice that there are some pixels with intermediate ε values, right at the boundary between two materials. This has a completely different source. Internally, Meep's simulation is performed on a [Yee grid](Yee_Lattice.md), in which every field component is stored on a slightly different grid which are offset from one another by half-pixels, and the ε values are also stored on this Yee grid. For output purposes, however, it is more user-friendly to output all fields etcetera on the same grid at the center of each pixel, so all quantities are interpolated onto this grid for output. Therefore, even though the internal ε values are indeed discontinuous when you disable subpixel averaging, the *output* file will still contain some "averaged" values at interfaces due to the interpolation from the Yee grid to the center-pixel grid.
