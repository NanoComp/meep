---
# Introduction
---

Meep implements the **finite-difference time-domain** (**FDTD**) method for computational electromagnetics. This is a widely used technique in which space is divided into a discrete grid and then the fields are evolved in time using discrete time steps &mdash; as the grid and the time steps are made finer and finer, this becomes a closer and closer approximation for the true continuous equations, and one can simulate many practical problems essentially exactly.

In this section, we introduce the equations and the electromagnetic units employed by Meep, the FDTD method, and Meep's approach to FDTD. Also, FDTD is only one of several useful methods in computational electromagnetics, each of which has their own special uses &mdash; we mention a few of the other methods, and try to give some hints as to which applications FDTD is well suited for and when you should consider a different method.

This introduction does *not* describe the user interface with which you can tell Meep to perform these tasks. Instead, we focus here on the *concepts* that are being simulated. The user interface is introduced in [Tutorial/Basics](Python_Tutorials/Basics).

[TOC]

Maxwell's Equations
-------------------

Meep simulates [Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell's_equations), which describe the interactions of electric (**E**) and magnetic (**H**) fields with one another and with matter and sources. In particular, the equations for the evolution of the fields are:

<center>

$\frac{d\mathbf{B}}{dt} = -\nabla\times\mathbf{E} - \mathbf{J}_B - \sigma_B \mathbf{B}$

$\mathbf{B} = \mu \mathbf{H}$

$\frac{d\mathbf{D}}{dt} = \nabla\times\mathbf{H} - \mathbf{J} - \sigma_D \mathbf{D}$

$\mathbf{D} = \varepsilon \mathbf{E}$

</center>

Where **D** is the displacement field, $\varepsilon$ is the dielectric constant, **J** is the current density (of electric charge), and **J**<sub>*B*</sub> is the *magnetic-charge* current density. Magnetic currents are a convenient computational fiction in some situations. **B** is the magnetic flux density (often called the magnetic field), $\mu$ is the magnetic permeability, and **H** is the magnetic field. The $\sigma_B$ and $\sigma_D$ terms correspond to (frequency-independent) magnetic and electric conductivities, respectively. The divergence equations are implicitly:

<center>

$\nabla \cdot \mathbf{B} = - \int^t \nabla \cdot (\mathbf{J}_B(t') + \sigma_B \mathbf{B}) dt'$

$\nabla \cdot \mathbf{D} = - \int^t \nabla \cdot (\mathbf{J}(t') + \sigma_D \mathbf{D})dt' \equiv \rho$

</center>

Most generally, $\varepsilon$ depends not only on position but also on frequency (material dispersion) and on the field **E** itself (nonlinearity), and may include loss or gain. These effects are supported in Meep and are described in [Materials](Materials.md).

Meep supports simulation in [Cylindrical Coordinates](Cylindrical_Coordinates.md).

### Units in Meep

You may have noticed the lack of annoying constants like $\varepsilon$<sub>0</sub>, $\mu$<sub>0</sub>, and [c](https://en.wikipedia.org/wiki/Speed_of_light) &mdash; that's because Meep uses **dimensionless** units where all these constants are unity. As a practical matter, almost everything you might want to compute (transmission spectra, frequencies, etcetera) is expressed as a ratio anyway, so the units end up cancelling.

In particular, because Maxwell's equations are scale invariant (multiplying the sizes of everything by 10 just divides the corresponding solution frequencies by 10), it is convenient in electromagnetic problems to choose **scale-invariant units**. See Chapter 2 of [Photonic Crystals: Molding the Flow of Light (second edition)](http://ab-initio.mit.edu/book). That means that we pick some characteristic lengthscale in the system, $a$, and use that as our unit of distance.

Moreover, since $c=1$ in Meep units, $a$ (or $a/c$) is our unit of *time* as well. In particular, the frequency *f* in Meep (corresponding to a time dependence $e^{-i 2\pi f t}$) is specified in units of $c/a$ (or equivalently $\omega$ is specified in units of $2\pi c/a$), which is equivalent to specifying *f* as $1/T$: the inverse of the optical period $T$ in units of $a/c$. This, in turn, is equivalent to specifying *f* as $a/\lambda$ where $\lambda$ is the vacuum wavelength. A similar scheme is used in [MPB](https://mpb.readthedocs.io).

For example, suppose we are describing some photonic structure at infrared frequencies, where it is convenient to specify distances in microns. Thus, we let $a$ = 1 &#956;m. Then, if we want to specify a source corresponding to $\lambda$ = 1.55 &#956;m, we specify the frequency *f* as 1/1.55 = 0.6452. If we want to run our simulation for 100 periods, we then run it for 155 time units (= 100 / *f*).

A transmission spectrum, for example, would be a ratio of transmitted to incident intensities, so the units of **E** are irrelevant unless there are [nonlinearities](Units_and_Nonlinearity/).

The Bloch wavevector (see below) **k** is specified in Cartesian coordinates in units of $2\pi/a$. This is *different* from MPB: it is equivalent to taking MPB's k-points and transforming them with `reciprocal->cartesian`.

Boundary Conditions and Symmetries
----------------------------------

Since only a finite region of space can be simulated, the simulation must always be terminated with some **boundary conditions**. Three basic types of terminations are supported in Meep: **Bloch-periodic boundaries**, **metallic walls**, and **PML absorbing layers**. Also, one can exploit **symmetries** of a problem to further reduce the computational requirements.

With ordinary periodic boundaries in a cell of size $L$, the field components satisfy $f(x+L) = f(x)$. **Bloch periodicity** is a generalization where $f(x+L) = e^{ik_x L} f(x)$ for some *Bloch wavevector* $\mathbf{k}$. This can be used to solve for the modes of waveguides, gratings, and so on, much like in [MPB](https://mpb.readthedocs.io). See Chapter 3 of [Photonic Crystals: Molding the Flow of Light (second edition)](http://ab-initio.mit.edu/book).

An even simpler boundary condition is a metallic wall, where the fields are simply forced to be zero on the boundaries, as if the cell were surrounded by a perfect metal (zero absorption, zero skin depth). More generally, you can place perfect metal materials anywhere you want in the computational cell, e.g. to simulate [metallic cavities](Python_Tutorials/Local_Density_of_States/) of an arbitrary shape.

To simulate open boundary conditions, one would like the boundaries to absorb all waves incident on them, with no reflections. This is implemented with something called **perfectly matched layers** (PML). PML is, strictly speaking, not a boundary condition &mdash; rather, it is a special absorbing material placed adjacent to the boundaries. PML is actually a fictitious (non-physical) material, designed to have zero reflections at its interface. Although PML is reflectionless in the theoretical continuous system, in the actual discretized system it has some small reflections which make it imperfect. For this reason, one always gives the PML some finite thickness in which the absorption gradually turns on. For more information, see [Perfectly Matched Layer](Perfectly_Matched_Layer.md).

Another way in which the computational cell is reduced in size is by **symmetry**. For example, if you know that your system has a mirror symmetry plane (both in the structure and in the current sources), then you can save a factor of two by only simulating half of the structure and obtaining the other half by mirror reflection. Meep can exploit several kinds of mirror and rotational symmetries — it is designed so that the symmetry is purely an optimization, and other than specifying the symmetry your computation is set up in exactly the same way. See [Exploiting Symmetry](Exploiting_Symmetry.md).

Finite-Difference Time-Domain Methods
-------------------------------------

FDTD methods divide space and time into a finite rectangular grid. As [described below](#the-illusion-of-continuity), Meep tries to hide this discreteness from the user as much as possible, but there are a few consequences of discretization that it is good to be familiar with.

Perhaps the most important thing you need to know is this: if the grid has some spatial resolution $\Delta x$, then our discrete time-step $\Delta t$ is given by $\Delta t = S \Delta x$, where $S$ is the Courant factor and must satisfy $S < n_\textrm{min} / \sqrt{\mathrm{\# dimensions}}$, where $n_\textrm{min}$ is the minimum refractive index (usually 1), in order for the method to be stable (not diverge). In Meep, $S$=0.5 by default (which is sufficient for 1 to 3 dimensions), but can be changed by the user. This means that **when you double the grid resolution, the number of time steps doubles as well** (for the same simulation period). Thus, in three dimensions, if you double the resolution, then the amount of memory increases by 8 and the amount of computational time increases by (at least) 16.

The second most important thing you should know is that, in order to discretize the equations with second-order accuracy, FDTD methods **store different field components at different grid locations**. This discretization is known as a [Yee lattice](Yee_Lattice.md). As a consequence, **Meep must interpolate the field components to a common point** whenever you want to combine, compare, or output the field components (e.g. in computing energy density or flux). Most of the time, you don't need to worry too much about this interpolation since it is automatic. However, because it is a simple linear interpolation, while **E** and **D** may be discontinuous across dielectric boundaries, it means that the interpolated **E** and **D** fields may be less accurate than you might expect right around dielectric interfaces.

Many references are available on FDTD methods for computational electromagnetics. See, for example:

- A. Taflove and S.C. Hagness, [Computational Electrodynamics: The Finite-Difference Time-Domain Method](https://www.amazon.com/Computational-Electrodynamics-Finite-Difference-Time-Domain-Method/dp/1580538320), Artech: Norwood, MA, 2005.

- A. Taflove, A. Oskooi, and S.G. Johnson, [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707), Artech: Norwood, MA, 2013.

### The Illusion of Continuity

Although FDTD inherently uses discretized space and time, as much as possible Meep attempts to maintain the illusion that you are using a continuous system. At the beginning of the simulation, you specify the spatial resolution, but from that point onwards you generally work in continuous coordinates in your chosen units. See [Units in Meep](Introduction.md#units-in-meep), above.

For example, you specify the dielectric function as a function $\varepsilon$(**x**) of continuous **x**, or as a set of solid objects like spheres, cylinders, etcetera, and Meep is responsible for figuring out how they are to be represented on a discrete grid. Or if you want to specify a point source, you simply specify the point **x** where you want the source to reside &mdash; Meep will figure out the closest grid points to **x** and add currents to those points, weighted according to their distance from **x**. If you change **x** continuously, the current in Meep will also change continuously by changing the weights. If you ask for the flux through a certain rectangle, then Meep will linearly interpolate the field values from the grid onto that rectangle.

In general, the philosophy of the Meep interface is **pervasive interpolation**, so that if you change any input continuously then the response of the Meep simulation will change continuously as well, so that it will converge as rapidly and as smoothly as possible to the continuous solution as you increase the spatial resolution.

For example, the $\varepsilon$ function used internally by Meep is not simply a discretely sampled version of the $\varepsilon$(**x**) specified by the user. Rather, each grid point is a kind of average of the $\varepsilon$ in the surrounding pixel. Meep's subpixel averaging, described in Chapter 6 ("Accurate FDTD Simulation of Discontinuous Materials by Subpixel Smoothing") of the book [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707), is specially designed in order to minimize the "staircasing" and other errors caused by sharp interfaces.

Other Numerical Methods in Computational Electromagnetics
---------------------------------------------------------

FDTD is, of course, not the only numerical method in computational electromagnetics, nor is it always the best one. In general, we advocate having several tools in your toolbox, and selecting the most convenient one for each task. See Appendix D of [Photonic Crystals: Molding the Flow of Light (second edition)](http://ab-initio.mit.edu/book).

For example, although FDTD can be used to compute electromagnetic eigenmodes (below), in lossless structures it is often quicker, easier, and more reliable to use a specialized eigenmode solver such as [MPB](http://mpb.readthedocs.io). See also the [frequency vs. time domain](http://mpb.readthedocs.io/en/latest/Introduction/) discussion in the MPB manual and the [resonant modes](Introduction.md#resonant-modes) discussion below.

For computing the field pattern or response of a structure at a *single frequency*, it may be more efficient to directly solve the corresponding linear equation rather than iterating in time. Indeed, we have an experimental implementation of this method directly in Meep (i.e. a finite-difference frequency-domain solver) — see the [frequency-domain solver](Scheme_User_Interface.md#frequency-domain-solver). However, especially in cases where there are large differences in scale (e.g. with metals with a shallow skin depth), it may be better to use a method that allows a variable resolution in different spatial regions, such as a finite-element or boundary-element method. Boundary-element methods are especially powerful when you have a large volume-to-surface ratio, such as for scattering calculations over small objects in a large (∞) volume.

A strength of time-domain methods is their ability to obtain the entire frequency spectrum of responses (or eigenfrequencies) in a single simulation, by Fourier-transforming the response to a short pulse or using more sophisticated signal-processing methods such as [Harminv](https://github.com/stevengj/harminv/blob/master/README.md). Finite-element methods can also be used for time-evolving fields, but they suffer a serious disadvantage compared to finite-difference methods: finite-element methods, for stability, must typically use some form of *implicit time-stepping*, where they must invert a matrix (solve a linear system) at every time step.

Finally, in systems that are composed of a small number of easily-analyzed pieces, such as a sequence of constant-cross-section waveguides, a collection of cylinders, or a multi-layer film, transfer-matrix/scattering-matrix methods may be especially attractive. These methods treat the individual simple elements in some analytic or semi-analytic fashion, enabling the entire structure to be simulated with great speed and accuracy. There are too many such techniques to easily summarize here, but one useful free tool that can handle a wide variety of structures is [CAMFR](http://camfr.sf.net/).

Applications of FDTD
--------------------

In this section, we sketch out a few of the basic ways in which FDTD can be used to analyze electromagnetic problems. Specific examples of how to use these techniques in Meep are described in [Tutorial/Basics](Python_Tutorials/Basics).

### Field Patterns and Green's Functions

The most obvious thing that you can do with a time-domain simulation, of course, is to simply get a picture of the field pattern resulting from a given source, or perhaps a movie showing the field evolution in time.

The field pattern from a given localized source at a particular frequency $\omega$ is a form of the **Green's function** of the system. More specifically, one typically writes the "dyadic" Green's function

$$G_{ij}(\omega; \mathbf{x}, \mathbf{x}')$$ which gives the $i$<sup>th</sup> component of (say) **E** at **x** from a point current source **J** at $\mathbf{x'}$, such that $\mathbf{J}(\mathbf{x})=\hat{\mathbf{e}_j} \cdot \exp(-i\omega t) \cdot \delta(\mathbf{x}-\mathbf{x}')$. To obtain this in FDTD, you simply place the requisite point source at $\mathbf{x'}$ and wait for a long enough time for all other frequency components to die out (noting that the mere act of turning on a current source at $t=0$ introduces a spectrum of frequencies). Alternatively, you can use Meep's frequency-domain solver to find the response directly (by solving the associated linear equation).

Given the Green's function, one can then compute a wide variety of useful things, from the radiated flux, to the local density of states (proportional to $\sum_i G_{ii}$), to Born approximations for small scatterers. Even more powerfully, one can compute many such quantities *for multiple frequencies simultaneously* using the Fourier transform of a short pulse as described below.

### Transmission/Reflection Spectra

Perhaps the most common task to which FDTD is applied is that of computing the transmission or scattering spectra from some finite structure, such as a resonant cavity, in response to some stimulus. One could, of course, compute the fields (and thus the transmitted flux) at each frequency $\omega$ separately, as described above. However, it is much more efficient to compute a broadband response via a single computation by Fourier-transforming the response to a short pulse.

For example, suppose we want the transmitted power through some structure. For fields at a given frequency $\omega$, this is the integral of the Poynting vector (in the normal $\hat{\mathbf{n}}$ direction) over a plane on the far side of the structure:

$$P(\omega) = \mathrm{Re}\, \hat{\mathbf{n}}\cdot \int \mathbf{E}_\omega(\mathbf{x})^* \times \mathbf{H}_\omega(\mathbf{x}) \, d^2\mathbf{x}$$ Now, if we input a short pulse, it is tempting to compute the integral $P(t)$ of the Poynting vector at each time, and then Fourier-transform this to find $P(\omega)$. That is **incorrect**, however, because what we want is the flux of the Fourier-transformed *fields* **E** and **H**, which is not the same as the transform of the flux (because the flux is not a linear function of the fields).

Instead, what one does is to accumulate the Fourier transforms $\mathbf{E}_\omega(\mathbf{x})$ and $\mathbf{H}_\omega(\mathbf{x})$ for every point in the flux plane via summation over the discrete time steps $n$:

$$\tilde{f}(\omega) = \frac{1}{\sqrt{2\pi}}  \sum_n e^{i\omega n \Delta t} f(n\Delta t) \Delta t \approx \frac{1}{\sqrt{2\pi}} \int e^{i\omega t} f(t) dt$$ and then, at the end of the time-stepping, computing $P(\omega)$ by the fluxes of these Fourier-transformed fields. Meep takes care of all of this for you automatically, of course—you simply specify the regions over which you want to integrate the flux, and the frequencies that you want to compute.

There are other possible methods of time-series analysis, of course. One method that is sometimes very effective is construct a [Padé approximant](https://en.wikipedia.org/wiki/Padé_approximant) of the time series of field values at some point, from which one can often [extrapolate a very accurate discrete-time Fourier transform](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=923035), including sharp peaks and other resonant features, from a relatively short time series. Meep does not provide a Padé computation for you, but of course you can output the fields at a point over time, ideally in a single-mode waveguide for transmission spectra via a single point, and compute the Padé approximant yourself by standard methods.

The power $P(\omega)$ by itself is not very useful — one needs to *normalize*, dividing by the incident power at each frequency, to get the transmission spectrum. Typically, this is done by running the simulation *twice*: once with only the incident wave and no scattering structure, and once with the scattering structure, where the first calculation is used for normalization.

It gets more tricky if one wants to compute the *reflection* spectrum as well as the transmission. You can't simply compute the flux in the backwards direction, because this would give you the sum of the reflected and the incident power. You also can't simply subtract the incident power from backwards flux to get the transmitted power, because in general there will be interference effects (between incident and reflected waves) that are not subtracted. Rather, you have to subtract the Fourier-transformed incident *fields* $\mathbf{E}_\omega^{(0)}(\mathbf{x})$ and $\mathbf{H}_\omega^{(0)}(\mathbf{x})$ to get the reflected/scattered power:

$$P_r(\omega) = \mathrm{Re}\,\hat{\mathbf{n}}\cdot\int \left[ \mathbf{E}_\omega(\mathbf{x}) - \mathbf{E}_\omega^{(0)}(\mathbf{x}) \right]^* \times \left[ \mathbf{H}_\omega(\mathbf{x}) - \mathbf{H}_\omega^{(0)}(\mathbf{x}) \right] \, d^2\mathbf{x}$$ Again, you can do this easily in practice by running the simulation twice, once without and once with the scatterer, and telling Meep to subtract the Fourier transforms in the reflected plane before computing the flux. And again, after computing the reflected power you will normalize by the incident power to get the reflection spectrum.

Meep is designed to make these kinds of calculations easy, as long as you have some idea of what is going on.

### Resonant Modes

Another common task in FDTD is to compute resonant modes or eigenmodes of a given structure. For example, suppose you have a photonic crystal (periodic dielectric structure) or a waveguide and you want to know its harmonic (definite-$\omega$) modes at a given wavevector **k**. Or, suppose you have a resonant cavity that traps light in a small region for a long time, and you want to know the resonant frequency $\omega$ and the decay lifetime (quality factor) *Q*. And, of course, you may want the field patterns of these modes as well.

In order to extract the frequencies and lifetimes (which may be infinite in a lossless system) with FDTD, the basic strategy is simple. You set up the structure with Bloch-periodic and/or absorbing boundaries, depending on whether it is a periodic or open system. Then you excite the mode(s) with a short pulse (broad bandwidth) from a current placed directly inside the cavity/waveguide/whatever. Finally, once the current source is turned off, you have some fields bouncing around inside the system, and you analyze them to extract the frequencies and decay rates.

The simplest form of harmonic analysis would be to compute the Fourier transform of the fields at some point — harmonic modes will yield sharp peaks in the spectrum. This method has serious drawbacks, however, in that high frequency resolution requires a very long running time, and moreover the problem of extracting the decay rates leads to a poorly-conditioned nonlinear fitting problem. Instead, Meep allows you to perform a more sophisticated signal processing algorithm borrowed from NMR spectroscopy — the algorithm is called *filter diagonalization* and is implemented by [Harminv](https://github.com/stevengj/harminv/blob/master/README.md) package. Harminv extracts all of the frequencies and their decay rates (and amplitudes) within a short time to high accuracy; for example, we have used it to find lifetimes *Q* of 10<sup>9</sup> periods in a computational run of only a few hundred periods.

Once you know the frequencies of the modes, if you want the field patterns you will need to run the simulation again with a *narrow*-bandwidth (long-time) pulse to excite only the mode in question. Unless you want the longest-lifetime mode, in which case you can just run long enough for the other modes to decay away. Given the field patterns, you can then perform other analyses (e.g. decomposing the *Q* into decay rates into different directions via flux computations, finding modal volumes, etcetera).

Why should you use Meep instead of [MPB](https://mpb.readthedocs.io) to compute the modes? Unlike MPB, Meep supports metallic and absorbing materials, can compute lossy resonant modes, can quickly compute large numbers of $\omega$'s at once by a single short pulse, and can efficiently extract modes in the interior of the spectrum (e.g. in a band gap). Why should you ever use MPB, then? MPB is quicker at computing the lowest-$\omega$ modes than Meep, gives you both the $\omega$ and the fields at once, and has no problem resolving closely-spaced frequencies. Moreover, computing modes in time domain is somewhat subtle and requires care &mdash; for example, one will occasionally miss a mode if the source happens to be nearly orthogonal to it or if it is too close to another mode; conversely, the signal processing will sometimes accidentally identify spurious peak frequencies. Also, studying periodic systems with non-rectangular unit cells is more subtle in Meep than in MPB. MPB is much more straightforward and reliable, albeit more limited in some ways.
