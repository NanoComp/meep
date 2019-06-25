---
# Scheme User Interface
---

This page is a listing of the functions exposed by the Scheme interface. For a gentler introduction, see [Tutorial/Basics](Scheme_Tutorials/Basics.md). This page does not document the Scheme language or the functions provided by [libctl](https://libctl.readthedocs.io). Also, note that this page is not a complete listing of all functions. In particular, because of the [SWIG wrappers](#swig-wrappers), every function in the C++ interface is accessible from Scheme, but not all of these functions are documented or intended for end users. See also the instructions for [parallel Meep](Parallel_Meep.md).

**Note:** The Scheme interface is being deprecated and has been replaced by the [Python interface](Python_User_Interface.md).

[TOC]

Input Variables
---------------

These are global variables that you can set to control various parameters of the Meep computation. In brackets after each variable is the type of value that it should hold. The classes, complex datatypes like `geometric-object`, are described in a later subsection. The basic datatypes, like `integer`, `boolean`, `cnumber`, and `vector3`, are defined by [libctl](https://libctl.readthedocs.io/en/latest/Libctl_User_Reference/).

**`geometry` [ list of `geometric-object` class ]**
—
Specifies the geometric objects making up the structure being simulated. When objects overlap, later objects in the list take precedence. Defaults to no objects (empty list).

**`geometry-center` [ `vector3` class ]**
—
Specifies the coordinates of the center of the cell. Defaults to (0, 0, 0), but changing this allows you to shift the coordinate system used in Meep (for example, to put the origin at the corner).

**`sources` [ list of `source` class ]**
—
Specifies the current sources to be present in the simulation. Defaults to none.

**`symmetries` [ list of `symmetry` class ]**
—
Specifies the spatial symmetries (mirror or rotation) to exploit in the simulation. Defaults to none. The symmetries must be obeyed by *both* the structure and the sources. See also [Exploiting Symmetry](Exploiting_Symmetry.md).

**`pml-layers` [ list of `pml` class ]**
—
Specifies the [PML](Perfectly_Matched_Layer.md) absorbing boundary layers to use. Defaults to none.

**`geometry-lattice` [`lattice` class ]**
—
Specifies the size of the unit cell which is centered on the origin of the coordinate system. Any sizes of `no-size` imply a reduced-dimensionality calculation. A 2d calculation is especially optimized. See `dimensions` below. Defaults to a cubic cell of unit size.

**`default-material` [`material-type` class ]**
—
Holds the default material that is used for points not in any object of the geometry list. Defaults to `air` (ε=1). See also `epsilon-input-file` below.

**`epsilon-input-file` [`string`]**
—
If this string is not empty (the default), then it should be the name of an HDF5 file whose first/only dataset defines a scalar, real-valued, frequency-independent dielectric function over some discrete grid. Alternatively, the dataset name can be specified explicitly if the string is in the form "filename:dataset". This dielectric function is then used in place of the ε property of `default-material` (i.e. where there are no `geometry` objects). The grid of the epsilon file dataset need *not* match the computational grid; it is scaled and/or linearly interpolated as needed to map the file onto the cell. The structure is warped if the proportions of the grids do not match. **Note:** the file contents only override the ε property of the `default-material`, whereas other properties (μ, susceptibilities, nonlinearities, etc.) of `default-material` are still used.

**`dimensions` [`integer`]**
—
Explicitly specifies the dimensionality of the simulation, if the value is less than 3. If the value is 3 (the default), then the dimensions are automatically reduced to 2 if possible when `geometry-lattice` size in the $z$ direction is `no-size`. If `dimensions` is the special value of `CYLINDRICAL`, then cylindrical coordinates are used and the $x$ and $z$ dimensions are interpreted as $r$ and $z$, respectively. If `dimensions` is 1, then the cell must be along the $z$ direction and only $E_x$ and $H_y$ field components are permitted. If `dimensions` is 2, then the cell must be in the $xy$ plane.

**`m` [`number`]**
—
For `CYLINDRICAL` simulations, specifies that the angular $\phi$ dependence of the fields is of the form $e^{im\phi}$ (default is `m=0`). If the simulation cell includes the origin $r=0$, then `m` must be an integer.

**`accurate-fields-near-cylorigin?` [`boolean`]**
—For `CYLINDRICAL` simulations with |*m*| &gt; 1, compute more accurate fields near the origin $r=0$ at the expense of requiring a smaller Courant factor. Empirically, when this option is set to `true`, a Courant factor of roughly $\min[0.5, 1 / (|m| + 0.5)]$ or smaller seems to be needed. Default is `false`, in which case the $D_r$, $D_z$, and $B_r$ fields within |*m*| pixels of the origin are forced to zero, which usually ensures stability with the default Courant factor of 0.5, at the expense of slowing convergence of the fields near $r=0$.

**`resolution` [`number`]**
—
Specifies the computational grid resolution in pixels per distance unit. Default is 10.

**`k-point` [`false` or `vector3`]**
—
If `false` (the default), then the boundaries are perfect metallic (zero electric field). If a `vector3`, then the boundaries are Bloch-periodic: the fields at one side are $\exp(i\mathbf{k}\cdot\mathbf{R})$ times the fields at the other side, separated by the lattice vector $\mathbf{R}$. A non-zero `vector3` will produce complex fields. The `k-point` vector is specified in Cartesian coordinates in units of 2π/distance. Note: this is *different* from [MPB](https://mpb.readthedocs.io), equivalent to taking MPB's `k-points` through its function `reciprocal->cartesian`.

**`ensure-periodicity` [`boolean`]**
—
If `true` (the default) *and* if the boundary conditions are periodic (`k-point` is not `false`), then the geometric objects are automatically repeated periodically according to the lattice vectors which define the size of the cell.

**`eps-averaging?` [`boolean`]**
—
If `true` (the default), then subpixel averaging is used when initializing the dielectric function. For details, see Section 3 ("Interpolation and the illusion of continuity") of [Computer Physics Communications, Vol. 181, pp. 687-702, 2010](http://ab-initio.mit.edu/~oskooi/papers/Oskooi10.pdf). The input variables `subpixel-maxeval` (default 10<sup>4</sup>) and `subpixel-tol` (default 10<sup>-4</sup>) specify the maximum number of function evaluations and the integration tolerance for subpixel averaging. Increasing/decreasing these, respectively, will cause a more accurate but slower computation of the average ε with diminishing returns for the actual FDTD error.

**`force-complex-fields?` [`boolean`]**
—
By default, Meep runs its simulations with purely real fields whenever possible. It uses complex fields which require twice the memory and computation if the `k-point` is non-zero or if `m` is non-zero. However, by setting `force-complex-fields?` to `true`, Meep will always use complex fields.

**`filename-prefix` [`string`]**
—
A string prepended to all output filenames. If empty (the default), then Meep uses the name of the current ctl file, with ".ctl" replaced by "-" (e.g. `foo.ctl` uses a `"foo-"` prefix). See also [Output File Names](Scheme_User_Interface.md#output-file-names).

**`Courant` [`number`]**
—
Specify the Courant factor $S$ which relates the time step size to the spatial discretization: $cΔ t = SΔ x$. Default is 0.5. For numerical stability, the Courant factor must be *at most* $n_\textrm{min}/\sqrt{\textrm{# dimensions}}$, where $n_\textrm{min}$ is the minimum refractive index (usually 1), and in practice $S$ should be slightly smaller.

**`output-volume` [`meep::volume*`]**
—
Specifies the default region of space that is output by the HDF5 output functions (below); see also the `(volume ...)` function to create `meep::volume*` objects. Default is `'()` (null), which means that the whole cell is output. Normally, you should use the `(in-volume ...)` function to modify the output volume instead of setting `output-volume` directly.

**`output-single-precision?` [`boolean`]**
—
Meep performs its computations in [double precision](https://en.wikipedia.org/wiki/double_precision), and by default its output HDF5 files are in the same format. However, by setting this variable to `true` (default is `false`) you can instead output in [single precision](https://en.wikipedia.org/wiki/single_precision) which saves a factor of two in space.

**`progress-interval` [`number`]**
—
Time interval (seconds) after which Meep prints a progress message. Default is 4 seconds.

**`extra-materials` [ list of `material-type` class ]**
—
By default, Meep turns off support for material dispersion (via susceptibilities or conductivity) or nonlinearities if none of the objects in `geometry` have materials with these properties &mdash; since they are not needed, it is faster to omit their calculation. This doesn't work, however, if you use a `material-function`: materials via a user-specified function of position instead of just geometric objects. If your material function only returns a nonlinear material, for example, Meep won't notice this unless you tell it explicitly via `extra-materials`. `extra-materials` is a list of materials that Meep should look for in the cell in addition to any materials that are specified by geometric objects. You should list any materials other than scalar dielectrics that are returned by `material-function` here.

The following require a bit more understanding of the inner workings of Meep to use. See also [SWIG Wrappers](#swig-wrappers).

**`structure` [`meep::structure*`]**
—
Pointer to the current structure being simulated; initialized by `(init-structure)` which is called automatically by `(init-fields)` which is called automatically by any of the [`(run)` functions](#run-functions).

**`fields` [`meep::fields*`]**
—
Pointer to the current fields being simulated; initialized by `(init-fields)` which is called automatically by any of the `(run)` functions.

**`num-chunks` [`integer`]**
—
Minimum number of "chunks" (subarrays) to divide the structure/fields into (default 0). Actual number is determined by number of processors, PML layers, etcetera. Mainly useful for debugging.

Predefined Variables
--------------------

**`air`, `vacuum` [`material-type` class ]**
—
Two aliases for a predefined material type with a dielectric constant of 1.

**`perfect-electric-conductor` or `metal` [`material-type` class ]**
—
A predefined material type corresponding to a perfect electric conductor at the boundary of which the parallel electric field is zero. Technically, $\varepsilon = -\infty$.

**`perfect-magnetic-conductor` [`material-type` class ]**
—
A predefined material type corresponding to a perfect magnetic conductor at the boundary of which the parallel magnetic field is zero. Technically, $\mu = -\infty$.

**`nothing` [`material-type` class ]**
—
A material that, effectively, punches a hole through other objects to the background (`default-material`).

**`infinity` [`number`]**
—
A big number (10<sup>20</sup>) to use for "infinite" dimensions of objects.

**`pi` [`number`]**
—
π (3.14159...).

Constants (Enumerated Types)
----------------------------

Several of the functions/classes in Meep ask you to specify e.g. a field component or a direction in the grid. These should be one of the following constants:

**`direction` constants**
—
Specify a direction in the grid. One of `X`, `Y`, `Z`, `R`, `P` for $x$, $y$, $z$, $r$, $\phi$, respectively.

**`side` constants**
—
Specify particular boundary in the positive `High` (e.g., +`X`) or negative `Low` (e.g., -`X`) direction.

**`component` constants**
—
Specify a particular field or other component. One of `Ex`, `Ey`, `Ez`, `Er`, `Ep`, `Hx`, `Hy`, `Hz`, `Hy`, `Hp`, `Hz`, `Bx`, `By`, `Bz`, `By`, `Bp`, `Bz`, `Dx`, `Dy`, `Dz`, `Dr`, `Dp`, `Dielectric`, `Permeability`, for $E_x$, $E_y$, $E_z$, $E_r$, $E_\phi$, $H_x$, $H_y$, $H_z$, $H_r$, $H_\phi$, $B_x$, $B_y$, $B_z$, $B_r$, $B_\phi$, $D_x$, $D_y$, $D_z$, $D_r$, $D_\phi$, ε, μ, respectively.

**`derived-component` constants**
—
These are additional components which are not actually stored by Meep but are computed as needed, mainly for use in output functions. One of `Sx`, `Sy`, `Sz`, `Sr`, `Sp`, `EnergyDensity`, `D-EnergyDensity`, `H-EnergyDensity` for $S_x$, $S_y$, $S_z$, $S_r$, $S_\phi$ (components of the Poynting vector $\mathrm{Re}\,\mathbf{E}^* \times \mathbf{H}$), $(\mathbf{E}^* \cdot \mathbf{D} + \mathbf{H}^* \cdot \mathbf{B})/2$, $\mathbf{E}^* \cdot \mathbf{D}/2$, $\mathbf{H}^* \cdot \mathbf{B}/2$, respectively.

Classes
-------

Classes are complex datatypes with various properties which may have default values. Classes can be "subclasses" of other classes. Subclasses inherit all the properties of their superclass and can be used in any place the superclass is expected. An object of a class is constructed with:

```scm
(make class (prop1 val1) (prop2 val2) ...)
```

See also the [libctl manual](https://libctl.readthedocs.io).

Meep defines several types of classes, the most numerous of which are the various geometric object classes which are the same as those used in [MPB](https://mpb.readthedocs.io). You can also get a list of the available classes, along with their property types and default values, at runtime with the `(help)` command.

### lattice

The `lattice` class is normally used only for the `geometry-lattice` variable, which sets the size of the cell. In [MPB](https://mpb.readthedocs.io), you can use this to specify a variety of affine lattice structures. In [Meep](index.md), only rectangular Cartesian cells are supported, so the only property of lattice that you should normally use is its `size`.

**`size` [`vector3`]**
—
The size of the cell. Defaults to unit lengths.

If any dimension has the special size `no-size`, then the dimensionality of the problem is essentially reduced by one. Strictly speaking, the dielectric function is taken to be uniform along that dimension.

Because Maxwell's equations are scale invariant, you can use any units of distance you want to specify the cell size: nanometers, microns, centimeters, etc. However, it is usually convenient to pick some characteristic lengthscale of your problem and set that length to 1. See also [Units](Introduction.md#units-in-meep).

### material

This class is used to specify the materials that geometric objects are made of. Currently, there are three subclasses, `dielectric`, `perfect-metal`, and `material-function`.

**`medium`**

An electromagnetic medium which is possibly nonlinear and/or dispersive. See also [Materials](Materials.md). For backwards compatibility, a synonym for `medium` is `dielectric`. It has several properties:

**`epsilon` [`number`]**
—The frequency-independent isotropic relative permittivity or dielectric constant. Default is 1. You can also use `(index n)` as a synonym for `(epsilon (* n n))`; note that this is not really the refractive index if you also specify μ, since the true index is $\sqrt{\mu\varepsilon}$.

Using `(epsilon ep)` is actually a synonym for `(epsilon-diag ep ep ep)`.

**`epsilon-diag` and `epsilon-offdiag` [`vector3`]**
—
These properties allow you to specify ε as an arbitrary real-symmetric tensor by giving the diagonal and offdiagonal parts. Specifying `(epsilon-diag a b c)` and/or `(epsilon-offdiag u v w)` corresponds to a relative permittivity ε tensor
\\begin{pmatrix} a & u & v \\\\ u & b & w \\\\ v & w & c \\end{pmatrix}

Default is the identity matrix ($a = b = c = 1$ and $u = v = w = 0$).

**`mu` [`number`]**
—
The frequency-independent isotropic relative permeability μ. Default is 1. Using `(mu pm)` is actually a synonym for `(mu-diag pm pm pm)`.

**`mu-diag` and `mu-offdiag` [`vector3`]**
—
These properties allow you to specify μ as an arbitrary real-symmetric tensor by giving the diagonal and offdiagonal parts exactly as for ε above. Default is the identity matrix.

**`D-conductivity` [`number`]**
—
The frequency-independent electric conductivity $σ_D$. Default is 0. You can also specify a diagonal anisotropic conductivity tensor by using the property `D-conductivity-diag` which takes three numbers or a `vector3` to give the $σ_D$ tensor diagonal. See also [Conductivity](Materials.md#conductivity-and-complex).

**`B-conductivity` [`number`]**
—
The frequency-independent magnetic conductivity $σ_B$. Default is 0. You can also specify a diagonal anisotropic conductivity tensor by using the property `B-conductivity-diag` which takes three numbers or a `vector3` to give the $σ_B$ tensor diagonal. See also [Conductivity](Materials.md#conductivity-and-complex).

**`chi2` [`number`]**
—
The nonlinear ([Pockels](https://en.wikipedia.org/wiki/Pockels_effect)) susceptibility $\chi^{(2)}$. Default is 0. See also [Nonlinearity](Materials.md#nonlinearity).

**`chi3` [`number`]**
—
The nonlinear ([Kerr](https://en.wikipedia.org/wiki/Kerr_effect)) susceptibility $\chi^{(3)}$. Default is 0. See also [Nonlinearity](Materials.md#nonlinearity).

**`E-susceptibilities` [ list of `susceptibility` class ]**
—
List of dispersive susceptibilities (see below) added to the dielectric constant ε in order to model material dispersion. Defaults to none. See also [Material Dispersion](Materials.md#material-dispersion). For backwards compatibility, synonyms of `E-susceptibilities` are `E-polarizations` and `polarizations`.

**`H-susceptibilities` [ list of `susceptibility` class ]**
—
List of dispersive susceptibilities (see below) added to the permeability μ in order to model material dispersion. Defaults to none. See also [Material Dispersion](Materials.md#material-dispersion).

**`perfect-metal`**

A perfectly-conducting metal. This class has no properties and you normally just use the predefined `metal` object, above. To model imperfect conductors, use a dispersive dielectric material. See also the [Predefined Variables](#predefined-variables): `metal`, `perfect-electric-conductor`, and `perfect-magnetic-conductor`.

**`material-function`**

This material type allows you to specify the material as an arbitrary function of position. It has one property:

**`material-func` [`function`]**
—
A function of one argument, the position `vector3`, that returns the material at that point. Note that the function you supply can return *any* material. It's even possible to return another `material-function` object which would then have its function invoked in turn.

Instead of `material-func`, you can use `epsilon-func`: give it a function of position that returns the dielectric constant at that point.

**Important:** If your material function returns nonlinear, dispersive (Lorentzian or conducting), or magnetic materials, you should also include a list of these materials in the `extra-materials` input variable (above) to let Meep know that it needs to support these material types in your simulation. For dispersive materials, you need to include a material with the *same* values of γ<sub>*n*</sub> and ω<sub>*n*</sub>, so you can only have a finite number of these, whereas σ<sub>*n*</sub> can vary continuously and a matching σ<sub>*n*</sub> need not be specified in `extra-materials`. For nonlinear or conductivity materials, your `extra-materials` list need not match the actual values of σ or χ returned by your material function, which can vary continuously.

**Complex ε and μ**: you cannot specify a frequency-independent complex ε or μ in Meep where the imaginary part is a frequency-independent loss but there is an alternative. That is because there are only two important physical situations. First, if you only care about the loss in a narrow bandwidth around some frequency, you can set the loss at that frequency via the [conductivity](Materials.md#conductivity-and-complex). Second, if you care about a broad bandwidth, then all physical materials have a frequency-dependent complex ε and/or μ, and you need to specify that frequency dependence by fitting to Lorentzian and/or Drude resonances via the `lorentzian-susceptibility` or `drude-susceptibility` classes below.

Dispersive dielectric and magnetic materials, above, are specified via a list of objects that are subclasses of type `susceptibility`.

### susceptibility

Parent class for various dispersive susceptibility terms, parameterized by an anisotropic amplitude σ. See [Material Dispersion](Materials.md#material-dispersion).

**`sigma` [`number`]**
—
The scale factor σ. You can also specify an anisotropic σ tensor by using the property `sigma-diag` which takes three numbers or a `vector3` to give the σ$_n$ tensor diagonal, and `sigma-offdiag` which specifies the offdiagonal elements (defaults to 0). That is, `(sigma-diag a b c)` and `(sigma-offdiag u v w)` corresponds to a σ tensor

\\begin{pmatrix} a & u & v \\\\ u & b & w \\\\ v & w & c \\end{pmatrix}

### lorentzian-susceptibility

Specifies a single dispersive susceptibility of Lorentzian (damped harmonic oscillator) form. See [Material Dispersion](Materials.md#material-dispersion), with the parameters (in addition to σ):

**`frequency` [`number`]**
—
The resonance frequency $f_n = \omega_n / 2\pi$.

**`gamma` [`number`]**
—
The resonance loss rate $γ_n / 2\pi$.

Note: multiple objects with identical values for the `frequency` and `gamma` but different `sigma` will appear as a *single* Lorentzian susceptibility term in the preliminary simulation info output.

### drude-susceptibility

Specifies a single dispersive susceptibility of Drude form. See [Material Dispersion](Materials.md#material-dispersion), with the parameters (in addition to σ):

**`frequency` [`number`]**
—
The frequency scale factor $f_n = \omega_n / 2\pi$ which multiplies σ (not a resonance frequency).

**`gamma` [`number`]**
—
The loss rate $γ_n / 2\pi$.

Meep also supports a somewhat unusual polarizable medium, a Lorentzian susceptibility with a random noise term added into the damped-oscillator equation at each point. This can be used to directly model thermal radiation in both the [far field](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.93.213905) and the [near field](http://math.mit.edu/~stevenj/papers/RodriguezIl11.pdf). Note, however that it is more efficient to compute far-field thermal radiation using Kirchhoff's law of radiation, which states that emissivity equals absorptivity. Near-field thermal radiation can usually be computed more efficiently using frequency-domain methods, e.g. via [SCUFF-EM](http://homerreid.dyndns.org/scuff-EM/).

### multilevel-atom

Specifies a multievel atomic susceptibility for modeling saturable gain and absorption. This is a subclass of `E-susceptibilities` which contains two objects: (1) `transitions`: a list of atomic `transition`s (defined below), and (2) `initial-populations`: a list of numbers defining the initial population of each atomic level. See [Materials/Saturable Gain and Absorption](Materials.md#saturable-gain-and-absorption).

#### transition

**`frequency` [`number`]**
—
The radiative transition frequency $f = \omega / 2\pi$.

**`gamma` [`number`]**
—
The loss rate $\gamma = \gamma / 2\pi$.

**`sigma` [`number`]**
—
The coupling strength $\sigma$.

**`from-level` [`number`]**
—
The atomic level from which the transition occurs.

**`to-level` [`number`]**
—
The atomic level to which the transition occurs.

**`transition-rate` [`number`]**
—
The non-radiative transition rate $f = \omega / 2\pi$. Default is 0.

**`pumping-rate` [`number`]**
—
The pumping rate $f = \omega / 2\pi$. Default is 0.

### noisy-lorentzian-susceptibility or noisy-drude-susceptibility

Specifies a single dispersive susceptibility of Lorentzian (damped harmonic oscillator) or Drude form. See [Material Dispersion](Materials.md#material-dispersion), with the same `sigma`, `frequency`, and `gamma` parameters, but with an additional Gaussian random noise term (uncorrelated in space and time, zero mean) added to the **P** damped-oscillator equation.

**`noise-amp` [`number`]**
—
The noise has root-mean square amplitude σ $\times$ `noise-amp`.

### gyrotropic-lorentzian-susceptibility or gyrotropic-drude-susceptibility

(**Experimental feature**) Specifies a single dispersive [gyrotropic susceptibility](Materials.md#gyrotropic-media) of [Lorentzian (damped harmonic oscillator) or Drude form](Materials.md#gyrotropic-drude-lorentz-model). Its parameters are `sigma`, `frequency`, and `gamma`, which have the [usual meanings](#susceptibility), and an additional 3-vector `bias`:

**`bias` [`vector3`]**
—
The gyrotropy vector.  Its direction determines the orientation of the gyrotropic response, and the magnitude is the precession frequency $|\mathbf{b}_n|/2\pi$.

### gyrotropic-saturated-susceptibility

(**Experimental feature**) Specifies a single dispersive [gyrotropic susceptibility](Materials.md#gyrotropic-media) governed by a [linearized Landau-Lifshitz-Gilbert equation](Materials.md#gyrotropic-saturated-dipole-linearized-landau-lifshitz-gilbert-model). This class takes parameters `sigma`, `frequency`, and `gamma`, whose meanings are different from the Lorentzian and Drude case. It also takes a 3-vector `bias` parameter and an `alpha` parameter:

**`sigma` [`number`]**
—
The coupling factor $\sigma_n / 2\pi$ between the polarization and the driving field. In magnetic ferrites, this is the Larmor precession frequency at the saturation field.

**`frequency` [`number`]**
—
The Larmor precession frequency, $f_n = \omega_n / 2\pi$.

**`gamma` [`number`]**
—
The loss rate $\gamma_n / 2\pi$ in the off-diagonal response.

**`alpha` [`number`]**
—
The loss factor $\alpha_n$ in the diagonal response. Note that this parameter is dimensionless and contains no 2π factor.

**`bias` [`vector3`]**
—
Vector specifying the orientation of the gyrotropic response. Unlike the similarly-named `bias` parameter for the [gyrotropic Lorentzian/Drude susceptibilities](#gyrotropiclorentziansusceptibility-or-gyrotropicdrudesusceptibility), the magnitude is ignored; instead, the relevant precession frequencies are determined by the `sigma` and `frequency` parameters.

### geometric-object

This class, and its descendants, are used to specify the solid geometric objects that form the dielectric structure being simulated. The base class is:

**`geometric-object`**

Properties:

**`material` [`material-type` class ]**
—
The material that the object is made of (usually some sort of dielectric). No default value (must be specified).

**`center` [`vector3`]**
—
Center point of the object. No default value.

One normally does not create objects of type `geometric-object` directly, however; instead, you use one of the following subclasses. Recall that subclasses inherit the properties of their superclass, so these subclasses automatically have the `material` and `center` properties which must be specified, since they have no default values.

In a 2d calculation, only the intersections of the objects with the $xy$ plane are considered.

### sphere

A sphere. Properties:

**`radius` [`number`]**
—
Radius of the sphere. No default value.

### cylinder

A cylinder, with circular cross-section and finite height. Properties:

**`radius` [`number`]**
—
Radius of the cylinder's cross-section. No default value.

**`height` [`number`]**
—
Length of the cylinder along its axis. No default value.

**`axis` [`vector3`]**
—
Direction of the cylinder's axis; the length of this vector is ignored. Defaults to point parallel to the $z$ axis.

### cone

A cone, or possibly a truncated cone. This is actually a subclass of `cylinder`, and inherits all of the same properties, with one additional property. The radius of the base of the cone is given by the `radius` property inherited from `cylinder`, while the radius of the tip is given by the new property, `radius2`. The `center` of a cone is halfway between the two circular ends.

**`radius2` [`number`]**
—
Radius of the tip of the cone (i.e. the end of the cone pointed to by the `axis` vector). Defaults to zero (a "sharp" cone).

### block

A parallelepiped (i.e., a brick, possibly with non-orthogonal axes).

**`size` [`vector3`]**
—
The lengths of the block edges along each of its three axes. Not really a 3-vector, but it has three components, each of which should be nonzero. No default value.

**`e1`, `e2`, `e3` [`vector3`]**
—
The directions of the axes of the block; the lengths of these vectors are ignored. Must be linearly independent. They default to the three lattice directions.

### ellipsoid

An ellipsoid. This is actually a subclass of `block`, and inherits all the same properties, but defines an ellipsoid inscribed inside the block.

### prism

Polygonal prism type.

**`vertices` [list of `vector3`]**
—
The vertices that define the polygonal *floor* of the prism; the vertices must be coplanar, and if `axis` is specified it must be normal to the plane of the vertices. Note that infinite prism lengths are not supported. To simulate infinite geometry, just extend the edge of the prism beyond the cell. The *ceiling* of the prism is just its floor polygon rigidly translated through the displacement vector `height*axis`.

**`height` [`number`]**
—
The prism thickness, extruded in the direction of `axis`. `infinity` can be used for infinite height.

**`axis` [`vector3`]**
—
(optional) specifies the extrusion axis, which must be normal to the plane of the vertices. If `axis` is not specified, the extrusion axis is taken to be the normal vector to the plane of the vertices, with sign determined by a right-hand rule with respect to the first two vertices: if your right-hand fingers point from vertex 1 to 2, your thumb points in the direction of `axis.` In vector language, `axis` is determined by computing a vector cross product and normalizing to unit magnitude:
$$ \mathbf{a}
  =(\mathbf{v}_1 - \overline{\mathbf v})
   \times
   (\mathbf{v}_2 - \overline{\mathbf{v}}),
   \quad
   \texttt{axis}\equiv \frac{\mathbf{a}}{|\mathbf{a}|}
$$
where $\mathbf{v}_{1,2}$ are the first and second `vertices` and $\overline{\mathbf{v}}\equiv\frac{1}{N}\sum_{n=1}^N \mathbf{v}_n$
is the *centroid* of the polygon (with $N\ge 3$ the length of the
`vertices` array).

There are two options for specifying the `center` of a prism. In contrast to the other types of `geometric-object`, the center of a prism does not need to be explicitly specified, because it may be calculated from `vertices`, `height`, and `axis.` (Specifically, we have `center = centroid + 0.5*height*axis,` where the `centroid` was defined above). To create a `prism` with the center computed automatically in this way, simply initialize the `center` field of the `prism` class (inherited from `geometric_object`) to the special initializer keyword `auto-center`. On the other hand, in some cases you may want to override this automatic calculation and instead specify your own `center` for a prism; this will have the effect of rigidly translating the entire prism so that it is centered at the point you specify. See below for examples of both possibilities.

These are some examples of geometric objects created using the above classes:

```scm
; A cylinder of infinite radius and height 0.25 pointing along the x axis,
; centered at the origin:

(make cylinder (center 0 0 0) (material (make dielectric (index 3.5))) 
               (radius infinity) (height 0.25) (axis 1 0 0))
```

```scm
; An ellipsoid with its long axis pointing along (1,1,1), centered on
; the origin (the other two axes are orthogonal and have equal semi-axis lengths)

(make ellipsoid (center 0 0 0) (material (make dielectric (epsilon 12.0)))
                (size 0.8 0.2 0.2)
                (e1 1 1 1)
                (e2 0 1 -1)
                (e3 -2 1 1))
```

```scm
; A unit cube of material metal with a spherical air hole of radius 0.2 at
; its center, the whole thing centered at (1,2,3):

(set! geometry (list
               (make block (center 1 2 3) (material metal) (size 1 1 1))
               (make sphere (center 1 2 3) (material air) (radius 0.2))))
```

```scm
; A hexagonal prism defined by six vertices centered on the origin
; and extruded in the z direction to a height of 1.5
; of material crystalline silicon (from the materials library)

(set! geometry
      (list
       (make prism
         (vertices
           (list
                 (vector3 -1 0 0)
                 (vector3 -0.5 (/ (sqrt 3) 2) 0)
                 (vector3 0.5 (/ (sqrt 3) 2) 0)
                 (vector3 1 0 0)
                 (vector3 0.5 (/ (sqrt 3) -2) 0)
                 (vector3 -0.5 (/ (sqrt 3) -2) 0)))
         (axis 0 0 1)
         (height 1.5)
         (center auto-center)
         (material cSi))))
```
Note the use of `(center auto-center)` to establish that the prism center will be computed automatically from the vertices, axes, and height &mdash; which, in this case, will put the center at $(0,0,0.75)$.

```scm
; The same hexagonal prism, but now rigidly displaced so that
; its center lies at (0.4, 0.8, -0.2):

(set! geometry
      (list
       (make prism
         (vertices
           (list
                 (vector3 -1 0 0)
                 (vector3 -0.5 (/ (sqrt 3) 2) 0)
                 (vector3 0.5 (/ (sqrt 3) 2) 0)
                 (vector3 1 0 0)
                 (vector3 0.5 (/ (sqrt 3) -2) 0)
                 (vector3 -0.5 (/ (sqrt 3) -2) 0)))
         (axis 0 0 1)
         (height 1.5)
         (center 0.4 0.8 -0.2)
         (material cSi))))
```

### symmetry

This class is used for the `symmetries` input variable to specify symmetries which must preserve both the structure *and* the sources. Any number of symmetries can be exploited simultaneously but there is no point in specifying redundant symmetries: the cell can be reduced by at most a factor of 4 in 2d and 8 in 3d. See also [Exploiting Symmetry](Exploiting_Symmetry.md).

**`symmetry`**

A single symmetry to exploit. This is the base class of the specific symmetries below, so normally you don't create it directly. However, it has two properties which are shared by all symmetries:

**`direction` [`direction` constant ]**
—
The direction of the symmetry (the normal to a mirror plane or the axis for a rotational symmetry). e.g. `X`, `Y`, `Z` (only Cartesian/grid directions are allowed). No default value.

**`phase` [`cnumber`]**
—
An additional phase to multiply the fields by when operating the symmetry on them. Default is +1, e.g. a phase of -1 for a mirror plane corresponds to an *odd* mirror. Technically, you are essentially specifying the representation of the symmetry group that your fields and sources transform under.

The specific symmetry sub-classes are:

**`mirror-sym`**
—
A mirror symmetry plane. `direction` is the direction *normal* to the mirror plane.

**`rotate2-sym`**
—
A 180° (twofold) rotational symmetry (a.k.a. $C_2$). `direction` is the axis of the rotation.

**`rotate4-sym`**
—
A 90° (fourfold) rotational symmetry (a.k.a. $C_4$). `direction` is the axis of the rotation.

### pml

This class is used for specifying the PML absorbing boundary layers around the cell, if any, via the `pml-layers` input variable. See also [Perfectly Matched Layers](Perfectly_Matched_Layer.md). `pml-layers` can be zero or more `pml` objects, with multiple objects allowing you to specify different PML layers on different boundaries.

**`pml`**

A single PML layer specification, which sets up one or more PML layers around the boundaries according to the following properties.

**`thickness` [`number`]**
—
The spatial thickness of the PML layer which extends from the boundary towards the *inside* of the cell. The thinner it is, the more numerical reflections become a problem. No default value.

**`direction` [`direction` constant ]**
—
Specify the direction of the boundaries to put the PML layers next to. e.g. if `X`, then specifies PML on the $\pm x$ boundaries (depending on the value of `side`, below). Default is the special value `ALL`, which puts PML layers on the boundaries in all directions.

**`side` [`boundary-side` constant ]**
—
Specify which side, `Low` or `High` of the boundary or boundaries to put PML on. e.g. if side is `Low` and direction is `X`, then a PML layer is added to the $-x$ boundary. Default is the special value `ALL`, which puts PML layers on both sides.

**`strength` [`number`]**
—
A strength (default is 1.0) to multiply the PML absorption coefficient by. A strength of 2.0 will *square* the theoretical asymptotic reflection coefficient of the PML (making it smaller), but will also increase numerical reflections. Alternatively, you can change `R-asymptotic`, below.

**`R-asymptotic` [`number`]**
—
The asymptotic reflection in the limit of infinite resolution or infinite PML thickness, for reflections from air (an upper bound for other media with index &gt; 1). For a finite resolution or thickness, the reflection will be *much larger*, due to the discretization of Maxwell's equation. Default value is 10<sup>−15</sup>, which should suffice for most purposes. You want to set this to be small enough so that waves propagating within the PML are attenuated sufficiently, but making `R-asymptotic` too small will increase the numerical reflection due to discretization.

**`pml-profile` [`function`]**
—
By default, Meep turns on the PML conductivity quadratically within the PML layer &mdash; one doesn't want to turn it on suddenly, because that exacerbates reflections due to the discretization. More generally, with `pml-profile` one can specify an arbitrary PML "profile" function $f(u)$ that determines the shape of the PML absorption profile up to an overall constant factor. *u* goes from 0 to 1 at the start and end of the PML, and the default is $f(u) = u^2$. In some cases where a very thick PML is required, such as in a periodic medium (where there is technically no such thing as a true PML, only a pseudo-PML), it can be advantageous to turn on the PML absorption more smoothly. See [Optics Express, Vol. 16, pp. 11376-92, 2008](http://www.opticsinfobase.org/abstract.cfm?URI=oe-16-15-11376). For example, one can use a cubic profile $f(u) = u^3$ by specifying `(pml-profile (lambda (u) (* u u u)))`.

#### `absorber`

Instead of a `pml` layer, there is an alternative class called `absorber` which is a **drop-in** replacement for `pml`. For example, you can do `(set! pml-layers (list (make absorber (thickness 2))))` instead of `(set! pml-layers (list (make pml (thickness 2))))`. All the parameters are the same as for `pml`, above. You can have a mix of `pml` on some boundaries and `absorber` on others.

The `absorber` class does *not* implement a perfectly matched layer (PML), however (except in 1d). Instead, it is simply a scalar electric **and** magnetic conductivity that turns on gradually within the layer according to the `pml-profile` (defaulting to quadratic). Such a scalar conductivity gradient is only reflectionless in the limit as the layer becomes sufficiently thick.

The main reason to use `absorber` is if you have **a case in which PML fails:**

-   No true PML exists for *periodic* media, and a scalar absorber is computationally less expensive and generally just as good. See [Optics Express, Vol. 16, pp. 11376-92, 2008](http://www.opticsinfobase.org/abstract.cfm?URI=oe-16-15-11376).
-   PML can lead to *divergent* fields for certain waveguides with "backward-wave" modes; this can readily occur in metals with surface plasmons, and a scalar absorber is your only choice. See [Physical Review E, Vol. 79, 065601, 2009](http://math.mit.edu/~stevenj/papers/LohOs09.pdf).
-   PML can fail if you have a waveguide hitting the edge of your cell *at an angle*. See [J. Computational Physics, Vol. 230, pp. 2369-77, 2011](http://math.mit.edu/~stevenj/papers/OskooiJo11.pdf).

### source

The `source` class is used to specify the current sources via the `sources` input variable. Note that all sources in Meep are separable in time and space, i.e. of the form $\mathbf{J}(\mathbf{x},t) = \mathbf{A}(\mathbf{x}) \cdot f(t)$ for some functions $\mathbf{A}$ and $f$. Non-separable sources can be simulated, however, by modifying the sources after each time step. When real fields are being used (which is the default in many cases; see the `force-complex-fields?` input variable), only the real part of the current source is used.

**Important note**: These are *current* sources (**J** terms in Maxwell's equations), even though they are labelled by electric/magnetic field components. They do *not* specify a particular electric/magnetic field which would be what is called a "hard" source in the FDTD literature. There is no fixed relationship between the current source and the resulting field amplitudes; it depends on the surrounding geometry, as described in the [FAQ](FAQ#how-does-the-current-amplitude-relate-to-the-resulting-field-amplitude) and in Section 4.4 ("Currents and Fields: The Local Density of States") in [Chapter 4](http://arxiv.org/abs/arXiv:1301.5366) ("Electromagnetic Wave Source Conditions") of the book [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707).

**`source`**

The source class has the following properties:

**`src` [`src-time` class ]**
—
Specify the time-dependence of the source (see below). No default.

**`component` [`component` constant ]**
—
Specify the direction and type of the current component: e.g. `Ex`, `Ey`, etcetera for an electric-charge current, and `Hx`, `Hy`, etcetera for a magnetic-charge current. Note that currents pointing in an arbitrary direction are specified simply as multiple current sources with the appropriate amplitudes for each component. No default.

**`center` [`vector3`]**
—
The location of the center of the current source in the cell. No default.

**`size` [`vector3`]**
—
The size of the current distribution along each direction of the cell. Default is (0,0,0): a point-dipole source.

**`amplitude` [`cnumber`]**
—
An overall complex amplitude multiplying the the current source. Default is 1.0.

**`amp-func` [`function`]**
—
A Scheme function of a single argument, that takes a vector3 giving a position and returns a complex current amplitude for that point. The position argument is *relative* to the `center` of the current source, so that you can move your current around without changing your function. Default is `'()` (null), meaning that a constant amplitude of 1.0 is used. Note that your amplitude function (if any) is *multiplied* by the `amplitude` property, so both properties can be used simultaneously.

As described in Section 4.2 ("Incident Fields and Equivalent Currents") in [Chapter 4](http://arxiv.org/abs/arXiv:1301.5366) ("Electromagnetic Wave Source Conditions") of the book [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707), it is also possible to supply a source that is designed to couple exclusively into a single waveguide mode (or other mode of some cross section or periodic region) at a single frequency, and which couples primarily into that mode as long as the bandwidth is not too broad. This is possible if you have [MPB](https://mpb.readthedocs.io) installed: Meep will call MPB to compute the field profile of the desired mode, and uses the field profile to produce an equivalent current source. Note: this feature does *not* work in cylindrical coordinates. To do this, instead of a `source` you should use an `eigenmode-source`:

### eigenmode-source

This is a subclass of `source` and has **all of the properties** of `source` above. However, you normally do not specify a `component`. Instead of `component`, the current source components and amplitude profile are computed by calling MPB to compute the modes, $\mathbf{u}_{n,\mathbf{k}}(\mathbf{r}) e^{i \mathbf{k} \cdot \mathbf{r}}$, of the dielectric profile in the region given by the `size` and `center` of the source, with the modes computed as if the *source region were repeated periodically in all directions*. If an `amplitude` and/or `amp-func` are supplied, they are *multiplied* by this current profile. The desired eigenmode and other features are specified by the following properties:

**`eig-band` [`integer`]**
—
The index *n* (1,2,3,...) of the desired band ω<sub>*n*</sub>(**k**) to compute in MPB where 1 denotes the lowest-frequency band at a given **k** point, and so on.

**`direction` [`X`, `Y`, or `Z;` default `AUTOMATIC`], `eig-match-freq?` [`boolean;` default `true`], `eig-kpoint` [`vector3`]**
—
By default (if `eig-match-freq?` is `true`), Meep tries to find a mode with the same frequency ω<sub>*n*</sub>(**k**) as the `src` property (above), by scanning **k** vectors in the given `direction` using MPB's `find-k` functionality. Alternatively, if `eig-kpoint` is supplied, it is used as an initial guess for **k**. By default, `direction` is the direction normal to the source region, assuming `size` is $d$–1 dimensional in a $d$-dimensional simulation (e.g. a plane in 3d). If `direction` is set to `NO-DIRECTION`, then `eig_kpoint` is not only the initial guess and the search direction of the **k** vectors, but is also taken to be the direction of the waveguide, allowing you to [launch modes in oblique ridge waveguides](Scheme_Tutorials/Eigenmode_Source.md#index-guided-modes-in-a-ridge-waveguide) (not perpendicular to the source plane). If `eig-match-freq?` is `false`, then the specific **k** vector of the desired mode is specified with `eig-kpoint` (in Meep units of 2π/(unit length)). By default, the **k** components in the plane of the source region are zero.  However, if the source region spans the *entire* cell in some directions, and the cell has Bloch-periodic boundary conditions via the `k-point` parameter, then the mode's **k** components in those directions will match `k-point` so that the mode satisfies the Meep boundary conditions, regardless of `eig-kpoint`. Note that once **k** is either found by MPB, or specified by `eig-kpoint`, the field profile used to create the current sources corresponds to the [Bloch mode](https://en.wikipedia.org/wiki/Bloch_wave), $\mathbf{u}_{n,\mathbf{k}}(\mathbf{r})$, multiplied by the appropriate exponential factor, $e^{i \mathbf{k} \cdot \mathbf{r}}$.

**`eig-parity` [`NO-PARITY` (default), `EVEN-Z`, `ODD-Z`, `EVEN-Y`, `ODD-Y`]**
—
The parity (= polarization in 2d) of the mode to calculate, assuming the structure has $z$ and/or $y$ mirror symmetry *in the source region*, with respect to the `center` of the source region.  (In particular, it does not matter if your simulation as a whole has that symmetry, only the cross section where you are introducing the source.) If the structure has both $y$ and $z$ mirror symmetry, you can combine more than one of these, e.g. `EVEN-Z + ODD-Y`. Default is `NO-PARITY`, in which case MPB computes all of the bands which will still be even or odd if the structure has mirror symmetry, of course. This is especially useful in 2d simulations to restrict yourself to a desired polarization.

**`eig-resolution` [`integer`, defaults to same as Meep resolution ]**
—
The spatial resolution to use in MPB for the eigenmode calculations. This defaults to the same resolution as Meep, but you can use a higher resolution in which case the structure is linearly interpolated from the Meep pixels.

**`eig-tolerance` [`number`, defaults to 10<sup>–12</sup> ]**
—
The tolerance to use in the MPB eigensolver. MPB terminates when the eigenvalues stop changing to less than this fractional tolerance.  (Note that this is the tolerance for the frequency eigenvalue ω; the tolerance for the mode profile is effectively the square root of this.)

**`component` [as above, but defaults to `ALL-COMPONENTS`]**
—
Once the MPB modes are computed, equivalent electric and magnetic sources are created within Meep. By default, these sources include magnetic and electric currents in *all* transverse directions within the source region, corresponding to the mode fields as described in Section 4.2 ("Incident Fields and Equivalent Currents") in [Chapter 4](http://arxiv.org/abs/arXiv:1301.5366) ("Electromagnetic Wave Source Conditions") of the book [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707). If you specify a `component` property, however, you can include only one component of these currents if you wish. Most users won't need this feature.

**`eig-lattice-size` [`vector3`], `eig-lattice-center` [`vector3`]**
—
Normally, the MPB computational unit cell is the same as the source volume given by the `size` and `center` parameters. However, occasionally you want the unit cell to be larger than the source volume. For example, to create an eigenmode source in a periodic medium, you need to pass MPB the entire unit cell of the periodic medium, but once the mode is computed then the actual current sources need only lie on a cross section of that medium. To accomplish this, you can specify the optional `eig-lattice-size` and `eig-lattice-center`, which define a volume (which must enclose `size` and `center`) that is used for the unit cell in MPB with the dielectric function ε taken from the corresponding region in the Meep simulation.

Note that Meep's MPB interface only supports dispersionless non-magnetic materials but it does support anisotropic ε. Any nonlinearities, magnetic responses μ, conductivities σ, or dispersive polarizations in your materials will be *ignored* when computing the eigenmode source. PML will also be ignored.

The `src-time` object, which specifies the time dependence of the source, can be one of the following three classes.

### continuous-src

A continuous-wave (CW) source is proportional to $\exp(-i\omega t)$, possibly with a smooth (exponential/tanh) turn-on/turn-off. In practice, the CW source [never produces an exact single-frequency response](FAQ.md#why-doesnt-the-continuous-wave-cw-source-produce-an-exact-single-frequency-response).

**`frequency` [`number`]**
—
The frequency *f* in units of $c$/distance or ω in units of 2π$c$/distance. See [Units](Introduction.md#units-in-meep). No default value. You can instead specify `(wavelength x)` or `(period x)`, which are both a synonym for `(frequency (/ 1 x))`; i.e. 1/ω in these units is the vacuum wavelength or the temporal period.

**`start-time` [`number`]**
—
The starting time for the source. Default is 0 (turn on at $t=0$).

**`end-time` [`number`]**
—
The end time for the source. Default is `infinity` (never turn off).

**`width` [`number`]**
—
Roughly, the temporal width of the smoothing (technically, the inverse of the exponential rate at which the current turns off and on). Default is 0 (no smoothing). You can instead specify `(fwidth x)`, which is a synonym for `(width (/ 1 x))` (i.e. the frequency width is proportional to the inverse of the temporal width).

**`slowness` [`number`]**
—
Controls how far into the exponential tail of the tanh function the source turns on. Default is 3.0. A larger value means that the source turns on more gradually at the beginning.

### gaussian-src

A Gaussian-pulse source roughly proportional to $\exp(-i\omega t - (t-t_0)^2/2w^2)$. Technically, the "Gaussian" sources in Meep are the (discrete-time) derivative of a Gaussian, i.e. they are $(-i\omega)^{-1} \frac{\partial}{\partial t} \exp(-i\omega t - (t-t_0)^2/2w^2)$, but the difference between this and a true Gaussian is usually irrelevant.

**`frequency` [`number`]**
—
The center frequency $f$ in units of $c$/distance (or ω in units of 2π$c$/distance). See [Units](Introduction.md#units-in-meep). No default value. You can instead specify `(wavelength x)` or `(period x)`, which are both a synonym for `(frequency (/ 1 x))`; i.e. 1/ω in these units is the vacuum wavelength or the temporal period.

**`width` [`number`]**
—
The width $w$ used in the Gaussian. No default value. You can instead specify `(fwidth x)`, which is a synonym for `(width (/ 1 x))` (i.e. the frequency width is proportional to the inverse of the temporal width).

**`start-time` [`number`]**
—
The starting time for the source. Default is 0 (turn on at $t=0$). This is not the time of the peak. See below.

**`cutoff` [`number`]**
—
How many `width`s the current decays for before it is cut off and set to zero &mdash; this applies for both turn-on and turn-off of the pulse. Default is 5.0. A larger value of `cutoff` will reduce the amount of high-frequency components that are introduced by the start/stop of the source, but will of course lead to longer simulation times. The peak of the Gaussian is reached at the time $t_0$=`start-time + cutoff*width`.

### custom-src

A user-specified source function $f(t)$. You can also specify start/end times at which point your current is set to zero whether or not your function is actually zero. These are optional, but you must specify an `end-time` explicitly if you want functions like `run-sources` to work, since they need to know when your source turns off. For a demonstration of a [linear-chirped pulse](FAQ.md#how-do-i-create-a-chirped-pulse), see [`examples/chirped-pulse.ctl`](https://github.com/NanoComp/meep/blob/master/scheme/examples/chirped-pulse.ctl).

**`src-func` [`function`]**
—
The function $f(t)$ specifying the time-dependence of the source. It should take one argument (the time in Meep units) and return a complex number.

**`start-time` [`number`]**
—
The starting time for the source. Default is `(-infinity)`: turn on at $t=-\infty$. Note, however, that the simulation normally starts at $t=0$ with zero fields as the initial condition, so there is implicitly a sharp turn-on at $t=0$ whether you specify it or not.

**`end-time` [`number`]**
—
The end time for the source. Default is `infinity` (never turn off).

### flux-region

A `flux-region` object is used with [`add-flux`](#flux-spectra) to specify a region in which Meep should accumulate the appropriate Fourier-transformed fields in order to compute a flux spectrum.

**`flux-region`**
—
A region (volume, plane, line, or point) in which to compute the integral of the Poynting vector of the Fourier-transformed fields.

**`center` [`vector3`]**
—The center of the flux region (no default).

**`size` [`vector3`]**
—The size of the flux region along each of the coordinate axes. Default is `(0,0,0)`; a single point.

**`direction` [`direction` constant ]**
—The direction in which to compute the flux (e.g. `X`, `Y`, etcetera). Default is `AUTOMATIC`, in which the direction is determined by taking the normal direction if the flux region is a plane (or a line, in 2d). If the normal direction is ambiguous (e.g. for a point or volume), then you *must* specify the `direction` explicitly (not doing so will lead to an error).

**`weight` [`cnumber`]**
—A weight factor to multiply the flux by when it is computed. Default is 1.0.

Note that the flux is always computed in the *positive* coordinate direction, although this can effectively be flipped by using a `weight` of -1.0. This is useful, for example, if you want to compute the outward flux through a box, so that the sides of the box add instead of subtract.

Miscellaneous Functions
-----------------------

Here, we describe a number of miscellaneous useful functions provided by Meep.

### Geometry Utilities

Some utility functions are provided to help you manipulate geometric objects:

**`(shift-geometric-object obj shift-vector)`**
—
Translate `obj` by the 3-vector `shift-vector`.

**`(geometric-object-duplicates shift-vector min-multiple max-multiple obj)`**
—
Return a list of duplicates of `obj`, shifted by various multiples of `shift-vector` from `min-multiple` to `max-multiple`, inclusive, in steps of 1.

**`(geometric-objects-duplicates shift-vector min-multiple max-multiple obj-list)`**
—
Same as `geometric-object-duplicates`, except operates on a list of objects, `obj-list`. If *A* appears before *B* in the input list, then all the duplicates of *A* appear before all the duplicates of *B* in the output list.

**`(geometric-objects-lattice-duplicates obj-list [ ux uy uz ])`**
—
Duplicates the objects in `obj-list` by multiples of the Cartesian basis vectors, making all possible shifts of the "primitive cell" (see below) that fit inside the lattice cell. The primitive cell to duplicate is `ux` by `uy` by `uz`, in units of the Cartesian basis vectors. These three parameters are optional; any that you do not specify are assumed to be `1`.

**`point_in_object(point, obj)`**
—
Returns whether or not the given 3-vector `point` is inside the geometric object `obj`.

**`(point-in-periodic-object? point obj)`**
—
As `point-in-object?`, but also checks translations of the given object by the lattice vectors.

**`(display-geometric-object-info indent-by obj)`**
—
Outputs some information about the given `obj`, indented by `indent-by` spaces.

### Output File Names

The output file names used by Meep, e.g. for HDF5 files, are automatically prefixed by the input variable `filename-prefix`. If `filename-prefix` is `""` (the default), however, then Meep constructs a default prefix based on the current ctl file name with `".ctl"` replaced by `"-"`: e.g. `test.ctl` implies a prefix of `"test-"`. You can get this prefix by running:

**`(get-filename-prefix)`**
—
Return the current prefix string that is prepended, by default, to all file names.

If you don't want to use any prefix, then you should set `filename-prefix` to `false`.

In addition to the filename prefix, you can also specify that all the output files be written into a newly-created directory (if it does not yet exist). This is done by running:

**`(use-output-directory [dirname])`**
—
Put output in a subdirectory, which is created if necessary. If the optional argument dirname is specified, that is the name of the directory. Otherwise, the directory name is the current ctl file name with `".ctl"` replaced by `"-out"`: e.g. `test.ctl` implies a directory of `"test-out"`.

### Output Volume

**`(volume (center ...) (size ...))`**
—
Many Meep functions require you to specify a volume in space, corresponding to the C++ type `meep::volume`. This function creates such a volume object, given the `center` and `size` properties (just like e.g. a `block` object). If the `size` is not specified, it defaults to `(0,0,0)`, i.e. a single point.

### Simulation Time

**`(meep-time)`**
—
Return the current simulation time in simulation time units (e.g. during a run function). This is not the wall-clock time.

Occasionally, e.g. for termination conditions of the form *time* &lt; *T*?, it is desirable to round the time to single precision in order to avoid small differences in roundoff error from making your results different by one timestep from machine to machine (a difference much bigger than roundoff error); in this case you can call `(meep-round-time)` instead, which returns the time rounded to single precision.

### Field Computations

Meep supports a large number of functions to perform computations on the fields. Most of them are accessed via the lower-level C++/SWIG interface. Some of them are based on the following simpler, higher-level versions.

**`(get-field-point c pt)`**
—
Given a `component` or `derived-component` constant `c` and a `vector3` `pt`, returns the value of that component at that point.

**`(get-epsilon-point pt)`**
—
Equivalent to `(get-field-point Dielectric pt)`.

**`(add-dft-fields cs freq-min freq-max nfreq [where])`**
—
Given a list of field components `cs`, compute the Fourier transform of these fields for `nfreq` equally spaced frequencies covering the frequency range `freq-min` to `freq-max` over the `volume` specified by `where` (default to the entire cell).

**`(flux-in-box dir box)`**
—
Given a `direction` constant, and a `meep::volume*`, returns the flux (the integral of $\Re [\mathbf{E}^* \times \mathbf{H}]$) in that volume. Most commonly, you specify a volume that is a plane or a line, and a direction perpendicular to it, e.g. `(flux-in-box `X (volume (center 0) (size 0 1 1)))`.

**`(electric-energy-in-box box)`**
—
Given a `meep::volume*`, returns the integral of the electric-field energy $\mathbf{E}^* \cdot \mathbf{D}/2$ in the given volume. If the volume has zero size along a dimension, a lower-dimensional integral is used.

**`(magnetic-energy-in-box box)`**
—
Given a `meep::volume*`, returns the integral of the magnetic-field energy $\mathbf{H}^* \cdot \mathbf{B}/2$ in the given volume. If the volume has zero size along a dimension, a lower-dimensional integral is used.

**`(field-energy-in-box box)`**
—
Given a `meep::volume*`, returns the integral of the electric- and magnetic-field energy $\mathbf{E}^* \cdot \mathbf{D}/2 + \mathbf{H}^* \cdot \mathbf{B}/2$in the given volume. If the volume has zero size along a dimension, a lower-dimensional integral is used.

Note that if you are at a fixed frequency and you use complex fields (via Bloch-periodic boundary conditions or `fields-complex?=true`), then one half of the flux or energy integrals above corresponds to the time average of the flux or energy for a simulation with real fields.

Often, you want the integration box to be the entire cell. A useful function to return this box, which you can then use for the `box` arguments above, is `(meep-fields-total-volume fields)`, where `fields` is the global variable (above) holding the current `meep::fields` object.

One versatile feature is that you can supply an arbitrary function $f(\mathbf{x},c_1,c_2,\ldots)$ of position $\mathbf{x}$ and various field components $c_1,\ldots$ and ask Meep to integrate it over a given volume, find its maximum, or output it (via `output-field-function`, described later). This is done via the functions:

**`(integrate-field-function cs func [where] [fields-var])`**
—
Returns the integral of the complex-valued function `func` over the `meep::volume` specified by `where` (defaults to entire cell) for the `meep::fields` specified by `fields-var` (defaults to `fields`). `func` is a function of position (a `vector3`, its first argument) and zero or more field components specified by `cs`: a list of `component` constants. `func` can be real- or complex-valued.

If any dimension of `where` is zero, that dimension is not integrated over. In this way you can specify 1d, 2d, or 3d integrals.

**`(max-abs-field-function cs func [where] [fields-var])`**
—
As `integrate-field-function`, but returns the maximum absolute value of `func` in the volume `where` instead of its integral.

The integration is performed by summing over the grid points with a simple trapezoidal rule, and the maximum is similarly over the grid points. See [Field Functions](Field_Functions.md) for examples of how to call `integrate-field-function` and `max-abs-field-function`. See [Synchronizing the Magnetic and Electric Fields](Synchronizing_the_Magnetic_and_Electric_Fields.md) if you want to do computations combining the electric and magnetic fields.

Occasionally, one wants to compute an integral that combines fields from two separate simulations (e.g. for nonlinear coupled-mode calculations). This functionality is supported in Meep, as long as the two simulations have the *same* cell, the same resolution, the same boundary conditions and symmetries (if any), and the same PML layers (if any).

**`(integrate2-field-function fields2 cs1 cs2 func [where] [fields-var])`**
—
Similar to `integrate-field-function`, but takes additional parameters `fields2` and `cs2`. `fields2` is a `meep::fields*` object similar to the global `fields` variable (see below) specifying the fields from another simulation. `cs1` is a list of components to integrate with from `fields-var` (defaults to `fields`), as for `integrate-field-function`, while `cs2` is a list of components to integrate from `fields2`. Similar to `integrate-field-function`, `func` is a function that returns an number given arguments consisting of: the position vector, followed by the values of the components specified by `cs1` (in order), followed by the values of the components specified by `cs2` (in order).

To get two fields in memory at once for `integrate2-field-function`, the easiest way is to run one simulation within a given Scheme (.ctl) file, then save the results in another fields variable, then run a second simulation. This would look something like:

```scm
...set up and run first simulation...
(define fields2 fields) ; save the fields in a variable
(set! fields '()) ; prevent the fields from getting deallocated by reset-meep
(reset-meep)
...set up and run second simulation...
```

It is also possible to timestep both fields simultaneously (e.g. doing one timestep of one simulation then one timestep of another simulation, and so on, but this requires you to call much lower-level functions like `(meep-fields-step fields)`.

### Reloading Parameters

Once the fields/simulation have been initialized, you can change the values of various parameters by using the following functions:

**`(reset-meep)`**
—
Reset all of Meep's parameters, deleting the fields, structures, etcetera, from memory as if you had not run any computations.

**`(restart-fields)`**
—
Restart the fields at time zero, with zero fields. Does *not* reset the Fourier transforms of the flux planes, which continue to be accumulated.

**`(change-k-point! k)`**
—
Change the `k-point` (the Bloch periodicity).

**`(change-sources! new-sources)`**
—
Change the `sources` input variable to `new-sources`, and changes the sources used for the current simulation.

### Flux Spectra

Given a bunch of [`flux-region`](#flux-region) objects, you can tell Meep to accumulate the Fourier transforms of the fields in those regions in order to compute the Poynting flux spectra. (Note: as a matter of convention, the "intensity" of the electromagnetic fields refers to the Poynting flux, *not* to the [energy density](#energy-density-spectra).) See also the [Introduction](Introduction.md#transmittancereflectance-spectra) and [Tutorial/Basics](Scheme_Tutorials/Basics.md#transmittance-spectrum-of-a-waveguide-bend). The most important function is:

**`(add-flux fcen df nfreq flux-regions...)`**
—
Add a bunch of `flux-region`s to the current simulation (initializing the fields if they have not yet been initialized), telling Meep to accumulate the appropriate field Fourier transforms for `nfreq` equally spaced frequencies covering the frequency range `fcen-df/2` to `fcen+df/2`. Return a *flux object*, which you can pass to the functions below to get the flux spectrum, etcetera.

As described in the tutorial, you normally use `add-flux` via statements like:

**`(define transmission (add-flux ...))`**
—
to store the flux object in a variable. `add-flux` initializes the fields if necessary, just like calling `run`, so you should only call it *after* setting up your `geometry`, `sources`, `pml-layers`, `k-point`, etcetera. You can create as many flux objects as you want, e.g. to look at powers flowing in different regions or in different frequency ranges. Note, however, that Meep has to store (and update at every time step) a number of Fourier components equal to the number of grid points intersecting the flux region multiplied by the number of electric and magnetic field components required to get the Poynting vector multiplied by `nfreq`, so this can get quite expensive (in both memory and time) if you want a lot of frequency points over large regions of space.

Once you have called `add-flux`, the Fourier transforms of the fields are accumulated automatically during time-stepping by the [run functions](#run-functions). At any time, you can ask for Meep to print out the current flux spectrum via:

**`(display-fluxes fluxes...)`**
—
Given a number of flux objects, this displays a comma-separated table of frequencies and flux spectra, prefixed by "flux1:" or similar (where the number is incremented after each run). All of the fluxes should be for the same `fcen`/`df`/`nfreq`. The first column are the frequencies, and subsequent columns are the flux spectra.

You might have to do something lower-level if you have multiple flux regions corresponding to *different* frequency ranges, or have other special needs. `(display-fluxes f1 f2 f3)` is actually equivalent to `(display-csv "flux" (get-flux-freqs f1) (get-fluxes f1) (get-fluxes f2) (get-fluxes f3))`, where `display-csv` takes a bunch of lists of numbers and prints them as a comma-separated table; this involves calling two lower-level functions:

**`(get-flux-freqs flux)`**
—
Given a flux object, returns a list of the frequencies that it is computing the spectrum for.

**`(get-fluxes flux)`**
—
Given a flux object, returns a list of the current flux spectrum that it has accumulated.

As described in [Tutorial/Basics](Scheme_Tutorials/Basics.md), for a reflection spectrum you often want to save the Fourier-transformed fields from a "normalization" run and then load them into another run to be subtracted. This can be done via:

**`(save-flux filename flux)`**
—
Save the Fourier-transformed fields corresponding to the given flux object in an HDF5 file of the given `filename` without the ".h5" suffix (the current filename-prefix is prepended automatically).

**`(load-flux filename flux)`**
—
Load the Fourier-transformed fields into the given flux object (replacing any values currently there) from an HDF5 file of the given `filename` without the ".h5" suffix (the current filename-prefix is prepended automatically). You must load from a file that was saved by `save-flux` in a simulation of the same dimensions (for both the cell and the flux regions) with the same number of processors.

**`(load-minus-flux filename flux)`**
—
As `load-flux`, but negates the Fourier-transformed fields after they are loaded. This means that they will be *subtracted* from any future field Fourier transforms that are accumulated.

**`(scale-flux-fields s flux)`**
—
Scale the Fourier-transformed fields in `flux` by the complex number `s`. e.g. `load-minus-flux` is equivalent to `load-flux` followed by `scale-flux-fields` with `s=-1`.

### Mode Decomposition

Given a structure, Meep can decompose the Fourier-transformed fields into a superposition of its harmonic modes. For a theoretical background, see [Mode Decomposition](Mode_Decomposition.md).

**`(get-eigenmode-coefficients flux bands eig-parity eig-vol eig-resolution eig-tolerance kpoint-func verbose=False)`**
—
Given a flux object and list of band indices, return a list with the following data:

+ `alpha`: the complex eigenmode coefficients as a 3d Guile [array](https://www.gnu.org/software/guile/manual/html_node/Arrays.html#Arrays) of size (`(length bands)`, `flux.Nfreq`, `2`). The last/third dimension refers to modes propagating in the forward (+) or backward (-) directions.
+ `vgrp`: the group velocity as a Guile array.
+ `kpoints`: a list of `vector3`s of the `kpoint` used in the mode calculation.
+ `kdom`: a list of `vector3`s of the mode's dominant wavevector.

Here is an example of calling `get-eigenmode-coefficients` overriding the default `eig-parity` with a keyword argument, and then printing the coefficient for first band, first frequency, and forward direction:

```scheme
(let ((result (get-eigenmode-coefficients flux (list 1) #:eig-parity (+ ODD-Z EVEN-Y))))
    (print (array-ref (list-ref result 0) 0 0 0)))
```
The flux object must be created using `add-mode-monitor` (an alias for `add-flux`). `eig-vol` is the volume passed to [MPB](https://mpb.readthedocs.io) for the eigenmode calculation (based on interpolating the discretized materials from the Yee grid); in most cases this will simply be the volume over which the frequency-domain fields are tabulated, which is the default (i.e. `(meep-dft-flux-where-get flux)`). `eig-parity` should be one of [`NO-PARITY` (default), `EVEN-Z`, `ODD-Z`, `EVEN-Y`, `ODD-Y`]. It is the parity (= polarization in 2d) of the mode to calculate, assuming the structure has $z$ and/or $y$ mirror symmetry *in the source region*, just as for `eigenmode-source` above. If the structure has both $y$ and $z$ mirror symmetry, you can combine more than one of these, e.g. `(+ EVEN-Z ODD-Y)`. Default is `NO-PARITY`, in which case MPB computes all of the bands which will still be even or odd if the structure has mirror symmetry, of course. This is especially useful in 2d simulations to restrict yourself to a desired polarization. `eig-resolution` is the spatial resolution to use in MPB for the eigenmode calculations. This defaults to the same resolution as Meep, but you can use a higher resolution in which case the structure is linearly interpolated from the Meep pixels. `eig-tolerance` is the tolerance to use in the MPB eigensolver. MPB terminates when the eigenvalues stop changing to less than this fractional tolerance. Defaults to `1e-12`.  (Note that this is the tolerance for the frequency eigenvalue ω; the tolerance for the mode profile is effectively the square root of this.)

Technically, MPB computes `ωₙ(k)` and then inverts it with Newton's method to find the wavevector `k` normal to `eig-vol` and mode for a given frequency; in rare cases (primarily waveguides with *nonmonotonic* dispersion relations, which doesn't usually happen in simple dielectric waveguides), MPB may need you to supply an initial "guess" for `k` in order for this Newton iteration to converge.  You can supply this initial guess with `kpoint-func`, which is a function `(kpoint-func f n)` that supplies a rough initial guess for the `k` of band number `n` at frequency `f = ω/2π`.  (By default, the **k** components in the plane of the `eig-vol` region are zero.  However, if this region spans the *entire* cell in some directions, and the cell has Bloch-periodic boundary conditions via the `k-point` parameter, then the mode's **k** components in those directions will match `k-point` so that the mode satisfies the Meep boundary conditions, regardless of `kpoint-func`.)

**Note:** for planewaves in homogeneous media, the `kpoints` may *not* necessarily be equivalent to the actual wavevector of the mode. This quantity is given by `kdom`.

**`(add_mode_monitor fcen df nfreq ModeRegions...)`**
—
Similar to `add-flux`, but for use with `get-eigenmode-coefficients`.

`add-mode-monitor` works properly with arbitrary symmetries, but may be suboptimal because the Fourier-transformed region does not exploit the symmetry.  As an optimization, if you have a mirror plane that bisects the mode monitor, you can instead use `add-flux` to gain a factor of two, but in that case you *must* also pass the corresponding `eig-parity` to `get-eigenmode-coefficients` in order to only compute eigenmodes with the corresponding mirror symmetry.

### Energy Density Spectra

Very similar to flux spectra, you can also compute **energy density spectra**: the energy density of the electromagnetic fields as a function of frequency, computed by Fourier transforming the fields and integrating the energy density:

$$ \frac{1}{2}ε|\mathbf{E}|^2 + \frac{1}{2}μ|\mathbf{H}|^2 $$

The usage is similar to the flux spectra: you define a set of `energy-region` objects telling Meep where it should compute the Fourier-transformed fields and energy densities, and call `add-energy` to add these regions to the current simulation over a specified frequency bandwidth, and then use `display-electric-energy`, `display-magnetic-energy`, or `display-total-energy` to display the energy density spectra at the end. There are also `save-energy`, `load-energy`, and `load-minus-energy` functions that you can use to subtract the fields from two simulation, e.g. in order to compute just the energy from scattered fields, similar to the flux spectra. These types and functions are defined as follows:

**`energy-region`**

A region (volume, plane, line, or point) in which to compute the integral of the energy density of the Fourier-transformed fields. Its properties are:

**`center` [`vector3`]**
—
The center of the energy region (no default).

**`size` [`vector3`]**
—
The size of the energy region along each of the coordinate axes. Default is (0,0,0): a single point.

**`weight` [`cnumber`]**
—
A weight factor to multiply the energy density by when it is computed. Default is 1.0.

**`(add-energy fcen df nfreq energy-regions...)`**
—
Add a bunch of `energy-region`s to the current simulation (initializing the fields if they have not yet been initialized), telling Meep to accumulate the appropriate field Fourier transforms for `nfreq` equally spaced frequencies covering the frequency range `fcen-df/2` to `fcen+df/2`. Return an *energy object*, which you can pass to the functions below to get the energy spectrum, etcetera.

As for energy regions, you normally use `add-energy` via statements like:

```scm
(define En (add-energy ...))
```

to store the energy object in a variable. `add-energy` initializes the fields if necessary, just like calling `run`, so you should only call it *after* setting up your `geometry`, `sources`, `pml-layers`, `k-point`, etcetera. You can create as many energy objects as you want, e.g. to look at the energy densities in different objects or in different frequency ranges. Note, however, that Meep has to store (and update at every time step) a number of Fourier components equal to the number of grid points intersecting the energy region multiplied by `nfreq`, so this can get quite expensive (in both memory and time) if you want a lot of frequency points over large regions of space.

Once you have called `add-energy`, the Fourier transforms of the fields are accumulated automatically during time-stepping by the `run` functions. At any time, you can ask for Meep to print out the current energy density spectrum via:

**`(display-electric-energy energy...)`, `(display-magnetic-energy energy...)`, `(display-total-energy energy...)` **
—
Given a number of energy objects, this displays a comma-separated table of frequencies and energy density spectra for the electric, magnetic and total fields, respectively prefixed by "electric-energy1:", "magnetic-energy1:," "total-energy1:," or similar (where the number is incremented after each run). All of the energy should be for the same `fcen`/`df`/`nfreq`. The first column are the frequencies, and subsequent columns are the energy density spectra.

You might have to do something lower-level if you have multiple energy regions corresponding to *different* frequency ranges, or have other special needs. `(display-electric-energy e1 e2 e3)` is actually equivalent to `(display-csv "electric-energy" (get-energy-freqs e1) (get-electric-energy e1) (get-electric-energy e2) (get-electric-energy e3))`, where `display-csv` takes a bunch of lists of numbers and prints them as a comma-separated table; this involves calling two lower-level functions:

**`(get-energy-freqs energy)`**
—
Given an energy object, returns a list of the frequencies that it is computing the spectrum for.

**`(get-electric-energy energy)`, `(get-magnetic-energy energy)`, `(get-total-energy energy)`**
—
Given an energy object, returns a list of the current energy density spectrum for the electric, magnetic, or total fields, respectively that it has accumulated.

As described in [Tutorial/Basics](Scheme_Tutorials/Basics.md), to compute the energy density from the scattered fields you often want to save the Fourier-transformed fields from a "normalization" run and then load them into another run to be subtracted. This can be done via:

**`(save-energy filename energy)`**
—
Save the Fourier-transformed fields corresponding to the given energy object in an HDF5 file of the given `filename` without the ".h5" suffix (the current filename-prefix is prepended automatically).

**`(load-energy filename energy)`**
—
Load the Fourier-transformed fields into the given energy object (replacing any values currently there) from an HDF5 file of the given `filename` without the ".h5" suffix (the current filename-prefix is prepended automatically). You must load from a file that was saved by `save-energy` in a simulation of the same dimensions for both the cell and the energy regions with the same number of processors.

**`(load-minus-energy filename energy)`**
—
As `load-energy`, but negates the Fourier-transformed fields after they are loaded. This means that they will be *subtracted* from any future field Fourier transforms that are accumulated.

### Force Spectra

Very similar to flux spectra, you can also compute **force spectra**: forces on an object as a function of frequency, computed by Fourier transforming the fields and integrating the vacuum [Maxwell stress tensor](https://en.wikipedia.org/wiki/Maxwell_stress_tensor):

$$σ_{ij} = E_i^*E_j + H_i^*H_j - \frac{1}{2} δ_{ij} \left( |\mathbf{E}|^2 + |\mathbf{H}|^2 \right)$$

over a surface $S$ via $\mathbf{F} = \int_S σ d\mathbf{A}$. You should normally **only evaluate the stress tensor over a surface lying in vacuum**, as the interpretation and definition of the stress tensor in arbitrary media is often problematic (the subject of extensive and controversial literature). It is fine if the surface *encloses* an object made of arbitrary materials, as long as the surface itself is in vacuum.

See also [Tutorial/Optical Forces](Scheme_Tutorials/Optical_Forces.md).

Most commonly, you will want to **normalize** the force spectrum in some way, just as for flux spectra. Most simply, you could divide two different force spectra to compute the ratio of forces on two objects. Often, you will divide a force spectrum by a flux spectrum, to divide the force $F$ by the incident power $P$ on an object, in order to compute the useful dimensionless ratio $Fc$/$P$ where $c=1$ in Meep units. For example, it is a simple exercise to show that the force $F$ on a perfectly reflecting mirror with normal-incident power $P$ satisfies $Fc$/$P=2$, and for a perfectly absorbing (black) surface $Fc$/$P=1$.

The usage is similar to the flux spectra: you define a set of `force-region` objects telling Meep where it should compute the Fourier-transformed fields and stress tensors, and call `add-force` to add these regions to the current simulation over a specified frequency bandwidth, and then use `display-forces` to display the force spectra at the end. There are also `save-force`, `load-force`, and `load-minus-force` functions that you can use to subtract the fields from two simulation, e.g. in order to compute just the force from scattered fields, similar to the flux spectra. These types and functions are defined as follows:

**`force-region`**

A region (volume, plane, line, or point) in which to compute the integral of the stress tensor of the Fourier-transformed fields. Its properties are:

**`center` [`vector3`]**
—
The center of the force region (no default).

**`size` [`vector3`]**
—
The size of the force region along each of the coordinate axes. Default is (0,0,0): a single point.

**`direction` [`direction constant`]**
—
The direction of the force that you wish to compute (e.g. `X`, `Y`, etcetera). Unlike `flux-region`, you must specify this explicitly, because there is not generally any relationship between the direction of the force and the orientation of the force region.

**`weight` [`cnumber`]**
—
A weight factor to multiply the force by when it is computed. Default is 1.0.

In most circumstances, you should define a set of `force-region`s whose union is a closed surface lying in vacuum and enclosing the object that is experiencing the force.

**`(add-force fcen df nfreq force-regions...)`**
—
Add a bunch of `force-region`s to the current simulation (initializing the fields if they have not yet been initialized), telling Meep to accumulate the appropriate field Fourier transforms for `nfreq` equally spaced frequencies covering the frequency range `fcen-df/2` to `fcen+df/2`. Return a *force object*, which you can pass to the functions below to get the force spectrum, etcetera.

As for force regions, you normally use `add-force` via statements like:

```scm
(define Fx (add-force ...))
```

to store the force object in a variable. `add-force` initializes the fields if necessary, just like calling `run`, so you should only call it *after* setting up your `geometry`, `sources`, `pml-layers`, etcetera. You can create as many force objects as you want, e.g. to look at forces on different objects, in different directions, or in different frequency ranges. Note, however, that Meep has to store (and update at every time step) a number of Fourier components equal to the number of grid points intersecting the force region, multiplied by the number of electric and magnetic field components required to get the stress vector, multiplied by `nfreq`, so this can get quite expensive (in both memory and time) if you want a lot of frequency points over large regions of space.

Once you have called `add-force`, the Fourier transforms of the fields are accumulated automatically during time-stepping by the `run` functions. At any time, you can ask for Meep to print out the current force spectrum via:

**`(display-forces forces...)`**
—
Given a number of force objects, this displays a comma-separated table of frequencies and force spectra, prefixed by "force1:" or similar (where the number is incremented after each run). All of the forces should be for the same `fcen`/`df`/`nfreq`. The first column are the frequencies, and subsequent columns are the force spectra.

You might have to do something lower-level if you have multiple force regions corresponding to *different* frequency ranges, or have other special needs. `(display-forces f1 f2 f3)` is actually equivalent to `(display-csv "force" (get-force-freqs f1) (get-forces f1) (get-forces f2) (get-forces f3))`, where `display-csv` takes a bunch of lists of numbers and prints them as a comma-separated table; this involves calling two lower-level functions:

**`(get-force-freqs force)`**
—
Given a force object, returns a list of the frequencies that it is computing the spectrum for.

**`(get-forces force)`**
—
Given a force object, returns a list of the current force spectrum that it has accumulated.

As described in [Tutorial/Basics](Scheme_Tutorials/Basics.md), to compute the force from scattered fields you often want to save the Fourier-transformed fields from a "normalization" run and then load them into another run to be subtracted. This can be done via:

**`(save-force filename force)`**
—
Save the Fourier-transformed fields corresponding to the given force object in an HDF5 file of the given `filename` without the ".h5" suffix (the current filename-prefix is prepended automatically).

**`(load-force filename force)`**
—
Load the Fourier-transformed fields into the given force object (replacing any values currently there) from an HDF5 file of the given `filename` without the ".h5" suffix (the current filename-prefix is prepended automatically). You must load from a file that was saved by `save-force` in a simulation of the same dimensions for both the cell and the force regions with the same number of processors.

**`(load-minus-force filename force)`**
—
As `load-force`, but negates the Fourier-transformed fields after they are loaded. This means that they will be *subtracted* from any future field Fourier transforms that are accumulated.

### LDOS spectra

Meep can also calculate the LDOS (local density of states) spectrum, as described in [Tutorial/Local Density of States](Scheme_Tutorials/Local_Density_of_States.md). To do this, you simply pass the following step function to your `run` command:

**`(dft-ldos fcen df nfreq)`**
—
Compute the power spectrum of the sources (usually a single point dipole source), normalized to correspond to the LDOS, in a frequency bandwidth `df` centered at `fcen`, at `nfreq` frequency points.

**`(get-ldos-freqs ldos)`**
—
Given an ldos object, returns a list of the frequencies that it is computing the spectrum for.

The resulting spectrum is outputted as comma-delimited text, prefixed by `ldos:,`, and is also stored in the `dft-ldos-data` global variable after the `run` is complete.

Analytically, the per-polarization LDOS is exactly proportional to the power radiated by an $\ell$-oriented point-dipole current, $p(t)$, at a given position in space. For a more mathematical treatment of the theory behind the LDOS, refer to the relevant discussion in Section 4.4 ("Currents and Fields: The Local Density of States") in [Chapter 4](http://arxiv.org/abs/arXiv:1301.5366) ("Electromagnetic Wave Source Conditions") of the book [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707), but for now it is defined as:

$$\operatorname{LDOS}_{\ell}(\vec{x}_0,\omega)=-\frac{2}{\pi}\varepsilon(\vec{x}_0)\frac{\operatorname{Re}[\hat{E}_{\ell}(\vec{x}_0,\omega)\hat{p}(\omega)^*]}{|\hat{p}(\omega)|^2}$$

where the $|\hat{p}(\omega)|^2$ normalization is necessary for obtaining the power exerted by a unit-amplitude dipole (assuming linear materials), and hats denote Fourier transforms. It is this quantity that is computed by the `dft-ldos` command for a single dipole source. For a volumetric source, the numerator and denominator are both integrated over the current volume, but "LDOS" computation is less meaningful in this case.

### Near-to-Far-Field Spectra

Meep can compute a near-to-far-field transformation in the frequency domain as described in [Tutorial/Near-to-Far Field Spectra](Scheme_Tutorials/Near_to_Far_Field_Spectra.md): given the fields on a "near" bounding surface inside the cell, it can compute the fields arbitrarily far away using an analytical transformation, assuming that the "near" surface and the "far" region lie in a single homogeneous non-periodic 2d or 3d region. That is, in a simulation *surrounded by PML* that absorbs outgoing waves, the near-to-far-field feature can compute the fields outside the cell as if the outgoing waves had not been absorbed (i.e. in the fictitious infinite open volume). Moreover, this operation is performed on the Fourier-transformed fields: like the flux and force spectra above, you specify a set of desired frequencies, Meep accumulates the Fourier transforms, and then Meep computes the fields at *each frequency* for the desired far-field points.

This is based on the principle of equivalence: given the Fourier-transformed tangential fields on the "near" surface, Meep computes equivalent currents and convolves them with the analytical Green's functions in order to compute the fields at any desired point in the "far" region. For details, see Section 4.2.1 ("The Principle of Equivalence") in [Chapter 4](http://arxiv.org/abs/arXiv:1301.5366) ("Electromagnetic Wave Source Conditions") of the book [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707).

Note: in order for the far-field results to be accurate, the [far region must be separated from the near region](https://en.wikipedia.org/wiki/Near_and_far_field) by *at least* 2D<sup>2</sup>/λ, the Fraunhofer distance, where D is the largest dimension of the radiator and λ is the vacuum wavelength.

There are three steps to using the near-to-far-field feature: first, define the "near" surface(s) as a set of surfaces capturing *all* outgoing radiation in the desired direction(s); second, run the simulation, typically with a pulsed source, to allow Meep to accumulate the Fourier transforms on the near surface(s); third, tell Meep to compute the far fields at any desired points (optionally saving the far fields from a grid of points to an HDF5 file). To define the near surfaces, use:

**`(add-near2far fcen df nfreq near2far-regions... nperiods)`**
—
Add a bunch of `near2far-region`s to the current simulation (initializing the fields if they have not yet been initialized), telling Meep to accumulate the appropriate field Fourier transforms for `nfreq` equally-spaced frequencies covering the frequency range `fcen-df/2` to `fcen+df/2`. Return a `near2far` object, which you can pass to the functions below to get the far fields. `nperiods` is a keyword argument that defaults to one, and can be passed after the list of `near2far-regions` like so: `(add-near2far fcen df nfreq region1 region2 region3 #:nperiods 2)`

Each `near2far-region` is identical to `flux-region` except for the name: in 3d, these give a set of planes (**important:** all these "near surfaces" must lie in a single *homogeneous* material with *isotropic* ε and μ &mdash; and they should *not* lie in the PML regions) surrounding the source(s) of outgoing radiation that you want to capture and convert to a far field. Ideally, these should form a closed surface, but in practice it is sufficient for the `near2far-region`s to capture all of the radiation in the direction of the far-field points. **Important:** as for flux computations, each `near2far-region` should be assigned a `weight` of &#177;1 indicating the direction of the outward normal relative to the +coordinate direction. So, for example, if you have six regions defining the six faces of a cube, i.e. the faces in the +x, -x, +y, -y, +z, and -z directions, then they should have weights +1, -1, +1, -1, +1, and -1 respectively. Note that, neglecting discretization errors, all near-field surfaces that enclose the same outgoing fields are equivalent and will yield the same far fields with a discretization-induced difference that vanishes with increasing resolution etc.

After the simulation run is complete, you can compute the far fields. This is usually for a pulsed source so that the fields have decayed away and the Fourier transforms have finished accumulating.

**`(get_farfield near2far x)`**
—
Given a `vector3` point `x` which can lie anywhere outside the near-field surface, including outside the cell and a near2far object, returns the computed (Fourier-transformed) "far" fields at `x` as list of length 6`nfreq`, consisting of fields (E<sub>x</sub><sup>1</sup>,E<sub>y</sub><sup>1</sup>,E<sub>z</sub><sup>1</sup>,H<sub>x</sub><sup>1</sup>,H<sub>y</sub><sup>1</sup>,H<sub>z</sub><sup>1</sup>,E<sub>x</sub><sup>2</sup>,E<sub>y</sub><sup>2</sup>,E<sub>z</sub><sup>2</sup>,H<sub>x</sub><sup>2</sup>,H<sub>y</sub><sup>2</sup>,H<sub>z</sub><sup>2</sup>,...) for the frequencies 1,2,…,`nfreq`.

**`(get-near2far-freqs near2far)`**
—
Given a `near2far` object, returns a list of the frequencies that it is computing the spectrum for.

**`(output-farfields near2far fname where resolution)`**
—
Given an HDF5 file name `fname` (does *not* include the `.h5` suffix), a `volume` given by `where` (may be 0d, 1d, 2d, or 3d), and a `resolution` (in grid points / distance unit), outputs the far fields in `where` (which may lie *outside* the cell) in a grid with the given resolution (which may differ from the FDTD grid resolution) to the HDF5 file as a set of twelve array datasets `ex.r`, `ex.i`, ..., `hz.r`, `hz.i`, giving the real and imaginary parts of the Fourier-transformed $E$ and $H$ fields on this grid. Each dataset is an nx&#215;ny&#215;nz&#215;nfreq 4d array of space&#215;frequency although dimensions that =1 are omitted.

Note that far fields have the same units and scaling as the *Fourier transforms* of the fields, and hence cannot be directly compared to time-domain fields. In practice, it is easiest to use the far fields in computations where overall scaling (units) cancel out or are irrelevant, e.g. to compute the fraction of the far fields in one region vs. another region.

(Multi-frequency `output-farfields` can be accelerated by
[compiling Meep](Build_From_Source.md#meep) with `--with-openmp` and using the
`OMP_NUM_THREADS` environment variable to specify multiple threads.)

For a scattered-field computation, you often want to separate the scattered and incident fields. Just as is described in [Tutorial/Basics/Transmittance Spectrum of a Waveguide Bend](Scheme_Tutorials/Basics.md#transmittance-spectrum-of-a-waveguide-bend) for flux computations, you can do this by saving the Fourier-transformed incident from a "normalization" run and then load them into another run to be subtracted. This can be done via:

**`(save-near2far filename near2far)`**
—
Save the Fourier-transformed fields corresponding to the given `near2far` object in an HDF5 file of the given `filename` (without the ".h5" suffix). The current filename-prefix is prepended automatically.

**`(load-near2far filename near2far)`**
—
Load the Fourier-transformed fields into the given `near2far` object replacing any values currently there from an HDF5 file of the given `filename` (without the ".h5" suffix) the current filename-prefix is prepended automatically. You must load from a file that was saved by `save-near2far` in a simulation of *the same dimensions* for both the cell and the near2far regions with the same number of processors.

**`(load-minus-near2far filename near2far)`**
—
As `load-near2far`, but negates the Fourier-transformed fields after they are loaded. This means that they will be *subtracted* from any future field Fourier transforms that are accumulated.

**`(scale-near2far-fields s near2far)`**
—
Scale the Fourier-transformed fields in `near2far` by the complex number `s`. e.g. `load-minus-near2far` is equivalent to `load-near2far` followed by `scale-near2far-fields` with `s=-1`.

**`(flux near2far direction where resolution)`**
—
Given a `volume` `where` (may be 0d, 1d, 2d, or 3d) and a `resolution` (in grid points / distance unit), compute the far fields in `where` (which may lie *outside* the cell) in a grid with the given resolution (which may differ from the FDTD solution) and return its Poynting flux in `direction` as a list. The dataset is a 1d array of `nfreq` dimensions.

### Load and Dump Structure

These functions dump the raw ε data to disk and load it back for doing multiple simulations with the same materials but different sources etc. The only prerequisite is that the dump/load simulations have the same [chunks](Chunks_and_Symmetry.md) (i.e. the same grid, number of processors, and PML). Currently only stores ε and μ, and not nonlinear coefficients or polarizability.

**`(meep-structure-dump structure fname)`**
—
Dumps the structure to the file `fname` using the global `structure` object
(which is initialized after you execute `run` or `init-structure`).

**`(meep-structure-load structure fname)`**
—
Loads a structure from the file `fname`.   This should be called after
`(init-structure)` so that the global `structure` object is initialized,
and you should generally `(set! geometry '())` to skip initializing the
geometry (since it will be overwritten by `meep-structure-load` anyway).

### Frequency-Domain Solver

Meep contains a frequency-domain solver that computes the fields produced in a geometry in response to a [continuous-wave (CW) source](https://en.wikipedia.org/wiki/Continuous_wave). This is based on an [iterative linear solver](https://en.wikipedia.org/wiki/Iterative_method) instead of time-stepping. For details, see Section 5.3 ("Frequency-domain solver") of [Computer Physics Communications, Vol. 181, pp. 687-702, 2010](http://ab-initio.mit.edu/~oskooi/papers/Oskooi10.pdf). Benchmarking results have shown that in many instances, such as cavities (e.g., ring resonators) with long-lived resonant modes, this solver converges much faster than simply running an equivalent time-domain simulation with a CW source, time-stepping until all transient effects from the source turn-on have disappeared, especially if the fields are desired to a very high accuracy. To use it, simply define a `continuous-src` with the desired frequency, [initialize the fields and geometry](#initializing-the-structure-and-fields) via `(init-fields)`, and then:

**`(meep-fields-solve-cw fields tol maxiters L)`**

After the `fields` variable (a global variable pointing to the `meep::fields*` object initialized by `init-fields`, see [Input Variables](Scheme_User_Interface.md#input-variables)), the next two parameters to the frequency-domain solver are the tolerance `tol` for the iterative solver (10<sup>−8</sup>, by default) and a maximum number of iterations `maxiters` (10<sup>4</sup>, by default). Finally, there is a parameter $L$ that determines a tradeoff between memory and work per step and convergence rate of the iterative algorithm, biconjugate gradient stabilized ([BiCGSTAB-L](https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method)), that is used; larger values of $L$ will often lead to faster convergence at the expense of more memory and more work per iteration. Default is $L=2$, and normally a value ≥ 2 should be used.

The frequency-domain solver supports arbitrary geometries, PML, boundary conditions, symmetries, parallelism, conductors, and arbitrary nondispersive materials. Lorentz-Drude dispersive materials are not currently supported in the frequency-domain solver, but since you are solving at a known fixed frequency rather than timestepping, you should be able to pick conductivities etcetera in order to obtain any desired complex ε and μ at that frequency.

The frequency-domain solver requires you to use complex-valued fields, via `(set! force-complex-fields? true)`.

After `meep-fields-solve-cw` completes, it should be as if you had just run the simulation for an infinite time with the source at that frequency. You can call the various field-output functions and so on as usual at this point.

Run and Step Functions
----------------------

The actual work in Meep is performed by `run` functions, which time-step the simulation for a given amount of time or until a given condition is satisfied.

The run functions, in turn, can be modified by use of [step functions](#predefined-step-functions): these are called at every time step and can perform any arbitrary computation on the fields, do outputs and I/O, or even modify the simulation. The step functions can be transformed by many [modifier functions](#step-function-modifiers), like `at-beginning`, `during-sources`, etcetera which cause them to only be called at certain times, etcetera, instead of at every time step.

A common point of confusion is described in [The Run Function Is Not A Loop](The_Run_Function_Is_Not_A_Loop.md). Read this article if you want to make Meep do some customized action on each time step, as many users make the same mistake. What you really want to in that case is to write a step function, as described below.

### Run Functions

The following run functions are available. You can also write your own, using the lower-level [C++/SWIG functions](#swig-wrappers), but these should suffice for most needs.

**`(run-until cond?/time step-functions...)`**
—
Run the simulation until a certain time or condition, calling the given step functions (if any) at each timestep. The first argument is *either* a number, in which case it is an additional time (in Meep units) to run for, *or* it is a function (of no arguments) which returns `true` when the simulation should stop.

**`(run-sources step-functions...)`**
—
Run the simulation until all sources have turned off, calling the given step functions (if any) at each timestep. Note that this does *not* mean that the fields will be zero at the end: in general, some fields will still be bouncing around that were excited by the sources.

**`(run-sources+ cond?/time step-functions...)`**
—
As `run-sources`, but with an additional first argument: either a number, in which case it is an *additional* time (in Meep units) to run for after the sources are off, *or* it is a function (of no arguments). In the latter case, the simulation runs until the sources are off *and* `(cond?)` returns `true`.

In particular, a useful first argument to `run-sources+` or `run-until` is often as shown below which is demonstrated in [Tutorial/Basics](Scheme_Tutorials/Basics.md):

**`(stop-when-fields-decayed dT c pt decay-by)`**
—
Return a `cond?` function, suitable for passing to `run-until`/`run-sources+`, that examines the component `c` (e.g. `Ex`, etc.) at the point `pt` (a `vector3`) and keeps running until its absolute value *squared* has decayed by at least `decay-by` from its maximum previous value. In particular, it keeps incrementing the run time by `dT` (in Meep units) and checks the maximum value over that time period &mdash; in this way, it won't be fooled just because the field happens to go through 0 at some instant.

Note that, if you make `decay-by` very small, you may need to increase the `cutoff` property of your source(s), to decrease the amplitude of the small high-frequency components that are excited when the source turns off. High frequencies near the [Nyquist frequency](https://en.wikipedia.org/wiki/Nyquist_frequency) of the grid have slow group velocities and are absorbed poorly by [PML](Perfectly_Matched_Layer.md).

Finally, another two run functions, useful for computing ω(**k**) band diagrams, are

**`(run-k-point T k)`**
—
Given a `vector3 k`, runs a simulation for each *k* point (i.e. specifying Bloch-periodic boundary conditions) and extracts the eigen-frequencies, and returns a list of the complex frequencies. In particular, you should have specified one or more Gaussian sources. It will run the simulation until the sources are turned off plus an additional $T$ time units. It will run [Harminv](#harminv) at the same point/component as the first Gaussian source and look for modes in the union of the frequency ranges for all sources.

**`(run-k-points T k-points)`**
—
Given a list `k-points` of *k* vectors, runs `run-k-point` for each one, and returns a list of lists of frequencies (one list of frequencies for each *k*). Also prints out a comma-delimited list of frequencies, prefixed by `freqs:`, and their imaginary parts, prefixed by `freqs-im:`. See [Tutorial/Resonant Modes and Transmission in a Waveguide Cavity](Scheme_Tutorials/Resonant_Modes_and_Transmission_in_a_Waveguide_Cavity.md).

### Predefined Step Functions

Several useful step functions are predefined by Meep.

#### Output Functions

The most common step function is an output function, which outputs some field component to an [HDF5](https://en.wikipedia.org/wiki/HDF5) file. Normally, you will want to modify this by one of the `at-*` functions, below, as outputting a field at *every* time step can get quite time- and storage-consuming.

Note that although the various field components are stored at different places in the [Yee lattice](Yee_Lattice.md), when they are outputted they are all linearly interpolated to the same grid: to the points at the *centers* of the Yee cells, i.e. $(i+0.5,j+0.5,k+0.5)\cdotΔ$ in 3d.

The predefined output functions are:

**`output-epsilon`**
—
Output the dielectric function (relative permittivity) ε. Note that this only outputs the real, frequency-independent part of ε (the $\omega\to\infty$ limit).

**`output-mu`**
—
Output the relative permeability function μ. Note that this only outputs the real, frequency-independent part of μ (the $\omega\to\infty$ limit).

**`(output-dft dft-fields fname [where])`**
—
Output the Fourier-transformed fields in `dft-fields` (created by `add-dft-fields`) to an HDF5 file with name `fname` (does *not* include the `.h5` suffix). The `volume` `where` defaults to the entire cell.

**`output-poynting`**
—
Output the Poynting flux $\mathrm{Re}\{\mathbf{E}^*\times\mathbf{H}\}$. Note that you might want to wrap this step function in `synchronized-magnetic` to compute it more accurately. See [Synchronizing the Magnetic and Electric Fields](Synchronizing_the_Magnetic_and_Electric_Fields.md).

**`output-hpwr`**
—
Output the magnetic-field energy density $\mathbf{H}^* \cdot \mathbf{B} / 2$

**`output-dpwr`**
—
Output the electric-field energy density $\mathbf{E}^* \cdot \mathbf{D} / 2$

**`output-tot-pwr`**
—
Output the total electric and magnetic energy density. Note that you might want to wrap this step function in `synchronized-magnetic` to compute it more accurately. See [Synchronizing the Magnetic and Electric Fields](Synchronizing_the_Magnetic_and_Electric_Fields.md).

**`output-Xfield-x, output-Xfield-y, output-Xfield-z, output-Xfield-r, output-Xfield-p`**
—
Output the $x$, $y$, $z$, $r$, or $\phi$ component respectively, of the field *X*, where *X* is either `h`, `b`, `e`, `d`, or `s` for the magnetic, electric, displacement, or Poynting flux, respectively. If the field is complex, outputs two datasets, e.g. `ex.r` and `ex.i`, within the same HDF5 file for the real and imaginary parts, respectively. Note that for outputting the Poynting flux, you might want to wrap the step function in `synchronized-magnetic` to compute it more accurately. See [Synchronizing the Magnetic and Electric Fields](Synchronizing_the_Magnetic_and_Electric_Fields.md).

**`output-Xfield`**
—
Outputs *all* the components of the field *X*, where *X* is either `h`, `b`, `e`, `d`, or `s` as above, to an HDF5 file. That is, the different components are stored as different datasets within the *same* file.

**`(output-png component h5topng-options)`**
—
Output the given field component (e.g. `Ex`, etc.) as a [PNG](https://en.wikipedia.org/wiki/PNG) image, by first outputting the HDF5 file, then converting to PNG via [h5topng](https://github.com/NanoComp/h5utils/blob/master/README.md), then deleting the HDF5 file. The second argument is a string giving options to pass to h5topng (e.g. `"-Zc bluered"`). See also [Tutorial/Basics](Scheme_Tutorials/Basics.md#output-tips-and-tricks).

It is often useful to use the `h5topng` `-C` or `-A` options to overlay the dielectric function when outputting fields. To do this, you need to know the name of the dielectric-function `.h5` file which must have been previously output by `output-epsilon`. To make this easier, a built-in shell variable `$EPS` is provided which refers to the last-output dielectric-function `.h5` file. So, for example `(output-png Ez "-C $EPS")` will output the $E_z$ field and overlay the dielectric contours.

**`(output-png+h5 component h5topng-options)`**
—
Like `output_png`, but also outputs the `.h5` file for the component. In contrast, `output_png` deletes the `.h5` when it is done.

More generally, it is possible to output an arbitrary function of position and zero or more field components, similar to the `integrate-field-function` described above. This is done by:

**`(output-field-function name cs func)`**
—
Output the field function `func` to an HDF5 file in the datasets named *`name`*`.r` and *`name`*`.i` for the real and imaginary parts. Similar to `integrate-field-function`, `func` is a function of position (a `vector3`) and the field components corresponding to `cs`: a list of `component` constants.

**`(output-real-field-function name cs func)`**
—
As `output-field-function`, but only outputs the real part of `func` to the dataset given by the string `name`.

See also [Field Functions](Field_Functions.md), and [Synchronizing the Magnetic and Electric Fields](Synchronizing_the_Magnetic_and_Electric_Fields.md) if you want to do computations combining the electric and magnetic fields.

#### Harminv

The following step function collects field data from a given point and runs [Harminv](https://github.com/NanoComp/harminv) on that data to extract the frequencies, decay rates, and other information.

**`(harminv c pt fcen df [maxbands])`**
—
Returns a step function that collects data from the field component `c` (e.g. $E_x$, etc.) at the given point `pt` (a `vector3`). Then, at the end of the run, it uses Harminv to look for modes in the given frequency range (center `fcen` and width `df`), printing the results to standard output (prefixed by `harminv:`) as comma-delimited text, and also storing them to the variable `harminv-results`. The optional argument `maxbands` is the maximum number of modes to search for. Defaults to 100.

**Important:** normally, you should only use `harminv` to analyze data *after the sources are off*. Wrapping it in `(after-sources (harminv ...))` is sufficient.

In particular, Harminv takes the time series $f(t)$ corresponding to the given field component as a function of time and decomposes it (within the specified bandwidth) as:

$$f(t) = \sum_n a_n e^{-i\omega_n t}$$

The results are stored in the list `harminv-results`, which is a list of tuples holding the frequency, amplitude, and error of the modes. Given one of these tuples, you can extract its various components with one of the accessor functions:

**`(harminv-freq result)`**
—
Return the complex frequency ω (in the usual Meep $2\pi c$ units).

**`(harminv-freq-re result)`**
—
Return the real part of the frequency ω.

**`(harminv-freq-im result)`**
—
Return the imaginary part of the frequency ω.

**`(harminv-Q result)`**
—
Return dimensionless lifetime, or quality factor, $Q$, defined as $-\mathrm{Re}\,\omega / 2 \mathrm{Im}\,\omega$.

**`(harminv-amp result)`**
—
Return the complex amplitude $a$.

**`(harminv-err result)`**
—
A crude measure of the error in the frequency (both real and imaginary)...if the error is much larger than the imaginary part, for example, then you can't trust the $Q$ to be accurate. **Note**: this error is only the uncertainty in the signal processing, and tells you nothing about the errors from finite resolution, finite cell size, and so on.

For example, `(map harminv-freq-re harminv-results)` gives a list of the real parts of the frequencies, using the Scheme built-in `map`.

### Step-Function Modifiers

Rather than writing a brand-new step function every time something a bit different is required, the following "modifier" functions take a bunch of step functions and produce *new* step functions with modified behavior. See also [Tutorial/Basics](Scheme_Tutorials/Basics.md) for examples.

#### Miscellaneous Step-Function Modifiers

**`(combine-step-funcs step-functions...)`**
—
Given zero or more step functions, return a new step function that on each step calls all of the passed step functions.

**`(synchronized-magnetic step-functions...)`**
—
Given zero or more step functions, return a new step function that on each step calls all of the passed step functions with the magnetic field synchronized in time with the electric field. See [Synchronizing the Magnetic and Electric Fields](Synchronizing_the_Magnetic_and_Electric_Fields.md).

#### Controlling When a Step Function Executes

**`(when-true cond? step-functions...)`**
—
Given zero or more step functions and a condition function `cond?` (a function of no arguments), evaluate the step functions whenever `(cond?)` returns `true`.

**`(when-false cond? step-functions...)`**
—
Given zero or more step functions and a condition function `cond?` (a function of no arguments), evaluate the step functions whenever `(cond?)` returns `false`.

**`(at-every dT step-functions...)`**
—
Given zero or more step functions, evaluates them at every time interval of $dT$ units (rounded up to the next time step).

**`(after-time T step-functions...)`**
—
Given zero or more step functions, evaluates them only for times after a $T$ time units have elapsed from the start of the run.

**`(before-time T step-functions...)`**
—
Given zero or more step functions, evaluates them only for times before a $T$ time units have elapsed from the start of the run.

**`(at-time T step-functions...)`**
—
Given zero or more step functions, evaluates them only once, after a $T$ time units have elapsed from the start of the run.

**`(after-sources step-functions...)`**
—
Given zero or more step functions, evaluates them only for times after all of the sources have turned off.

**`(after-sources+ T step-functions...)`**
—
Given zero or more step functions, evaluates them only for times after all of the sources have turned off, plus an additional $T$ time units have elapsed.

**`(during-sources step-functions...)`**
—
Given zero or more step functions, evaluates them only for times *before* all of the sources have turned off.

**`(at-beginning step-functions...)`**
—
Given zero or more step functions, evaluates them only once, at the beginning of the run.

**`(at-end step-functions...)`**
—
Given zero or more step functions, evaluates them only once, at the end of the run.

#### Modifying HDF5 Output

**`(in-volume v step-functions...)`**
—
Given zero or more step functions, modifies any output functions among them to only output a subset (or a superset) of the cell, corresponding to the `meep::volume* v` (created by the `volume` function).

**`(in-point pt step-functions...)`**
—
Given zero or more step functions, modifies any output functions among them to only output a single *point* of data, at `pt` (a `vector3`).

**`(to-appended filename step-functions...)`**
—
Given zero or more step functions, modifies any output functions among them to *append* their data to datasets in a single newly-created file named `filename` (plus an `.h5` suffix and the current filename prefix). They append by adding an *extra dimension* to their datasets, corresponding to time.

**`(with-prefix prefix step-functions...)`**
—
Given zero or more step functions, modifies any output functions among them to prepend the string `prefix` to the file names (much like `filename-prefix`, above).

### Writing Your Own Step Functions

A step function can take two forms. The simplest is just a function of no arguments, which is called at every time step (unless modified by one of the modifier functions above). e.g.

**`(define (my-step) (print "Hello world!\n"))`**

If one then does `(run-until 100 my-step)`, Meep will run for 100 time units and print "Hello world!" at every time step.

This suffices for most purposes. However, sometimes you need a step function that opens a file, or accumulates some computation, and you need to clean up (e.g. close the file or print the results) at the end of the run. For this case, you can write a step function of one argument: that argument will either be `'step` when it is called during time-stepping, or `'finish` when it is called at the end of the run.

Low-Level Functions
-------------------

By default, Meep reads input functions like `sources` and `geometry` and creates *global* variables `structure` and `fields` to store the corresponding C++ objects. Given these, you can then call essentially *any* function in the C++ interface, because all of the C++ functions are automatically made accessible to Scheme by the wrapper-generator program [SWIG](https://en.wikipedia.org/wiki/SWIG).

### Initializing the Structure and Fields

The `structure` and `fields` variables are automatically initialized when any of the run functions is called, or by various other functions such as `add-flux`. To initialize them separately, you can call `(init-fields)` manually, or `(init-structure k-point)` to just initialize the structure.

If you want to time step more than one field simultaneously, the easiest way is probably to do something like:

```scm
(init-fields)
(define my-fields fields)
(set! fields '())
(reset-meep)
```

and then change the geometry etc. and re-run `(init-fields)`. Then you'll have two field objects in memory.

### SWIG Wrappers

If you look at a function in the C++ interface, then there are a few simple rules to infer the name of the corresponding Scheme function.

-   First, all functions in the `meep::` namespace are prefixed with `meep-` in the Scheme interface.
-   Second, any method of a class is prefixed with the name of the class and a hyphen. For example, `meep::fields::step`, which is the function that performs a time-step, is exposed to Scheme as `meep-fields-step`. Moreover, you pass the object as the first argument in the Scheme wrapper. e.g. `f.step()` becomes `(meep-fields-step f)`.
-   To call the C++ constructor for a type, you use `new-*`. e.g. `(new-meep-fields ...)` returns a new `meep::fields` object. Conversely, to call the destructor and deallocate an object, you use `delete-*`; most of the time, this is not necessary because objects are automatically garbage collected.

Some argument type conversion is performed automatically, e.g. types like complex numbers are converted to `complex<double>`, etcetera. `vector3` vectors are converted to `meep::vec`, but to do this it is necessary to know the dimensionality of the problem in C++. The problem dimensions are automatically initialized by `init-structure`, but if you want to pass vector arguments to C++ before that time you should call `(require-dimensions!)`, which infers the dimensions from the `geometry-lattice`, `k-point`, and `dimensions` variables.
