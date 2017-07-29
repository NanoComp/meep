---
# Scheme User Interface
---

The Scheme user interface is documented in this page. We do not document the Scheme language or the functions provided by [libctl](https://libctl.readthedocs.io). See also the [libctl User Reference](https://libctl.readthedocs.io/en/latest/Libctl_User_Reference/) section of the [libctl manual](https://libctl.readthedocs.io).

This page is simply a compact listing of the functions exposed by the interface. For a gentler introduction, see the [tutorial](Scheme_Tutorial.md). Also, we note that this page is not a complete listing of all functions. In particular, because of the [SWIG wrappers](#swig-wrappers), every function in the C++ interface is accessible from Scheme, but not all of these functions are documented or intended for end users.

See also the instructions for [parallel Meep](Parallel_Meep.md) for MPI machines.

[TOC]

Input Variables
---------------

These are global variables that you can set to control various parameters of the Meep computation. In brackets after each variable is the type of value that it should hold. The classes, complex datatypes like `geometric-object`, are described in a later subsection. The basic datatypes, like `integer`, `boolean`, `cnumber`, and `vector3`, are defined by [libctl](https://libctl.readthedocs.io/en/latest/Libctl_User_Reference/).

**`geometry` [ list of `geometric-object` class ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Specifies the geometric objects making up the structure being simulated. When objects overlap, later objects in the list take precedence. Defaults to no objects (empty list).

**`sources` [ list of `source` class ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Specifies the current sources to be present in the simulation. Defaults to none.

**`symmetries` [ list of `symmetry` class ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Specifies the spatial (mirror/rotation) symmetries to exploit in the simulation. Defaults to none. The symmetries *must* be obeyed by *both* the structure *and* by the sources. See also: [Exploiting Symmetry](Exploiting_Symmetry.md).

**`pml-layers` [ list of `pml` class ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Specifies the absorbing [PML](Perfectly_Matched_Layer.md) boundary layers to use; defaults to none.

**`geometry-lattice` [`lattice` class ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Specifies the size of the unit cell which is centered on the origin of the coordinate system. Any sizes of `no-size` imply a reduced-dimensionality calculation. A 2d calculation is especially optimized. See `dimensions` below. Defaults to a cubic cell of unit size.

**`default-material` [`material-type` class ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Holds the default material that is used for points not in any object of the geometry list. Defaults to `air` (ε of 1). See also `epsilon-input-file` below.

**`epsilon-input-file` [`string`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
If this string is not empty (the default), then it should be the name of an HDF5 file whose first/only dataset defines a scalar dielectric function over some discrete grid. Alternatively, the dataset name can be specified explicitly if the string is in the form "filename:dataset". This dielectric function is then used in place of the ε property of `default-material` (*i.e.* where there are no `geometry` objects). The grid of the epsilon file dataset need not match Meep's computational grid; it is scaled and/or linearly interpolated as needed to map the file onto the computational cell (which warps the structure if the proportions of the grids do not match, however). **Note:** the file contents *only* override the ε property of the `default-material`, whereas other properties (μ, susceptibilities, nonlinearities, etc.) of `default-material` are still used.

**`dimensions` [`integer`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Explicitly specifies the dimensionality of the simulation, if the value is less than 3. If the value is 3 (the default), then the dimensions are automatically reduced to 2 if possible when `geometry-lattice` size in the $z$ direction is `no-size`. If `dimensions` is the special value of `CYLINDRICAL`, then cylindrical coordinates are used and the $x$ and $z$ dimensions are interpreted as $r$ and $z$, respectively. If `dimensions` is 1, then the cell must be along the $z$ direction and only $E_x$ and $H_y$ field components are permitted. If `dimensions` is 2, then the cell must be in the $xy$ plane.

**`m` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
For `CYLINDRICAL` simulations, specifies that the angular φ dependence of the fields is of the form $e^{im\phi}$ (default is `m=0`). If the simulation cell includes the origin $r=0$, then `m` must be an integer.

**`accurate-fields-near-cylorigin?` [`boolean`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;For `CYLINDRICAL` simulations with |*m*| &gt; 1, compute more accurate fields near the origin $r=0$ at the expense of requiring a smaller Courant factor. Empirically, when this option is set to `true`, a Courant factor of roughly $\min[0.5, 1 / (|m| + 0.5)]$ or smaller seems to be needed. The default is `false`, in which case the $D_r$, $D_z$, and $B_r$ fields within |*m*| pixels of the origin are forced to zero, which usually ensures stability with the default Courant factor of 0.5, at the expense of slowing convergence of the fields near $r=0$.

**`resolution` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Specifies the computational grid resolution in pixels per distance unit. Defaults to 10.

**`k-point` [`false` or `vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
If `false` (the default), then the boundaries are perfect metallic (zero electric field). If a vector, then the boundaries are Bloch-periodic: the fields at one side are $\exp(i\mathbf{k}\cdot\mathbf{R})$ times the fields at the other side, separated by the lattice vector $\mathbf{R}$. The `k-point` vector is specified in *Cartesian* coordinates in units of 2π/distance. Note: this is *different* from [MPB](http://mpb.readthedocs.io), equivalent to taking MPB's `k-points` through the function `reciprocal->cartesian`.

**`ensure-periodicity` [`boolean`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
If `true` (the default) *and* if the boundary conditions are periodic (`k-point` is not `false`), then the geometric objects are automatically repeated periodically according to the lattice vectors which define the size of the computational cell.

**`eps-averaging?` [`boolean`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
If `true` (the default), then subpixel averaging is used when initializing the dielectric function. See the [reference publication](License_and_Copyright.md#referencing). The input variables `subpixel-maxeval` (default 100000) and `subpixel-tol` (default 1.0e-4) specify the maximum number of function evaluations and the integration tolerance for subpixel averaging. Increasing/decreasing these, respectively, will cause a more accurate but slower computation of the average ε with diminishing returns for the actual FDTD error.

**`force-complex-fields?` [`boolean`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
By default, Meep runs its simulations with purely real fields whenever possible. It uses complex fields which require twice the memory and computation if the `k-point` is non-zero or if `m` is non-zero. However, by setting `force-complex-fields?` to `true`, Meep will always use complex fields.

**`filename-prefix` [`string`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
A string prepended to all output filenames. If empty (the default), then Meep uses the name of the current ctl file, with ".ctl" replaced by "-" (e.g. `foo.ctl` uses a `"foo-"` prefix). See also: [Output File Names](Scheme_User_Interface.md#output-file-names).

**`Courant` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Specify the Courant factor $S$ which relates the time step size to the spatial discretization: $c\Delta t = S\Delta x$. Default is `0.5`. For numerical stability, the Courant factor must be *at most* $n_\textrm{min}/\sqrt{\textrm{# dimensions}}$, where $n_\textrm{min}$ is the minimum refractive index (usually 1), and in practice $S$ should be slightly smaller.

**`output-volume` [`meep::geometric_volume*`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Specifies the default region of space that is output by the HDF5 output functions (below); see also the `(volume` `...)` function to create `meep::geometric_volume*` objects. The default is `'()` (null), which means that the whole computational cell is output. Normally, you should use the `(in-volume` `...)` function to modify the output volume instead of setting `output-volume` directly.

**`output-single-precision?` [`boolean`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Meep performs its computations in [double precision](https://en.wikipedia.org/wiki/double_precision), and by default its output HDF5 files are in the same format. However, by setting this variable to `true` (default is `false`) you can instead output in [single precision](https://en.wikipedia.org/wiki/single_precision) which saves a factor of two in space.

**`progress-interval` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Time interval (seconds) after which Meep prints a progress message. Default is 4 seconds.

**`extra-materials` [ list of `material-type` class ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
By default, Meep turns off support for material dispersion, nonlinearities, and similar properties if none of the objects in the `geometry` have materials with these properties—since they are not needed, it is faster to omit their calculation. This doesn't work however if you use a `material-function`: materials via a user-specified function of position instead of just geometric objects. If only your material function returns a nonlinear material, for example, Meep won't notice this unless you tell it explicitly via `extra-materials`. `extra-materials` is a list of materials that Meep should look for in the computational cell in addition to any materials that are specified by geometric objects. You should list any materials other than scalar dielectrics that are returned by `material-function`s here.

The following require a bit more understanding of the inner workings of Meep to use. See also: [SWIG wrappers](#swig-wrappers).

**`structure` [`meep::structure*`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Pointer to the current structure being simulated; initialized by `(init-structure)` which is called automatically by `(init-fields)` which is called automatically by any of the `(run)` functions.

**`fields` [`meep::fields*`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Pointer to the current fields being simulated; initialized by `(init-fields)` which is called automatically by any of the `(run)` functions.

**`num-chunks` [`integer`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Minimum number of "chunks" (subarrays) to divide the structure/fields into (default 0). Actual number is determined by number of processors, PML layers, etcetera. Mainly useful for debugging.

Predefined Variables
--------------------

**`air`, `vacuum` [`material-type` class ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Two aliases for a predefined material type with a dielectric constant of 1.

**`perfect-electric-conductor` or `metal` [`material-type` class ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
A predefined material type corresponding to a perfect electric conductor at the boundary of which the parallel electric field is zero. Technically, $\varepsilon = -\infty$.

**`perfect-magnetic-conductor` [`material-type` class ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
A predefined material type corresponding to a perfect magnetic conductor at the boundary of which the parallel magnetic field is zero. Technically, $\mu = -\infty$.

**`nothing` [`material-type` class ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
A material that, effectively, punches a hole through other objects to the background (`default-material`).

**`infinity` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A big number (1.0e20) to use for "infinite" dimensions of objects.

**`pi` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;π (3.14159...).

Constants (Enumerated Types)
----------------------------

Several of the functions/classes in Meep ask you to specify e.g. a field component or a direction in the grid. These should be one of the following constants:

**`direction` constants**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Specify a direction in the grid. One of `X`, `Y`, `Z`, `R`, `P` for $x$, $y$, $z$, $r$, $\phi$, respectively.

**`boundary-side` constants**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Specify particular boundary in some direction (e.g. $+x$ or $-x$). One of `Low` or `High`.

**`component` constants**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Specify a particular field or other component. One of `Ex`, `Ey`, `Ez`, `Er`, `Ep`, `Hx`, `Hy`, `Hz`, `Hy`, `Hp`, `Hz`, `Bx`, `By`, `Bz`, `By`, `Bp`, `Bz`, `Dx`, `Dy`, `Dz`, `Dr`, `Dp`, `Dielectric`, `Permeability`, for $E_x$, $E_y$, $E_z$, $E_r$, $E_\phi$, $H_x$, $H_y$, $H_z$, $H_r$, $H_\phi$, $B_x$, $B_y$, $B_z$, $B_r$, $B_\phi$, $D_x$, $D_y$, $D_z$, $D_r$, $D_\phi$, $\varepsilon$, $\mu$, respectively.

**`derived-component` constants**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
These are additional components which are not actually stored by Meep but are computed as needed, mainly for use in output functions. One of `Sx`, `Sy`, `Sz`, `Sr`, `Sp`, `EnergyDensity`, `D-EnergyDensity`, `H-EnergyDensity` for $S_x$, $S_y$, $S_z$, $S_r$, $S_\phi$ (components of the Poynting vector $\mathrm{Re}\,\mathbf{E}^* \times \mathbf{H}$), $(\mathbf{E}^* \cdot \mathbf{D} + \mathbf{H}^* \cdot \mathbf{B})/2$, $\mathbf{E}^* \cdot \mathbf{D}/2$, $\mathbf{H}^* \cdot \mathbf{B}/2$, respectively.

Classes
-------

Classes are complex datatypes with various "properties" which may have default values. Classes can be "subclasses" of other classes. Subclasses inherit all the properties of their superclass and can be used in any place the superclass is expected. An object of a class is constructed with:

```
(make class (prop1 val1) (prop2 val2) ...)
```

See also the [libctl manual](https://libctl.readthedocs.io).

Meep defines several types of classes, the most numerous of which are the various geometric object classes which are the same as those used in [MPB](http://mpb.readthedocs.io). You can also get a list of the available classes, along with their property types and default values, at runtime with the `(help)` command.

### lattice

The `lattice` class is normally used only for the `geometry-lattice` variable, which sets the size of the computational cell. In [MPB](http://mpb.readthedocs.io), you can use this to specify a variety of affine lattice structures. In [Meep](index.md), only rectangular Cartesian computational cells are supported, so the only property of lattice that you should normally use is its `size`.

**`size` [`vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The size of the computational cell. Defaults to unit lengths.

If any dimension has the special size `no-size`, then the dimensionality of the problem is essentially reduced by one. Strictly speaking, the dielectric function is taken to be uniform along that dimension.

Because Maxwell's equations are scale-invariant, you can use any units of distance you want to specify the cell size: nanometers, inches, parsecs, whatever. However, it is usually convenient to pick some characteristic lengthscale of your problem and set that length to 1. See also: [Units](Introduction.md#units-in-meep).

### material

This class is used to specify the materials that geometric objects are made of. Currently, there are three subclasses, `dielectric`, `perfect-metal`, and `material-function`.

**`medium`**  

An electromagnetic medium which is possibly nonlinear and/or dispersive. See also: [Materials](Materials.md). For backwards compatibility, a synonym for `medium` is `dielectric`. It has several properties:

**`epsilon` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The frequency-independent isotropic relative permittivity or dielectric constant. Default is 1. You can also use `(index n)` as a synonym for `(epsilon (* n n))`; note that this is not really the refractive index if you also specify μ, since the true index is $\sqrt{\mu\epsilon}$.

Using `(epsilon ε)` is actually a synonym for `(epsilon-diag ε ε ε)`.

**`epsilon-diag` and `epsilon-offdiag` [`vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
These properties allow you to specify ε as an arbitrary real-symmetric tensor by giving the diagonal and offdiagonal parts. Specifying `(epsilon-diag` `a` `b` `c)` and/or `(epsilon-offdiag` `u` `v` `w)` corresponds to a relative permittivity $\varepsilon$ tensor
\\begin{pmatrix} a & u & v \\\\ u & b & w \\\\ v & w & c \\end{pmatrix}

The default is the identity matrix ($a = b = c = 1$ and $u = v = w = 0$).

**`mu` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The frequency-independent isotropic relative permeability μ. Default is 1. Using `(mu μ)` is actually a synonym for `(mu-diag μ μ μ)`.

**`mu-diag` and `mu-offdiag` [`vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
These properties allow you to specify μ as an arbitrary real-symmetric tensor by giving the diagonal and offdiagonal parts exactly as for ε above. Default is the identity matrix.

**`D-conductivity` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The frequency-independent electric conductivity $\sigma_D$. Default is 0. You can also specify an diagonal anisotropic conductivity tensor by using the property `D-conductivity-diag` [vector3], which takes three numbers or a [vector3] to give the $\sigma_D$ tensor diagonal. See also [Conductivity](Materials.md#conductivity-and-complex).

**`B-conductivity` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The (frequency-independent) magnetic conductivity $\sigma_B$. Default is 0. You can also specify an diagonal anisotropic conductivity tensor by using the property `B-conductivity-diag` [vector3], which takes three numbers or a [vector3] to give the $\sigma_B$ tensor diagonal. See also [Conductivity](Materials.md#conductivity-and-complex).

**`chi2` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The nonlinear (Pockels) susceptibility $\chi^{(2)}$. Default is 0. See also [Nonlinearity](Materials.md#nonlinearity).

**`chi3` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The nonlinear (Kerr) susceptibility $\chi^{(3)}$. Default is 0. See also [Nonlinearity](Materials.md#nonlinearity).

**`E-susceptibilities` [ list of `susceptibility` class ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
List of dispersive susceptibilities (see below) added to the dielectric constant ε in order to model material dispersion. Defaults to none. See also: [Material Dispersion](Materials.md#material-dispersion). For backwards compatibility, synonyms of `E-susceptibilities` are `E-polarizations` and `polarizations`.

**`H-susceptibilities` [ list of `susceptibility` class ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
List of dispersive susceptibilities (see below) added to the permeability μ in order to model material dispersion. Defaults to none. See also: [Material Dispersion](Materials.md#material-dispersion). For backwards compatibility, a synonym of `H-susceptibilities` is `perfect-metal`.

**`perfect-metal`**  

A perfectly conducting metal. This class has no properties and you normally just use the predefined `metal` object, above. To model imperfect conductors, use a dispersive dielectric material. See also the [predefined variables](#Predefined_Variables) `metal`, `perfect-electric-conductor`, and `perfect-magnetic-conductor` above.

**`material-function`**  

This material type allows you to specify the material as an arbitrary function of position. It has one property:

**`material-func` [`function`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
A function of one argument, the position `vector3`, that returns the material at that point. Note that the function you supply can return *any* material. It's even possible to return another `material-function` object which would then have its function invoked in turn.

Instead of `material-func`, you can use `epsilon-func`: give it a function of position that returns the dielectric constant at that point.

**Important:** If your material function returns nonlinear, dispersive (Lorentzian or conducting), or magnetic materials, you should also include a list of these materials in the `extra-materials` input variable (above) to let Meep know that it needs to support these material types in your simulation. For dispersive materials, you need to include a material with the *same* γ<sub>*n*</sub> and ω<sub>*n*</sub> values, so you can only have a finite number of these, whereas σ<sub>*n*</sub> can vary continuously if you want and a matching σ<sub>*n*</sub> need not be specified in `extra-materials`. For nonlinear or conductivity materials, your `extra-materials` list need not match the actual values of σ or χ returned by your material function, which can vary continuously if you want.

**Complex ε and μ**: you cannot specify a frequency-independent complex ε or μ in Meep where the imaginary part is a frequency-independent loss but there is an alternative. That is because there are only two important physical situations. First, if you only care about the loss in a narrow bandwidth around some frequency, you can set the loss at that frequency via the conductivity (see [Conductivity](Materials.md#conductivity)). Second, if you care about a broad bandwidth, then all physical materials have a frequency-dependent complex ε and/or μ, and you need to specify that frequency dependence by fitting to Lorentzian and/or Drude resonances via the `lorentzian-susceptibility` or `drude-susceptibility` classes below.

Dispersive dielectric and magnetic materials, above, are specified via a list of objects that are subclasses of type `susceptibility`.

**`susceptibility`**  

Parent class for various dispersive susceptibility terms, parameterized by an anisotropic amplitude σ. See [material dispersion](Materials.md#material-dispersion).

**`sigma` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The scale factor σ. You can also specify an anisotropic σ tensor by using the property `sigma-diag` [vector3], which takes three numbers or a [vector3] to give the $\sigma_n$ tensor diagonal, and `sigma-offdiag` [vector3] which specifies the offdiagonal elements (defaults to 0). That is, `(sigma-diag a b c)` and `(sigma-offdiag u v w)` corresponds to a σ tensor

\\begin{pmatrix} a & u & v \\\\ u & b & w \\\\ v & w & c \\end{pmatrix}

**`lorentzian-susceptibility`**  

Specifies a single dispersive susceptibility of Lorentzian (damped harmonic oscillator) form. See [material dispersion](Materials.md#material-dispersion), with the parameters (in addition to σ):

**`frequency` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The resonance frequency $f_n = \omega_n / 2\pi$.

**`gamma` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The resonance loss rate $\gamma_n / 2\pi$.

**`drude-susceptibility`**  

Specifies a single dispersive susceptibility of Drude form. See [material dispersion](Materials.md#material-dispersion), with the parameters (in addition to σ):

**`frequency` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The frequency scale factor $f_n = \omega_n / 2\pi$ which multiplies σ (not a resonance frequency).

**`gamma` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The loss rate $\gamma_n / 2\pi$.

Meep also supports a somewhat unusual polarizable medium, a Lorentzian susceptibility with a random noise term added into the damped-oscillator equation at each point. This can be used to directly model thermal radiation in both the [far field](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.93.213905) and the [near field](http://math.mit.edu/~stevenj/papers/RodriguezIl11.pdf). Note, however that it is more efficient to compute far-field thermal radiation using Kirchhoff's law of radiation, which states that emissivity equals absorptivity. Near-field thermal radiation can usually be [computed more efficiently](http://math.mit.edu/~stevenj/papers/RodriguezRe13-heat.pdf) using frequency domain techniques, e.g. via our [SCUFF-EM](http://homerreid.dyndns.org/scuff-EM/) code.

**`noisy-lorentzian-susceptibility` or `noisy-drude-susceptibility`**  

Specifies a single dispersive susceptibility of Lorentzian (damped harmonic oscillator) or Drude form. See [material dispersion](Materials.md#material-dispersion), with the same σ, `frequency`, and `gamma`) parameters, but with an additional Gaussian random noise term (uncorrelated in space and time, zero mean) added to the **P** damped-oscillator equation.

**`noise-amp` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The noise has root-mean square amplitude σ⋅`noise-amp`.

### geometric-object

This class, and its descendants, are used to specify the solid geometric objects that form the dielectric structure being simulated. The base class is:

**`geometric-object`**  

Properties:

**`material` [`material-type` class ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The material that the object is made of (usually some sort of dielectric). No default value (must be specified).

**`center` [`vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Center point of the object. No default value.

One normally does not create objects of type `geometric-object` directly, however; instead, you use one of the following subclasses. Recall that subclasses inherit the properties of their superclass, so these subclasses automatically have the `material` and `center` properties (which must be specified, since they have no default values).

In a two-dimensional calculation, only the intersections of the objects with the x-y plane are considered.

**`sphere`**

A sphere. Properties:

**`radius` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Radius of the sphere. No default value.

**`cylinder`**

A cylinder, with circular cross-section and finite height. Properties:

**`radius` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Radius of the cylinder's cross-section. No default value.

**`height` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Length of the cylinder along its axis. No default value.

**`axis` [`vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Direction of the cylinder's axis; the length of this vector is ignored. Defaults to point parallel to the $z$ axis.

**`cone`**

A cone, or possibly a truncated cone. This is actually a subclass of `cylinder`, and inherits all of the same properties, with one additional property. The radius of the base of the cone is given by the `radius` property inherited from `cylinder`, while the radius of the tip is given by the new property, `radius2`. The `center` of a cone is halfway between the two circular ends.

**`radius2` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Radius of the tip of the cone (i.e. the end of the cone pointed to by the `axis` vector). Defaults to zero (a "sharp" cone).

**`block`**

A parallelepiped (i.e., a brick, possibly with non-orthogonal axes). Properties:

**`size` [`vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The lengths of the block edges along each of its three axes. Not really a 3-vector, but it has three components, each of which should be nonzero. No default value.

**`e1`, `e2`, `e3` [`vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The directions of the axes of the block; the lengths of these vectors are ignored. Must be linearly independent. They default to the three lattice directions.

**`ellipsoid`**

An ellipsoid. This is actually a subclass of `block`, and inherits all the same properties, but defines an ellipsoid inscribed inside the block.

Here are some examples of geometric objects created using the above classes, assuming `mat` is some material we have defined:

```
; A cylinder of infinite radius and height 0.25 pointing along the x axis,
; centered at the origin:
(make cylinder (center 0 0 0) (material mat) 
               (radius infinity) (height 0.25) (axis 1 0 0))
```


```
; An ellipsoid with its long axis pointing along (1,1,1), centered on
; the origin (the other two axes are orthogonal and have equal
; semi-axis lengths):
(make ellipsoid (center 0 0 0) (material mat)
                (size 0.8 0.2 0.2)
               (e1 1 1 1)
               (e2 0 1 -1)
               (e3 -2 1 1))
```


```
; A unit cube of material m with a spherical air hole of radius 0.2 at
; its center, the whole thing centered at (1,2,3):
(set! geometry (list
               (make block (center 1 2 3) (material mat) (size 1 1 1))
               (make sphere (center 1 2 3) (material air) (radius 0.2))))
```

### symmetry

This class is used for the `symmetries` input variable to specify symmetries which must preserve both the structure *and* the sources for Meep to exploit. Any number of symmetries can be exploited simultaneously but there is no point in specifying redundant symmetries: the computational cell can be reduced by at most a factor of 4 in 2d and 8 in 3d. See also: [Exploiting Symmetry](Exploiting_Symmetry.md).

**`symmetry`**  

A single symmetry to exploit. This is the base class of the specific symmetries below, so normally you don't create it directly. However, it has two properties which are shared by all symmetries:

**`direction` [`direction` constant ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The direction of the symmetry (the normal to a mirror plane or the axis for a rotational symmetry). e.g. `X`, `Y`, ... (only Cartesian/grid directions are allowed). No default value.

**`phase` [`cnumber`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
An additional phase to multiply the fields by when operating the symmetry on them; defaults to `1.0`. e.g. a phase of `-1` for a mirror plane corresponds to an *odd* mirror. Technically, you are essentially specifying the representation of the symmetry group that your fields and sources transform under.

The specific symmetry sub-classes are:

**`mirror-sym`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
A mirror symmetry plane. Here, the `direction` is the direction *normal* to the mirror plane.

**`rotate2-sym`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
A 180° (twofold) rotational symmetry (a.k.a. $C_2$). Here, the `direction` is the axis of the rotation.

**`rotate4-sym`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
A 90° (fourfold) rotational symmetry (a.k.a. $C_4$). Here, the `direction` is the axis of the rotation.

### pml

This class is used for specifying the PML absorbing boundary layers around the cell, if any, via the `pml-layers` input variable. See also [Perfectly Matched Layers](Perfectly_Matched_Layer.md). `pml-layers` can be zero or more `pml` objects, with multiple objects allowing you to specify different PML layers on different boundaries.

**`pml`**  

A single PML layer specification, which sets up one or more PML layers around the boundaries according to the following properties.

**`thickness` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The spatial thickness of the PML layer (which extends from the boundary towards the *inside* of the computational cell). The thinner it is, the more numerical reflections become a problem. No default value.

**`direction` [`direction` constant ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Specify the direction of the boundaries to put the PML layers next to. e.g. if `X`, then specifies PML on the $\pm x$ boundaries (depending on the value of `side`, below). Default is the special value `ALL`, which puts PML layers on the boundaries in all directions.

**`side` [`boundary-side` constant ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Specify which side, `Low` or `High` of the boundary or boundaries to put PML on. e.g. if side is `Low` and direction is `X`, then a PML layer is added to the $-x$ boundary. Default is the special value `ALL`, which puts PML layers on both sides.

**`strength` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
A strength (default is `1.0`) to multiply the PML absorption coefficient by. A strength of `2.0` will *square* the theoretical asymptotic reflection coefficient of the PML (making it smaller), but will also increase numerical reflections. Alternatively, you can change `R-asymptotic`, below.

**`R-asymptotic` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The asymptotic reflection in the limit of infinite resolution or infinite PML thickness, for refections from air (an upper bound for other media with index &gt; 1). (For a finite resolution or thickness, the reflection will be *much larger*, due to the discretization of Maxwell's equation.) The default value is 10<sup>−15</sup>, which should suffice for most purposes. (You want to set this to be small enough so that waves propagating within the PML are attenuated sufficiently, but making `R-asymptotic` too small will increase the numerical reflection due to discretization.)

**`pml-profile` [`function`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
By default, Meep turns on the PML conductivity quadratically within the PML layer—one doesn't want to turn it on suddenly, because that exacerbates reflections due to the discretization. More generally, with `pml-profile` one can specify an arbitrary PML "profile" function *f*(*u*) that determines the shape of the PML absorption profile up to an overall constant factor. *u* goes from 0 to 1 at the start and end of the PML, and the default is *f*(*u*)=*u*<sup>2</sup>. In some cases where a very thick PML is required, such as in a periodic medium (where there is technically no such thing as a true PML, only a pseudo-PML), it can be advantageous to turn on the PML absorption more smoothly (see [Oskooi et al., 2008](http://www.opticsinfobase.org/abstract.cfm?URI=oe-16-15-11376)). For example, one can use a cubic profile *f*(*u*)=*u*<sup>3</sup> by specifying `(pml-profile` `(lambda` `(u)` `(*` `u` `u` `u)))`.

#### `absorber`

Instead of a `pml` layer, there is an alternative class called `absorber` which is a **drop-in** replacement for `pml`. For example, you can do `(set!` `pml-layers` `(list` `(make` `absorber` `(thickness` `2))))` instead of `(set!` `pml-layers` `(list` `(make` `pml` `(thickness` `2))))`. All the parameters are the same as for `pml`, above. You can have a mix of `pml` on some boundaries and `absorber` on others.

The `absorber` class does *not* implement a perfectly matched layer (PML), however (except in 1d). Instead, it is simply a scalar electric **and** magnetic conductivity that turns on gradually within the layer according to the `pml-profile` (defaulting to quadratic). Such a scalar conductivity gradient is only reflectionless in the limit as the layer becomes sufficiently thick.

The main reason to use `absorber` is if you have **a case in which PML fails:**

-   No true PML exists for *periodic* media, and a scalar absorber is cheaper and is generally just as good. See [Oskooi et al. (2008)](http://math.mit.edu/~stevenj/papers/papers_abstracts.html#OskooiZh08).
-   PML can lead to *divergent* fields for certain waveguides with "backward-wave" modes; this can easily happen in metallic with surface plasmons, and a scalar absorber is your only choice. See [Loh et al. (2009)](http://math.mit.edu/~stevenj/papers/papers_abstracts.html#LohOs09).
-   PML can fail if you have a waveguide hitting the edge of your computational cell *at an angle*. See [Oskooi et. al. (2011).](http://math.mit.edu/~stevenj/papers/papers_abstracts.html#OskooiJo11)

### source

The `source` class is used to specify the current sources (via the `sources` input variable). Note that all sources in Meep are separable in time and space, i.e. of the form $\mathbf{J}(\mathbf{x},t) = \mathbf{A}(\mathbf{x}) \cdot f(t)$ for some functions $\mathbf{A}$ and $f$. (Non-separable sources can be simulated, however, by modifying the sources after each time step.) When real fields are being used (which is the default in many cases...see the `force-complex-fields?` input variable), only the real part of the current source is used by Meep.

**Important note**: These are *current* sources (**J** terms in Maxwell's equations), even though they are labelled by electric/magnetic field components. They do *not* specify a particular electric/magnetic field (which would be what is called a "hard" source in the FDTD literature). There is no fixed relationship between the current source and the resulting field amplitudes; it depends on the surrounding geometry, as described in the [Meep FAQ](Meep_FAQ#How_does_the_current_amplitude_relate_to_the_resulting_field_amplitude?.md) and in [our book chapter online](http://arxiv.org/abs/arXiv:1301.5366).

**`source`**  

The source class has the following properties:

**`src` [`src-time` class ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Specify the time-dependence of the source (see below). No default.

**`component` [`component` constant ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Specify the direction and type of the current component: e.g. `Ex`, `Ey`, etcetera for an electric-charge current, and `Hx`, `Hy`, etcetera for a magnetic-charge current. Note that currents pointing in an arbitrary direction are specified simply as multiple current sources with the appropriate amplitudes for each component. No default.

**`center` [`vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The location of the center of the current source in the computational cell; no default.

**`size` [`vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The size of the current distribution along each direction of the computational cell. The default is (0,0,0): a point-dipole source.

**`amplitude` [`cnumber`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
An overall (complex) amplitude multiplying the the current source. Default is `1.0`.

**`amp-func` [`function`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
A Scheme function of a single argument, that takes a vector3 giving a position and returns a (complex) current amplitude for that point. The position argument is *relative* to the `center` of the current source, so that you can move your current around without changing your function. The default is `'()` (null), meaning that a constant amplitude of 1.0 is used. Note that your amplitude function (if any) is *multiplied* by the `amplitude` property, so both properties can be used simultaneously.

As described in section 4.2 of [our book chapter online](http://arxiv.org/abs/arXiv:1301.5366), it is also possible to supply a source that is designed to couple exclusively into a single waveguide mode (or other mode of some cross section or periodic region) at a single frequency, and which couples primarily into that mode as long as the bandwidth is not too broad. This is possible if you have [MPB](http://mpb.readthedocs.io) installed: Meep will call MPB to compute the field profile of the desired mode, and uses the field profile to produce an equivalent current source. (Note: this feature does *not* work in cylindrical coordinates.) To do this, instead of a `source` you should use an `eigenmode-source`:

**`eigenmode-source`**  

This is a subclass of `source` and has **all of the properties** of `source` above. However, you normally do not specify a `component`. Instead of `component`, the current source components and amplitude profile are computed by calling MPB to compute the modes of the dielectric profile in the region given by the `size` and `center` of the source, with the modes computed as if the *source region were repeated periodically in all directions*. If an `amplitude` and/or `amp-func` are supplied, they are *multiplied* by this current profile. The desired eigenmode and other features are specified by the following properties:

**`eig-band` [`integer`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The index *n* (1,2,3,...) of the desired band ω<sub>*n*</sub>(**k**) to compute in MPB (1 denotes the lowest-frequency band at a given **k** point, and so on).

**`direction` [`X`, `Y`, or `Z;` default `AUTOMATIC`], `eig-match-freq?` [`boolean;` default `true`], `eig-kpoint` [`vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
By default (if `eig-match-freq?` is `true`), Meep tries to find a mode with the same frequency ω<sub>*n*</sub>(**k**) as the `src` property (above), by scanning **k** vectors in the given `direction` using MPB's find-k functionality. Alternatively, if `eig-kpoint` is supplied, it is used as an initial guess and direction for **k**. By default, `direction` is the direction normal to the source region, assuming `size` is *d*–1 dimensional in a *d*-dimensional simulation (e.g. a plane in 3d). Alternatively (if `eig-match-freq?` is `false`), you can specify a particular **k** vector of the desired mode with `eig-kpoint` (in Meep units of 2π/*a*).

**`eig-parity` [`NO-PARITY` (= default), `EVEN-Z` (= `TE`), `ODD-Z` (= `TM`), `EVEN-Y`, `ODD-Y`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The parity (= polarization in 2d) of the mode to calculate, assuming the structure has *z* and/or *y* mirror symmetry *in the source region*. If the structure has both *y* and *z* mirror symmetry, you can combine more than one of these, e.g. `EVEN-Z` `+` `ODD-Y`. The default is `NO-PARITY`, in which case MPB computes all of the bands which will still be even or odd if the structure has mirror symmetry, of course. This is especially useful in 2d simulations to restrict yourself to a desired polarization.

**`eig-resolution` [`integer`, defaults to same as Meep resolution ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The spatial resolution to use in MPB for the eigenmode calculations. This defaults to the same resolution as Meep, but you can use a higher resolution in which case the structure is linearly interpolated from the Meep pixels.

**`eig-tolerance` [`number`, defaults to 10<sup>–7</sup> ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The tolerance to use in the MPB eigensolver. MPB terminates when the eigenvalues stop changing to less than this fractional tolerance.

**`component` [as above, but defaults to `ALL-COMPONENTS`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Once the MPB modes are computed, equivalent electric and magnetic sources are created within Meep. By default, these sources include magnetic and electric currents in *all* transverse directions within the source region, corresponding to the mode fields as described in our book chapter. If you specify a `component` property, however, you can include only one component of these currents if you wish. Most people won't need this feature.

**`eig-lattice-size` [`vector3`], `eig-lattice-center` [`vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Normally, the MPB computational unit cell is the same as the source volume (given by the `size` and `center` parameters). However, occasionally you want the unit cell to be larger than the source volume. For example, to create an eigenmode source in a periodic medium (photonic crystal), you need to pass MPB the entire unit cell of the periodic medium, but once the mode is computed then the actual current sources need only lie on a cross section of that medium. To accomplish this, you can specify the optional `eig-lattice-size` and `eig-lattice-center`, which define a volume (which must enclose `size` and `center`) that is used for the unit cell in MPB with the dielectric function ε taken from the corresponding region in the Meep simulation.

Note that MPB only supports dispersionless non-magnetic materials but it does support anisotropic ε. Any nonlinearities, magnetic responses µ, conductivities σ, or dispersive polarizations in your materials will be *ignored* when computing the eigenmode source. PML will also be ignored.

The `src-time` object, which specifies the time dependence of the source, can be one of the following three classes.

**`continuous-src`**  

A continuous-wave source proportional to $\exp(-i\omega t)$, possibly with a smooth (exponential/tanh) turn-on/turn-off. It has the properties:

**`frequency` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The frequency *f* in units of c/distance (or ω in units of 2πc/distance). See: [Units](Meep_Introduction#Units_in_Meep.md). No default value. You can instead specify `(wavelength` `x)` or `(period` `x)`, which are both a synonym for `(frequency` `(/` `1` `x))`; i.e. 1/ω in these units is the vacuum wavelength or the temporal period.

**`start-time` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The starting time for the source; default is `0` (turn on at $t=0$).

**`end-time` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The end time for the source; default is `infinity` (never turn off).

**`width` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Roughly, the temporal width of the smoothing (technically, the inverse of the exponential rate at which the current turns off and on). Default is `0` (no smoothing). You can instead specify `(fwidth` `x)`, which is a synonym for `(width` `(/` `1` `x))` (i.e. the frequency width is proportional to the inverse of the temporal width).

**`cutoff` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
How many `width`s the current decays for before we cut it off and set it to zero; the default is `3.0`. A larger value of `cutoff` will reduce the amount of high-frequency components that are introduced by the start/stop of the source, but will of course lead to longer simulation times.

**`gaussian-src`**  

A Gaussian-pulse source roughly proportional to $\exp(-i\omega t - (t-t_0)^2/2w^2)$. Technically, the "Gaussian" sources in Meep are the (discrete-time) derivative of a Gaussian, i.e. they are $(-i\omega)^{-1} \frac{\partial}{\partial t} \exp(-i\omega t - (t-t_0)^2/2w^2)$, but the difference between this and a true Gaussian is usually irrelevant. It has the properties:

**`frequency` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The center frequency *f* in units of c/distance (or ω in units of 2πc/distance). See [Units](Introduction.md#units-in-meep). No default value. You can instead specify `(wavelength` `x)` or `(period` `x)`, which are both a synonym for `(frequency` `(/` `1` `x))`; i.e. 1/ω in these units is the vacuum wavelength or the temporal period.

**`width` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The width $w$ used in the Gaussian. No default value. You can instead specify `(fwidth` `x)`, which is a synonym for `(width` `(/` `1` `x))` (i.e. the frequency width is proportional to the inverse of the temporal width).

**`start-time` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The starting time for the source; default is `0` (turn on at $t=0$). (Not the time of the peak! See below.)

**`cutoff` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
How many `width`s the current decays for before we cut it off and set it to zero—this applies for both turn-on and turn-off of the pulse. The default is `5.0`. A larger value of `cutoff` will reduce the amount of high-frequency components that are introduced by the start/stop of the source, but will of course lead to longer simulation times. The peak of the gaussian is reached at the time *t*<sub>0</sub>= `start-time` + `cutoff`\*`width`.

**`custom-src`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
A user-specified source function $f(t)$. You can also specify start/end times (at which point your current is set to zero whether or not your function is actually zero). These are optional, but you *must specify* an `end-time` explicitly if you want functions like `run-sources` to work, since they need to know when your source turns off.

**`src-func` [`function`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The function $f(t)$ specifying the time-dependence of the source. It should take one argument (the time in Meep units) and return a complex number.

**`start-time` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The starting time for the source. Default is `(-infinity)`: turn on at $t=-\infty$). Note, however, that the simulation normally starts at $t=0$ with zero fields as the initial condition, so there is implicitly a sharp turn-on at $t=0$ whether you specify it or not.

**`end-time` [`number`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The end time for the source; default is `infinity` (never turn off).

### flux-region

A `flux-region` object is used with `add-flux` [below](#Flux_spectra.md) to specify a region in which Meep should accumulate the appropriate Fourier-transformed fields in order to compute a flux spectrum.

**`flux-region`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
A region (volume, plane, line, or point) in which to compute the integral of the Poynting vector of the Fourier-transformed fields. Its properties are:

**`center` [`vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The center of the flux region (no default).

**`size` [`vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The size of the flux region along each of the coordinate axes; default is (0,0,0) (a single point).

**`direction` [`direction` constant ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The direction in which to compute the flux (e.g. `X`, `Y`, etcetera). The default is `AUTOMATIC`, in which the direction is determined by taking the normal direction if the flux region is a plane (or a line, in 2d). If the normal direction is ambiguous (e.g. for a point or volume), then you *must* specify the `direction` explicitly (not doing so will lead to an error).

**`weight` [`cnumber`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A weight factor to multiply the flux by when it is computed; default is `1.0`.

Note that the flux is always computed in the *positive* coordinate direction, although this can effectively be flipped by using a `weight` of `-1.0`. (This is useful, for example, if you want to compute the outward flux through a box, so that the sides of the box add instead of subtract!)

Miscellaneous Functions
-----------------------

Here, we describe a number of miscellaneous useful functions provided by Meep.

See also the [reference section](https://libctl.readthedocs.io/en/latest/Libctl_User_Reference/) of the libctl manual, which describes a number of useful functions defined by libctl.

### Geometry Utilities

Some utility functions are provided to help you manipulate geometric objects:

**`(shift-geometric-object obj shift-vector)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Translate `obj` by the 3-vector `shift-vector`.

**`(geometric-object-duplicates shift-vector min-multiple max-multiple obj)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Return a list of duplicates of `obj`, shifted by various multiples of `shift-vector` from `min-multiple` to `max-multiple`, inclusive, in steps of 1.

**`(geometric-objects-duplicates shift-vector min-multiple max-multiple obj-list)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Same as `geometric-object-duplicates`, except operates on a list of objects, `obj-list`. If *A* appears before *B* in the input list, then all the duplicates of *A* appear before all the duplicates of *B* in the output list.

**`(geometric-objects-lattice-duplicates obj-list [ ux uy uz ])`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Duplicates the objects in `obj-list` by multiples of the Cartesian basis vectors, making all possible shifts of the "primitive cell" (see below) that fit inside the lattice cell. The primitive cell to duplicate is `ux` by `uy` by `uz`, in units of the Cartesian basis vectors. These three parameters are optional; any that you do not specify are assumed to be `1`.

**`(point-in-object? point obj)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Returns whether or not the given 3-vector `point` is inside the geometric object `obj`.

**`(point-in-periodic-object? point obj)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
As `point-in-object?`, but also checks translations of the given object by the lattice vectors.

**`(display-geometric-object-info indent-by obj)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Outputs some information about the given `obj`, indented by `indent-by` spaces.

### Output File Names

The output file names used by Meep, e.g. for HDF5 files, are automatically prefixed by the input variable `filename-prefix`. If `filename-prefix` is `""` (the default), however, then Meep constructs a default prefix based on the current ctl file name with `".ctl"` replaced by `"-"`: e.g. `tst.ctl` implies a prefix of `"tst-"`. You can get this prefix by running:

**`(get-filename-prefix)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Return the current prefix string that is prepended, by default, to all file names.

If you don't want to use any prefix, then you should set `filename-prefix` to `false`.

In addition to the filename prefix, you can also specify that all the output files be written into a newly-created directory (if it does not yet exist). This is done by running:

**`(use-output-directory [dirname])`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Put output in a subdirectory, which is created if necessary. If the optional argument dirname is specified, that is the name of the directory. Otherwise, the directory name is the current ctl file name with `".ctl"` replaced by `"-out"`: e.g. `tst.ctl` implies a directory of `"tst-out"`.

### Misc.

**`(volume (center ...) (size ...))`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Many Meep functions require you to specify a volume in space, corresponding to the C++ type `meep::geometric_volume`. This function creates such a volume object, given the `center` and `size` properties (just like e.g. a `block` object). If the `size` is not specified, it defaults to (0,0,0), i.e. a single point.

**`(meep-time)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Return the current simulation time in simulation time units, not wall-clock time!. (e.g. during a run function.)

Occasionally, e.g. for termination conditions of the form *time* &lt; *T*?, it is desirable to round the time to single precision in order to avoid small differences in roundoff error from making your results different by one timestep from machine to machine (a difference much bigger than roundoff error); in this case you can call `(meep-round-time)` instead, which returns the time rounded to single precision.

### Field Computations

Meep supports a large number of functions to perform computations on the fields. Currently, most of them are accessed via the lower-level C++/SWIG interface, but we are slowly adding simpler, higher-level versions of them here.

**`(get-field-point c pt)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a `component` or `derived-component` constant `c` and a `vector3` `pt`, returns the value of that component at that point.

**`(get-epsilon-point pt)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Equivalent to `(get-field-point Dielectric pt)`.

**`(flux-in-box dir box)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a `direction` constant, and a `meep::volume*`, returns the flux (the integral of $\Re [\mathbf{E}^* \times \mathbf{H}]$) in that volume. Most commonly, you specify a volume that is a plane or a line, and a direction perpendicular to it, e.g. `(flux-in-box` `X` `(volume` `(center` `0)` `(size` `0` `1` `1)))`.

**`(electric-energy-in-box box)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a `meep::volume*`, returns the integral of the electric-field energy $\mathbf{E}^* \cdot \mathbf{D}/2$ in the given volume. (If the volume has zero size along a dimension, a lower-dimensional integral is used.)

**`(magnetic-energy-in-box box)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a `meep::volume*`, returns the integral of the magnetic-field energy $\mathbf{H}^* \cdot \mathbf{B}/2$in the given volume. (If the volume has zero size along a dimension, a lower-dimensional integral is used.)

**`(field-energy-in-box box)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a `meep::volume*`, returns the integral of the electric+magnetic-field energy $\mathbf{E}^* \cdot \mathbf{D}/2 + \mathbf{H}^* \cdot \mathbf{B}/2$in the given volume. (If the volume has zero size along a dimension, a lower-dimensional integral is used.)

Note that if you are at a fixed frequency and you use complex fields (Bloch-periodic boundary conditions or `fields-complex?=true`), then one half of the flux or energy integrals above corresponds to the time-average of the flux or energy for a simulation with real fields.

Often, you want the integration box to be the entire computational cell. A useful function to return this box, which you can then use for the `box` arguments above, is `(meep-fields-total-volume` `fields)`, where `fields` is the global variable (above) holding the current `meep::fields` object.

One powerful feature is that you can supply an arbitrary function $f(\mathbf{x},c_1,c_2,\ldots)$ of position $\mathbf{x}$ and various field components $c_1,\ldots$ and ask Meep to integrate it over a given volume, find its maximum, or output it (via `output-field-function`, described later). This is done via the functions:

**`(integrate-field-function cs func [where] [fields-var])`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Returns the integral of the complex-valued function `func` over the `meep::geometric_volume` specified by `where` (defaults to entire computational cell) for the `meep::fields` specified by `fields-var` (defaults to `fields`). `func` is a function of position (a `vector3`, its first argument) and zero or more field components specified by `cs`: a list of `component` constants. `func` can be real- or complex-valued.

If any dimension of `where` is zero, that dimension is not integrated over. In this way you can specify one-, two-, or three-dimensional integrals.

**`(max-abs-field-function cs func [where] [fields-var])`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
As `integrate-field-function`, but returns the maximum absolute value of `func` in the volume `where` instead of its integral.

The integration is performed by summing over the grid points with a simple trapezoidal rule, and the maximum is similarly over the grid points. See also [Field Function Examples](Field_Function_Examples.md) for illustrations of how to call `integrate-field-function` and `max-abs-field-function`. See also [Synchronizing the Magnetic and Electric Fields](Synchronizing_the_Magnetic_and_Electric_Fields.md) if you want to do computations combining the electric and magnetic fields.

Occasionally, one wants to compute an integral that combines fields from two separate simulations (e.g. for nonlinear coupled-mode calculations). This functionality is supported in Meep, as long as the two simulations have the *same* computational cell, the same resolution, the same boundary conditions and symmetries (if any), and the same PML layers (if any).

**`(integrate2-field-function fields2 cs1 cs2 func [where] [fields-var])`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Similar to `integrate-field-function`, but takes additional parameters `fields2` and `cs2`. `fields2` is a `meep::fields*` object similar to the global `fields` variable (see below) specifying the fields from another simulation. `cs1` is a list of components to integrate with from `fields-var` (defaults to `fields`), as for `integrate-field-function`, while `cs2` is a list of components to integrate from `fields2`. Similar to `integrate-field-function`, `func` is a function that returns an number given arguments consisting of: the position vector, followed by the values of the components specified by `cs1` (in order), followed by the values of the components specified by `cs2` (in order).

To get two fields in memory at once for `integrate2-field-function`, the easiest way is to run one simulation within a given .ctl file, then save the results in another fields variable, then run a second simulation. This would look something like:

```
...set up and run first simulation...
(define fields2 fields) ; save the fields in a variable
(set! fields '()) ; prevent the fields from getting deallocated by reset-meep
(reset-meep)
...set up and run second simulation...
```

It is also possible to timestep both fields simultaneously (e.g. doing one timestep of one simulation then one timestep of another simulation, and so on, but this requires you to call much lower-level functions like `(meep-fields-step` `fields)`.

### Reloading Parameters

Once the fields/simulation have been initialized, you can change the values of various parameters by using the following functions:

**`(reset-meep)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Reset all of Meep's parameters, deleting the fields, structures, etcetera, from memory as if you had not run any computations.

**`(restart-fields)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Restart the fields at time zero, with zero fields. (Does *not* reset the Fourier transforms of the flux planes, which continue to be accumulated.)

**`(change-k-point! k)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Change the `k-point` (the Bloch periodicity).

**`(change-sources! new-sources)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Change the `sources` input variable to `new-sources`, and changes the sources used for the current simulation.

### Flux Spectra

Given a bunch of `flux-region` objects (see above), you can tell Meep to accumulate the Fourier transforms of the fields in those regions in order to compute flux spectra. See also the [transmission/reflection spectra introduction](Introduction.md#transmissionreflection-spectra) and the [tutorial](Scheme_Tutorial.md). The most important function is:

**`(add-flux fcen df nfreq flux-regions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Add a bunch of `flux-region`s to the current simulation (initializing the fields if they have not yet been initialized), telling Meep to accumulate the appropriate field Fourier transforms for `nfreq` equally spaced frequencies covering the frequency range `fcen-df/2` to `fcen+df/2`. Return a *flux object*, which you can pass to the functions below to get the flux spectrum, etcetera.

As described in the tutorial, you normally use `add-flux` via statements like:

**`(define transmission (add-flux ...))`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
to store the flux object in a variable. `add-flux` initializes the fields if necessary, just like calling `run`, so you should only call it *after* setting up your `geometry`, `sources`, `pml-layers`, etcetera. You can create as many flux objects as you want, e.g. to look at powers flowing in different regions or in different frequency ranges. Note, however, that Meep has to store (and update at every time step) a number of Fourier components equal to the number of grid points intersecting the flux region multiplied by the number of electric and magnetic field components required to get the Poynting vector multiplied by `nfreq`, so this can get quite expensive (in both memory and time) if you want a lot of frequency points over large regions of space.

Once you have called `add-flux`, the Fourier transforms of the fields are accumulated automatically during time-stepping by the `run` functions. At any time, you can ask for Meep to print out the current flux spectrum via:

**`(display-fluxes fluxes...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a number of flux objects, this displays a comma-separated table of frequencies and flux spectra, prefixed by "flux1:" or similar (where the number is incremented after each run). All of the fluxes should be for the same `fcen`/`df`/`nfreq`. The first column are the frequencies, and subsequent columns are the flux spectra.

You might have to do something lower-level if you have multiple flux regions corresponding to *different* frequency ranges, or have other special needs. `(display-fluxes` `f1` `f2` `f3)` is actually equivalent to `(display-csv` `"flux"` `(get-flux-freqs` `f1)` `(get-fluxes` `f1)` `(get-fluxes` `f2)` `(get-fluxes` `f3)`, where `display-csv` takes a bunch of lists of numbers and prints them as a comma-separated table, and we are calling two lower-level functions:

**`(get-flux-freqs flux)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a flux object, returns a list of the frequencies that it is computing the spectrum for.

**`(get-fluxes flux)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a flux object, returns a list of the current flux spectrum that it has accumulated.

As described in the [tutorial](Scheme_Tutorial.md), for a reflection spectrum you often want to save the Fourier-transformed fields from a "normalization" run and then load them into another run to be subtracted. This can be done via:

**`(save-flux filename flux)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Save the Fourier-transformed fields corresponding to the given flux object in an HDF5 file of the given name (without the ".h5" suffix) (the current filename-prefix is prepended automatically).

**`(load-flux filename flux)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Load the Fourier-transformed fields into the given flux object (replacing any values currently there) from an HDF5 file of the given name (without the ".h5" suffix) (the current filename-prefix is prepended automatically). You must load from a file that was saved by `save-flux` in a simulation of *the same dimensions* (for both the computational cell and the flux regions) with the *same number of processors*.

**`(load-minus-flux filename flux)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
As `load-flux`, but negates the Fourier-transformed fields after they are loaded. This means that they will be *subtracted* from any future field Fourier transforms that are accumulated.

**`(scale-flux-fields s flux)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Scale the Fourier-transformed fields in `flux` by the complex number `s`. e.g. `load-minus-flux` is equivalent to `load-flux` followed by `scale-flux-fields` with `s=-1`.

### Force Spectra

Very similar to flux spectra, you can also compute **force spectra**: forces on an object as a function of frequency, computed by Fourier transforming the fields and integrating the vacuum [Maxwell stress tensor](https://en.wikipedia.org/wiki/Maxwell_stress_tensor)

$$\sigma_{ij} = E_i^*E_j + H_i^*H_j - \frac{1}{2} \delta_{ij} \left( |\mathbf{E}|^2 + |\mathbf{H}|^2 \right)$$

over a surface *S* via $\mathbf{F} = \int_S \sigma d\mathbf{A}$. We recommend that you normally **only evaluate the stress tensor over a surface lying in vacuum**, as the interpretation and definition of the stress tensor in arbitrary media is often problematic (the subject of extensive and controversial literature). It is fine if the surface *encloses* an object made of arbitrary materials, as long as the surface itself is in vacuum.

See also: [Optical Forces Tutorial](Scheme_Tutorials/Optical_Forces.md).

Most commonly, you will want to **normalize** the force spectrum in some way, just as for flux spectra. Most simply, you could divide two different force spectra to compute the ratio of forces on two objects. Often, you will divide a force spectrum by a flux spectrum, to divide the force *F* by the incident power *P* on an object, in order to compute the useful dimensionless ratio *Fc*/*P* where *c*=1 in Meep units. For example, it is a simple exercise to show that the force *F* on a perfectly reflecting mirror with normal-incident power *P* satisfies *Fc*/*P*=2, and for a perfectly absorbing (black) surface *Fc*/*P*=1.

The usage is similar to the flux spectra: you define a set of `force-region` objects telling Meep where it should compute the Fourier-transformed fields and stress tensors, and call `add-force` to add these regions to the current simulation over a specified frequency bandwidth, and then use `display-forces` to display the force spectra at the end. There are also `save-force`, `load-force`, and `load-minus-force` functions that you can use to subtract the fields from two simulation, e.g. in order to compute just the force from scattered fields, similar to the flux spectra. These types and functions are defined as follows:

**`force-region`**  

A region (volume, plane, line, or point) in which to compute the integral of the stress tensor of the Fourier-transformed fields. Its properties are:

**`center [ vector3 ]`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The center of the force region (no default).

**`size [ vector3 ]`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The size of the force region along each of the coordinate axes; default is (0,0,0) (a single point).

**`direction [ direction constant ]`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The direction of the force that you wish to compute (e.g. `X`, `Y`, etcetera). Unlike `flux-region`, you must specify this explicitly, because there is not generally any relationship between the direction of the force and the orientation of the force region.

**`weight [ cnumber ]`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
A weight factor to multiply the force by when it is computed; default is `1.0`.

In most circumstances, you should define a set of `force-region`s whose union is a closed surface (lying in vacuum and enclosing the object that is experiencing the force).

**`(add-force fcen df nfreq force-regions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Add a bunch of `force-region`s to the current simulation (initializing the fields if they have not yet been initialized), telling Meep to accumulate the appropriate field Fourier transforms for `nfreq` equally spaced frequencies covering the frequency range `fcen-df/2` to `fcen+df/2`. Return a *force object*, which you can pass to the functions below to get the force spectrum, etcetera.

As for flux regions, you normally use `add-force` via statements like:

```
(define Fx (add-force ...))
```

to store the flux object in a variable. `add-force` initializes the fields if necessary, just like calling `run`, so you should only call it *after* setting up your `geometry`, `sources`, `pml-layers`, etcetera. You can create as many force objects as you want, e.g. to look at forces on different objects, in different directions, or in different frequency ranges. Note, however, that Meep has to store (and update at every time step) a number of Fourier components equal to the number of grid points intersecting the force region, multiplied by the number of electric and magnetic field components required to get the stress vector, multiplied by `nfreq`, so this can get quite expensive (in both memory and time) if you want a lot of frequency points over large regions of space.

Once you have called `add-force`, the Fourier transforms of the fields are accumulated automatically during time-stepping by the `run` functions. At any time, you can ask for Meep to print out the current force spectrum via:

**`(display-forces forces...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a number of force objects, this displays a comma-separated table of frequencies and force spectra, prefixed by "force1:" or similar (where the number is incremented after each run). All of the forces should be for the same `fcen`/`df`/`nfreq`. The first column are the frequencies, and subsequent columns are the force spectra.

You might have to do something lower-level if you have multiple force regions corresponding to *different* frequency ranges, or have other special needs. `(display-forces` `f1` `f2` `f3)` is actually equivalent to `(display-csv` `"force"` `(get-force-freqs` `f1)` `(get-forces` `f1)` `(get-forces` `f2)` `(get-forces` `f3)`, where `display-csv` takes a bunch of lists of numbers and prints them as a comma-separated table, and we are calling two lower-level functions:

**`(get-force-freqs flux)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a force object, returns a list of the frequencies that it is computing the spectrum for.

**`(get-forces flux)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a force object, returns a list of the current force spectrum that it has accumulated.

As described in the [tutorial](Scheme_Tutorial.md), to compute the force from scattered fields often want to save the Fourier-transformed fields from a "normalization" run and then load them into another run to be subtracted. This can be done via:

**`(save-force filename force)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Save the Fourier-transformed fields corresponding to the given force object in an HDF5 file of the given name (without the ".h5" suffix) (the current filename-prefix is prepended automatically).

**`(load-force filename force)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Load the Fourier-transformed fields into the given force object (replacing any values currently there) from an HDF5 file of the given name (without the ".h5" suffix) (the current filename-prefix is prepended automatically). You must load from a file that was saved by `save-force` in a simulation of *the same dimensions* (for both the computational cell and the force regions) with the *same number of processors*.

**`(load-minus-force filename force)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
As `load-force`, but negates the Fourier-transformed fields after they are loaded. This means that they will be *subtracted* from any future field Fourier transforms that are accumulated.

### LDOS spectra

Meep can also calculate the LDOS (local density of states) spectrum, as described in the [LDOS tutorial](Scheme_Tutorials/Local_Density_of_States.md). To do this, you simply pass the following step function to your `run` command:

**`(dft-ldos fcen df nfreq)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Compute the power spectrum of the sources (usually a single point dipole source), normalized to correspond to the LDOS, in a frequency bandwith `df` centered at `fcen`, at `nfreq` frequency points.

The resulting spectrum is outputted as comma-delimited text, prefixed by `ldos:,`, and is also stored in the `dft-ldos-data` global variable after the `run` is complete.

Analytically, the per-polarization LDOS is exactly proportional to the power radiated by an $\ell$-oriented point-dipole current, $p(t)$, at a given position in space. For a more mathematical treatment of the theory behind the LDOS, we refer you to the relevant discussion in [chapter 4](http://arxiv.org/abs/1301.5366) of our [book](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707), but for now we simply give the definition:

$$\operatorname{LDOS}_{\ell}(\vec{x}_0,\omega)=-\frac{2}{\pi}\varepsilon(\vec{x}_0)\frac{\operatorname{Re}[\hat{E}_{\ell}(\vec{x}_0,\omega)\hat{p}(\omega)^*]}{|\hat{p}(\omega)|^2}$$

where the $|\hat{p}(\omega)|^2$ normalization is necessary for obtaining the power exerted by a unit-amplitude dipole (assuming linear materials), and hats denote Fourier transforms. It is this quantity that is computed by the `dft-ldos` command for a single dipole source. For a volumetric source, the numerator and denominator are both integrated over the current volume, but "LDOS" computation is less meaningful in this case.

### Near-to-Far-Field Spectra

Meep can compute a "near-to-far-field transformation" in the frequency domain as described in this [tutorial](Scheme_Tutorials/Near_to_Far_Field_Spectra.md): given the fields on a "near" bounding surface inside the computational cell, it can compute the fields arbitrarily far away using an analytical transformation, assuming that the "near" surface and the "far" region lie in a single homogeneous non-periodic 2d or 3d region. That is, in a simulation *surrounded by PML* that absorbs outgoing waves, the near-to-far-field feature can compute the fields outside the computational cell *as if* the outgoing waves had not been absorbed (i.e. in the fictitious infinite open volume). Moreover, this operation is performed on the Fourier-transformed fields: like the flux and force spectra above, you specify a set of desired frequencies, Meep accumulates the Fourier transforms, and then Meep computes the fields at *each frequency* for the desired far-field points.

This is based on the [principle of equivalence](http://arxiv.org/abs/1301.5366) — given the Fourier-transformed tangential fields on the "near" surface, Meep computes equivalent currents and convolves them with the analytical Green's functions in order to compute the fields at any desired point in the "far" region.

There are three steps to using the near-to-far-field feature: first, define the "near" surface(s) as a set of surfaces capturing *all* outgoing radiation in the desired direction(s); second, run the simulation, typically with a pulsed source, to allow Meep to accumulate the Fourier transforms on the near surface(s); third, tell Meep to compute the far fields at any desired points (optionally saving the far fields from a grid of points to an HDF5 file). To define the near surfaces, use:

**`(add-near2far fcen df nfreq near2far-regions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Add a bunch of `near2far-region`s to the current simulation (initializing the fields if they have not yet been initialized), telling Meep to accumulate the appropriate field Fourier transforms for `nfreq` equally spaced frequencies covering the frequency range `fcen-df/2` to `fcen+df/2`. Return a *near2far object*, which you can pass to the functions below to get the far fields.

Each `near2far-region` is identical to `flux-region` except for the name: in 3d, these give a set of planes (**important:** all these "near surfaces" must lie in a single *homogeneous* material with *isotropic* ε and μ — and they should *not* lie in the PML regions) surrounding the source(s) of outgoing radiation that you want to capture and convert to a far field. Ideally, these should form a closed surface, but in practice it is sufficient for the `near2far-region`s to capture all of the radiation in the direction of the far-field points. **Important:** as for flux computations, each `near2far-region` should be assigned a `weight` of ±1 indicating the direction of the outward normal relative to the +coordinate direction. So, for example, if you have six regions defining the six faces of a cube, i.e. the faces in the +x, -x, +y, -y, +z, and -z directions, then they should have weights +1, -1, +1, -1, +1, and -1 respectively. Note that, neglecting discretization errors, all near-field surfaces that enclose the same outgoing fields are equivalent and will yield the same far fields (with a discretization-induced difference that vanishes with increasing resolution etc.).

*After* the simulation `run` is complete, you can compute the far fields. This is usually for a pulsed source .so that the fields have decayed away and the Fourier transforms have finished accumulating.

**`(get-farfield near2far x)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a `vector3` point `x` (which can lie anywhere outside the near-field surface, including outside the computational cell) and a near2far object, returns the computed (Fourier-transformed) "far" fields at x as list of length 6`nfreq`, consisting of fields (Ex1,Ey1,Ez1,Hx1,Hy1,Hz1,Ex2,Ey2,Ez2,Hx2,Hy2,Hz2,...) for the frequencies 1,2,…,`nfreq`.

**`(output-farfields near2far fname where resolution)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given an HDF5 file name `fname` (does *not* include the `.h5` suffix), a `volume` given by `where` (may be 0d, 1d, 2d, or 3d), and a `resolution` (in grids points / distance *a*), outputs the far fields in `where` (which may lie *outside* the computational cell) in a grid with the given resolution (which may differ from the FDTD grid resolution) to the HDF5 file as a set of twelve array datasets `ex.r`, `ex.i`, ..., `hz.r`, `hz.i`, giving the real and imaginary parts of the Fourier-transformed E and H fields on this grid. Each dataset is an nx×ny×nz×nfreq 4-dimensional array of space×frequency (although dimensions that =1 are omitted).

Note that far fields have the same units and scaling as the *Fourier transforms* of the fields, and hence cannot be directly compared to time-domain fields. In practice, it is easiest to use the far fields in computations where overall scaling (units) cancel out or are irrelevant, e.g. to compute the fraction of the far fields in one region vs. another region.

For a scattered-field computation, you often want to separate the scattered and incident fields. Just as is described in the [tutorial](Scheme_Tutorial.md) for flux computations, you can do this by saving the Fourier-transformed incident from a "normalization" run and then load them into another run to be subtracted. This can be done via:

**`(save-near2far filename near2far)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Save the Fourier-transformed fields corresponding to the given near2far object in an HDF5 file of the given name (without the ".h5" suffix) (the current filename-prefix is prepended automatically).

**`(load-near2far filename near2far)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Load the Fourier-transformed fields into the given near2far object (replacing any values currently there) from an HDF5 file of the given name (without the ".h5" suffix) (the current filename-prefix is prepended automatically). You must load from a file that was saved by `save-near2far` in a simulation of *the same dimensions* (for both the computational cell and the near2far regions) with the *same number of processors*.

**`(load-minus-near2far filename near2far)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
As `load-near2far`, but negates the Fourier-transformed fields after they are loaded. This means that they will be *subtracted* from any future field Fourier transforms that are accumulated.

**`(scale-near2far-fields s near2far)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Scale the Fourier-transformed fields in `near2far` by the complex number `s`. e.g. `load-minus-near2far` is equivalent to `load-near2far` followed by `scale-near2far-fields` with `s=-1`.

### Frequency-Domain Solver

Meep contains a frequency-domain solver that directly computes the fields produced in a geometry in response to a constant-frequency source, using an [iterative linear solver](http://www.netlib.org/linalg/html_templates/Templates.html) instead of timestepping. Preliminary tests have shown that in many instances, this solver converges much faster than simply running an equivalent time domain simulation with a continuous wave source, timestepping until all transient effects from the source turn-on have disappeared, especially if the fields are desired to a very high accuracy. To use it, simply define a `continuous-src` with the desired frequency, [initialize the fields and geometry](#Initializing_the_Structure_and_Fields) via `(init-fields)`, and then:

**`(meep-fields-solve-cw fields tol maxiters L)`**

After the `fields` variable (a global variable pointing to the `meep::fields*` object initialized by `init-fields`, [see above](Scheme_User_Interface.md#input-variables)), the next two parameters to the frequency-domaine solver are the tolerance `tol` for the iterative solver (10<sup>−8</sup>, by default) and a maximum number of iterations `maxiters` (10<sup>4</sup>, by default). Finally, there is a parameter `L` that determines a tradeoff between memory and work per step and convergence rate of the iterative algorithm [biCGSTAB-(L)](http://www.math.uu.nl/people/sleijpen/CGSTAB_software/CGSTAB.html) that is used; larger values of `L` of will often lead to faster convergence at the expense of more memory and more work per iteration. The default is `L`=2, and normally a value ≥ 2 should be used.

The frequency-domain solver supports arbitrary geometries, PML, boundary conditions, symmetries, parallelism, conductors, and arbitrary nondispersive materials. Lorentz–Drude dispersive materials are not currently supported in the frequency-domain solver, but since you are solving at a known fixed frequency rather than timestepping, you should be able to pick conductivities etcetera in order to obtain any desired complex ε and μ at that frequency.

The frequency-domain solver requires you to use complex-valued fields, via `(set! force-complex-fields? true)`.

After `meep-fields-solve-cw` completes, it should be as if you had just run the simulation for an infinite time with the source at that frequency. You can call the various field-output functions and so on as usual at this point.

Run and Step Functions
----------------------

The actual work in Meep is performed by *run* functions, which time-step the simulation for a given amount of time or until a given condition is satisfied.

The run functions, in turn, can be modified by use of *step functions*: these are called at every time step and can perform any arbitrary computation on the fields, do outputs and I/O, or even modify the simulation. The step functions can be transformed by many *modifier functions*, like *at-beginning*, *during-sources*, etcetera which cause them to only be called at certain times, etcetera, instead of at every time step.

A common point of confusion is described in the article: [The Run Function Is Not A Loop](The_Run_Function_Is_Not_A_Loop.md). Please read this article if you want to make Meep do some customized action on each time step, as many users make the same mistake. What you really want to in that case is to write a step function, as described below.

### Run Functions

The following run functions are available. (You can also write your own, using the lower-level [C++/SWIG functions](#swig-wrappers), but these should suffice for most needs.)

**`(run-until cond?/time step-functions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Run the simulation until a certain time or condition, calling the given step functions (if any) at each timestep. The first argument is *either* a number, in which case it is an additional time (in Meep units) to run for, *or* it is a function (of no arguments) which returns `true` when the simulation should stop.

**`(run-sources step-functions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Run the simulation until all sources have turned off, calling the given step functions (if any) at each timestep. Note that this does *not* mean that the fields will be zero at the end: in general, some fields will still be bouncing around that were excited by the sources.

**`(run-sources+ cond?/time step-functions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
As `run-sources`, but with an additional first argument: either a number, in which case it is an *additional* time (in Meep units) to run for after the sources are off, *or* it is a function (of no arguments). In the latter case, the simulation runs until the sources are off *and* `(cond?)` returns `true`.

In particular, a useful first argument to `run-sources+` or `run-until` is often given by as in the [tutorial](Scheme_Tutorial.md):

**`(stop-when-fields-decayed dT c pt decay-by)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Return a `cond?` function, suitable for passing to `run-until`/`run-sources+`, that examines the component `c` (e.g. `Ex`, etc.) at the point `pt` (a `vector3`) and keeps running until its absolute value *squared* has decayed by at least `decay-by` from its maximum previous value. In particular, it keeps incrementing the run time by `dT` (in Meep units) and checks the maximum value over that time period—in this way, it won't be fooled just because the field happens to go through 0 at some instant.

Note that, if you make `decay-by` very small, you may need to increase the `cutoff` property of your source(s), to decrease the amplitude of the small high-frequency components that are excited when the source turns off. High frequencies near the [Nyquist frequency](https://en.wikipedia.org/wiki/Nyquist_frequency) of the grid have slow group velocities and are absorbed poorly by [PML](Perfectly_Matched_Layer.md).

Finally, another two run functions, useful for computing ω(**k**) band diagrams, are

**`(run-k-point T k)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a `vector3 k`, runs a simulation for each *k* point (i.e. specifying Bloch-periodic boundary conditions) and extracts the eigen-frequencies, and returns a list of the (complex) frequencies. In particular, you should have specified one or more Gaussian sources. It will run the simulation until the sources are turned off plus an additional `T` time units. It will run `harminv` (see below) at the same point/component as the first Gaussian source and look for modes in the union of the frequency ranges for all sources.

**`(run-k-points T k-points)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a list `k-points` of *k* vectors, runs `run-k-point` for each one, and returns a list of lists of frequencies (one list of frequencies for each *k*). Also prints out a comma-delimited list of frequencies, prefixed by `freqs:`, and their imaginary parts, prefixed by `freqs-im:`. See e.g. this [band diagram tutorial](Scheme_Tutorials/Resonant_Modes_and_Transmission_in_a_Waveguide_Cavity.md).

### Predefined Step Functions

Several useful step functions are predefined for you by Meep.

#### Output Functions

The most common step function is an output function, which outputs some field component to an [HDF5](https://en.wikipedia.org/wiki/HDF5) file. Normally, you will want to modify this by one of the `at-` functions, below, as outputting a field at *every* time step can get quite time- and storage-consuming.

Note that although the various field components are stored at different places in the [Yee lattice](Yee_Lattice.md), when they are outputted they are all linearly interpolated to the same grid: to the points at the *centers* of the Yee cells, i.e. $(i+0.5,j+0.5,k+0.5)\cdot\Delta$ in 3d.

The predefined output functions are:

**`output-epsilon`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Output the dielectric function (relative permittivity) ε. Note that this only outputs the frequency-independent part of ε (the $\omega\to\infty$ limit).

**`output-mu`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Output the relative permeability function μ. Note that this only outputs the frequency-independent part of μ (the $\omega\to\infty$ limit).

**`output-hpwr`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Output the magnetic-field energy density $\mathbf{H}^* \cdot \mathbf{B} / 2$

**`output-dpwr`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Output the electric-field energy density $\mathbf{E}^* \cdot \mathbf{D} / 2$

**`output-tot-pwr`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Output the total electric and magnetic energy density. Note that you might want to wrap this step function in `synchronized-magnetic` to compute it more accurately. See [Synchronizing the Magnetic and Electric Fields](Synchronizing_the_Magnetic_and_Electric_Fields.md).

**`output-Xfield-x, output-Xfield-y, output-Xfield-z, output-Xfield-r, output-Xfield-p`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Output the $x$, $y$, $z$, $r$, or $\phi$ component respectively, of the field *X*, where *X* is either `h`, `b`, `e`, `d`, or `s` for the magnetic, electric, displacement, or Poynting field, respectively. If the field is complex, outputs two datasets, e.g. `ex.r` and `ex.i`, within the same HDF5 file for the real and imaginary parts, respectively. Note that for outputting the Poynting field, you might want to wrap the step function in `synchronized-magnetic` to compute it more accurately; see [Synchronizing the magnetic and electric fields](Synchronizing_the_Magnetic_and_Electric_Fields.md).

**`output-Xfield`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Outputs *all* the components of the field *X*, where *X* is either `h`, `b`, `e`, `d`, or `s` as above, to an HDF5 file. (That is, the different components are stored as different datasets within the *same* file.)

**`(output-png component h5topng-options)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Output the given field component (e.g. `Ex`, etc.) as a [PNG](https://en.wikipedia.org/wiki/PNG) image, by first outputting the HDF5 file, then converting to PNG via [h5topng](http://ab-initio.mit.edu/wiki/index.php/H5utils), then deleting the HDF5 file. The second argument is a string giving options to pass to h5topng (e.g. `"-Zc` `bluered"`). See also the [tutorial](Scheme_Tutorial.md#output-tips-and-tricks).

It is often useful to use the `h5topng` `-C` or `-A` options to overlay the dielectric function when outputting fields. To do this, you need to know the name of the dielectric-function `.h5` file (which must have been previously output by output-epsilon). To make this easier, a built-in shell variable `$EPS` is provided which refers to the last-output dielectric-function `.h5` file. So, for example `(output-png` `Ez` `"-C` `$EPS")` will output the $E_z$ field and overlay the dielectric contours.

**`(output-png+h5 component h5topng-options)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Like `output-png`, but also outputs the `.h5` file for the component. (In contrast, `output-png` deletes the `.h5` when it is done.)

More generally, it is possible to output an arbitrary function of position and zero or more field components, similar to the `integrate-field-function` described above. This is done by:

**`(output-field-function name cs func)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Output the field function `func` to an HDF5 file in the datasets named *`name`*`.r` and *`name`*`.i` (for the real and imaginary parts). Similar to `integrate-field-function`, `func` is a function of position (a `vector3`) and the field components corresponding to `cs`: a list of `component` constants.

**`(output-real-field-function name cs func)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
As `output-field-function`, but only outputs the real part of `func` to the dataset given by the string `name`.

See also [Field Function Examples](Field_Function_Examples.md). See also [Synchronizing the Magnetic and Electric Fields](Synchronizing_the_Magnetic_and_Electric_Fields.md) if you want to do computations combining the electric and magnetic fields.

#### Harminv

The following step function collects field data from a given point and runs [Harminv](http://ab-initio.mit.edu/wiki/index.php/harminv) on that data to extract the frequencies, decay rates, and other information.

**`(harminv c pt fcen df [maxbands])`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Returns a step function that collects data from the field component `c` (e.g. `Ex`, etc.) at the given point `pt` (a `vector3`). Then, at the end of the run, it uses Harminv to look for modes in the given frequency range (center `fcen` and width `df`), printing the results to standard output (prefixed by `harminv:`) as comma-delimited text, and also storing them to the variable `harminv-results`. The optional argument `maxbands` is the maximum number of modes to search for; defaults to `100`.

**Important:** normally, you should only use `harminv` to analyze data *after the sources are off*. Wrapping it in `(after-sources` `(harminv` `...))` is sufficient.

In particular, Harminv takes the time series $f(t)$ corresponding to the given field component as a function of time and decomposes it (within the specified bandwidth) as:

$$f(t) = \sum_n a_n e^{-i\omega_n t}$$

The results are stored in the list `harminv-results`, which is a list of tuples holding the frequency, amplitude, and error of the modes. Given one of these tuples, you can extract its various components with one of the accessor functions:

**`(harminv-freq result)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Return the complex frequency ω (in the usual Meep 2πc units).

**`(harminv-freq-re result)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Return the real part of the frequency ω.

**`(harminv-freq-im result)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Return the imaginary part of the frequency ω.

**`(harminv-Q result)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Return dimensionless lifetime, or "quality factor", $Q$, defined as $-\mathrm{Re}\,\omega / 2 \mathrm{Im}\,\omega$.

**`(harminv-amp result)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Return the complex amplitude $a$.

**`(harminv-err result)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
A crude measure of the error in the frequency (both real and imaginary)...if the error is much larger than the imaginary part, for example, then you can't trust the $Q$ to be accurate. **Note**: *this error is only the uncertainty in the signal processing*, and tells you nothing about the errors from finite resolution, finite cell size, and so on!

For example, `(map harminv-freq-re harminv-results)` gives you a list of the real parts of the frequencies, using the Scheme built-in `map`.

### Step-Function Modifiers

Rather than writing a brand-new step function every time we want to do something a bit different, the following "modifier" functions take a bunch of step functions and produce *new* step functions with modified behavior. See also the [tutorial](Scheme_Tutorial.md) for examples.

#### Miscellaneous Step-Function Modifiers

**`(combine-step-funcs step-functions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given zero or more step functions, return a new step function that (on each step) calls all of the passed step functions.

**`(synchronized-magnetic step-functions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given zero or more step functions, return a new step function that (on each step) calls all of the passed step functions with the magnetic field synchronized in time with the electric field. See [Synchronizing the Magnetic and Electric Fields](Synchronizing_the_Magnetic_and_Electric_Fields.md).

#### Controlling When a Step Function Executes

**`(when-true cond? step-functions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given zero or more step functions and a condition function `cond?` ( a function of no arguments), evaluate the step functions whenever `(cond?)` returns `true`.

**`(when-false cond? step-functions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given zero or more step functions and a condition function `cond?` ( a function of no arguments), evaluate the step functions whenever `(cond?)` returns `false`.

**`(at-every dT step-functions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given zero or more step functions, evaluates them at every time interval of `dT` units (rounded up to the next time step).

**`(after-time T step-functions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given zero or more step functions, evaluates them only for times after a `T` time units have elapsed from the start of the run.

**`(before-time T step-functions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given zero or more step functions, evaluates them only for times before a `T` time units have elapsed from the start of the run.

**`(at-time T step-functions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given zero or more step functions, evaluates them only once, after a `T` time units have elapsed from the start of the run.

**`(after-sources step-functions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given zero or more step functions, evaluates them only for times after all of the sources have turned off.

**`(after-sources+ T step-functions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given zero or more step functions, evaluates them only for times after all of the sources have turned off, plus an additional `T` time units have elapsed.

**`(during-sources step-functions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given zero or more step functions, evaluates them only for times *before* all of the sources have turned off.

**`(at-beginning step-functions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given zero or more step functions, evaluates them only once, at the beginning of the run.

**`(at-end step-functions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given zero or more step functions, evaluates them only once, at the end of the run.

#### Modifying HDF5 Output

**`(in-volume v step-functions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given zero or more step functions, modifies any output functions among them to only output a subset (or a superset) of the computational cell, corresponding to the `meep::geometric_volume*` `v` (created by the `volume` function).

**`(in-point pt step-functions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given zero or more step functions, modifies any output functions among them to only output a single *point* of data, at `pt` (a `vector3`).

**`(to-appended filename step-functions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given zero or more step functions, modifies any output functions among them to *append* their data to datasets in a single newly-created file named `filename` (plus an `.h5` suffix and the current filename prefix). They append by adding an *extra dimension* to their datasets, corresponding to time.

**`(with-prefix prefix step-functions...)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given zero or more step functions, modifies any output functions among them to prepend the string `prefix` to the file names (much like `filename-prefix`, above).

### Writing Your Own Step Functions

A step function can take two forms. The simplest is just a function of no arguments, which is called at every time step (unless modified by one of the modifier functions above). e.g.

**`(define (my-step) (print "Hello world!\n"))`**  

If one then does `(run-until 100 my-step)`, Meep will run for 100 time units and print "Hello world!" at every time step.

This suffices for most purposes. However, sometimes you need a step function that opens a file, or accumulates some computation, and you need to clean up (e.g. close the file or print the results) at the end of the run. For this case, you can write a step function of one argument: that argument will either be `'step` when it is called during time-stepping, or `'finish` when it is called at the end of the run.

Low-Level Functions
-------------------

By default, Meep reads input functions like `sources` and `geometry` and creates *global* variables `structure` and `fields` to store the corresponding C++ objects. Given these, you can then call essentially *any* function in the C++ interface, because all of the C++ functions are automatically made accessible to Scheme by a wrapper-generator program called [SWIG](https://en.wikipedia.org/wiki/SWIG).

### Initializing the Structure and Fields

The `structure` and `fields` variables are automatically initialized when any of the run functions is called, or by various other functions such as `add-flux`. To initialize them separately, you can call `(init-fields)` manually, or `(init-structure k-point)` to just initialize the structure.

If you want to time step more than one field simultaneously, the easiest way is probably to do something like:

```
(init-fields)
(define my-fields fields)
(set! fields '())
(reset-meep)
```

and then change the geometry etc. and re-run `(init-fields)`. Then you'll have two field objects in memory.

### SWIG Wrappers

If you look at a function in the C++ interface, then there are a few simple rules to infer the name of the corresponding Scheme function.

-   First, all Meep functions (in the `meep::` namespace) are prefixed with `meep-` in the Scheme interface.
-   Second, any method of a class is prefixed with the name of the class and a hyphen. For example, `meep::fields::step`, which is the function that performs a time-step, is exposed to Scheme as `meep-fields-step`. Moreover, you pass the object as the first argument in the Scheme wrapper. e.g. `f.step()` becomes `(meep-fields-step` `f)`.
-   To call the C++ constructor for a type, you use `new-`*`type`*. e.g. `(new-meep-fields` `...arguments...)` returns a new `meep::fields` object. Conversely, to call the destructor and deallocate an object, you use `delete-`*`type`*; most of the time, this is not necessary because objects are automatically garbage-collected.

Some argument type conversion is performed automatically, e.g. types like complex numbers are converted to `complex<double>`, etcetera. `vector3` vectors are converted to `meep::vec`, but to do this we need to know the dimensionality of the problem in C++. The problem dimensions are automatically initialized by `init-structure`, but if you want to pass vector arguments to C++ before that time you should call `(require-dimensions!)`, which infers the dimensions from the `geometry-lattice`, `k-point`, and `dimensions` variables.
