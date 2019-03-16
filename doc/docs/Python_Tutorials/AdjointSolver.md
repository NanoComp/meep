<!----------------------------------------------------->
<!- LaTeX macros used only for adjoint documentation  ->
<!----------------------------------------------------->
$$
   \newcommand{\vb}{\mathbf}
   \newcommand{\wt}{\widetilde}
   \newcommand{\mc}{\mathcal}
   \newcommand{\bmc}[1]{\boldsymbol{\mathcal{#1}}}
   \newcommand{\sup}[1]{^{\text{#1}}}
   \newcommand{\pard}[2]{\frac{\partial #1}{\partial #2}}
   \newcommand{\VMV}[3]{ \Big\langle #1 \Big| #2 \Big| #3 \Big\rangle}
$$

---
# Automated design optimization via adjoint-based sensitivity analysis
---

This tutorial introduces <span class=SC>meep</span>'s support
for [*adjoint-based sensitivity analysis*](https://en.wikipedia.org/wiki/Adjoint_state_method)
to facilitate automated design optimization using derivative-based numerical optimizers.

## Overview
--------------------------------

### Adjoint-based optimization

A common task in electromagnetic engineering is to custom-tune the design
of some component of a system---a waveguide taper, a power splitter,
an input coupler, an antenna, etc.---to optimize the performance of the system
as defined by some problem-specific metric. For our purposes,
a "design" will consist of a specification of the spatially-varying
scalar permittivity $\epsilon(\mathbf x)$ in some subregion
of a <span class=SC>meep</span> geometry, and the performance metric
will be a physical quantity computed from frequency-domain
fields---a 
[power flux](../Python_User_Interface.md#get_fluxes),
an [energy density](../Python_User_Interface.md#dft_energy),
an
[eigenmode expansion coefficient](../Python_User_Interface.md#get_eigenmode_coefficients),
or perhaps some mathematical function of one or more of these
quantities. We will shortly present a smorgasbord of examples; for now,
perhaps a good one to have in mind is the hole cloak problem, in which a
chunk of material has been removed from an otherwise perfect waveguide
section, ruining the otherwise perfectly unidirectional (no scattering or reflection)
flow of power from a source at one end of the guide to a sink at the other;
our task is to tweak the permittivity in an annular region
surrounding the defect so as to restore as much as possible the reflectionless 
transfer of power across the waveguide---thus "hiding" or "cloaking"
the defect.

Now, given a candidate design
$\epsilon\sup{trial}(\mathbf{x})$, it's easy enough to see
how we can use <span class=SC>meep</span> to evaluate
our objective function---just define a
<span class=SC>meep</span> geometry with $\epsilon\sup{trial}$ as a
[spatially-varying permittivity function](../Python_User_Interface.md#eps_func)
in the design region,
add [DFT cells](../Python_User_Interface.md#FluxSpectra)
to tabulate frequency-domain fields in the regions of interest,
[timestep](../Python_User_Interface.md#RunStepFunctions) until 
the DFTs converge, and use post-processing routines like
[`get_fluxes()`](../Python_User_Interface.md#get_fluxes)
or perhaps
[`get_eigenmode_coefficients()`](../Python_User_Interface.md#get_eigenmode_coefficients)
to get the quantities needed to evaluate the performance of the device.
Thus, for the cost of one full <span class=SC>meep</span> timestepping
run we obtain the value of our objective function at one point
in the parameter space of possible inputs.

But *now* what do we do?! The difficulty is that the computation
just described furnishes only the *value* of the objective function
for a given input, not its *derivatives* with respect to the
design variables---and thus yields zero insight into how we should
tweak the design to improve performance.
In simple cases we might hope to get somewhere with just
engineering design intuition---like in the old days!---while
for small problems with just a few parameters we might try our luck with a
[derivative-free optimization algorithm](https://en.wikipedia.org/wiki/Derivative-free_optimization);
however, both of these approaches will run out of steam long before
we scale up to 
the full complexity of a practical problem with thousands
of degrees of freedom.
Alternatively, we could get approximate derivative information by brute-force
finite-differencing---slightly tweaking one design variable, repeating 
the full timestepping run, and asking how the results changed---but 
proceeding this way to compute derivatives with respect to all $D$ 
design variables would require fully $D$ separate timestepping runs;
for the problem sizes we have in mind, this would make calculating the 
objective-function gradient
*several thousand times* more costly than calculating its value.
So we face a dilemma: How can we obtain the derivative information
necessary for effective optimization in a reasonable amount of time?
This is where adjoints come to the rescue.

The *adjoint method* of sensitivity analysis is a technique in which
we exploit certain facts about the physics of a problem and the
consequent mathematical structure---specifically, in this case, the
linearity and reciprocity of Maxwell's equations---to rearrange the
calculation of derivatives in a way that yields an *enormous* speedup
over the brute-force finite-difference approach. More specifically,
after we have computed the objective-function value by doing
the full <span class=SC>meep</span> timestepping run mentioned
above---the "forward" run in adjoint-method parlance---we can magically
compute its derivatives with respect to *all* design variables by doing
just *one* additional timestepping run with a funny-looking choice
of sources and outputs (the "adjoint" run).
Thus, whereas gradient computation via finite-differencing is at least $D$
times more expensive than computing the objective function value,
with adjoints we get both value and gradient for roughly just *twice* the
cost of the value alone. Such a bargain! At this modest cost, derivative-based 
optimization becomes entirely feasible.

??? note "**More general materials**"
    <small>
    Although for simplicity we focus here on
    the case of isotropic, non-magnetic materials,
    the adjoint solver is also capable of optimizing
    geometries involving permeable ($\mu\ne 1$)
    and anisotropic
    (tensor-valued $\boldsymbol{\epsilon},\boldsymbol{\mu}$)
    media.
    </small>


### Examples of optimization problems

Throughout this tutorial we will refer to a running collection of simple optimization
problems to illustrate the mechanics of optimization in <span class=SC>meep</span>.



### Common elements of optimization geometries: Objective regions, objective functions, design regions, basis sets

The examples above, distinct though they all are, illustrate
some common features that will be present in every
<span class=SC>meep</span> optimization problem:

+   One or more [regions over which to tabulate frequency-domain fields (DFT cells)](../Python_User_Interface.md#dft_obj) for use in computing power fluxes, mode-expansion coefficients, and other frequency-domain quantities used in characterizing device performance.  Because these regions are used to evaluate objective functions, we refer to them as *objective regions.* 

??? note "**Objective regions may or may not have zero thickness**"
    In the examples above, it happens that all objective regions are one-dimensional
    (zero-thickness) flux monitors, indicated by magenta lines; in a 3D geometry they
    would be two-dimensional flux planes, still of zero thickness in the normal 
    direction.  However, objective regions may also be of nonzero thickness, as for
    instance if the objective function involves the [field energy in a box-shaped
    subregion of a geometry.](../Python_User_Interface.md#energy)

+    A specification of which quantities (power fluxes, mode coefficients, energies, etc.) 
     are to be computed for each objective region, and of how those quantities are to be
     crunched mathematically to yield a single number measuring device performance. We
     refer to the individual quantities as *objective quantities*, while the overall
     function that inputs multiple objective quantities and outputs a single numerical
     score is the *objective function.*
     For example,

+   A specification of the region over which the material design is to be
    optimized, i.e. the region in which the permittivity is given by the
    design quantity $\epsilon\sup{des}(\mathbf x)$.
    We refer to this as the *design region* $\mathcal{V}\sup{des}$.
    (In some problems the design region may be a union of two or more
    disconnected subregions, i.e. if we are trying to optimize the
    $\epsilon$ distribution in two or more separated subregions of 
    the geometry.)

+   Because the design variable $\epsilon\sup{des}(\mathbf x)$
    is a continuous function defined throughout a finite volume of space,
    technically it involves infinitely many degrees of freedom.
    To yield a finite-dimensional optimization problem, it is convenient
    to approximate $\epsilon\sup{des}$ as a finite expansion in some
    convenient set of basis functions, i.e.
    $$ \epsilon(\mathbf x) \equiv \sum_{n=1}^N \beta_n \mathcal{b}_n(\mathbf x),
       \qquad \mathbf x\in \mathcal{V}\sup{des},
    $$
    where $\{\mathcal{b}_n(\mathbf x)\} is a set of $N$ scalar-valued
    basis functions defined for $\mathbf x\in\mathcal{V}\sup{des}$.
    For adjoint optimization in <span class=SC>meep</span>, the
    basis set is chosen by the user (either from among a predefined collection of
    common basis sets, or as an arbitrary user-defined set),
    and the task of the optimizer becomes to determine
    numerical values for the $N$-vector of coefficients 
    $\boldsymbol{\beta}=\{\beta_n\},n=1,\cdots,N.$
    
### Mechanics of <span class=SC>meep</span> design optimization: A choice of two tracks
is 

### Objective regions and design regions
--------------------------------

The methods we'll illustrate are applicable to a broad class of
optimization problems in which we seek to optimize Poynting fluxes
and/or mode-expansion coefficients in one or more regions of a
geometry by tuning the structure---that is, the spatial
permittivity distribution $\epsilon(\mathbf x)$---of some portion
of the geometry.

More specifically, we'll consider problems with the following characteristics:

+   The objective function is of the form
    $$F\Big( \big\{S_n\big\}, \big\{|\alpha_{nm}|^2\big\}\Big)$$
    where

    + $\{S_n\}=\{S_1, \cdots, S_N\}$ are the Poynting fluxes
        at some user-specified set of $N$ 
        [flux regions](https://meep.readthedocs.io/en/latest/Python_User_Interface/#fluxregion)
        $\{\mathcal{V}^o_1, \cdots, \mathcal{V}^o_N\}$, which we call the *objective regions* for
        the optimization problem  

    + $\alpha_{nm}$ is the
      [eigenmode-expansion coefficient](https://meep.readthedocs.io/en/latest/Mode_Decomposition)
      for the $m$th eigenmode in the $n$th objective region    

    + $F$ is an arbitrary user-specified function of its inputs    

+   The quantity to be optimized is the permittivity distribution $\epsilon(\mathbf x)$
    within some user-specified [volume](https://meep.readthedocs.io/en/latest/Python_User_Interface/#volume)
    $\mathcal{V}^d$, which we call the *design region*.

+   The permittivity in the design region is expressed as an expansion
    in some (arbitrary, user-specified) set of basis functions:

    Here $\{\psi_p(\mathbf x)\}, p=1,\cdots,P$ is a set of $P$ scalar-valued basis functions,
    defined for $\mathbf x \in \mathcal{V}^d$ the objective region, and the expansion coefficients $\{a_p\}$
    are the variable parameters in our optimization.

Thus, problems for adjoint-based solvers in meep, you will specify

+ a list of objective regions $\{\mathcal{V}^o_n\}$

+ a design region $\mathcal{V}^d$    

+ a basis set $\{b_d(\mathbf{x})\}$    

+ an objective function $F$    

### Sample geometries

+ Warmup: circular defect in waveguide

A simple warmup problem is a uniform waveguide with a
circular inclusion, in which we want to optimize the 
permittivity distribution in the circle to maximize
power flux through the waveguide:


In this case we have only a single objective region,
and the objective function is simply

$$ F\Big\(\mathcal{V}^o_1\Big)=S_1.$$

+ Bidirectional coupler

## Basis functions and parameterized geometries

### Specifying basis sets

The basis of expansion functions for the permittivity in the design region
is described by a single routine that inputs the normalized coordinates of a point
$\mathbf{x}$ in the design region and returns a vector of length $P$
containing $[\psi_1(\mathbf{x}), \cdots, \psi_P(\mathbf{x})]$.
Here "normalized coordinates" refers to a shifted and scaled coordinate system
in which the design region is centered at the origin and has length 1 in
all directions; so that the normalized coordinates $\overline{\mathbf{x}}$
for points in the design region satisfy $-0.5<\overline{x}_i<0.5$ for all $i$.
then every component of $\overline{\mathbf{x}}$ lies in the range $[-0.5,0.5]$

Here are some examples of basis sets for 2D design regions:

+ A simple polynomial basis set with 4 functions $\{1,x,y,xy\}$:

```python
  def basis(pbar):
      x=pbar[0]
      y=pbar[1]
      return [1,x,y,x*y]
```

+ A basis set for the circular-hole example above
  that contains a single function that vanishes
  outside a circle of (normalized) radius 0.25
  and is 1 inside that circle:

```python
  def basis(pbar):
      x=pbar[0]
      y=pbar[1]
      return [1.0] if (x*x+y*y)<=0.25*0.25 else [0.0]
```

   (Note that normalized radius 0.5 would correspond to the
    maximal inscribed circle in the design region.)

The file `adjoint.py` defines some convenient predefined basis sets:

+ `fourier_basis(kxMax, kyMax)`

    2D Fourier basis. `kxMax, kyMax` are integers specifying the
    maximum spatial frequency in each direction as a multiple of the
    base frequency. The size of the basis is `(2*kxMax+1)*(2*kyMax+1).`
    (The return values of e.g. `fourier_basis(4,3)` is a function
    that may be passed as the `basis` parameter to functions like
    `custom_dielectric` below.)

+ `annular_basis(NR=2, kMax=2, rho_max=0.5, rho_min=0.0)`

    A basis for a disc or annular region, consisting of a
    sinusoids in the angular direction paired with Legendre polynomials
    in the radial direction.

### Parameterized geometries

Having defined a set of basis functions, a design region
with parametrized permittivity described by a basis-set function `basis`
and a vector of basis-coefficients `coeffs`
may be added to a meep geometry by including the following block among its objects: 

```python
 mp.Block(center=design_center, size=design_size,
                        material=custom_dielectric(design_center, design_size, basis, coeffs)
         )
```

Here `design_center, design_size` are the center and size of the design region
and `custom_dielectric` is a function defined in `adjoint.py` that defines
a meep material function given by the weighted basis set.

## Computing objective functions and gradients

Having defined an optimization problem, we compute the value of the
objective function, and its gradient with respect to the expansion coefficients,
by calling

```
  f,gradf = get_objective_and_gradient(sim, forward_sources,
                                       objective_regions, design_region, basis,
                                       objective)
```

where the inputs are

+ `sim`: `simulation` containing your geometry (but no sources, DFT regions, etc.)
+ `forward_sources`: sources for the fields considered by your objective function.
+ `objective_regions`: List of DFT regions on which arguments to your objective function are defined.
+ `design_region`: Design region of your problem
+ `basis`: Function describing basis set, as defined above
+ `objective`: Function describing objective function, as defined above

## Defining basis sets

### Predefined basis sets

**Nomenclature for sinusoids**

In what follows it will be convenient to introduce some
shorthand notation for sinusoids.
For $n=0,1,2,\cdots$ we define
$$ s_n(\theta)\equiv
   \begin{cases}
   cos \nu\theta , \qquad &n=2\nu   \text{even} \\
   sin \nu\theta , \qquad &n=2\nu-1 \text{odd}
   \end{cases}
$$
In general we will specify the size of a Fourier basis
by stating its nonnegative-integer-valued maximum frequency
$\nu\sup{max}\ge 0$, in which case the set contains
$(2\nu\sup{max}+1)$ basis functions
$\{s_0, \cdots, s_{2\nu\sup{max}}\}.
(The $n$ index ranges from 0 to $2\nu\sup{max}$ inclusive.)
The first few sinusoids are
$$ s_0(\theta)=1,
   \qquad
   s_1(\theta)=\sin \theta,
   \qquad
   s_2(\theta)=\cos \theta,
   \qquad
   s_3(\theta)=\sin 2\theta,
   \qquad
   s_4(\theta)=\cos 2\theta, \cdots
$$

#### Rectangular design regions: Plane-wave basis

#### Disc-shaped or annular design regions:
     Legendre-Fourier basis

```
    disc_basis(ell_max=2, nu_max=2, rhobar_max=0.5, rhobar_min=0.0)
```

Returns a basis suitable for a disc-shaped or annular region defined by
$$ \overline{\rho}\sup{min}\le  \overline{\rho} \overline{\rho}\sup{max},
   \qquad 0\le \theta < 2\pi
$$
Here $\overline{\rho}$ is a normalized radial coordinate running from 0
(at the center of the design region) to 1 (on the inscribed circle
of the design region).

The basis functions are products of Legendre polynomials in
$\overline{\rho$} with sinusoids in $\theta$:
$$ b_{\ell n}(\overline{\rho},\theta)\equiv P_\ell(u)s_n(\theta),
   \qquad 
   0\le \ell \le \ell\sup{max}, \qquad
   0\le n \le 2\nu\sup{max}
$$
where $u\equiv 2\overline{\rho}-1$ runs over the usual range
$[-1,1]$ of the Legendre polynomials as $\overline{\rho}$
runs from 0 to 1.$ The total number of basis functions
is $(\ell\sup{max}+1)(2\nu\sup{max}+1)$, and the index 
of $b_{\ell n}$ within the full basis set is
$ (2\nu\sup{max}+1)\ell + n.$ For example, with $\nu\sup{max}=2$ the
first few functions are
$$ \begin{array}{cccccccccc}
   b_0(\overline{\rho},\theta})&=&1,
   b_1(\overline{\rho},\theta})&=&\sin \theta,
   b_2(\overline{\rho},\theta})&=&\cos \theta,
   b_3(\overline{\rho},\theta})&=&\sin 2\theta,
   b_4(\overline{\rho},\theta})&=&\cos 2\theta,
\\
   b_5(\overline{\rho},\theta})&=&u,
   b_6(\overline{\rho},\theta})&=&u\sin \theta,
   b_7(\overline{\rho},\theta})&=&u\cos \theta,
   b_8(\overline{\rho},\theta})&=&u\sin 2\theta,
   b_9(\overline{\rho},\theta})&=&u\cos 2\theta.
   \end{array}
$$
