--8<-- "AdjointSolver/AdjointDocumentationStyleHeader.md"

---
# `meep.adjoint:` Adjoint sensitivity analysis for automated design optimization via 
---

This section of the <span class=SC>meep</span> documentation
covers `meep.adjoint,` a submodule of the <span class=SC>meep</span> python module
that implements an [*adjoint-based sensitivity solver*](https://en.wikipedia.org/wiki/Adjoint_state_method)
to facilitate automated design optimization via derivative-based numerical optimizers.

The `meep.adjoint` documentation is divided into a number of subsections:
 
+ This **Overview** page reviews some basic facts about adjoints and optimizers,
  outlines the steps needed to prepare a <span class=SC>meep</span>
  geometry for optimization, and sketches the mechanics of
  the `meep.adjoint` design process.
  (This page is designed to be a gentle introduction for the
  adjoint neophyte; experts may want only to skim it before
  skipping to the next section.)

+ The [**Reference Manual**](ReferenceManual.md) fills in the details of 
  the topics oulined on this page, spelling out exactly how to
  write the python script that will drive your optimization
  session.

+ The [**Example Gallery**](ExampleGallery.md) presents a number
  of worked examples that illustrate how `meep.adjoint` tackles
  practical problems in various settings.

+ The [**Implementation Notes**](ImplementationNotes.md) page
  offers a glimpse of what's behind the hood---the physical 
   and mathematical basis of the adjoint method and how they
  are implemented by `meep.adjoint.` An understanding of this 
  content is not strictly necessary to use the solver, but may
  help you get more out of the process.

+ Although logically independent of the adjoint-solver
  implementation, the `Visualization` package bundled
  with the `meep.adjoint` module offers several general-purpose
  utilities for convenient visualization of various aspects
  of <span class=SC>meep</span>calculations, which are
  useful in *any* meep calculation whether adjoint-related
  or not.


## Overview: Adjoint-based optimization

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
perhaps a good one to have in mind is the
[hole cloak ](#HoleCloak) discussed below, in which a
chunk of material has been removed from an otherwise perfect waveguide
section, ruining the otherwise perfectly unidirectional (no scattering or reflection)
flow of power from a source at one end of the guide to a sink at the other;
our task is to tweak the permittivity in an annular region
surrounding the defect (the *cloak*) so as to restore 
as much as possible the reflectionless transfer of power 
across the waveguide---thus hiding or "cloaking"
the presence of defect from external detection.

Now, given a candidate design
$\epsilon\sup{trial}(\mathbf{x})$, it's easy enough to see
how we can use <span class=SC>meep</span> to evaluate
its performance---just create a
<span class=SC>meep</span> geometry with $\epsilon\sup{trial}$ as a
[spatially-varying permittivity function](../Python_User_Interface.md#eps_func)
in the design region,
add [DFT cells](../Python_User_Interface.md#FluxSpectra)
to tabulate the frequency-domain Poynting flux entering and departing
the cloak region,
[timestep](../Python_User_Interface.md#RunStepFunctions) until 
the DFTs converge, then use post-processing routines like
[`#!py3 get_fluxes('hello')`](../Python_User_Interface.md#get_fluxes)
or perhaps
[`#!py3 get_eigenmode_coefficients()`](../Python_User_Interface.md#get_eigenmode_coefficients)
to get the quantities needed to evaluate the performance of the device.
Thus, for the cost of one full <span class=SC>meep</span> timestepping
run we obtain the value of our objective function at one point
in the parameter space of possible inputs. 

But *now* what do we do?! The difficulty is that the computation
just described furnishes only the *value* of the objective function
for a given input, not its *derivatives* with respect to the
design variables---and thus yields zero insight into how we should
tweak the design to improve performance.
In simple cases we might hope to proceed on the basis of physical
intuition, while
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


## Examples of optimization problems

Throughout this tutorial we will refer to a running collection of simple optimization
problems to illustrate the mechanics of optimization in <span class=SC>meep</span>.

| ![images/
![](https://example.com/qqq.img){: style="height:150px;width:150px"}



| ![blue_pc_icon](https://cloud.githubusercontent.com/assets/10035308/22727034/7d475af2-ed8b-11e6-87bf-0872a2cd006f.png)|

## Common elements of optimization geometries: Objective regions, objective functions, design regions, basis sets

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
    $$ \epsilon(\mathbf x) \equiv \sum_{d=1}^N \beta_d \mathcal{b}_d(\mathbf x),
       \qquad \mathbf x\in \mathcal{V}\sup{des},
    $$
    where $\{\mathcal{b}_n(\mathbf x)\}$ is a set of $D$ scalar-valued
    basis functions defined for $\mathbf x\in\mathcal{V}\sup{des}$.
    For adjoint optimization in <span class=SC>meep</span>, the
    basis set is chosen by the user (either from among a predefined collection of
    common basis sets, or as an arbitrary user-defined set),
  and    and the task of the optimizer becomes to determine
    numerical values for the $N$-vector of coefficients 
    $\boldsymbol{\beta}=\{\beta_n\},n=1,\cdots,N.$
    
## Mechanics of <span class=SC>meep</span> design optimization

With that by way of overview, we're in a position to sketch
the 

The first step is to communicate to the
You will do this by writing


```python
import foo.bar
import meep as mp
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

n = 3.4
w = 1
r = 1
```
