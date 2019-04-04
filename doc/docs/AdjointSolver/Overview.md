--8<-- "doc/docs/AdjointSolver/AdjointDocumentationStyleHeader.md"

---
# Adjoint sensitivity analysis for automated design optimization
---

This section of the Meep documentation
covers `meep.adjoint`, a submodule of the Meep python module
that implements an [*adjoint-based sensitivity solver*](https://en.wikipedia.org/wiki/Adjoint_state_method)
to facilitate automated design optimization via derivative-based numerical optimizers.

> :bookmark:{.tocfloat .summary} **`table of contents`**
>
> The `meep.adjoint` documentation is divided into a number of subsections:
> 
> + This **Overview** page reviews some basic facts about adjoints and optimizers,
>   outlines the steps needed to prepare a Meep
>   geometry for optimization, and sketches the mechanics of
>   the `meep.adjoint` design process.
>   (This page is designed to be a gentle introduction for the
>   adjoint neophyte; experts may want only to skim it before
>   skipping to the next section.)
> 
> + The [**Reference Manual**](ReferenceManual.md) fills in the details of 
>   the topics outlined on this page, spelling out exactly how to
>   write the python script that will drive your optimization
>   session.
> 
> + The [**Example Gallery**](ExampleGallery.md) presents a number
>   of worked examples that illustrate how `meep.adjoint` tackles
>   practical problems in various settings.
> 
> + The [**Implementation Notes**](AdjointImplementationNotes.md) page
>   offers a glimpse of what's behind the hood---the physical 
>    and mathematical basis of the adjoint method and how they
>   are implemented by `meep.adjoint.` An understanding of this 
>   content is not strictly necessary to use the solver, but may
>   help you get more out of the process.
> 
> + Although logically independent of the adjoint solver,
>   the [**Visualization**](Visualization.md) package bundled
>  with the `meep.adjoint` module offers several general-purpose
>   utilities for convenient visualization of various aspects
>  of Meep calculations, which are
>  useful in *any* meep calculation whether adjoint-related
>  or not.


## Overview: Adjoint-based optimization

A common task in electromagnetic engineering is to custom-tune the design
of some component of a system---a waveguide taper, a power splitter,
an input coupler, an antenna, etc.---to optimize the performance of the system
as defined by some problem-specific metric. For our purposes,
a "design" will consist of a specification of the spatially-varying
scalar permittivity $\epsilon(\mathbf x)$ in some subregion
of a Meep geometry, and the performance metric
will be a physical quantity computed from frequency-domain
fields---a [power flux][GetFluxes],
an [energy density][DFTEnergy],
an
[eigenmode expansion coefficient][EigenCoefficients],
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
how we can use Meep to evaluate
its performance---just create a
Meep geometry with $\epsilon\sup{trial}$ as a
[spatially-varying permittivity function][EpsFunc],
in the design region,
add [DFT cells][FluxSpectra]
to tabulate the frequency-domain Poynting flux entering and departing
the cloak region,
[timestep][RunStepFunctions] until
the DFTs converge, then use post-processing routines like
[`get_fluxes()`][GetFluxes]
or perhaps
[`get_eigenmode_coefficients()`][EigenCoefficients]
to get the quantities needed to evaluate the performance of the device.
Thus, for the cost of one full Meep timestepping
run we obtain the value of our objective function at one point
in the parameter space of possible inputs. 

But *now* what do we do?! The difficulty is that the computation
izinust described furnishes only the *value* of the objective function
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
the full Meep timestepping run mentioned
above---the "forward" run in adjoint-method parlance---we can magically
compute its derivatives with respect to *all* design variables by doing
just *one* additional timestepping run with a funny-looking choice
of sources and outputs (the "adjoint" run).
Thus, whereas gradient computation via finite-differencing is at least $D$
times more expensive than computing the objective function value,
with adjoints we get both value and gradient for roughly just *twice* the
cost of the value alone. Such a bargain! At this modest cost, derivative-based 
optimization becomes entirely feasible.

> :memo:{.center .pct80} **More general materials**
>     Although for simplicity we focus here on case of isotropic,
>     non-magnetic materials, the adjoint solver is also capable
>     of optimizing geometries involving permeable ($\mu\ne 1$) and
>     anisotropic (tensor-valued $\boldsymbol{\epsilon},\boldsymbol{\mu}$) media.


---

## Examples of optimization problems

Throughout the `meep.adjoint` documentation we will refer to a running collection of 
simple optimization problems to illustrate the mechanics of optimization,
among which are the following; click the geometry images to view 
in higher resolution.   
    

--------------------------------------------------


### The Holey Waveguide

By way of warm-up, a useful toy version of an optimization problem
is an otherwise pristine length of dielectric slab waveguide in
which a careless technician has torn a circular `hole` of variable
permittivity $\epsilon\sup{hole}$.     

    

    

> :bookmark:{.center}
>
> ![zoomify](images/HoleyWaveguideGeometry.png)


 

Incident power from an
[eigenmode source][EigenModeSource] (cyan line in figure)
travels leftward through the waveguide, but is partially 
reflected by the hole, resulting in less than 100% power
the waveguide output (as may be 
characterized in Meep
by observing power flux and/or
eigenmode expansion coefficients at the two 
flux monitors, labeled `east` and `west`).
Our objective is to tweak the value of
$\epsilon\sup{hole}$ to maximize transmission
as assessed by one of these metrics.
The simplicity of this model makes it a useful
initial warm-up and sanity check for making sure we know
what we are doing in design optimization; for example, 
[in this worked example][AdjointVsFDTest]
we use it to confirm the numerical accuracy of
adjoint-based gradients computed by `mp.adjoint`


--------------------------------------------------

### The Hole Cloak

We obtain a more challenging variant of the holey-waveguide problem
be supposing that the material in the hole region is *not* a
tunable design parameter---it is fixed at vacuum, say, or 
perfect metal---but that we *are* allowed to vary the permittivity
in an annular region surrounding the hole in such a way
as to mimic the effect of filling in the hole, i.e. of hiding
or "cloaking" the hole  as much as  possible from external 
 detection.

> :bookmark:{.center}
>
> ![zoomify](images/HoleCloakBGeometry.png)

For the hole-cloak optimization, the objective function
will most likely the same as that considered above---namely,
to maximize the Poynting flux through the flux monitor
labeled `east` (a quantity we label $S\subs{east}$)
 or perhaps to maximize the overlap coefficient
between the actual fields passing through monitor
``east`` and the fields of (say)
the $n$th forward- or backward-traveling eigenmode
of the waveguide (which we label $\{P,M\}_{n,\text{east}}$
with $P,M$ standing for "plus and minus.")
On the other hand, the design space here is more 
complicated than for the simple hole, consisting
of all possible scalar functions $\epsilon(r,\theta)$ 
defined on the annular cloak region.


--------------------------------------------------

### The cross-router

A different flavor of waveguide-optimization problem arises when we
consider the *routing* of signals from given inputs to 
given destinations. One example is the *cross-router*, involving
an intersection between $x-$directed and $y-$directed waveguides,
with center region of variable permittivity that we may
tweak to control the routing of power through it.

> :bookmark:{.center}
>
> ![zoomify](images/RouterGeometry_Iter0.png)

Whereas in the previous examples there was more or less
only one reasonable design objective one might realistically
want to optimize,
for a problem like this there are many possibilities.
For example, given fixed input power supplied by an eigenmode
source on the "western" branch (cyan line),
we might be less interested in the absolute output
power at any port and more concerned with 
achieving maximal *equality* of output 
power among the north, south, and east outputs,
whereupon we might minimize an objective function of
the form
$$f\sub{obj}  =
   \Big( S\sub{north} - S\sub{south}\Big)^2
  +\Big( S\sub{north} - S\sub{east}\Big)^2
 + \Big( S\sub{east} - S\sub{south}\Big)^2
$$
(or a similar functional form involving eigenmode 
coefficients).
Alternatively, perhaps we don't care what happens in
the southern branch, but we really want the fields 
traveling past the `north` monitor 
to have twice as much
overlap with the forward-traveling 3rd eigenmode of that
waveguide 
as the `east` fields have with their backward-traveling
2nd eigenmode:

$$ f\sub{obj} \equiv \Big( P\sub{3,north} - 2M\sub{2,east}\Big)^2$$

The point is that the definition of an optimization problem
involves not only a set of physical quantities  (power fluxes, eigenmode coefficients,
etc.) that we compute from Meep calculations,
but also a rule (the objective function $f$) for crunching those 
numbers in some specific way to define a single scalar figure of merit. 

In  `mp.adjoint` we use the collective term *objective quantities*
for the power fluxes, eigenmode coefficients, and other physical quantities
needed to compute the objective function.
Similarly, the special geometric subregions of 
Meep geometries with
which objective quantities are associated---the
cross-sectional flux planes of `DFTFlux` cells or 
field-energy boxes of `DFTField` cells----are known as *objective regions.*

The [Example Gallery][ExampleGallery.md] includes a worked example
of a full successful iterative optimization in which
`mp.adjoint` begins with the design shown above and thoroughly rejiggers
it over the course of 50 iterations to yield a device
that efficiently routs power around a 90&degree; bend
from the eigenmode source (cyan line above)
to the 'north' output port.
 
--------------------------------------------------



### The asymmetric splitter

A `splitter` seeks to divide incoming power from one source
in some specific way among two or more destinations.,
We will consider an asymmetric splitter in which power
arriving from a single incoming waveguide is to be routed
into two outgoing waveguides by varying the design of the 
central coupler region:

> :bookmark:{.center}
>
> ![zoomify](images/SplitterGeometry.png)


--------------------------------------------------


## Common elements of optimization geometries: Objective regions, objective functions, design regions, basis sets

The examples above, distinct though they all are, illustrate
the common defining features that are present in every
Meep optimization problem:

+ **Objective regions:** One or more [regions over which to tabulate frequency-domain fields (DFT cells)][DFTObj]
  for use in computing power fluxes, mode-expansion coefficients, and other frequency-domain
   quantities used in characterizing device performance.  Because these regions are used to evaluate
   objective functions, we refer to them as *objective regions.*

>:memo:{.center .pct80} **Objective regions may or may not have zero thickness**
>     In the examples above, it happens that all objective regions are one-dimensional
>     (zero-thickness) flux monitors, indicated by magenta lines; in a 3D geometry they
>     would be two-dimensional flux planes, still of zero thickness in the normal 
>     direction.  However, objective regions may also be of nonzero thickness, as for
>     instance if the objective function involves the [field energy in a box-shaped
>     subregion of a geometry.][Energy]
 
 + **Objective quantities and the objective function:** 
      A specification of which quantities (power fluxes, mode coefficients,
      energies, etc.) are to be computed for each objective region, and how
      those quantities are to be crunched mathematically to yield a single number 
      measuring device performance. We refer to the individual quantities as 
      *objective quantities*, while the overall function that inputs one more more
      objective quantities and outputs a single numerical score is the 
      *objective function.*

+ **Design region:** A specification of the region over which the material design is to be
    optimized, i.e. the region in which the permittivity is given by the
    design quantity $\epsilon\sup{des}(\mathbf x)$.
    We refer to this as the *design region* $\mathcal{V}\sup{des}$.

+ **Basis:** Because the design variable $\epsilon\sup{des}(\mathbf x)$
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
    The task of the optimizer then becomes to determine
    numerical values for the $N$-vector of coefficients 
    $\boldsymbol{\beta}=\{\beta_n\},n=1,\cdots,N.$

    For adjoint optimization in Meep, the
    basis set is chosen by the user, either from among a predefined collection of
    common basis sets, or as an arbitrary user-defined basis set specified by
    subclassing an abstract base class in `mp.adjoint.`
    
## Mechanics of Meep design optimization

With all that by way of background, here's a quick rundown of the 
process you'll follow to optimize a geometry in `meep.adjoint.`

1. You write a python script that implements a subclass of
   `OptimizationProblem` (an abstract base class in `meep.adjoint`)
   to describe your specific problem. In particular, your 
   class must override the following two pure virtual methods
   in `OptimizationProblem:`

??? summary "`init_problem`: One-time initialization"
    Inputs an `args` structure describing command-line options
    and returns a 5-tuple
    ```py
       fstr, objective_regions, extra_regions, design_region, basis
    ```
    defining your objective function, the objective regions on which 
    its inputs (the objective variables) are defined, the design region,
    and an expansion basis.


??? summary "`create_sim`: Instantiation of design-dependent geometries"
    Inputs a vector of expansion coefficients `beta_vector` and 
    returns a `meep.simulation` describing a geometry with the 
    corresponding spatially-varying permittivity.


2. You run computations on your geometry either by executing your
   script from the shell with command-line options:

```bash
  % python HoleyWaveguide.py --beta 0 2.3 --eval_gradient
```

or equivalently from a python script or console by 
calling its `run()` method:

```py
  from HoleyWaveguide import HoleyWaveguide

  HW=HoleyWaveguide(cmdline='--beta 0 2.3 --eval_gradient')
  HW.run()
```

The actual calculations that may be run in this way
range from a single non-iterative computation of the objective
function and (optionally) its gradient at a given set of design-variable
values, to [full-blown iterative design optimization][CrossRouterExample].

Here, in their entirety, are the python scripts implementing the 4 examples
described above. (These may also be found in the `python/examples/adjoint_optimization`
subdirectory of your Meep installation.)


??? example "`HoleyWaveguide.py`"
    ```py
    
    import sys
    import argparse
    import numpy as np
    import meep as mp
    
    from meep.adjoint import (OptimizationProblem, DFTCell, adjoint_options,
                              xHat, yHat, zHat, origin, FluxLine,
                              ParameterizedDielectric, FourierLegendreBasis)
    
    ##################################################
    ##################################################
    ##################################################
    class HoleyWaveguide(OptimizationProblem):
    
        ##################################################
        ##################################################
        ##################################################
        def add_args(self, parser):
    
            # add new problem-specific arguments
            parser.add_argument('--dair',        type=float, default=-1.0, help='')
            parser.add_argument('--w_wvg',       type=float, default=3.0,  help='')
            parser.add_argument('--eps_wvg',     type=float, default=6.0,  help='')
            parser.add_argument('--r_disc',      type=float, default=0.5,  help='')
            parser.add_argument('--nr_max',      type=int,   default=0,    help='')
            parser.add_argument('--kphi_max',    type=int,   default=0,    help='')
    
            # set problem-specific defaults for existing (general) arguments
            parser.set_defaults(fcen=0.5)
            parser.set_defaults(df=0.2)
            parser.set_defaults(dpml=1.0)
    
        ##################################################
        ##################################################
        ##################################################
        def init_problem(self, args):
    
            #----------------------------------------
            # size of computational cell
            #----------------------------------------
            lcen       = 1.0/args.fcen
            dpml       = 0.5*lcen if args.dpml==-1.0 else args.dpml
            dair       = 0.5*args.w_wvg if args.dair==-1.0 else args.dair
            L          = 3.0*lcen
            Lmin       = 6.0*dpml + 2.0*args.r_disc
            L          = max(L,Lmin)
            sx         = dpml+L+dpml
            sy         = dpml+dair+args.w_wvg+dair+dpml
            cell_size  = mp.Vector3(sx,sy)
    
            #----------------------------------------
            #- design region
            #----------------------------------------
            design_center = origin
            design_size   = mp.Vector3(2.0*args.r_disc, 2.0*args.r_disc)
            design_region = mp.Volume(center=design_center, size=design_size)
    
            #----------------------------------------
            #- objective regions
            #----------------------------------------
            x0_east       =  args.r_disc + dpml
            x0_west       = -args.r_disc - dpml
            y0            = 0.0
            flux_length   = 2.0*args.w_wvg
            east          = FluxLine(x0_east,y0,flux_length,mp.X,'east')
            west          = FluxLine(x0_west,y0,flux_length,mp.X,'west')
    
            objective_regions  = [east, west]
    
            #----------------------------------------
            #- optional extra regions for visualization
            #----------------------------------------
            extra_regions      = [mp.Volume(center=origin, size=cell_size)] if args.full_dfts else []
    
            #----------------------------------------
            # basis set
            #----------------------------------------
            basis = FourierLegendreBasis(radius=args.r_disc, nr_max=args.nr_max, kphi_max=args.kphi_max)
    
            #----------------------------------------
            #- source location
            #----------------------------------------
            source_center    = (x0_west - dpml)*xHat
            source_size      = flux_length*yHat
    
            #----------------------------------------
            #- objective function
            #----------------------------------------
            fstr='Abs(P1_east)**2+0.0*(P2_east+P1_west+P2_west+M1_east+M2_east+M1_west+M2_west+S_east+S_west)'
    
            #----------------------------------------
            #- internal storage for variables needed later
            #----------------------------------------
            self.args            = args
            self.dpml            = dpml
            self.cell_size       = cell_size
            self.basis           = basis
            self.design_center   = design_center
            self.source_center   = source_center
            self.source_size     = source_size
    
            return fstr, objective_regions, extra_regions, design_region, basis
    
        ##############################################################
        ##############################################################
        ##############################################################
        def create_sim(self, beta_vector, vacuum=False):
    
            args=self.args
            sx=self.cell_size.x
    
            wvg=mp.Block(center=origin, material=mp.Medium(epsilon=args.eps_wvg),
                         size=mp.Vector3(self.cell_size.x,args.w_wvg))
            disc=mp.Cylinder(center=self.design_center, radius=args.r_disc,
                             epsilon_func=ParameterizedDielectric(self.design_center,
                                                                  self.basis,
                                                                  beta_vector))
    
            geometry=[wvg] if vacuum else [wvg, disc]
    
            envelope = mp.GaussianSource(args.fcen,fwidth=args.df)
            amp=1.0
            if callable(getattr(envelope, "fourier_transform", None)):
                amp /= envelope.fourier_transform(args.fcen)
            sources=[mp.EigenModeSource(src=envelope,
                                        center=self.source_center,
                                        size=self.source_size,
                                        eig_band=self.args.source_mode,
                                        amplitude=amp
                                       )
                    ]
    
            sim=mp.Simulation(resolution=args.res, cell_size=self.cell_size,
                              boundary_layers=[mp.PML(args.dpml)], geometry=geometry,
                              sources=sources)
    
            if args.complex_fields:
                sim.force_complex_fields=True
    
            return sim
    
    ######################################################################
    # if executed as a script, we look at our own filename to figure out
    # the name of the class above, create an instance of this class called
    # opt_prob, and call its run() method.
    ######################################################################
    if __name__ == '__main__':
        opt_prob=globals()[__file__.split('/')[-1].split('.')[0]]()
        opt_prob.run()

    ```

--------------------------------------------------


??? example "`HoleCloak.py`"
    ```py
    
    import sys
    import argparse
    import numpy as np
    import meep as mp
    
    from meep.adjoint import (OptimizationProblem, DFTCell, adjoint_options,
                              xHat, yHat, zHat, origin, FluxLine,
                              ParameterizedDielectric, FourierLegendreBasis)
    
    ##################################################
    ##################################################
    ##################################################
    class HoleCloak(OptimizationProblem):
    
        ##################################################
        ##################################################
        ##################################################
        def add_args(self, parser):
    
            # add new problem-specific arguments
            parser.add_argument('--dair',        type=float, default=-1.0, help='')
            parser.add_argument('--w_wvg',       type=float, default=4.0,  help='')
            parser.add_argument('--eps_wvg',     type=float, default=6.0,  help='')
            parser.add_argument('--r_disc',      type=float, default=0.5,  help='')
            parser.add_argument('--r_cloak',     type=float, default=1.5,  help='')
            parser.add_argument('--nr_max',      type=int,   default=3,    help='')
            parser.add_argument('--kphi_max',    type=int,   default=2,    help='')
            parser.add_argument('--eps_disc',    type=float, default=1.0,  help='permittivity in hole region (0.0 for PEC)')
    
            # set problem-specific defaults for existing (general) arguments
            parser.set_defaults(fcen=0.5)
            parser.set_defaults(df=0.2)
            parser.set_defaults(dpml=1.0)
    
        ##################################################
        ##################################################
        ##################################################
        def init_problem(self, args):
    
            #----------------------------------------
            # size of computational cell
            #----------------------------------------
            lcen       = 1.0/args.fcen
            dpml       = 0.5*lcen if args.dpml==-1.0 else args.dpml
            dair       = 0.5*args.w_wvg if args.dair==-1.0 else args.dair
            L          = 3.0*lcen
            Lmin       = 6.0*dpml + 2.0*args.r_cloak
            L          = max(L,Lmin)
            sx         = dpml+L+dpml
            sy         = dpml+dair+args.w_wvg+dair+dpml
            cell_size  = mp.Vector3(sx, sy, 0.0)
    
            #----------------------------------------
            #- design region
            #----------------------------------------
            design_center = origin
            design_size   = mp.Vector3(2.0*args.r_cloak, 2.0*args.r_cloak)
            design_region = mp.Volume(center=design_center, size=design_size)
    
            #----------------------------------------
            #- objective regions
            #----------------------------------------
            fluxW_center  =  (+args.r_cloak+ dpml)*xHat
            fluxE_center  =  (-args.r_cloak- dpml)*xHat
            flux_size     =  2.0*args.w_wvg*yHat
    
            #fluxW_region  = mp.FluxRegion(center=fluxW_center, size=flux_size, direction=mp.X)
            #fluxE_region  = mp.FluxRegion(center=fluxE_center, size=flux_size, direction=mp.X)
            x0_east       =  args.r_cloak + dpml
            x0_west       = -args.r_cloak - dpml
            y0            = 0.0
            flux_length   = 2.0*args.w_wvg
            east          = FluxLine(x0_east,y0,flux_length,mp.X,'east')
            west          = FluxLine(x0_west,y0,flux_length,mp.X,'west')
    
            objective_regions  = [east, west]
    
            #----------------------------------------
            #- optional extra regions for visualization
            #----------------------------------------
            extra_regions      = [mp.Volume(center=origin, size=cell_size)] if args.full_dfts else []
    
            #----------------------------------------
            # basis set
            #----------------------------------------
            basis = FourierLegendreBasis(outer_radius=args.r_cloak, inner_radius=args.r_disc,
                                         nr_max=args.nr_max, kphi_max=args.kphi_max)
    
            #----------------------------------------
            #- source location
            #----------------------------------------
            source_center    = (x0_west-dpml)*xHat
            source_size      = flux_length*yHat
    
            #----------------------------------------
            #- objective function
            #----------------------------------------
            fstr='Abs(P1_east)**2+0.0*(P1_west + M1_east + M1_west + S_west + S_east)'
    
            #----------------------------------------
            #- internal storage for variables needed later
            #----------------------------------------
            self.args            = args
            self.dpml            = dpml
            self.cell_size       = cell_size
            self.basis           = basis
            self.design_center   = design_center
            self.source_center   = source_center
            self.source_size     = source_size
    
            return fstr, objective_regions, extra_regions, design_region, basis
    
        ##############################################################
        ##############################################################
        ##############################################################
        def create_sim(self, beta_vector, vacuum=False):
    
            args=self.args
            sx=self.cell_size.x
    
            wvg=mp.Block(center=origin, material=mp.Medium(epsilon=args.eps_wvg),
                         size=mp.Vector3(self.cell_size.x,args.w_wvg))
            cloak=mp.Cylinder(center=self.design_center, radius=args.r_cloak,
                              epsilon_func=ParameterizedDielectric(self.design_center,
                                                                   self.basis,
                                                                   beta_vector))
            disc=mp.Cylinder(center=self.design_center, radius=args.r_disc,
                             material=(mp.metal if args.eps_disc==0 else
                                       mp.Medium(epsilon=args.eps_disc)))
    
            geometry=[wvg] if vacuum else [wvg, cloak, disc]
    
            envelope = mp.GaussianSource(args.fcen,fwidth=args.df)
            amp=1.0
            if callable(getattr(envelope, "fourier_transform", None)):
                amp /= envelope.fourier_transform(args.fcen)
            sources=[mp.EigenModeSource(src=envelope,
                                        center=self.source_center,
                                        size=self.source_size,
                                        eig_band=self.args.source_mode,
                                        amplitude=amp
                                       )
                    ]
    
            sim=mp.Simulation(resolution=args.res, cell_size=self.cell_size,
                              boundary_layers=[mp.PML(args.dpml)], geometry=geometry,
                              sources=sources)
        
            if args.complex_fields:
                sim.force_complex_fields=True
    
            return sim
    
    ######################################################################
    # if executed as a script, we look at our own filename to figure out
    # the name of the class above, create an instance of this class called
    # opt_prob, and call its run() method.
    ######################################################################
    if __name__ == '__main__':
        opt_prob=globals()[__file__.split('/')[-1].split('.')[0]]()
        opt_prob.run()
    ```
    
??? example "`CrossRouter.py`"

    ```py
    
    import numpy as np
    import meep as mp
    
    from meep.adjoint import (OptimizationProblem, FluxLine,
                              xHat, yHat, zHat, origin,
                              ParameterizedDielectric, FiniteElementBasis)
    
    ##################################################
    ##################################################
    ##################################################
    class CrossRouter(OptimizationProblem):
    
        ##################################################
        ##################################################
        ##################################################
        def add_args(self, parser):
    
            # add new problem-specific arguments
            parser.add_argument('--wh',       type=float, default=1.5,  help='width of horizontal waveguide')
            parser.add_argument('--wv',       type=float, default=1.5,  help='width of vertical waveguide')
            parser.add_argument('--l_stub',   type=float, default=3.0,  help='waveguide input/output stub length')
            parser.add_argument('--eps',      type=float, default=6.0,  help='waveguide permittivity')
            parser.add_argument('--r_design', type=float, default=0.0,  help='design region radius')
            parser.add_argument('--l_design', type=float, default=4.0,  help='design region side length')
            parser.add_argument('--nfe',      type=int,   default=2,    help='number of finite elements per unit length')
            parser.add_argument('--n_weight', type=float, default=1.00, help='')
            parser.add_argument('--s_weight', type=float, default=0.00, help='')
            parser.add_argument('--e_weight', type=float, default=0.00, help='')
    
            # set problem-specific defaults for existing (general) arguments
            parser.set_defaults(fcen=0.5)
            parser.set_defaults(df=0.2)
            parser.set_defaults(dpml=1.0)
            parser.set_defaults(epsilon_design=6.0)
    
        ##################################################
        ##################################################
        ##################################################
        def init_problem(self, args):
    
            #----------------------------------------
            # size of computational cell
            #----------------------------------------
            lcen          = 1.0/args.fcen
            dpml          = 0.5*lcen if args.dpml == -1.0 else args.dpml
            design_length = 2.0*args.r_design if args.r_design > 0.0 else args.l_design
            sx = sy       = dpml + args.l_stub + design_length + args.l_stub + dpml
            cell_size     = mp.Vector3(sx, sy, 0.0)
    
            #----------------------------------------
            #- design region bounding box
            #----------------------------------------
            design_center = origin
            design_size   = mp.Vector3(design_length, design_length)
            design_region = mp.Volume(center=design_center, size=design_size)
    
            #----------------------------------------
            #- objective and source regions
            #----------------------------------------
            gap            =  args.l_stub/6.0                    # gap between source region and flux monitor
            d_flux         =  0.5*(design_length + args.l_stub)  # distance from origin to NSEW flux monitors
            d_source       =  d_flux + gap                       # distance from origin to source
            d_flx2         =  d_flux + 2.0*gap
            l_flux_NS      =  2.0*args.wv
            l_flux_EW      =  2.0*args.wh
            north          =  FluxLine(0.0, +d_flux, l_flux_NS, mp.Y, 'north')
            south          =  FluxLine(0.0, -d_flux, l_flux_NS, mp.Y, 'south')
            east           =  FluxLine(+d_flux, 0.0, l_flux_EW, mp.X, 'east')
            west1          =  FluxLine(-d_flux, 0.0, l_flux_EW, mp.X, 'west1')
            west2          =  FluxLine(-d_flx2, 0.0, l_flux_EW, mp.X, 'west2')
    
            objective_regions  = [north, south, east, west1, west2]
    
            source_center  =  mp.Vector3(-d_source, 0.0)
            source_size    =  mp.Vector3(0.0,l_flux_EW)
    
            #----------------------------------------
            #- optional extra regions for visualization
            #----------------------------------------
            extra_regions  = [mp.Volume(center=origin, size=cell_size)] if args.full_dfts else []
    
            #----------------------------------------
            # basis set
            #----------------------------------------
            basis = FiniteElementBasis(lx=args.l_design, ly=args.l_design, density=args.nfe)
    
            #----------------------------------------
            #- objective function
            #----------------------------------------
            fstr=(   '   {:s}*Abs(P1_north)**2'.format('0.0' if args.n_weight==0.0 else '{}'.format(args.n_weight))
                   + ' + {:s}*Abs(M1_south)**2'.format('0.0' if args.s_weight==0.0 else '{}'.format(args.s_weight))
                   + ' + {:s}*Abs(P1_east)**2'.format('0.0'  if args.e_weight==0.0 else '{}'.format(args.e_weight))
                   + ' + 0.0*(P1_north + M1_south + P1_east + P1_west1 + P1_west2)'
                   + ' + 0.0*(M1_north + M1_south + M1_east + M1_west1 + M1_west2)'
                   + ' + 0.0*(S_north + S_south + S_east + S_west1 + S_west2)'
                 )
    
            #----------------------------------------
            #- internal storage for variables needed later
            #----------------------------------------
            self.args            = args
            self.dpml            = dpml
            self.cell_size       = cell_size
            self.basis           = basis
            self.design_center   = design_center
            self.design_size     = design_size
            self.source_center   = source_center
            self.source_size     = source_size
    
            if args.eps_design is None:
                args.eps_design = args.eps
    
            return fstr, objective_regions, extra_regions, design_region, basis
    
        ##############################################################
        ##############################################################
        ##############################################################
        def create_sim(self, beta_vector, vacuum=False):
    
            args=self.args
    
            hwvg=mp.Block(center=origin, material=mp.Medium(epsilon=args.eps),
                          size=mp.Vector3(self.cell_size.x,args.wh))
            vwvg=mp.Block(center=origin, material=mp.Medium(epsilon=args.eps),
                          size=mp.Vector3(args.wv,self.cell_size.y))
    
            if args.r_design>0.0:
                router=mp.Cylinder(center=self.design_center, radius=args.r_design,
                                   epsilon_func=ParameterizedDielectric(self.design_center,
                                                                        self.basis,
                                                                        beta_vector))
            else:
                router=mp.Block(center=self.design_center, size=self.design_size,
                                epsilon_func=ParameterizedDielectric(self.design_center,
                                                                     self.basis,
                                                                     beta_vector))
            geometry=[hwvg, vwvg, router]
    
            envelope = mp.GaussianSource(args.fcen,fwidth=args.df)
            amp=1.0
            if callable(getattr(envelope, "fourier_transform", None)):
                amp /= envelope.fourier_transform(args.fcen)
            sources=[mp.EigenModeSource(src=envelope,
                                        center=self.source_center,
                                        size=self.source_size,
                                        eig_band=args.source_mode,
                                        amplitude=amp
                                       )
                    ]
    
            sim=mp.Simulation(resolution=args.res, cell_size=self.cell_size,
                              boundary_layers=[mp.PML(self.dpml)], geometry=geometry,
                              sources=sources)
    
            if args.complex_fields:
                sim.force_complex_fields=True
    
            return sim
    
    ######################################################################
    # if executed as a script, we look at our own filename to figure out
    # the name of the class above, create an instance of this class called
    # op, and call its run() method.
    ######################################################################
    if __name__ == '__main__':
        op=globals()[__file__.split('/')[-1].split('.')[0]]()
        op.run()
    
    ```
    
??? example "`AsymmetricSplitter.py`"
   ```py
   import sys
   import argparse
   import numpy as np
   import meep as mp
   
   from meep.adjoint import (OptimizationProblem, DFTCell, adjoint_options,
                             xHat, yHat, zHat, origin, FluxLine,
                             ParameterizedDielectric, FiniteElementBasis)
   
   ##################################################
   ##################################################
   ##################################################
   class AsymmetricSplitter(OptimizationProblem):
   
       ##################################################
       ##################################################
       ##################################################
       def add_args(self, parser):
   
           # add new problem-specific arguments
           parser.add_argument('--dair',        type=float, default=-1.0, help='')
           parser.add_argument('--w_in',        type=float, default=1.0,  help='width of input waveguide')
           parser.add_argument('--w_out1',      type=float, default=0.5,  help='width of output waveguide 1')
           parser.add_argument('--w_out2',      type=float, default=0.5,  help='width of output waveguide 2')
           parser.add_argument('--l_stub',      type=float, default=3.0,  help='length of waveguide input/output stub')
           parser.add_argument('--l_design',    type=float, default=2.0,  help='length of design region')
           parser.add_argument('--h_design',    type=float, default=6.0,  help='height of design region')
           parser.add_argument('--eps_in',      type=float, default=6.0,  help='input waveguide permittivity')
           parser.add_argument('--eps_out1',    type=float, default=2.0,  help='output waveguide 1 permittivity')
           parser.add_argument('--eps_out2',    type=float, default=12.0, help='output waveguide 2 permittivity')
           parser.add_argument('--nfe',         type=int,   default=2,    help='number of finite elements per unit length')
   
           # set problem-specific defaults for existing (general) arguments
           parser.set_defaults(fcen=0.5)
           parser.set_defaults(df=0.2)
           parser.set_defaults(dpml=1.0)
   
       ##################################################
       ##################################################
       ##################################################
       def init_problem(self, args):
   
           #----------------------------------------
           # size of computational cell
           #----------------------------------------
           lcen       = 1.0/args.fcen
           dpml       = 0.5*lcen if args.dpml==-1.0 else args.dpml
           dair       = 0.5*args.w_in if args.dair==-1.0 else args.dair
           sx         = dpml + args.l_stub + args.l_design + args.l_stub + dpml
           sy         = dpml + dair + args.h_design + dair + dpml
           cell_size  = mp.Vector3(sx, sy, 0.0)
   
           #----------------------------------------
           #- design region
           #----------------------------------------
           design_center = origin
           design_size   = mp.Vector3(args.l_design, args.h_design, 0.0)
           design_region = mp.Volume(center=design_center, size=design_size)
   
           #----------------------------------------
           #- objective regions
           #----------------------------------------
           x_in          =  -0.5*(args.l_design + args.l_stub)
           x_out         =  +0.5*(args.l_design + args.l_stub)
           y_out1        =  +0.25*args.h_design
           y_out2        =  -0.25*args.h_design
   
           flux_in       =  FluxLine(x_in,     0.0, 2.0*args.w_in,   mp.X, 'in')
           flux_out1     =  FluxLine(x_out, y_out1, 2.0*args.w_out1, mp.X, 'out1')
           flux_out2     =  FluxLine(x_out, y_out2, 2.0*args.w_out2, mp.X, 'out2')
   
           objective_regions  = [flux_in, flux_out1, flux_out2]
   
           #----------------------------------------
           #- optional extra regions for visualization if the --full-dfts options is present.
           #----------------------------------------
           extra_regions      = [mp.Volume(center=origin, size=cell_size)] if args.full_dfts else []
   
           #----------------------------------------
           # forward source region
           #----------------------------------------
           source_center    =  (x_in - 0.25*args.l_stub)*xHat
           source_size      =  2.0*args.w_in*yHat
   
           #----------------------------------------
           # basis set
           #----------------------------------------
           basis = FiniteElementBasis(args.l_design, args.h_design, args.nfe)
   
           #----------------------------------------
           #- objective function
           #----------------------------------------
           fstr = (   'Abs(P1_out1)**2'
                    + '+0.0*(P1_out1 + M1_out1)'
                    + '+0.0*(P1_out2 + M1_out2)'
                    + '+0.0*(P1_in   + M1_in + S_out1 + S_out2 + S_in)'
                  )
   
           #----------------------------------------
           #- internal storage for variables needed later
           #----------------------------------------
           self.args            = args
           self.dpml            = dpml
           self.cell_size       = cell_size
           self.basis           = basis
           self.design_center   = origin
           self.design_size     = design_size
           self.source_center   = source_center
           self.source_size     = source_size
   
           return fstr, objective_regions, extra_regions, design_region, basis
   
       ##############################################################
       ##############################################################
       ##############################################################
       def create_sim(self, beta_vector, vacuum=False):
   
           args=self.args
           sx=self.cell_size.x
   
           x_in   = -0.5*(args.l_design + args.l_stub)
           x_out  = +0.5*(args.l_design + args.l_stub)
           y_out1 = +0.25*args.h_design
           y_out2 = -0.25*args.h_design
   
           wvg_in = mp.Block( center=mp.Vector3(x_in,0.0),
                              size=mp.Vector3(args.l_stub,args.w_in),
                              material=mp.Medium(epsilon=args.eps_in))
           wvg_out1 = mp.Block( center=mp.Vector3(x_out,y_out1),
                                size=mp.Vector3(args.l_stub,args.w_out1),
                                material=mp.Medium(epsilon=args.eps_out1))
           wvg_out2 = mp.Block( center=mp.Vector3(x_out,y_out2),
                                size=mp.Vector3(args.l_stub,args.w_out2),
                                material=mp.Medium(epsilon=args.eps_out2))
           design   = mp.Block( center=origin,
                                size=mp.Vector3(args.l_design,args.h_design),
                                epsilon_func=ParameterizedDielectric(self.design_center,
                                                                     self.basis,
                                                                     beta_vector)
                              )
   
           geometry=[wvg_in, wvg_out1, wvg_out2, design]
   
           envelope = mp.GaussianSource(args.fcen,fwidth=args.df)
           amp=1.0
           if callable(getattr(envelope, "fourier_transform", None)):
               amp /= envelope.fourier_transform(args.fcen)
           sources=[mp.EigenModeSource(src=envelope,
                                       center=self.source_center,
                                       size=self.source_size,
                                       eig_band=self.args.source_mode,
                                       amplitude=amp
                                      )
                   ]
   
           sim=mp.Simulation(resolution=args.res, cell_size=self.cell_size,
                             boundary_layers=[mp.PML(args.dpml)], geometry=geometry,
                             sources=sources)
   
           if args.complex_fields:
               sim.force_complex_fields=True
   
           return sim
   
   ######################################################################
   # if executed as a script, we look at our own filename to figure out
   # the name of the class above, create an instance of this class called
   # opt_prob, and call its run() method.
   ######################################################################
   if __name__ == '__main__':
       opt_prob=globals()[__file__.split('/')[-1].split('.')[0]]()
       opt_prob.run()
   
   ```
   
   
   
--8<-- "doc/docs/AdjointSolver/AdjointLinks.md"
