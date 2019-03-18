--8<-- "AdjointSolver/AdjointDocumentationStyleHeader.md"

---
# `meep.adjoint.Visualization`: Visualization of
---


The `meep.adjoint` python module that implements the
<span class=SC>meep</span> adjoint solver includes
an extensive set of tools for graphical visualization
of <span class=SC>meep</span> sessions via
[<span class=SC>matplotlib</span>](http://matplotlib.org).
Although these tools are packaged with and used by
`meep.adjoint`, they are independent of the adjoint solver
and can be used for convenient visualization of *any* python-driven
<span class=SC>meep</span> session---indeed, it is a specific
goal of the module that the visualization routines can
be invoked, with no arguments,
on any <span class=SC>meep</span> geometry
and will do something useful and standardized.

## Motivation and philosophy

### The usefulness of visual sanity checks in a text-based simulator

The fact that every aspect of a <span class=SC>meep</span>
calculation---from the
[material geometry](/Python_User_Interface.md#geometricobject),
to the [absorbing boundaries](/Python_User_Interface.md#pml),
to the placement and orientation
of the [sources](/Python_User_Interface.md#source),
to the positions and sizes of the
[flux monitors](/Python_User_Interface.md#flux-spectra)---is
defined by lines of text in script files offers a welcome
change of pace from the immersive GUI experience of commercial
CAD tools, as well as a natural platform for complex,
0ustomized calculations that would be unwieldy or impossible
to launch from a point-and-click pull-down interface.
Nonetheless, as convenient as it is to specify geometries
non-graphically, it's *also* quite useful---arguably even
*essential*---to review a graphical representation
of the geometry you input as interpreted by
<span class=SC>meep</span>, both to confirm that your input
was processed as you expected and to catch the sorts of
inadvertent errors---a flux monitor in the wrong place,
a source pointing the wrong direction---that are
all too easy to overlook in your python script,
but instantly obvious upon visual inspection.

In addition to the usefulness of visual sanity checks
on your input geometry, it can be helpful to look at
graphical visualizations of the electric and magnetic 
fields computed by <span class=SC>meep</span>---both
the time-domain and the frequency-domain fields---and
in particular to review how the field distributions 
evolve in space in the presence of your material geometry.
This is true even though, in many cases, the ultimate
desired output of your calculation will be non-graphical,
sucha a number---a power flux, an eigenmode coefficient, 
etc---whose computation, again, will be described by
lines of text you write in a script file. If the calculation
works correctly as you wrote it and you get 
a reasonable numerical output on your first try,
well, bully for you! You have no need of graphical
visualization aids. For the rest of us---who can't
figure out why our flux is 10 orders of magnitude 
too small until we look at the time-domain power flow and
realize we set things up backward, so that all the 
excitation power flows promptly *out of*, not into,
the computational cell---visualization tools can be 
a lifesaver.

### Isn't this a solved problem?

So we're agreed that visualization---of both inputs (geometries)
and outputs (fields)---is a good thing. But isn't this a 
solved problem? After all, <span class=SC>meep</span> offers
plenty of routines for retrieving as much raw data on
geometries and fields as any user could want---and, once one
has the raw data, it's just a matter of choosing from among the 
infinitude of available tools for plotting and visualization.
Indeed, already within the <span class=SC>meep</span> documentation
itself one can find many different types of visualization generated 
by many different types of tool. What more is left to say?

### Wanted: Quick, easy, *standard* visualization paradigms that *just work* on any geometry

Despite the arguably non-issue status of the situation, I would argue
that the current situation is suboptimal in at least two ways.

**(a)** The absence of a single, canonical solution means that, as a <span class=SC>meep</span> user,
every time you feel the urge to visualize something you have to spend some time and effort
figuring out how you are going to do it---and then going through the hassle of setting that up.
In my experience, this was especially true when it came to making movies of the evolution
time-domain fields---although the problem had been solved with documentation a couple of times
in the manual, still somehow it felt like a major chore to set this up every time
I wanted it.

**(b)** On a different note, the absence of generally-accepted visualization protocols
means that everybody's visualizations look different. This is not necessarily tragic, and
we wouldn't want to enforce sterile conformity, but it might be nice if there were *some*
notion of "canonical <span class=SC>meep</span> visualization format" to serve as a common
language.

### Desiderata motivating `meep.adjoint.Visualization`

With all of that by way of backstory, here were the guidelines that motivated
the design of these routines.

1. There are three basic scenarios calling for a canonical visualization format:

    + *Simulation geometry before timestepping.* A static image of the simulation geometry at the begininng of the calculation, 
      before timestepping has begun. The important targets for visualization are
      the following items and their relative positioning vis-a-vis one another:
      **(a)** the material geometry
      **(b)** the absorbers
      **(c)** the location and nature of sources,
      **(d)** the position and size of any DFT cells in which frequency-domain fields 
              are to be computed.

    + *Time-domain fields during timestepping.* An animation or movie showing the time-domain fields evolving throughout 
      the timestepping. Ideally, it would be nice to see the fields together
      with the geometry in some sort of superposed visualization complex.

    + *Frequency-domain fields after timestepping.* A static image showing the amplitudes of frequency-domain fields
      in regions of DFT cells, with the spatial distributiono of Poynting flux plotted for DFT flux cells.

   The module should provide functionality to generate each of these three types of visualization.

2. The visualization routines should

   **Zero-effort calling model for lazy forgetful slackers **

   For users who can't be bothered to memorize or look up library function
   prototypes, it should be possible to call the visualization routines
   on just the `simulation` entity alone *with no configuration arguments whatsoever,*
   with the module making reasonable choices to produce helpful graphics 
   in a standardized format.

2.  **Minimal-effort calling model for the slightly more motivated**
    For users who fall basically into the previous bucket but are
    motivated to change one or another option that otherwise
    obstructs the visualization (like, label font size too large),
    there should be an easy way to configure global settings once,
    after which the user can proceed lazily to call all visualization
    routines with no customization arguments.

3. **Full-customization calling model for diligent nitpickers:**
    Finally, for users who *do* want to invest time and effort to 
    configure plotting options to their liking, the visualization
    module should provide copious customization options to allow
    as much as possible of the visualization to be customized.
    
## Before timestepping: Inspecting and sanity-checking simulation geometries

The first visualization task we will consider is that of double-checking the various geometric inputs
we supply to <span class=SC>meep</span> when initiating a simulation: the computational cell
and material geometry, the absorbing boundaries, the exciting sources, and any DFT cells for which
we requested tabulation of frequency-domain fields. All of this information is presented in graphical
form

## During timestepping: Real-time visualization of time-domain field evolution

To support this, `meep.adjoint.Visualization` defines a new
[step function](/Python_User_Interface.md#run-and-step-functions)
that you can pass to the
[run function](/Python_User_Interface.md#run-functions) you use
to drive timestepping.
 
## After timestepping: Visualizing frequency-domain Poynting fluxes and energy densities

## 
