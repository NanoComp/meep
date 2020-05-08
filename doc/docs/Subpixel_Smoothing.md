---
# Subpixel Smoothing
---

Meep uses a [second-order accurate finite-difference scheme](https://en.wikipedia.org/wiki/Finite_difference_method#Accuracy_and_order) for discretizing [Maxwell's equations](Introduction.md#maxwells-equations). This means that the results from Meep converge to the "exact" result from the non-discretized (i.e., continuous) system quadratically with the resolution Δx. However, this second-order error O(Δx<sup>2</sup>) is generally spoiled to first-order error O(Δx) if the discretization involves a *discontinuous* material boundary (similar to [Gibbs phenomenon](https://en.wikipedia.org/wiki/Gibbs_phenomenon) in signal processing). Moreover, directly discretizing a discontinuity in ε or μ leads to "stairstepped" interfaces that can only be varied in discrete jumps of one pixel. Meep solves both of these problems by smoothing ε and μ: before discretizing, discontinuities are smoothed into continuous transitions over a distance of one pixel Δx, using a second-order accurate averaging procedure summarized [below](#smoothed-permittivity-tensor-via-perturbation-theory). This subpixel smoothing enables the discretized solution to converge as quickly as possible to the exact solution as the `resolution` increases.

<center>
![](images/subpixel_smoothing.png)
</center>

The subpixel smoothing has four limitations: (1) it only applies to frequency-independent, lossless dielectrics (i.e., silicon at λ=1.55 μm); dispersive materials are [not supported](FAQ.md#can-subpixel-averaging-be-applied-to-dispersive-materials), (2) it can be efficiently applied to [`GeometricObject`](Python_User_Interface.md#geometricobject)s (i.e. `Block`, `Prism`, `Sphere`, etc.) but [*not* to a user-defined `material_function`](FAQ.md#can-subpixel-averaging-be-applied-to-a-user-defined-material-function) which is [disabled by default](#enabling-averaging-for-material-function), (3) objects with sharp corners or edges are associated with field singularities which introduce an unavoidable error intermediate between first- and second-order, and (4) the fields directly *on* the interface are still at best first-order accurate. The improved accuracy from smoothing is therefore obtained for fields evaluated off of the interface as in the [scattered Poynting flux](Python_Tutorials/Basics.md#transmittance-spectrum-of-a-waveguide-bend) integrated over a surface away from the interface, for nonlocal properties such as [resonant frequencies](Python_Tutorials/Resonant_Modes_and_Transmission_in_a_Waveguide_Cavity.md#resonant-modes), and for overall integrals of fields and energies to which the interface contributes only O(Δx) of the integration domain.

[TOC]

Smoothed Permittivity Tensor via Perturbation Theory
----------------------------------------------------

Any scheme for smoothing the interface perturbs the problem you are solving, as shown in the figure above, and a second-order accurate smoothing scheme must mean that the perturbation's effect is zero to first order in the smoothing diameter (the resolution). This turns out to require that the smoothing scheme be anisotropic. Even if the initial interface is between isotropic materials, one obtains an effective tensor $\tilde{ε}$ (or $\tilde{μ}$) which uses the mean ε for fields parallel to the interface and the harmonic mean (inverse of mean of ε<sup>-1</sup>) for fields perpendicular to the interface:

<center>

$$ \tilde{ε}^{-1} = \textbf{P}\langleε^{-1}\rangle + \big(1-\textbf{P}\big)\langleε\rangle^{-1} $$

</center>

where $\textbf{P}$ is the projection matrix $P_{ij}=n_{i}n_{j}$ onto the normal $\vec{n}$. The $\langle\cdots\rangle$ denotes an average over the voxel $sΔx\times sΔy\times sΔz$ surrounding the grid point in question where $s$ is a smoothing diameter in grid units equal to 1/`resolution`. If the initial materials are anisotropic (via `epsilon_diag` and `epsilon_offdiag`), a more complicated formula is used. They key point is that, even if the structure consists entirely of isotropic materials, the discretized structure will use anisotropic materials. For interface pixels, Meep computes the effective permittivity tensor automatically at the start of the simulation prior to time stepping via analytic expressions for the filling fraction and local normal vector. For details involving derivation of the effective permittivity tensor and its implementation in Meep/FDTD, see [Optics Letters, Vol. 36, pp. 2972-4, 2006](https://www.osapublishing.org/ol/abstract.cfm?uri=ol-31-20-2972) and [Optics Letters, Vol. 35, pp. 2778-80, 2009](https://www.osapublishing.org/abstract.cfm?uri=ol-34-18-2778).

Continuously Varying Shapes and Results
---------------------------------------

A key feature of Meep's subpixel smoothing, particularly relevant for shape optimization (i.e., [Applied Physics Letters, Vol. 104, 091121, 2014](https://aip.scitation.org/doi/abs/10.1063/1.4867892) ([pdf](http://ab-initio.mit.edu/~oskooi/papers/Oskooi14_tandem.pdf))), is that continuously varying the `geometry` yields continuously varying results. This is demonstrated for a ring resonator: as the radius increases, the frequency of a resonant $H_z$-polarized mode decreases. Note: unlike the example in [Tutorial/Basics/Modes of a Ring Resonator](Python_Tutorials/Basics.md#modes-of-a-ring-resonator) involving $E_z$-polarized modes where the electric fields are always continuous (i.e., parallel to the interface), this example involves discontinuous fields. Also, the ring geometry contains no sharp corners/edges which tend to produce field singularities that degrade the error. The simulation script is shown below. The inner ring radius is varied from 1.8 to 2.0 μm in gradations of 0.005 μm. The ring width is constant (1 μm). The resolution is 10 pixels/μm. The gradations are therefore well below subpixel dimensions.

```py
import meep as mp
import numpy as np

resolution = 10         # pixels/μm

n = 3.4                 # index of waveguide
w = 1                   # width of waveguide
pad = 4                 # padding between waveguide and edge of PML
dpml = 2                # thickness of PML

for rad in np.arange(1.800,2.001,0.005):
    sxy = 2*(rad+w+pad+dpml)  # cell size

    # pulse center frequency (from third-order polynomial fit)
    fcen = -0.018765*rad**3 + 0.137685*rad**2 -0.393918*rad + 0.636202
    # pulse frequency width
    df = 0.02*fcen

    src = [mp.Source(mp.GaussianSource(fcen, fwidth=df),
                     component=mp.Hz,
                     center=mp.Vector3(rad+0.1*w)),
           mp.Source(mp.GaussianSource(fcen, fwidth=df),
                     component=mp.Hz,
                     center=mp.Vector3(-(rad+0.1*w)),
                     amplitude=-1)]

    symmetries = [mp.Mirror(mp.X,phase=+1),
                  mp.Mirror(mp.Y,phase=-1)]

    geometry = [mp.Cylinder(material=mp.Medium(index=n),
                            radius=rad+w,
                            height=mp.inf,
                            center=mp.Vector3()),
                mp.Cylinder(material=mp.vacuum,
                            radius=rad,
                            height=mp.inf,
                            center=mp.Vector3())]
    
    sim = mp.Simulation(cell_size=mp.Vector3(sxy,sxy),
                        geometry=geometry,
                        eps_averaging=True,
                        sources=src,
                        resolution=resolution,
                        symmetries=symmetries,
                        boundary_layers=[mp.PML(dpml)])

    sim.run(mp.after_sources(mp.Harminv(mp.Hz, mp.Vector3(rad+0.1), fcen, df)),
            until_after_sources=300)

    sim.reset_meep()
```

A plot of the resonant frequency versus the ring radius is shown below for subpixel smoothing (red) and no smoothing (blue). Included for reference is the "exact" (black) computed using *no smoothing* at a resolution of 60 pixels/μm. The no-smoothing result shows "staircasing" effects which are artifacts of the discretization. The subpixel-smoothing result varies continuously with the ring radius similar to the high-resolution result which is at a resolution six times larger. The inset shows the scalar $H_z$ field profile of the resonant mode for a structure with inner radius of 1.9 μm.

This particular resonant mode has a [quality (Q) factor](https://en.wikipedia.org/wiki/Q_factor) of ~10<sup>7</sup> at a frequency of 0.25 and radius of 2.0 μm. This means that roughly 4x10<sup>7</sup> optical periods are required to accurately resolve the field decay due to the Fourier uncertainty relation. Instead, [`Harminv`](Python_User_Interface.md#harminv) can resolve the $Q$ using just ~1000 periods. This is nearly a four orders of magnitude reduction in the run time.

<center>
![](images/ring_vary_radius.png)
</center>

To compare the convergence rate of the discretization error, the following plot shows the error in the resonant frequency (relative to the "exact" result at a resolution of 300 pixels/μm) as a function of the grid resolution for a ring geometry with a fixed radius of 2.0 μm. The no-smoothing results have a linear error due to the stairstepped interface discontinuities. The subpixel-smoothing results have roughly second-order convergence.

<center>
![](images/ring_subpixel_smoothing_rate.png)
</center>

Enabling Averaging for Material Function
----------------------------------------

By default, subpixel smoothing is automatically applied to any `GeometricObject` in the cell as `eps_averaging=True` in the [`Simulation`](Python_User_Interface.md#the-simulation-class) constructor. For a `material_function` however, subpixel smoothing [tends to be slow](FAQ.md#why-does-subpixel-averaging-take-so-long) due to an adaptive numerical integration method that involves callbacks from the low-level C++ routines and the Python-defined material functions. Because of this poor performance, subpixel smoothing is disabled by default for material functions (even though subpixel smoothing is still applied to other `GeometricObjects` which do not contain a `material_function`).

Subpixel smoothing can be enabled for a `material_function` by setting its `do_averaging` property to `True` as demonstrated in the following example.

```py

def ring_resonator(p):
    rr = (p.x**2+p.y**2)**0.5
    if (rr > rad) and (rr < rad+w):
        return mp.Medium(index=n)
    return mp.air

ring_resonator.do_averaging = True

geometry = [mp.Block(center=mp.Vector3(),
                     size=mp.Vector3(sxy,sxy),
                     material=ring_resonator)]

sim = mp.Simulation(cell_size=mp.Vector3(sxy,sxy),
                    geometry=geometry,
                    subpixel_tol=1e-4,
                    subpixel_maxeval=1000,
                    sources=src,
                    resolution=resolution,
                    symmetries=symmetries,
                    boundary_layers=[mp.PML(dpml)])
```

The adaptive numerical integration used for subpixel smoothing of material functions tends to be significantly slower than the analytic approach. To speed this up at the expense of reduced accuracy, the values for its two convergence parameters `subpixel_tol` (tolerance) and `subpixel_maxeval` (maximum number of function evaluations) can be increased/lowered.

What Happens When Subpixel Smoothing is Disabled?
-------------------------------------------------

When subpixel smoothing is disabled by either (1) setting `eps_averaging=False` in the [`Simulation`](Python_User_Interface.md#the-simulation-class) constructor or (2) using a [material_function](Subpixel_Smoothing.md#enabling-averaging-for-material-function) (as is typical in the [adjoint solver](Python_Tutorials/AdjointSolver.md)), each electric field component ($E_x$, $E_y$, $E_z$) in a given voxel is individually assigned a scalar permittivity (for isotropic materials) based on whatever the value of the permittivity is at that position in the [Yee grid](Yee_Lattice.md). This results in [staircasing artifacts](Subpixel_Smoothing.md) due to the discontinuous material interfaces as well as the staggered nature of the Yee grid points. Any change in the resolution which shifts the location of the Yee grid points relative to the material interfaces will result in unpredictable changes to any computed quantities. The coordinates the Yee grid points can be obtained using a [field function](Field_Functions.md#coordinates-of-the-yee-grid) which can be useful for debugging.

Subpixel Smoothing vs. Bilinear Interpolation
---------------------------------------------

In certain cases, using subpixel smoothing may be impractical given the poor runtime performance of the adaptive numerical integration method as discussed previously. A workaround to ensure that Meep responds continuously to changes in the simulation parameters is to *interpolate* any discontinuous structure onto the Yee grid. Otherwise, tiny changes in Meep's Yee grid due to e.g. small changes in the resolution could cause discontinuous jumps in ε.

As a demonstration of this effect, consider a ring resonator (inner radius: 2.0 μm, width: 1.0 μm) in which the "same" ring geometry can be represented using five different methods:

1. Two overlapping `Cylinder` objects (subpixel smoothing).
2. Two overlapping `Prism` objects, each with 40 vertices (subpixel smoothing). Could also be imported as a [GDSII file](Python_Tutorials/GDSII_Import/).
3. `material_function` (no smoothing).
4. Grid via `epsilon_input_file` (no smoothing; bilinear interpolation onto Yee grid).
5. Grid via `epsilon_input_file` (no smoothing; no interpolation).

Of these five methods, (3) and (5) produce discontinuous structures.

The HDF5 file used as the `epsilon_input_file` in (4) and (5) is generated by the function `output_epsilon` when using the `material_function` from (3) at a resolution of 80.

The `Prism` ring geometry in (2) is generated using:

```py
N = 40
phis = np.linspace(0,2*np.pi,N+1)
vertices_outer = []
vertices_inner = []
for phi in phis[:-1]:
    vertices_outer.append((rad+w)*mp.Vector3(np.cos(phi),np.sin(phi),0))
    vertices_inner.append(rad*mp.Vector3(np.cos(phi),np.sin(phi),0))
geometry = [mp.Prism(vertices_outer, height=mp.inf, material=mp.Medium(index=n)),
            mp.Prism(vertices_inner, height=mp.inf, material=mp.vacuum)]
```

The following convergence plot shows the frequency for the resonant mode with $H_z$ polarization and $Q$ of ~10<sup>7</sup> as a function of resolution.

<center>
![](images/ring_freq_vs_resolution.png)
</center>

There are three important items to note. (1) The grid and `Prism` representations are each converging to a *different* frequency than the `material_function` and `Cylinder`. This is because in the limit of infinite resolution, they are *different* structures than the cylinders. (2) The `material_function` is the same structure as the `Cylinder` with no smoothing. In the limit of infinite resolution, the `material_function` and `Cylinder` converge to the same frequency. The only difference is the *rate* of convergence: the `Cylinder` is second order (due to subpixel smoothing) whereas the `material_function` is first order. See the convergence plot in the third figure from the top. (3) The grid with no interpolation contains large irregularities when compared with the grid with interpolation. This is expected because one structure is discontinuous but the other is not. Also, because they are different structures the two grids converge to *different* frequencies. To see this trend clearly requires reducing the "jumpiness" in the results for the grid with no interpolation: the Meep resolution needs to be increased beyond 200 which is already ~3X the grid resolution.

Since the grid with interpolation has already been smoothed to a *continuous* $\varepsilon(x,y)$ function, subpixel smoothing (which is not supported for `epsilon_input_file`) is not really necessary once the Yee grid resolution exceeds the input image resolution. This can be seen in the above plot: for Meep resolutions of 80 (the grid resolution) and above, the changes in the results are much smaller than those at lower resolutions. Also, higher-order interpolation schemes are not necessary because the Yee discretization is already essentially equivalent to linear interpolation.

As a practical matter, increasing the Meep resolution beyond the resolution of a grid with no interpolation is not physically meaningful because this is trying to resolve the individual pixels of an imported image. In the case of a grid imported via `epsilon_input_file`, this is not an issue because the bilinear interpolation is performed automatically by default. However, no interpolation is applied to a `material_function` unless this is somehow incorporated into the `material_function` by the user (i.e., convolving an analytic function with a smoothing kernel).

The key point is that to reduce discretization errors, any discontinuous structure simulated using Meep must be made continuous somehow. This should be done preferably using subpixel smoothing or, as a fallback in cases where subpixel smoothing is inapplicable, bilinear interpolation. If the initial structure is already continuous, no additional preprocessing is necessary.