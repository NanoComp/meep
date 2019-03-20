---
# Subpixel Smoothing
---

Meep uses a second-order accurate finite-difference scheme for discretizing [Maxwell's equations](Introduction.md#maxwells-equations). This means that the results from Meep converge to the "exact" result from the non-discretized system quadratically with the resolution Δx. However, this second-order error O(Δx<sup>2</sup>) is generally spoiled to first-order error O(Δx) if one discretizes a *discontinuous* material boundary. Moreover, directly discretizing a discontinuity in ε or μ leads to "stairstepped" interfaces that can only be varied in discrete jumps of one pixel per iteration. Meep solves both of these problems by smoothing ε and μ: before discretizing, discontinuities are smoothed into continuous transitions over a distance of one pixel Δx, using a second-order accurate averaging procedure which is summarized [below](#smoothed-permittivity-tensor-via-perturbation-theory). This feature enables the discretized solution to converge as quickly as possible to the exact solution as you increase the `resolution`.

The subpixel smoothing has four limitations: (1) it only applies to frequency-independent, lossless dielectrics; dispersive materials are [not supported](FAQ.md#can-subpixel-averaging-be-applied-to-dispersive-materials), (2) it only applies to [`GeometricObject`](Python_User_Interface.md#geometricobject)s (i.e. `Block`, `Prism`, `Sphere`, etc.); user-defined `material_function` and `epsilon_input_file` are [not supported](FAQ.md#can-subpixel-averaging-be-applied-to-a-user-defined-material-function), (3) objects with sharp corners or edges are associated with field singularities which introduce an error intermediate between first- and second-order, and (4) the fields directly *on* the interface are still at best first-order accurate. These localized errors are equivalent to currents that radiate zero power to first order; their overall effect is second order. The improved accuracy from smoothing is therefore obtained for fields evaluated off of the interface as in the [scattered Poynting flux](Python_Tutorials/Basics.md#transmittance-spectrum-of-a-waveguide-bend) integrated over a surface away from the interface, for nonlocal properties such as resonant frequencies, and for overall integrals of fields and energies to which the interface contributes only O(Δx) of the integration domain.

[TOC]

Smoothed Permittivity Tensor via Perturbation Theory
----------------------------------------------------

A second-order accurate smoothing scheme must yield zero first-order perturbation to the interface. This requires the smoothing scheme to be anisotropic. Even if the initial interface is between isotropic materials, one obtains a tensor ε (or μ) which uses the mean ε for fields parallel to the interface and the harmonic mean (inverse of mean of ε<sup>-1</sup>) for fields perpendicular to the interface. If the initial materials are anisotropic, a more complicated formula is used. They key point is that, even if the structure consists entirely of isotropic materials, the discretized structure will use anisotropic materials. For interface pixels, Meep computes the effective permittivity tensor automatically at the start of the simulation prior to time stepping via analytic expressions for the filling fraction and local normal vector. For details involving derivation of the effective permittivity tensor and its implementation in Meep/FDTD, see [Optics Letters, Vol. 36, pp. 2972-4, 2006](https://www.osapublishing.org/ol/abstract.cfm?uri=ol-31-20-2972) and [Optics Letters, Vol. 35, pp. 2778-80, 2009](https://www.osapublishing.org/abstract.cfm?uri=ol-34-18-2778).

Continuously-Varying Shapes and Results
---------------------------------------

A key feature of Meep's subpixel smoothing, particularly relevant for shape optimization (e.g. [Applied Physics Letters, Vol. 104, 091121, 2014](https://aip.scitation.org/doi/abs/10.1063/1.4867892)), is that continuously-varying the `geometry` yields continuously-varying results. This is demonstrated for a [ring resonator](Python_Tutorials/Basics.md#modes-of-a-ring-resonator): as the radius increases, the frequency of a resonant H<sub>z</sub>-polarized mode decreases. (Note that unlike the ring-resonator example in [Tutorial/Basics](Python_Tutorials/Basics.md#modes-of-a-ring-resonator) involving E<sub>z</sub>-polarized modes where the fields are always parallel to the interface, this example involves discontinuous fields.) The simulation script is shown below. The resolution is 10 pixels/μm. The inner ring radius is varied from 1.8 to 2.0 μm in gradations of 0.005 μm.

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
    fcen = -0.022*rad**3 + 0.127*rad**2 - 0.269*rad + 0.373  
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

    geometry = [mp.Cylinder(radius=rad+w,
                            material=mp.Medium(index=n)),
                mp.Cylinder(radius=rad)]
    
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

A plot of the resonant frequency versus the ring radius is shown below for subpixel smoothing (red) and no smoothing (blue). Included for reference is the "exact" result (black) computed using *no smoothing* at a resolution of 50 pixels/μm. The no smoothing result shows a stairstepped discontinuous eigenfrequency. The subpixel-smoothing result varies continuously with the ring radius and is similar to the exact result which is at a resolution five times larger. Finally, the inset shows the scalar H<sub>z</sub> field profile of the resonant mode for a structure with inner radius of 1.9 μm.

<center>
![](../images/ring_vary_radius.png)
</center>

