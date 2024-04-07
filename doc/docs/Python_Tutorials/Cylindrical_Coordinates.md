---
# Cylindrical Coordinates
---

Meep supports the simulation of Maxwell's equations in [cylindrical coordinates](https://en.wikipedia.org/wiki/Cylindrical_coordinate_system) for structures that have [continuous rotational symmetry around the *z* axis](../Exploiting_Symmetry.md#cylindrical-symmetry). This reduces problems in 3d to 2d, and 2d to 1d, if there is sufficient symmetry.

[TOC]

Modes of a Ring Resonator
-------------------------

In [Tutorial/Basics/Modes of a Ring Resonator](Basics.md#modes-of-a-ring-resonator), the modes of a ring resonator were computed by performing a 2d simulation. This example involves simulating the *same* structure while [exploiting](../Exploiting_Symmetry.md) the fact that the system has *continuous* rotational symmetry, by performing the simulation in cylindrical coordinates. The simulation script is in [examples/ring-cyl.py](https://github.com/NanoComp/meep/blob/master/python/examples/ring-cyl.py).

As always, the starting point is to import the `meep` and other library modules:

```py
import meep as mp
import argparse

def main(args):
```

The parameters of the problem are defined with exactly the same values as in the 2d simulation:

```py
n = 3.4     # index of waveguide
w = 1       # width of waveguide
r = 1       # inner radius of ring
pad = 4     # padding between waveguide and edge of PML
dpml = 2    # thickness of PML
```

The dimensions and size of the computational cell are defined:

```py
sr = r + w + pad + dpml         # radial size (cell is from 0 to sr)
dimensions = mp.CYLINDRICAL
cell = mp.Vector3(sr, 0, 0)
```

The key thing is to set the `dimensions` parameter to `CYLINDRICAL`. This means that all vectors represent ($r$,$\phi$,$z$) coordinates instead of ($x$,$y$,$z$). The computational cell in the $r$ direction is of size `sr = r + w + pad + dpml`, and runs from `0` to `sr` (by default) rather than from `-sr/2` to `sr/2` as it would for any other dimension. Note that the $z$ size is 0 because it is in 2d. The $\phi$ size is also 0, corresponding to the continuous rotational symmetry. A finite $\phi$ size might correspond to discrete rotational symmetry, but this is not currently supported.

In particular, in systems with continuous rotational symmetry, by an analogue of Bloch's theorem, the angular dependence of the fields can always be chosen in the form $\exp(i m \phi)$ for some integer $m$. Meep uses this fact to treat the angular dependence analytically, with $m$ given by the input variable `m` which is set to a command-line argument that is 3 by default.

```py
m = args.m
```

This is essentially a 1d calculation where Meep must discretize the $r$ direction only. For this reason, it will be much faster than the previous 2d calculation.

The geometry is now specified by a single `Block` object &mdash; remember that this is a block in cylindrical coordinates, so that it really specifies an annular ring:

```py
geometry = [mp.Block(center=mp.Vector3(r + (w / 2)),
                     size=mp.Vector3(w, 1e20, 1e20),
                     material=mp.Medium(index=n))]

pml_layers = [mp.PML(dpml)]
resolution = 10
```

PMLs are on "all" sides. The $z$ direction has no thickness and therefore it is automatically periodic with no PML. PML is also omitted from the boundary at $r$=0 which is handled by the analytical reflection symmetry.

The remaining inputs are almost exactly the same as in the previous 2d simulation. A single Gaussian point source is added in the $z$ direction to excite $E_z$-polarized modes, with some center frequency and width:

```py
fcen = args.fcen  # pulse center frequency
df = args.df      # pulse width (in frequency)
sources = [mp.Source(src=mp.GaussianSource(fcen, fwidth=df),
                     component=mp.Ez,
                     center=mp.Vector3(r + 0.1))]
```

Note that this isn't really a point source, however, because of the cylindrical symmetry &mdash; it is really a ring source with $\phi$ dependence $\exp(i m \phi)$. Finally, as before, the fields are timestepped until the source has turned off, plus 200 additional time units during which [Harminv](../Python_User_Interface.md#harminv) is used to analyze the $E_z$ field at a given point to extract the frequencies and decay rates of the modes.

```py
sim = mp.Simulation(cell_size=cell,
                    geometry=geometry,
                    boundary_layers=pml_layers,
                    resolution=resolution,
                    sources=sources,
                    dimensions=dimensions,
                    m=m)

sim.run(mp.after_sources(mp.Harminv(mp.Ez, mp.Vector3(r + 0.1), fcen, df)),
        until_after_sources=200)
```

At the very end, one period of the fields is output to create an animation. A single field output would be a 1d dataset along the $r$ direction, so to make things more interesting `to_appended` is used to append these datasets to a single HDF5 file to obtain an $r \times t$ 2d dataset. Also `in_volume` is used to specify a larger output volume than just the computational cell: in particular, the output is from `-sr` to `sr` in the $r$ direction, where the $-r$ field values are automatically inferred from the reflection symmetry.

```py
sim.run(mp.in_volume(mp.Volume(center=mp.Vector3(), size=mp.Vector3(2 * sr)),
                     mp.at_beginning(mp.output_epsilon),
                     mp.to_appended("ez", mp.at_every(1 / fcen / 20, mp.output_efield_z))),
        until=1 / fcen)
```

The last component of the script involves defining the three command-line arguments and their default values:

```py
if __name__ == '__main__':
   parser = argparse.ArgumentParser()
   parser.add_argument('-fcen', type=float, default=0.15, help='pulse center frequency')
   parser.add_argument('-df', type=float, default=0.1, help='pulse frequency width')
   parser.add_argument('-m', type=int, default=3, help='phi (angular) dependence of the fields given by exp(i m phi)')
   args = parser.parse_args()
   main(args)
```

The simulation is ready to be run. Recall that, in the 2d calculation, three modes were obtained in this frequency range: (1) $\omega$=0.11785 with $Q$=77 and an $m$=3 field pattern, (2) $\omega$=0.14687 with $Q$=351 and an $m$=4 field pattern, and (3) $\omega$=0.17501 with $Q$=1630 and an $m$=5 field pattern. To verify the correctness of this script, the *same* modes should be obtained with some differences due to the finite resolution, except now *three* calculations are necessary, a separate one for each value of $m$. It will still be much faster than the 2d simulation because these simulations are 1d.

In particular, the three calculations are:

```sh
unix% python ring-cyl.py -m 3 | grep harminv
unix% python ring-cyl.py -m 4 | grep harminv
unix% python ring-cyl.py -m 5 | grep harminv
```

giving the combined output:

```
harminv0:, frequency, imag. freq., Q, |amp|, amplitude, error
harminv0:, 0.11835455441250631, -0.0006907792691647415, 85.66741917111612, 0.02570190626349302, (-0.02402703883357199-0.00912630212448642j), (5.286949731053267e-10+0j)
harminv0:, 0.1475578747705309, -0.0001938438860632441, 380.61008208014414, 0.19361245519715206, (0.1447225471614173+0.12861246887677943j), (5.889273063545974e-11+0j)
harminv0:, 0.1759448592380757, -4.900590034953583e-05, 1795.1395442502285, 0.0452479314013276, (-0.014395016792255884-0.042897072017212545j), (1.6343462235932872e-10+0j)
```

This is indeed very close to the 2d simulations: the frequencies are within 1% of the previous values. The $Q$ values (lifetimes) differ by a larger amount although they are still reasonably close.

Which is more accurate, the 2d or the cylindrical simulation? This question can be answered by increasing the resolution in both cases and seeing what they converge towards. In particular, let's focus on the $m$=4 mode. In the cylindrical case, if the resolution is doubled to 20, the mode is $\omega$=0.14748 and $Q$=384. In the 2d case, if the resolution is doubled to 20 the mode is $\omega$=0.14733 and $Q$=321. It looks like the frequencies are clearly converging together and that the cylindrical simulation is more accurate (as you might expect since it describes the $\phi$ direction analytically). But the $Q$ values seem to be getting *farther* apart &mdash; what's going on?

The problem is twofold. First, there is some signal-processing error in determining $Q$ in the 2d case, as indicated by the "error" column of the `harminv` output which is only 4e-7 for the 2d simulation vs. 6e-11 for the cylindrical case. This error can be reduced by running with a narrower bandwidth source, which excites just one mode and gives a cleaner signal, or by analyzing over a longer time than 200. Doing the former, we find that the 2d value of $Q$ at a resolution of 20 should really be $Q$=343. Second, [PML](../Perfectly_Matched_Layer.md) absorbing layers are really designed to absorb planewaves incident on flat interfaces, but here we have a *cylindrical* PML layer. Because of this, there are larger numerical reflections from the PML in the cylindrical simulation, which we can rectify by pushing the PML out to a larger radius (i.e. using a larger value of `pad`) and/or increasing the PML thickness (increasing `dpml`) so that it turns on more adiabatically. In the cylindrical simulation for `resolution = 20`, if the PML thickness is increased to `dpml = 16`, the $Q$ is 343, which is in much better agreement with the 2d calculation and if the PML thickness is increased to `dpml = 32` the $Q$ is the same 343, so it seems to be converged.

This illustrates the general principle that you need to [check several parameters to ensure that results are converged in time-domain simulations](../FAQ.md#checking-convergence): the resolution, the run time, the PML thickness, etcetera.

Finally, the field images are obtained. Since one mode per `m` is being excited here anyway, according to `harminv`, there is no real need for a narrow-band source. This will be used anyway just to remind you of the general procedure, however, e.g. for the $\omega$=0.118, $m$=3 mode:

```sh
unix% python ring-cyl.py -m 3 -fcen 0.118 -df 0.01
unix% h5topng -S 2 -Zc dkbluered -C ring-cyl-eps-001200.00.h5 ring-cyl-ez.h5
```

Note that, because of the `to_appended` command, the `ring-cyl-ez.h5` file is a 160$\times$18 dataset corresponding to an $r \times t$ slice. Repeating this for all three modes results in the images:

$E_z$ for $\omega$=0.118 $m$=3 mode:

![](../images/Ring-cyl-ez-0.118.png#center)


$E_z$ for $\omega$=0.148 $m$=4 mode:

![](../images/Ring-cyl-ez-0.148.png#center)


$E_z$ for $\omega$=0.176 $m$=5 mode:

![](../images/Ring-cyl-ez-0.176.png#center)


Because only the $\phi$=0 slice is used, the visual distinction between $m$ values is much less than with the 2d simulation. What is apparent is that, as the frequency increases, the mode becomes more localized in the waveguide and the radiating field (seen in the $r \times t$ slice as curved waves extending outward) becomes less, as expected.

Sensitivity Analysis via Perturbation Theory
--------------------------------------------

For a given mode of the ring resonator, it is often useful to know how sensitive the resonant frequency $\omega$ is to small changes in the ring radius $r$ by computing its derivative $\partial\omega/\partial r$. The gradient is also a useful quantity for shape optimization because it can be paired with fast iterative methods such as [BFGS](https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm) to find local optima. The "brute-force" approach for computing the gradient is via a finite-difference approximation requiring *two* simulations of the (1) unperturbed [$\omega(r)$] and (2) perturbed [$\omega(r+\Delta r)$] structures. Since each simulation is potentially costly, an alternative approach based on semi analytics is to use [perturbation theory](https://en.wikipedia.org/wiki/Perturbation_theory) to obtain the gradient from the fields of the unperturbed structure. This involves a single simulation and is often more computationally efficient than the brute-force approach although some care is required to set up the calculation properly.  (More generally, [adjoint methods](https://math.mit.edu/~stevenj/18.336/adjoint.pdf) can be used to compute any number of derivatives with a single additional simulation.)

[Perturbation theory for Maxwell equations involving high index-contrast dielectric interfaces](http://math.mit.edu/~stevenj/papers/KottkeFa08.pdf) is reviewed in Chapter 2 of [Photonics Crystals: Molding the Flow of Light, 2nd Edition (2008)](http://ab-initio.mit.edu/book/). The formula (equation 30 on p.19) for the frequency shift $\Delta \omega$ resulting from the displacement of a block of $\varepsilon_1$-material towards $\varepsilon_2$-material by a distance $\Delta h$ (perpendicular to the boundary) is:


$$ \Delta\omega = -\frac{\omega}{2} \frac{ \iint d^2 \vec{r} \big[ (\varepsilon_1 - \varepsilon_2) |\vec{E}_{\parallel}(\vec{r})|^2 - \big(\frac{1}{\varepsilon_1} - \frac{1}{\varepsilon_2}\big)|\varepsilon\vec{E}_{\perp}|^2\big] \Delta h}{\int d^3\vec{r} \varepsilon(\vec{r})|\vec{E}(\vec{r})|^2} + O(\Delta h^2) $$


In this expression, $\vec{E}_{\parallel}(\vec{r})$ is the component of $\vec{E}$ that is parallel to the surface, and $\varepsilon\vec{E}_{\perp}$ is the component of $\varepsilon\vec{E}$ that is perpendicular to the surface. These two components are guaranteed to be continuous across an interface between two isotropic dielectric materials. In this demonstration, $\partial\omega/\partial r$ is computed using this formula and the results are validated by comparing with the finite-difference approximation: $[\omega(r+\Delta r)-\omega(r)]/\Delta r$.

There are three parts to the calculation: (1) find the resonant frequency of the unperturbed geometry using a broadband Gaussian source, (2) find the resonant mode profile of the unperturbed geometry using a narrowband source and from these fields compute the gradient via the perturbation-theory formula, and (3) find the resonant frequency of the perturbed geometry and from this compute the gradient using the finite-difference approximation. The perturbation is applied only to the inner and outer ring radii. The ring width is constant. There are two types of modes which are computed in separate simulations using different source polarizations: parallel ($E_z$) and perpendicular ($H_z$) to the interface.

The simulation script is in [examples/perturbation_theory.py](https://github.com/NanoComp/meep/blob/master/python/examples/perturbation_theory.py). The notebook is [examples/perturbation_theory.ipynb](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/perturbation_theory.ipynb).

```py
import meep as mp
import numpy as np
import argparse


def main(args):
    if args.perpendicular:
        src_cmpt = mp.Hz
        fcen = 0.21         # pulse center frequency
    else:
        src_cmpt = mp.Ez
        fcen = 0.17         # pulse center frequency

    n = 3.4                 # index of waveguide
    w = 1                   # ring width
    r = 1                   # inner radius of ring
    pad = 4                 # padding between waveguide and edge of PML
    dpml = 2                # thickness of PML
    m = 5                   # angular dependence

    pml_layers = [mp.PML(dpml)]

    sr = r + w + pad + dpml        # radial size (cell is from 0 to sr)
    dimensions = mp.CYLINDRICAL    # coordinate system is (r,phi,z) instead of (x,y,z)
    cell = mp.Vector3(sr)

    geometry = [mp.Block(center=mp.Vector3(r + (w / 2)),
                         size=mp.Vector3(w, mp.inf, mp.inf),
                         material=mp.Medium(index=n))]

    # find resonant frequency of unperturbed geometry using broadband source

    df = 0.2*fcen           # pulse width (in frequency)

    sources = [mp.Source(mp.GaussianSource(fcen,fwidth=df),
                         component=src_cmpt,
                         center=mp.Vector3(r+0.1))]

    sim = mp.Simulation(cell_size=cell,
                        geometry=geometry,
                        boundary_layers=pml_layers,
                        resolution=args.res,
                        sources=sources,
                        dimensions=dimensions,
                        m=m)

    h = mp.Harminv(src_cmpt, mp.Vector3(r+0.1), fcen, df)
    sim.run(mp.after_sources(h),
            until_after_sources=100)

    frq_unperturbed = h.modes[0].freq

    sim.reset_meep()

    # unperturbed geometry with narrowband source centered at resonant frequency

    fcen = frq_unperturbed
    df = 0.05*fcen

    sources = [mp.Source(mp.GaussianSource(fcen,fwidth=df),
                         component=src_cmpt,
                         center=mp.Vector3(r+0.1))]

    sim = mp.Simulation(cell_size=cell,
                        geometry=geometry,
                        boundary_layers=pml_layers,
                        resolution=args.res,
                        sources=sources,
                        dimensions=dimensions,
                        m=m)

    sim.run(until_after_sources=100)

    deps = 1 - n**2
    deps_inv = 1 - 1/n**2

    if args.perpendicular:
        para_integral = deps*2*np.pi*(r*abs(sim.get_field_point(mp.Ep, mp.Vector3(r)))**2 - (r+w)*abs(sim.get_field_point(mp.Ep, mp.Vector3(r+w)))**2)
        perp_integral = deps_inv*2*np.pi*(-r*abs(sim.get_field_point(mp.Dr, mp.Vector3(r)))**2 + (r+w)*abs(sim.get_field_point(mp.Dr, mp.Vector3(r+w)))**2)
        numerator_integral = para_integral + perp_integral
    else:
        numerator_integral = deps*2*np.pi*(r*abs(sim.get_field_point(mp.Ez, mp.Vector3(r)))**2 - (r+w)*abs(sim.get_field_point(mp.Ez, mp.Vector3(r+w)))**2)

    denominator_integral = sim.electric_energy_in_box(center=mp.Vector3(0.5*sr), size=mp.Vector3(sr))
    perturb_theory_dw_dR = -frq_unperturbed * numerator_integral / (4 * denominator_integral)

    sim.reset_meep()

    # perturbed geometry with narrowband source

    dr = 0.01

    sources = [mp.Source(mp.GaussianSource(fcen,fwidth=df),
                         component=src_cmpt,
                         center=mp.Vector3(r + dr + 0.1))]

    geometry = [mp.Block(center=mp.Vector3(r + dr + (w / 2)),
                         size=mp.Vector3(w, mp.inf, mp.inf),
                         material=mp.Medium(index=n))]

    sim = mp.Simulation(cell_size=cell,
                        geometry=geometry,
                        boundary_layers=pml_layers,
                        resolution=args.res,
                        sources=sources,
                        dimensions=dimensions,
                        m=m)

    h = mp.Harminv(src_cmpt, mp.Vector3(r+dr+0.1), fcen, df)
    sim.run(mp.after_sources(h),
            until_after_sources=100)

    frq_perturbed = h.modes[0].freq

    finite_diff_dw_dR = (frq_perturbed - frq_unperturbed) / dr

    print("dwdR:, {} (pert. theory), {} (finite diff.)".format(perturb_theory_dw_dR,finite_diff_dw_dR))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-perpendicular', action='store_true', help='use perpendicular field source (default: parallel field source)')
    parser.add_argument('-res', type=int, default=100, help='resolution (default: 100 pixels/um)')
    args = parser.parse_args()
    main(args)
```

There are three things to note. First, there is a built-in function `electric_energy_in_box` which calculates the integral of $\vec{E}\cdot\vec{D}/2 = \varepsilon|E|^2/2$. This is exactly the integral in the denominator, divided by 2. In cylindrical coordinates $(r,\phi,z)$, the integrand is [multiplied](https://en.wikipedia.org/wiki/Cylindrical_coordinate_system#Line_and_volume_elements) by the circumference $2\pi r$, or equivalently the integral is over an annular volume. Second, for the case involving the $H_z$ source, both parallel ($E_{\parallel}=E_{\phi}$) and perpendicular ($\varepsilon E_{\perp}=D_r$) fields are present which must be included in the numerator as separate terms. Field values anywhere in the grid obtained with `get_field_point` are [automatically interpolated](../Introduction.md#the-illusion-of-continuity); i.e., no additional post-processing is necessary. Third, when comparing the perturbation-theory result to the finite-difference approximation, there are *two* convergence parameters: the resolution and $\Delta r$. In principle, to demonstrate agreement with perturbation theory, the limit of the resolution should be taken to infinity and *then* the limit of $\Delta r$ to 0. In practice, this can be obtained by doubling the resolution at a given $\Delta r$ until it is sufficiently converged, then halving $\Delta r$ and repeating.

For an $E_z$ source (parallel to the interface) and `resolution = 100` the results are:
```
dwdR:, -0.08544696397218979 (pert. theory), -0.08521249090736038 (finite diff.)
```

Doubling the resolution to 200, the results are:
```
dwdR:, -0.08544607322081005 (pert. theory), -0.08521153501551137 (finite diff.)
```

Both results have converged to at least five digits. The relative error at resolution 200 is 0.3%. The mode has a $\omega$ of 0.175 and $Q$ of 1800.

For an $H_z$ source (perpendicular to the interface) and `resolution = 100` the results are:
```
dwdR:, -0.0805038571770864 (pert. theory), -0.07980873307536773 (finite diff.)
```
Doubling the resolution to 200, the results are:
```
dwdR:, -0.08020283464036788 (pert. theory), -0.07980880151594316 (finite diff.)
```
Both results have converged to at least three digits. The relative error at resolution 200 is 0.5%. The error is larger in this case due to the presence of the [discontinuous fields at the dielectric interface](../Subpixel_Smoothing.md). The mode has a $\omega$ of 0.208 and $Q$ of 1200.

Finally, as reference, the same calculation can be set up in Cartesian coordinates as a 2d simulation. The simulation script is in [examples/perturbation_theory_2d.py](https://github.com/NanoComp/meep/blob/master/python/examples/perturbation_theory_2d.py). There is one major difference in the 2d calculation: the mode produced by a point source in 2d is actually the $\cos(m\phi)$ mode, *not* $\exp(im\phi)$, or equivalently it is the superposition of the $\pm m$ modes. This means that computing the numerator integral does not involve just multiplying the field at a single point on the surface by $2\pi r$ &mdash; rather, it is the integral of $\cos^2(m\phi)$ which gives a factor of 1/2. (For non-circular shapes in 2d, the surface integral must be computed numerically.) The results are comparable to the cylindrical coordinate case (a 1d calculation) but the 2d simulation is much slower and less accurate at the same grid resolution.

Scattering Cross Section of a Finite Dielectric Cylinder
--------------------------------------------------------

As an alternative to the "ring" sources of the previous examples, it is also possible to launch planewaves in cylindrical coordinates. This is demonstrated in this example which involves computing the scattering cross section of a finite-height dielectric cylinder. The results for the 2d simulation involving the cylindrical ($r$, $z$) or ($\rho$, $z$) cell are validated by comparing to the same simulation in 3d Cartesian ($x$, $y$, $z$) coordinates which tends to be much slower and less accurate at the same grid resolution.

The calculation of the scattering cross section is described in [Tutorial/Basics/Mie Scattering of a Lossless Dielectric Sphere](Basics.md#mie-scattering-of-a-lossless-dielectric-sphere) which is modified for this example. A linearly-polarized ($x$) planewave is normally incident on a $z$-oriented cylinder which is enclosed by a DFT flux box. Expressed in cylindrical coordinates, an $x$-polarized planewave propagating in the $z$ direction is the sum of two circularly-polarized planewaves of opposite chirality:

$$ \hat{E}_x = \frac{1}{2} \left[e^{i\phi}(\hat{E}_\rho + i\hat{E}_\phi) + e^{-i\phi}(\hat{E}_\rho - i\hat{E}_\phi)\right] $$

A $y$-polarized planewave involves subtracting rather than adding the two terms in parentheses:

$$ \hat{E}_y = \frac{1}{2} \left[e^{i\phi}(\hat{E}_\rho + i\hat{E}_\phi) - e^{-i\phi}(\hat{E}_\rho - i\hat{E}_\phi)\right] $$
(Note, however, that for axisymmetric problems the $\hat{E}_y$ solution is merely a 90° rotation of the $\hat{E}_x$ solution.)

In principle, this involves performing *two* separate simulations for $m=\pm 1$. The scattered power from each simulation is then simply summed since the cross term in the total Poynting flux cancels for the different $m$ values when integrated over the $\phi$ direction. As a simplification, in the case of a material with isotropic permittivity (and/or real permittivity), only one of the two simulations is necessary: the scattered power is the same for $m=\pm 1$ due to the mirror (and/or conjugate) symmetry of the structure.

If one has a gyromagnetic material (which breaks mirror symmetry, conjugate symmetry, and reciprocity), then ±m simulations are generally inequivalent and one may require two separate simulations. For a given linearly-polarized planewave, the solution is computed by combining the fields from the two current sources of opposite chirality in separate runs (and subsequently computing Poynting flux or other desired quantities).

Note that a linearly-polarized planewave is *not* $m=0$, which corresponds to a field pattern that is *invariant* under rotations similar to [TE<sub>01</sub>/TM<sub>01</sub> modes](https://en.wikipedia.org/wiki/Transverse_mode). A linear polarization is the superposition of left and right circularly-polarized waves ($m=\pm 1$) and is *not* rotationally invariant; it flips sign if it is rotated by 180°.

The simulation script is in [examples/cylinder_cross_section.py](https://github.com/NanoComp/meep/blob/master/python/examples/cylinder_cross_section.py). The notebook is [examples/cylinder_cross_section.ipynb](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/cylinder_cross_section.ipynb).

```py
import meep as mp
import numpy as np
import matplotlib.pyplot as plt

r = 0.7  # radius of cylinder
h = 2.3  # height of cylinder

wvl_min = 2*np.pi*r/10
wvl_max = 2*np.pi*r/2

frq_min = 1/wvl_max
frq_max = 1/wvl_min
frq_cen = 0.5*(frq_min+frq_max)
dfrq = frq_max-frq_min
nfrq = 100

## at least 8 pixels per smallest wavelength, i.e. np.floor(8/wvl_min)
resolution = 25

dpml = 0.5*wvl_max
dair = 1.0*wvl_max

pml_layers = [mp.PML(thickness=dpml)]

sr = r+dair+dpml
sz = dpml+dair+h+dair+dpml
cell_size = mp.Vector3(sr,0,sz)

sources = [mp.Source(mp.GaussianSource(frq_cen,fwidth=dfrq,is_integrated=True),
                     component=mp.Er,
                     center=mp.Vector3(0.5*sr,0,-0.5*sz+dpml),
                     size=mp.Vector3(sr)),
           mp.Source(mp.GaussianSource(frq_cen,fwidth=dfrq,is_integrated=True),
                     component=mp.Ep,
                     center=mp.Vector3(0.5*sr,0,-0.5*sz+dpml),
                     size=mp.Vector3(sr),
                     amplitude=-1j)]

sim = mp.Simulation(cell_size=cell_size,
                    boundary_layers=pml_layers,
                    resolution=resolution,
                    sources=sources,
                    dimensions=mp.CYLINDRICAL,
                    m=-1)

box_z1 = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(0.5*r,0,-0.5*h),size=mp.Vector3(r)))
box_z2 = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(0.5*r,0,+0.5*h),size=mp.Vector3(r)))
box_r = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(r),size=mp.Vector3(z=h)))

sim.run(until_after_sources=10)

freqs = mp.get_flux_freqs(box_z1)
box_z1_data = sim.get_flux_data(box_z1)
box_z2_data = sim.get_flux_data(box_z2)
box_r_data = sim.get_flux_data(box_r)

box_z1_flux0 = mp.get_fluxes(box_z1)

sim.reset_meep()

n_cyl = 2.0
geometry = [mp.Block(material=mp.Medium(index=n_cyl),
                     center=mp.Vector3(0.5*r),
                     size=mp.Vector3(r,0,h))]

sim = mp.Simulation(cell_size=cell_size,
                    geometry=geometry,
                    boundary_layers=pml_layers,
                    resolution=resolution,
                    sources=sources,
                    dimensions=mp.CYLINDRICAL,
                    m=-1)

box_z1 = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(0.5*r,0,-0.5*h),size=mp.Vector3(r)))
box_z2 = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(0.5*r,0,+0.5*h),size=mp.Vector3(r)))
box_r  = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(r),size=mp.Vector3(z=h)))

sim.load_minus_flux_data(box_z1, box_z1_data)
sim.load_minus_flux_data(box_z2, box_z2_data)
sim.load_minus_flux_data(box_r, box_r_data)

sim.run(until_after_sources=100)

box_z1_flux = mp.get_fluxes(box_z1)
box_z2_flux = mp.get_fluxes(box_z2)
box_r_flux = mp.get_fluxes(box_r)

scatt_flux = np.asarray(box_z1_flux)-np.asarray(box_z2_flux)-np.asarray(box_r_flux)
intensity = np.asarray(box_z1_flux0)/(np.pi*r**2)
scatt_cross_section = np.divide(-scatt_flux,intensity)

if mp.am_master():
    plt.figure(dpi=150)
    plt.loglog(2*np.pi*r*np.asarray(freqs),scatt_cross_section,'bo-')
    plt.grid(True,which="both",ls="-")
    plt.xlabel('(cylinder circumference)/wavelength, 2πr/λ')
    plt.ylabel('scattering cross section, σ')
    plt.title('Scattering Cross Section of a Lossless Dielectric Cylinder')
    plt.tight_layout()
    plt.savefig("cylinder_cross_section.png")
```

Note that the "closed" DFT flux box is comprised of just three flux objects: two along $z$ and one in the radial $r$ direction. The function `get_fluxes` which computes the integral of the Poynting vector does so over the annular volume in cylindrical coordinates. There is no need for additional post-processing of the flux values.

As shown below, the results for the scattering cross section computed using cylindrical coordinates agree well with the 3d Cartesian simulation. However, there is a large discrepancy in performance: for a single Intel Xeon 4.2GHz processor, the runtime of the cylindrical simulation is nearly 90 times shorter than the 3d simulation.

![](../images/cylinder_cross_section.png#center)


Scattering of Sphere with Oblique Planewave
-------------------------------------------

It is also possible to launch an oblique incident planewave in cylindrical coordinate by decomposing the planewave $A_xe^{ik_xx+ik_yy}\hat{x} + A_ye^{ik_xx+ik_yy}\hat{y}$ into $\sum_m (J_r(r, m)\hat{r} + J_\phi(r, m)\hat{\phi})e^{im\phi}$ through [Jacobi-Anger expansion](https://en.wikipedia.org/wiki/Jacobi%E2%80%93Anger_expansion). The exact expressions of $J_r(r,m)$ and $J_\phi(r,m)$ are given [here](http://github.com/zlin-opt/axisym_meta3d_inverse_design/blob/master/Implementation_of_FDFD_with_Cylindrical_Coordinates.pdf) by Zin Lin. In the simplest case of normal incidence, $J_r(r,m)$ and $J_\phi(r,m)$ are nonzero only when $m = \pm 1$, as shown in the [previous tutorial](https://meep.readthedocs.io/en/latest/Python_Tutorials/Cylindrical_Coordinates/#scattering-cross-section-of-a-finite-dielectric-cylinder).

Given the decomposition of planewave into the sum of different current sources at each $m$, we can run individual simulations at each $m$ with their corresponding source amplitudes and record the relevant physical quantities. For some quantities such as fields, linearity implies that we can simply sum the results from each simulations; for some other quantities such as flux, orthogonality implies cross terms will be zero, and we can again simply sum the results. Moreover, simulations
at each $m$ values are embarrassingly parallel so they can be run simultaneously.

We present an example below that calculates the scattered flux of a sphere. Because of the spherical symmetry, incidence at different angle should have identical results. We can thus use this feature to check our approach. Note that because of the axial symmetry in the cylindrical coordinates, we cannot distinguish different azimuthal angles but we can distinguish different polar angles. We thus simply choose our incidence to be of form $E_ye^{ik_xx}$, and we can vary the angle of incidence by varying $k_x$.

On the other hand, because the source amplitudes $J_r(r,m)$ and $J_\phi(r,m)$ are generally not constant and extend to infinity, we used the principle of equivalence (for reference, see [Electromagnetic wave source condition](https://arxiv.org/pdf/1301.5366.pdf)) to create equivalent sources that are of finite sizes. Specifically, with the chosen incidence, the E fields in space are $E_ye^{ik_xx+ik_zz}$, and thus H fields can be computed by taking the curl; then Jacobi-Anger expansion can express the dependencies in $x$ and $y$ in terms of $m$ and $r$; afterwards, we created a box of sources surrounding the geometry and specify sources of amplitude $J = n \times H$ and $K = - n \times E$.

Empirically, we found that the Courant factor has to scale as $1/(|m|+0.5)$ in cylindrical coordinate to maintain numerical stability. By default, Meep uses the same Courant factor but instead zeros out fields near axis for $|m| > 1 $. In this tutorial, we choose to scale the Courant factor accordingly and force Meep to use the actual fields near axis via `accurate_fields_near_cylorigin=True`.

```py
import numpy as np
from scipy import special
import meep as mp
mp.verbosity(0)
r = 0.6  # size of flux box
cyl_r = 0.5 # radius of sphere
h = 2 * r  # height/diameter of sphere

wvl = 2 * np.pi * cyl_r / 4
frq_cen = 1 / wvl
dfrq = 0.2
nfrq = 1
resolution, mrange = 50, 5
dpml = 0.5 * wvl
dair = 1.0 * wvl
pml_layers = [mp.PML(thickness=dpml)]
sr = r + dair + dpml
sz = dpml + dair + h + dair + dpml
cell_size = mp.Vector3(sr, 0, sz)
n_cyl = 2.0
geometry = [mp.Sphere(material=mp.Medium(index=n_cyl), center=mp.Vector3(), radius=cyl_r)]

k_cen = 2 * np.pi * frq_cen
alpha_list = [0, np.pi/36, np.pi/24, np.pi/18, np.pi/12]
alpha_range = len(alpha_list)


src_size_tb = 2*r
src_size_side = 3*r
src_center_top = mp.Vector3(src_size_tb/2, 0, src_size_side/2)
src_center_bottom = mp.Vector3(src_size_tb/2, 0, -src_size_side/2)
src_center_side = mp.Vector3(src_size_tb, 0, 0)

scatt_flux_m = np.zeros((alpha_range, mrange+1))
for alpha_i in range(alpha_range):
    alpha = alpha_list[alpha_i]
    kxy, kz = k_cen*np.sin(alpha), k_cen * np.cos(alpha)
    amp_side = lambda v3: np.exp(1j * kz*(v3.z+src_size_side/2))
    phase_top = amp_side(src_center_top)

    for cur_m in range(0, mrange+1):
        if alpha!=0 or cur_m == 1:
            coeff_p1 = 0.5 * (1j)**(cur_m+1)
            coeff_m1 = 0.5 * (1j)**(cur_m-1)

            src_cen = src_size_tb/2
            Jpm = lambda v3: coeff_p1 * special.jv(cur_m+1, kxy * (v3.x+src_cen)) + coeff_m1 * special.jv(cur_m-1, kxy * (v3.x+src_cen))
            Jrm = lambda v3: 1j * coeff_p1 * special.jv(cur_m+1, kxy * (v3.x+src_cen)) - 1j * coeff_m1 * special.jv(cur_m-1, kxy * (v3.x+src_cen))
            Jside = (1j)**cur_m * special.jv(cur_m, kxy*src_size_tb) * kxy/k_cen

            src_t  = mp.GaussianSource(frq_cen, fwidth=dfrq)
            sourcesp = [
                mp.Source(src_t,component=mp.Er, center=src_center_bottom,size=mp.Vector3(src_size_tb), amplitude = -kz/k_cen, amp_func = Jrm),
                mp.Source(src_t,component=mp.Ep, center=src_center_bottom,size=mp.Vector3(src_size_tb), amplitude = -kz/k_cen, amp_func = Jpm),
                mp.Source(src_t,component=mp.Hr, center=src_center_bottom,size=mp.Vector3(src_size_tb), amp_func = Jpm),
                mp.Source(src_t,component=mp.Hp, center=src_center_bottom,size=mp.Vector3(src_size_tb), amplitude = -1, amp_func = Jrm),
                mp.Source(src_t,component=mp.Er, center=src_center_top,size=mp.Vector3(src_size_tb), amplitude = phase_top*kz/k_cen, amp_func = Jrm),
                mp.Source(src_t,component=mp.Ep, center=src_center_top,size=mp.Vector3(src_size_tb), amplitude = phase_top*kz/k_cen, amp_func = Jpm),
                mp.Source(src_t,component=mp.Hr, center=src_center_top,size=mp.Vector3(src_size_tb), amplitude = -phase_top, amp_func = Jpm),
                mp.Source(src_t,component=mp.Hp, center=src_center_top,size=mp.Vector3(src_size_tb), amplitude = phase_top, amp_func = Jrm),
                mp.Source(src_t,component=mp.Ez, center=src_center_side,size=mp.Vector3(z=src_size_side), amplitude = -Jrm(src_center_top)*kz/k_cen, amp_func = amp_side),
                mp.Source(src_t,component=mp.Hz, center=src_center_side,size=mp.Vector3(z=src_size_side), amplitude = Jpm(src_center_top), amp_func = amp_side),
                mp.Source(src_t,component=mp.Ep, center=src_center_side,size=mp.Vector3(z=src_size_side), amplitude = Jside, amp_func = amp_side),
            ]


            sim = mp.Simulation(
                cell_size=cell_size,
                boundary_layers=pml_layers,
                resolution=resolution,
                sources=sourcesp,
                dimensions=mp.CYLINDRICAL,
                m=cur_m,
                force_complex_fields = True,
                accurate_fields_near_cylorigin=True,
                Courant=min(0.5, 1/(abs(cur_m)+0.5)))

            box_z1 = sim.add_flux(frq_cen, dfrq, nfrq,
                mp.FluxRegion(center=mp.Vector3(0.5 * r, 0, -0.5 * h), size=mp.Vector3(r)))
            box_z2 = sim.add_flux(frq_cen, dfrq, nfrq,
                mp.FluxRegion(center=mp.Vector3(0.5 * r, 0, +0.5 * h), size=mp.Vector3(r)))
            box_r = sim.add_flux(frq_cen, dfrq, nfrq,
                mp.FluxRegion(center=mp.Vector3(r), size=mp.Vector3(z=h)))


            sim.run(until_after_sources=10)

            freqs = mp.get_flux_freqs(box_z1)
            box_z1_data = sim.get_flux_data(box_z1)
            box_z2_data = sim.get_flux_data(box_z2)
            box_r_data = sim.get_flux_data(box_r)
            box_z1_flux0 = mp.get_fluxes(box_z1)


            sim.reset_meep()

            sim = mp.Simulation(
                cell_size=cell_size,
                geometry=geometry,
                boundary_layers=pml_layers,
                resolution=resolution,
                sources=sourcesp,
                dimensions=mp.CYLINDRICAL,
                m=cur_m,
                force_complex_fields = True,
                accurate_fields_near_cylorigin=True,
                Courant=min(0.5, 1/(abs(cur_m)+0.5)))

            box_z1 = sim.add_flux(frq_cen, dfrq, nfrq,
                mp.FluxRegion(center=mp.Vector3(0.5 * r, 0, -0.5 * h), size=mp.Vector3(r)))
            box_z2 = sim.add_flux(frq_cen, dfrq, nfrq,
                mp.FluxRegion(center=mp.Vector3(0.5 * r, 0, +0.5 * h), size=mp.Vector3(r)))
            box_r = sim.add_flux(frq_cen, dfrq, nfrq,
                mp.FluxRegion(center=mp.Vector3(r), size=mp.Vector3(z=h)))


            sim.load_minus_flux_data(box_z1, box_z1_data)
            sim.load_minus_flux_data(box_z2, box_z2_data)
            sim.load_minus_flux_data(box_r, box_r_data)


            sim.run(until_after_sources=100)

            box_z1_flux = mp.get_fluxes(box_z1)
            box_z2_flux = mp.get_fluxes(box_z2)
            box_r_flux = mp.get_fluxes(box_r)

            scatt_flux_m[alpha_i, cur_m] = box_z1_flux[0] - box_z2_flux[0] - box_r_flux[0]
            sim.reset_meep()

scatt_power_m = np.zeros((alpha_range, mrange+1))
for i in range(mrange+1):
    scatt_power_m[:,i] = - 2*np.sum(scatt_flux_m[:,0:(i+1)], axis=1) + scatt_flux_m[:,0]

print(scatt_power_m)

```

The resulting `scatt_power_m` array is a table where each row `scatt_power_m[j,:]` corresponds to one angle, and
the columns `scatt_power_m[j,M]` is the sum of power contributions for `|m| ≤ M`.  For `M` sufficiently large,
these sums approach the same value, because the scattering from a sphere is angle independent.  For `M=5` there are slight (≈2%)
discrepancies between angles due to primarily discretization errors (doubling the resolution more than halves this
error).

Focusing Properties of a Binary-Phase Zone Plate
------------------------------------------------

It is also possible to compute a [near-to-far field transformation](../Python_User_Interface.md#near-to-far-field-spectra) in cylindrical coordinates. This is demonstrated in this example for a binary-phase [zone plate](https://en.wikipedia.org/wiki/Zone_plate) which is a rotationally-symmetric diffractive lens used to focus a normally-incident planewave to a single spot.

Using [scalar theory](https://en.wikipedia.org/wiki/Zone_plate#Design_and_manufacture), the radius of the $n$<sup>th</sup> zone can be computed as:

$$ r_n^2 = n\lambda (f+\frac{n\lambda}{4})$$

where $n$ is the zone index (1,2,3,...,$N$), $f$ is the focal length, and $\lambda$ is the operating wavelength. The main design variable is the number of zones $N$. The design specifications of the zone plate are similar to the binary-phase grating in [Tutorial/Mode Decomposition/Diffraction Spectrum of a Binary Grating](Mode_Decomposition.md#diffraction-spectrum-of-a-binary-grating) with refractive index of 1.5 (glass), $\lambda$ of 0.5 μm, and height of 0.5 μm. The focusing property of the zone plate is verified by the concentration of the electric-field energy density at the focal length of 0.2 mm (which lies *outside* the cell). The planewave is incident from within a glass substrate and spans the entire length of the cell in the radial direction. The cell is surrounded on all sides by PML. A schematic of the simulation geometry for a design with 25 zones and flat-surface termination is shown below. The near-field monitor is positioned at the edge of the PML and captures the scattered fields in *all* directions.


![](../images/zone_plate_schematic.png#center)


The simulation script is in [examples/zone_plate.py](https://github.com/NanoComp/meep/blob/master/python/examples/zone_plate.py). The notebook is [examples/zone_plate.ipynb](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/zone_plate.ipynb).

```py
import math

import matplotlib.pyplot as plt
import meep as mp
import numpy as np


resolution_um = 25

pml_um = 1.0
substrate_um = 2.0
padding_um = 2.0
height_um = 0.5
focal_length_um = 200
scan_length_z_um = 100
farfield_resolution_um = 10

pml_layers = [mp.PML(thickness=pml_um)]

wavelength_um = 0.5
frequency = 1 / wavelength_um
frequench_width = 0.2 * frequency

# The number of zones in the zone plate.
# Odd-numbered zones impart a π phase shift and
# even-numbered zones impart no phase shift.
num_zones = 25

# Specify the radius of each zone using the equation
# from https://en.wikipedia.org/wiki/Zone_plate.
zone_radius_um = np.zeros(num_zones)
for n in range(1, num_zones + 1):
    zone_radius_um[n-1] = math.sqrt(
        n * wavelength_um *
        (focal_length_um + n * wavelength_um / 4)
    )

size_r_um = zone_radius_um[-1] + padding_um + pml_um
size_z_um = pml_um + substrate_um + height_um + padding_um + pml_um
cell_size = mp.Vector3(size_r_um, 0, size_z_um)

# Specify a (linearly polarized) planewave at normal incidence.
sources = [
    mp.Source(
        mp.GaussianSource(
            frequency,
            fwidth=frequench_width,
            is_integrated=True
        ),
        component=mp.Er,
        center=mp.Vector3(0.5 * size_r_um, 0, -0.5 * size_z_um + pml_um),
        size=mp.Vector3(size_r_um),
    ),
    mp.Source(
        mp.GaussianSource(
            frequency,
            fwidth=frequench_width,
            is_integrated=True
        ),
        component=mp.Ep,
        center=mp.Vector3(0.5 * size_r_um, 0, -0.5 * size_z_um + pml_um),
        size=mp.Vector3(size_r_um),
        amplitude=-1j,
    ),
]

glass = mp.Medium(index=1.5)

# Add the substrate.
geometry = [
    mp.Block(
        material=glass,
        size=mp.Vector3(size_r_um, 0, pml_um + substrate_um),
        center=mp.Vector3(
            0.5 * size_r_um,
            0,
            -0.5 * size_z_um + 0.5 * (pml_um + substrate_um)
        ),
    )
]

# Add the zone plates starting with the ones with largest radius.
for n in range(num_zones - 1, -1, -1):
    geometry.append(
        mp.Block(
            material=glass if n % 2 == 0 else mp.vacuum,
            size=mp.Vector3(zone_radius_um[n], 0, height_um),
            center=mp.Vector3(
                0.5 * zone_radius_um[n],
                0,
                -0.5 * size_z_um + pml_um + substrate_um + 0.5 * height_um
            ),
        )
    )

sim = mp.Simulation(
    cell_size=cell_size,
    boundary_layers=pml_layers,
    resolution=resolution_um,
    sources=sources,
    geometry=geometry,
    dimensions=mp.CYLINDRICAL,
    m=-1,
)

# Add the near-field monitor (must be entirely in air).
n2f_monitor = sim.add_near2far(
    frequency,
    0,
    1,
    mp.Near2FarRegion(
        center=mp.Vector3(
            0.5 * (size_r_um - pml_um),
            0,
            0.5 * size_z_um - pml_um
        ),
        size=mp.Vector3(size_r_um - pml_um, 0, 0),
    ),
    mp.Near2FarRegion(
        center=mp.Vector3(
            size_r_um - pml_um,
            0,
            0.5 * size_z_um - pml_um - 0.5 * (height_um + padding_um)
        ),
        size=mp.Vector3(0, 0, height_um + padding_um),
    ),
)

fig, ax = plt.subplots()
sim.plot2D(ax=ax)
if mp.am_master():
    fig.savefig("zone_plate_layout.png", bbox_inches="tight", dpi=150)

# Timestep the fields until they have sufficiently decayed away.
sim.run(
    until_after_sources=mp.stop_when_fields_decayed(
        50.0,
        mp.Er,
        mp.Vector3(0.5 * size_r_um, 0, 0),
        1e-6
    )
)

farfields_r = sim.get_farfields(
    n2f_monitor,
    farfield_resolution_um,
    center=mp.Vector3(
        0.5 * (size_r_um - pml_um),
        0,
        -0.5 * size_z_um + pml_um + substrate_um + height_um + focal_length_um
    ),
    size=mp.Vector3(size_r_um - pml_um, 0, 0),
)

farfields_z = sim.get_farfields(
    n2f_monitor,
    farfield_resolution_um,
    center=mp.Vector3(
        0,
        0,
        -0.5 * size_z_um + pml_um + substrate_um + height_um + focal_length_um
    ),
    size=mp.Vector3(0, 0, scan_length_z_um),
)

intensity_r = (
    np.absolute(farfields_r["Ex"]) ** 2
    + np.absolute(farfields_r["Ey"]) ** 2
    + np.absolute(farfields_r["Ez"]) ** 2
)
intensity_z = (
    np.absolute(farfields_z["Ex"]) ** 2
    + np.absolute(farfields_z["Ey"]) ** 2
    + np.absolute(farfields_z["Ez"]) ** 2
)

# Plot the intensity data and save the result to disk.
fig, ax = plt.subplots(ncols=2)

ax[0].semilogy(
    np.linspace(0, size_r_um - pml_um, intensity_r.size),
    intensity_r,
    "bo-"
)
ax[0].set_xlim(-2, 20)
ax[0].set_xticks(np.arange(0, 25, 5))
ax[0].grid(True, axis="y", which="both", ls="-")
ax[0].set_xlabel(r"$r$ coordinate (μm)")
ax[0].set_ylabel(r"energy density of far fields, |E|$^2$")

ax[1].semilogy(
    np.linspace(
        focal_length_um - 0.5 * scan_length_z_um,
        focal_length_um + 0.5 * scan_length_z_um,
        intensity_z.size,
    ),
    intensity_z,
    "bo-",
)
ax[1].grid(True, axis="y", which="both", ls="-")
ax[1].set_xlabel(r"$z$ coordinate (μm)")
ax[1].set_ylabel(r"energy density of far fields, |E|$^2$")

fig.suptitle(
    f"binary-phase zone plate with focal length $z$ = {focal_length_um} μm"
)

if mp.am_master():
    fig.savefig("zone_plate_farfields.png", dpi=200, bbox_inches="tight")
```

Note that the volume specified in `get_farfields` via `center` and `size` is in cylindrical coordinates. These points must therefore lie in the $\phi = 0$ ($rz = xz$) plane. The fields $E$ and $H$ returned by `get_farfields` can be thought of as either cylindrical ($r$,$\phi$,$z$) or Cartesian ($x$,$y$,$z$) coordinates since these are the same in the $\phi = 0$ plane (i.e., $E_r=E_x$ and $E_\phi=E_y$). Also, `get_farfields` tends to gradually *slow down* as the far-field point gets closer to the near-field monitor. This performance degradation is unavoidable and is due to the larger number of $\phi$ integration points required for accurate convergence of the integral involving the Green's function which diverges as the evaluation point approaches the source point.

Shown below is the far-field energy-density profile around the focal length for both the *r* and *z* coordinate directions for three lens designs with $N$ of 25, 50, and 100. The focus becomes sharper with increasing $N$ due to the enhanced constructive interference of the diffracted beam. As the number of zones $N$ increases, the size of the focal spot (full width at half maximum) at $z = 200$ μm decreases as $1/\sqrt{N}$ (see eq. 17 of the [reference](http://zoneplate.lbl.gov/theory)). This means that doubling the resolution (halving the spot width) requires quadrupling the number of zones.

![](../images/zone_plate_farfield.png#center)

Nonaxisymmetric Dipole Sources
-------------------------------

In [Tutorial/Local Density of States/Extraction Efficiency of a Light-Emitting Diode (LED)](Local_Density_of_States.md#extraction-efficiency-of-a-light-emitting-diode-led), the extraction efficiency of an LED was computed using an axisymmetric point-dipole source at $r = 0$. This involved a single simulation with $m = \pm 1$. Simulating a point-dipole source at $r > 0$ (as shown in the schematic below) is more challenging because it is nonaxisymmetric whereas any point source at $r > 0$ is equivalent to an axisymmetric ring source.

![](../images/cyl_nonaxisymmetric_source_layout.png#center)

A point-dipole source at $r_0 > 0$ can be represented as a Dirac delta function in space: $\delta(r - r_0)\delta(\phi)\delta(z) / r_0$. (The $r_0$ factor in the denominator is [necessary to ensure correct normalization](https://math.stackexchange.com/questions/398777/dirac-delta-in-polar-coordinates).) In order to set up such a source using only axisymmetric simulations, it is necessary to expand the $\delta(\phi)$ term as a Fourier series of $\phi$: $\delta(\phi) = \frac{1}{2\pi} \sum_m e^{im\phi}$. (The Fourier transform of a Dirac delta function is a [constant](https://en.wikipedia.org/wiki/Dirac_delta_function#Fourier_transform). Each spectral component has equal weighting in its Fourier-series expansion.)

Simulating a point-dipole source involves two parts: (1) perform a series of simulations for $m = 0, 1, 2, ..., M$ for some cutoff $M$ of the Fourier-series expansion (the solutions for $\pm m$ are simply complex conjugates), and (2) because of power orthogonality, sum the results from each $m$-simulation in post processing, where the $m > 0$ terms are multiplied by two to account for the $-m$ solutions. This procedure is described in more detail below.

Physically, the *total* field $E(x,y,z)$ is a sum of $E_m(r,z)e^{im\phi}$ terms, one for the solution at each $m$ (similarly for $H$). Computing the total Poynting flux, however, involves integrating $\Re [E \times H^*]$ over a surface that includes an integral over $\phi$ in the range $[0,2\pi]$. The key point is that the cross terms $E_mH^*_ne^{i(m-n)\phi}$ integrate to zero due to Fourier orthogonality. **The total Poynting flux is therefore a sum of the Poynting fluxes calculated separately for each $m$.**

A note regarding the source polarization at $r > 0$. The $\hat{x}$ polarization in 3d (the "in-plane" polarization) corresponds to the $\hat{r}$ polarization in cylindrical coordinates. An $\hat{r}$-polarized point-dipole source involves $\hat{r}$-polarized point sources in the $m$-simulations. Even though $\hat{r}$ is in fact $\phi$-dependent, $\hat{r}$ is only evaluated at $\phi = 0$ because of $\delta(\phi)$. $\hat{r}$ is therefore equivalent to $\hat{x}$. This property does not hold for an $\hat{x}$-polarized point source at $r = 0$ (where $\delta(\phi)$ is replaced by $1/2\pi$): in that case, we write $\hat{x} = \hat{r}\cos(\phi) - \hat{\phi}\sin(\phi)$, and the $\sin$ and $\cos$ terms yield simulations for $m = \pm 1$. See also [Tutorial/Scattering Cross Section of a Finite Dielectric Cylinder](#scattering-cross-section-of-a-finite-dielectric-cylinder) which demonstrates setting up a linearly polarized planewave using a similar approach. However, in practice, a single $\hat{r}$-polarized point source at $r = 0$ is necessary for $m = \pm 1$, because that gives a circularly polarized source that emits the same power as a linearly polarized source.

Two features of this method may provide a significant speedup compared to an identical 3d simulation:

1. Convergence of the Fourier series may require only a small number ($M + 1$) of simulations. For a given source position $r$, $M$ can be estimated analytically as $M \approx k r$ where $k = n\omega/c$ is the wavenumber of the source within the source medium. This comes from the fact that a source $\sim e^{im\phi}$ at $r$ oscillates in the angular direction with a spatial frequency $m/r$, but $m/r > k$ waves are evanescent, so for $m \gtrsim kr$ the radiated power tends to drop exponentially. As an example, a point-dipole source with wavelength of $1.0 \mu m$ at a radial position of $r = 1.0 \mu m$ within a medium of $n = 2.4$ would require roughly $M = 16$ simulations. (In practice, however, we can usually truncate the Fourier-series expansion earlier without significantly degrading accuracy whenever the radiated flux at some $m$ has dropped to some small fraction of its maximum value in the summation.) The plot below shows the radiated flux vs. $m$ for three different source positions used in this tutorial example. Generally, the farther the point source is from $r = 0$, the more simulations are required for the Fourier-series summation to converge.

2. Each $m$-simulation in the Fourier-series expansion is independent of the others. The simulations can therefore be executed simultaneously using an [embarrassingly parallel](https://meep.readthedocs.io/en/latest/Parallel_Meep/#different-forms-of-parallelization) approach.

![](../images/cyl_nonaxisymmetric_source_flux_vs_m.png#center)

Note: in a simulation with `m = 0`, the real and imaginary parts of the fields are decoupled. As a runtime optimization, Meep simulates only the *real* part of the fields for this case which roughly halves the number of floating-point operations during timestepping. However, using purely real fields effectively *halves* the current source. Combining the results of the different $m$-simulations correctly using the Fourier-series expansion of the fields requires either setting `force_complex_fields=True` or multiplying the power from the `m = 0` run by four. This tutorial uses the former approach since the cost for using complex fields for only a single run among many is usually insignificant.

As a demonstration, we compute the [extraction efficiency of an LED](https://meep.readthedocs.io/en/latest/Python_Tutorials/Local_Density_of_States/#extraction-efficiency-of-a-light-emitting-diode-led) from a point dipole at $r = 0$ and three different locations at $r > 0$. The test involves verifying that the extraction efficiency is independent of the dipole location. The results are compared to an [identical calculation in 3d](https://github.com/NanoComp/meep/blob/1fe38999997f1825054fc978e473327c77169671/python/examples/extraction_eff_ldos.py#L100-L187) for which the extraction efficiency is 0.333718.

Results are shown in the table below. At this resolution, the relative error is at most ~4% even when $M + 1$ is relatively large (141). The error decreases with increasing resolution.

| `rpos` | **extraction efficiency** | **relative error** |  $M + 1$ |
|:------:|:-------------------------:|:------------------:|:--------:|
|    0   |          0.319556         |        0.042       |     1    |
|   3.5  |          0.319939         |        0.041       |    56    |
|   6.7  |          0.321860         |        0.036       |    101   |
|   9.5  |          0.324270         |        0.028       |    141   |

The extraction efficiency computed thus far is for *all* angles. To compute the extraction efficiency within an angular cone (i.e., as part of an overall calculation of the [radiation pattern](Near_to_Far_Field_Spectra.md#radiation-pattern-of-an-antenna)), we would need to surround the emitting structure with a closed box of near-field monitors. However, because the LED slab is infinitely extended a non-closed box must be used.  This will introduce [truncation errors](Near_to_Far_Field_Spectra.md#truncation-errors-from-a-non-closed-near-field-surface) which are unavoidable.

In principle, computing extraction efficiency first involves computing the radiation pattern $P(\theta, \phi)$ (the power as a function of [spherical angles](https://en.wikipedia.org/wiki/Spherical_coordinate_system)), and then computing the fraction of this power (integrated over the azimuthal angle $\phi$) that lies within a given angular cone $\theta \in [0,\theta_0]$.   By convention, $\theta = 0$ is in the $+z$ direction (the "pole") and $\theta = \pi / 2$ is $+r$ (the "equator"). It turns out that there is a simplification because we can compute the azimuthal $P(\theta) = \int P(\theta, \phi) d\phi$ more efficiently without first computing $P(\theta, \phi)$.   However, it is instructive to explain how to compute both $P(\theta, \phi)$ and the extraction efficiency.

To compute the radiation pattern $P(\theta, \phi)$ requires three steps:

1. For each simulation in the Fourier-series expansion ($m = 0, 1, ..., M$), compute the far fields $\vec{E}_m$, $\vec{H}_m$ for the desired $\theta$ points in the $rz$ ($\phi = 0$) plane, at an "infinite" radius (i.e., $R \gg \lambda$) using a [near-to-far field transformation](../Python_User_Interface.md#near-to-far-field-spectra).
2. Obtain the *total* far fields at these points, for a given $\phi$ by summing the far fields from (1): $\vec{E}_{tot}(\theta, \phi) = \vec{E}_{m=0}(\theta)e^{im\phi} + 2\sum_{m=1}^M \vec{E}_m(\theta)e^{im\phi}$ and $\vec{H}_{tot}(\theta, \phi) = \vec{H}_{m=0}(\theta)e^{im\phi} + 2\sum_{m=1}^M \vec{H}_m(\theta)e^{im\phi}$.  Note that $\vec{E}_m$ and $\vec{H}_m$ are generally complex, and are conjugates for $\pm m$.
3. Compute the radial Poynting flux $P_i(\theta_i, \phi)$ for each of $N$ points $i = 0, 1, ..., N - 1$ on the circumference using $\Re\left[\left[\vec{E}_{tot}(\theta_i, \phi) \times \vec{H}^*_{tot}(\theta_i, \phi)\right]\cdot\hat{r}\right]$.

However, if you want to compute the extraction efficiency within an angular cone given $P(\theta) = \int P(\theta, \phi) d\phi$, the calculations simplify because the cross terms in $\vec{E}_{tot} \times \vec{H}^*_{tot}$ between different $m$'s integrate to zero when integrated over $\phi$ from $0$ to $2\pi$.  Thus, one can replace step (2) with a direct computation of the powers $P(\theta)$ rather than summing the fields.  As a result, the procedure for computing the extraction efficiency within an angular cone for a dipole source at $r > 0$ involves three steps:

1. For each simulation in the Fourier-series expansion ($m = 0, 1, ..., M$), compute the far fields $\vec{E}_m$, $\vec{H}_m$ for the desired $\theta$ points in the $rz$ ($\phi = 0$) plane, at an "infinite" radius (i.e., $R \gg \lambda$) using a near-to-far field transformation.
2. Obtain the powers $P(\theta)$ from these far fields by summing: $P(\theta) = \int_{0}^{2\pi} \left[ P_{m=0}(\theta) + 2\sum_{m=1}^{M} P_m(\theta)  \right] d\phi =  2\pi \Re\left[ \left[\vec{E}_{m=0}(\theta) \times \vec{H}^*_{m=0}(\theta)\right]\cdot\hat{r} + 2\sum_{m=1}^M \left[\vec{E}_{m}(\theta) \times \vec{H}^*_{m}(\theta)\right]\cdot\hat{r} \right]$.
3. Compute the extraction efficiency within an angular cone $\left[ \int_0^\theta P(\theta') \sin(\theta') d\theta' \right] / \left[ \int_0^{\pi/2} P(\theta') \sin(\theta') d\theta' \right]$ by some discretized integral, e.g. a [trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule). See [Tutorial/Cylindrical Coordinates/Radiation Pattern of a Disc in Cylindrical Coordinates](Near_to_Far_Field_Spectra.md#radiation-pattern-of-a-disc-in-cylindrical-coordinates).

The simulation script is in [examples/point_dipole_cyl.py](https://github.com/NanoComp/meep/blob/master/python/examples/point_dipole_cyl.py).

```py
from typing import Tuple

import meep as mp
import numpy as np


RESOLUTION_UM = 50
WAVELENGTH_UM = 1.0
N_SLAB = 2.4
SLAB_THICKNESS_UM = 0.7 * WAVELENGTH_UM / N_SLAB


def dipole_in_slab(zpos: float, rpos_um: float, m: int) -> Tuple[float, float]:
    """Computes the flux from a dipole in a slab.

    Args:
      zpos: position of dipole as a fraction of layer thickness.
      rpos_um: position of source in radial direction.
      m: angular φ dependence of the fields exp(imφ).

    Returns:
      A 2-tuple of the radiated and total flux.
    """
    pml_um = 1.0  # thickness of PML
    padding_um = 1.0  # thickness of air padding
    r_um = 20.0  # length of cell in r

    frequency = 1 / WAVELENGTH_UM  # center frequency of source/monitor

    # runtime termination criteria
    flux_decay_threshold = 1e-4

    size_r = r_um + pml_um
    size_z = SLAB_THICKNESS_UM + padding_um + pml_um
    cell_size = mp.Vector3(size_r, 0, size_z)

    boundary_layers = [
        mp.PML(pml_um, direction=mp.R),
        mp.PML(pml_um, direction=mp.Z, side=mp.High),
    ]

    src_pt = mp.Vector3(rpos_um, 0, -0.5 * size_z + zpos * SLAB_THICKNESS_UM)
    sources = [
        mp.Source(
            src=mp.GaussianSource(frequency, fwidth=0.05 * frequency),
            component=mp.Er,
            center=src_pt,
        ),
    ]

    geometry = [
        mp.Block(
            material=mp.Medium(index=N_SLAB),
            center=mp.Vector3(0, 0, -0.5 * size_z + 0.5 * SLAB_THICKNESS_UM),
            size=mp.Vector3(mp.inf, mp.inf, SLAB_THICKNESS_UM),
        )
    ]

    sim = mp.Simulation(
        resolution=RESOLUTION_UM,
        cell_size=cell_size,
        dimensions=mp.CYLINDRICAL,
        m=m,
        boundary_layers=boundary_layers,
        sources=sources,
        geometry=geometry,
        force_complex_fields=True
    )

    flux_mon = sim.add_flux(
        frequency,
        0,
        1,
        mp.FluxRegion(
            center=mp.Vector3(0.5 * r_um, 0, 0.5 * size_z - pml_um),
            size=mp.Vector3(r_um, 0, 0),
        ),
        mp.FluxRegion(
            center=mp.Vector3(r_um, 0, 0.5 * size_z - pml_um - 0.5 * padding_um),
            size=mp.Vector3(0, 0, padding_um),
        ),
    )

    sim.run(
        mp.dft_ldos(frequency, 0, 1),
        until_after_sources=mp.stop_when_dft_decayed(
            tol=flux_decay_threshold
        ),
    )

    radiated_flux = mp.get_fluxes(flux_mon)[0]

    # volume of the ring current source
    delta_vol = 2 * np.pi * rpos_um / (RESOLUTION_UM**2)

    # total flux from point source via LDOS
    source_flux = (-np.real(sim.ldos_Fdata[0] * np.conj(sim.ldos_Jdata[0])) *
                   delta_vol)

    print(f"flux-cyl:, {rpos_um:.2f}, {m:3d}, "
          f"{source_flux:.6f}, {radiated_flux:.6f}")

    return radiated_flux, source_flux


if __name__ == "__main__":
    dipole_height = 0.5

    # An Er source at r = 0 needs to be slightly offset.
    # https://github.com/NanoComp/meep/issues/2704
    dipole_rpos_um = 1.5 / RESOLUTION_UM

    # Er source at r = 0 requires a single simulation with m = ±1.
    m = 1
    radiated_flux, source_flux = dipole_in_slab(
        dipole_height,
        dipole_rpos_um,
        m,
    )
    extraction_efficiency = radiated_flux / source_flux
    print(f"exteff:, {dipole_rpos_um}, {extraction_efficiency:.6f}")

    # Er source at r > 0 requires Fourier-series expansion of φ.

    # Threshold flux to determine when to truncate expansion.
    flux_decay_threshold = 1e-2

    dipole_rpos_um = [3.5, 6.7, 9.5]
    for rpos_um in dipole_rpos_um:
        source_flux_total = 0
        radiated_flux_total = 0
        radiated_flux_max = 0
        m = 0
        while True:
            radiated_flux, source_flux = dipole_in_slab(
                dipole_height,
                rpos_um,
                m,
            )
            radiated_flux_total += radiated_flux * (1 if m == 0 else 2)
            source_flux_total += source_flux * (1 if m == 0 else 2)

            if radiated_flux > radiated_flux_max:
                radiated_flux_max = radiated_flux

            if (m > 0 and
                (radiated_flux / radiated_flux_max) < flux_decay_threshold):
                break
            else:
                m += 1

        extraction_efficiency = radiated_flux_total / source_flux_total
        print(f"exteff:, {rpos_um}, {extraction_efficiency:.6f}")
```
