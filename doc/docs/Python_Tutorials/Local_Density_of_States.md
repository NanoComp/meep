---
# Local Density of States
---

In this example, we will demonstrate the local [density of states](https://en.wikipedia.org/wiki/Density_of_states) (LDOS) feature by investigating the Purcell enhancement phenomena in two different metallic cavities. The LDOS, in general, has many important uses for understanding classical dipole sources, but also in many physical phenomena that can be understood semiclassically in terms of dipole currents &mdash; for example, the [spontaneous emission](https://en.wikipedia.org/wiki/Spontaneous_emission) rate of atoms (key to fluorescence and lasing phenomena) is proportional to the LDOS.

The LDOS is equivalent to the power radiated by a unit dipole, $P=\frac{1}{2}\operatorname{Re}[\mathbf{E}^*\cdot\mathbf{J}]$, which, alternatively, is really just a measure of how much the harmonic modes of a system overlap with the source point. Also, the LDOS is proportional to the [radiation resistance](https://en.wikipedia.org/wiki/Radiation_resistance) of a dipole antenna. It is a useful quantity in electromagnetism due to the fact that the *same* current radiates a *different* amount of power depending on the surrounding environment. Analytically, the per-polarization LDOS is exactly proportional to the power radiated by an $\ell$-oriented point-dipole current, $p(t)$, at a given position in space. For a more mathematical treatment of the theory behind the LDOS, see Section 4.4 ("Currents and Fields: The Local Density of States") in [Chapter 4](http://arxiv.org/abs/arXiv:1301.5366) ("Electromagnetic Wave Source Conditions") of the book [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707), but for now we simply give the result:

$$\operatorname{LDOS}_{\ell}(\vec{x}_0,\omega)=-\frac{2}{\pi}\varepsilon(\vec{x}_0)\frac{\operatorname{Re}[\hat{E}_{\ell}(\vec{x}_0,\omega)\hat{p}(\omega)^*]}{|\hat{p}(\omega)|^2}$$

where the $|\hat{p}(\omega)|^2$ normalization is necessary for obtaining the power exerted by a unit-amplitude dipole assuming linear materials. In FDTD, computing the LDOS is straightforward: excite a point dipole source and accumulate the Fourier transforms of the field at a given point in space to obtain the entire LDOS spectrum in a single calculation. This is implemented in the `dft_ldos` feature which is the subject of this tutorial.

[TOC]

Planar Cavity with Lossless Metallic Walls
------------------------------------------

The spontaneous-emission rate of a point-dipole emitter in a planar cavity bounded by a lossless metallic mirror can be tuned using the thickness of the cavity. A schematic of the cavity structure is shown in the figure below. The planar cavity consists of two mirrors separated in the $z$ direction by a distance $L$. The cavity consists of a homogeneous dielectric with $n$=2.4. The dipole wavelength ($\lambda$) is 1.0 μm with polarization *parallel* to the mirrors. The Purcell enhancement factor can be computed analytically as a function of the cavity thickness expressed in units of the wavelength in the cavity medium ($nL/\lambda$) using equation 7 of [IEEE J. Quantum Electronics, Vol. 34, pp. 71-76 (1998)](https://ieeexplore.ieee.org/abstract/document/655009). We will validate the simulated results using the analytic theory. Since this system has cylindrical symmetry, we can perform an identical simulation in [cylindrical coordinates](Cylindrical_Coordinates.md) which is significantly faster because it is a 2D calculation.

In this demonstration, the cavity thickness is swept over a range of 0.5 to 2.5. Below a thickness of 0.5 there are no guided modes and thus the Purcell enhancement factor is zero.

Two types of simulations are necessary for computing the Purcell enhancement factor: (1) bulk cavity medium and (2) cavity. The `dft_ldos` feature is used to compute the LDOS in each case. The Purcell enhancement factor is computed as the ratio of the LDOS measured in (2) to that from (1).

In 3D, each simulation uses three [mirror symmetries](../Exploiting_Symmetry.md#rotations-and-reflections) to reduce the size of the computation by a factor of eight. The dipole is polarized in the $x$ direction. The cavity is infinitely extended in the $xy$ plane and the cell is therefore terminated using PMLs in $x$ and $y$. Because Meep uses a default boundary condition of a perfect-electric conductor, there is no need to explicitly define the boundaries in the $z$ direction. The fields are timestepped until they have decayed away sufficiently due to absorption by the PMLs.

In cylindrical coordinates, the dipole is polarized in the $r$ direction. Setting up a linearly polarized source in cylindrical coordinates is demonstrated in [Tutorial/Cylindrical Coordinates/Scattering Cross Section of a Finite Dielectric Cylinder](Cylindrical_Coordinates.md#scattering-cross-section-of-a-finite-dielectric-cylinder). However, all that is necessary in this example which involves a single point dipole rather than a planewave is one simulation involving an $E_r$ source at $r=0$ with $m=-1$. This is actually a circularly polarized source but this is sufficient because the $m=+1$ simulation produces an identical result to the $m=-1$ simulation. It is therefore not necessary to perform two separate simulations for $m=\pm 1$ in order to average the results from the left- and right-circularly polarized sources.

One important parameter when setting up this calculation is the grid resolution.

A key feature of the LDOS in this geometry is that it experiences discontinuities, called  [Van Hove singularities](https://en.wikipedia.org/wiki/Van_Hove_singularity), any time the cavity thickness/λ passes through the cutoff for a waveguide mode, which occurs for cavity-thickness/λ values of 0.5, 1.5, 2.5, etc.   (Mathematically, Van Hove singularities depend strongly on the dimensionality — it is a discontinuity in this case because the waves are propagating along two dimensions, i.e. each cutoff is a minimum in the 2d dispersion relation $\omega(k_x,k_y)$.)  This discontinuity also means that the LDOS *exactly at* the cutoff thickness/λ is ill-defined and convergence with discretization can be problematic at this point.  (In consequence, the LDOS *exactly* at the Van Hove discontinuity can behave erratically with resolution, and should be viewed with caution.)

As shown in the plot below, the results from Meep for both coordinate systems agree well with the analytic theory over the entire range of values of the cavity thickness.


![](../images/planar_cavity_purcell_enhancement.png#center)


The simulation script is [examples/planar_cavity_ldos.py](https://github.com/NanoComp/meep/blob/master/python/examples/planar_cavity_ldos.py).

```py
import meep as mp
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt


# important note:
# Meep may round the cell dimensions to an integer number
# of pixels which could modify the cavity structure.
resolution = 70  # pixels/μm


dpml = 0.5       # thickness of PML
L = 6.0          # length of non-PML region
n = 2.4          # refractive index of surrounding medium
wvl = 1.0        # wavelength (in vacuum)

fcen = 1/wvl


def bulk_ldos_cyl():
    sr = L+dpml
    sz = L+2*dpml
    cell_size = mp.Vector3(sr,0,sz)

    pml_layers = [mp.PML(dpml)]

    sources = [mp.Source(src=mp.GaussianSource(fcen,fwidth=0.2*fcen),
                         component=mp.Er,
                         center=mp.Vector3())]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=pml_layers,
                        sources=sources,
                        dimensions=mp.CYLINDRICAL,
                        m=-1,
                        default_material=mp.Medium(index=n))

    sim.run(mp.dft_ldos(fcen,0,1),
            until_after_sources=mp.stop_when_fields_decayed(20,
                                                            mp.Er,
                                                            mp.Vector3(),
                                                            1e-6))

    return sim.ldos_data[0]


def cavity_ldos_cyl(sz):
    sr = L+dpml
    cell_size = mp.Vector3(sr,0,sz)

    pml_layers = [mp.PML(dpml,direction=mp.R)]

    sources = [mp.Source(src=mp.GaussianSource(fcen,fwidth=0.2*fcen),
                         component=mp.Er,
                         center=mp.Vector3())]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=pml_layers,
                        sources=sources,
                        dimensions=mp.CYLINDRICAL,
                        m=-1,
                        default_material=mp.Medium(index=n))

    sim.run(mp.dft_ldos(fcen,0,1),
            until_after_sources=mp.stop_when_fields_decayed(20,
                                                            mp.Er,
                                                            mp.Vector3(),
                                                            1e-6))

    return sim.ldos_data[0]


def bulk_ldos_3D():
    s = L+2*dpml
    cell_size = mp.Vector3(s,s,s)

    pml_layers = [mp.PML(dpml)]

    sources = [mp.Source(src=mp.GaussianSource(fcen,fwidth=0.2*fcen),
                         component=mp.Ex,
                         center=mp.Vector3())]

    symmetries = [mp.Mirror(direction=mp.X,phase=-1),
                  mp.Mirror(direction=mp.Y),
                  mp.Mirror(direction=mp.Z)]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=pml_layers,
                        sources=sources,
                        symmetries=symmetries,
                        default_material=mp.Medium(index=n))

    sim.run(mp.dft_ldos(fcen,0,1),
            until_after_sources=mp.stop_when_fields_decayed(20,
                                                            mp.Ex,
                                                            mp.Vector3(),
                                                            1e-6))

    return sim.ldos_data[0]


def cavity_ldos_3D(sz):
    sxy = L+2*dpml
    cell_size = mp.Vector3(sxy,sxy,sz)

    boundary_layers = [mp.PML(dpml,direction=mp.X),
                       mp.PML(dpml,direction=mp.Y)]

    sources = [mp.Source(src=mp.GaussianSource(fcen,fwidth=0.2*fcen),
                         component=mp.Ex,
                         center=mp.Vector3())]

    symmetries = [mp.Mirror(direction=mp.X,phase=-1),
                  mp.Mirror(direction=mp.Y),
                  mp.Mirror(direction=mp.Z)]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=boundary_layers,
                        sources=sources,
                        symmetries=symmetries,
                        default_material=mp.Medium(index=n))

    sim.run(mp.dft_ldos(fcen,0,1),
            until_after_sources=mp.stop_when_fields_decayed(20,
                                                            mp.Ex,
                                                            mp.Vector3(),
                                                            1e-6))

    return sim.ldos_data[0]


if __name__ == '__main__':
    ldos_bulk_cyl = bulk_ldos_cyl()
    ldos_bulk_3D = bulk_ldos_3D()

    # units of wavelength in cavity medium
    cavity_thickness = np.arange(0.50,2.55,0.05)

    gap = cavity_thickness*wvl/n

    ldos_cavity_cyl = np.zeros(len(cavity_thickness))
    ldos_cavity_3D = np.zeros(len(cavity_thickness))
    for idx,g in enumerate(gap):
        ldos_cavity_cyl[idx] = cavity_ldos_cyl(g)
        ldos_cavity_3D[idx] = cavity_ldos_3D(g)
        print("purcell-enh:, {:.3f}, "
              "{:.6f} (cyl.), {:.6f} (3D)".format(cavity_thickness[idx],
                                                  ldos_cavity_cyl[idx]/ldos_bulk_cyl,
                                                  ldos_cavity_3D[idx]/ldos_bulk_3D))

    # Purcell enhancement factor (relative to bulk medium)
    pe_meep_cyl = ldos_cavity_cyl / ldos_bulk_cyl
    pe_meep_3D = ldos_cavity_3D / ldos_bulk_3D

    # equation 7 of reference
    pe_theory = (3*np.fix(cavity_thickness+0.5)/(4*cavity_thickness) +
                 (4*np.power(np.fix(cavity_thickness+0.5),3) -
                  np.fix(cavity_thickness+0.5)) /
                 (16*np.power(cavity_thickness,3)))

    if mp.am_master():
        plt.plot(cavity_thickness,pe_meep_3D,'b-',label='Meep (3D)')
        plt.plot(cavity_thickness,pe_meep_cyl,'r-',label='Meep (cylin.)')
        plt.plot(cavity_thickness,pe_theory,'g-',label='theory')
        plt.plot(cavity_thickness,np.ones(len(cavity_thickness)),'k--')
        plt.xlabel('cavity thickness, $nL/\lambda$')
        plt.ylabel('Purcell enhancement factor')
        plt.title("parallel point dipole at λ=1.0 μm in a planar cavity\n"
                  "with n=2.4 and lossless metallic walls")
        plt.axis([0.5,2.5,0.4,3.1])
        plt.legend()
        plt.savefig('cavity_purcell_factor_vs_thickness.png',
                    bbox_inches='tight')
```

Square Box with a Small Opening
-------------------------------

We consider the simple example of a 2D perfect-metal $a$x$a$ cavity of finite thickness 0.1$a$, with a small notch of width $w$ on one side that allows the modes to escape. The nice thing about this example is that in the absence of the notch, the lowest-frequency $E_z$-polarized mode is known analytically to be $E_z^{(1)}=\frac{4}{a^2}\sin(\pi x/a)\sin(\pi \gamma/a)$, with a frequency $\omega^{(1)}=\sqrt{2}\pi c/a$ and modal volume $V^{(1)}=a^2/4$. The notch slightly perturbs this solution, but more importantly the opening allows the confined mode to radiate out into the surrounding air, yielding a finite $Q$. For $w \ll a$, this radiative escape occurs via an evanescent (sub-cutoff) mode of the channel waveguide formed by the notch, and it follows from inspection of the evanescent decay rate $\sqrt{(\pi/\omega)^2-(\omega^{(1)})^2}/c$ that the lifetime scales asymptotically as $Q^{(1)} \sim e^{\#/\omega}$ for some coefficient \#.

We will validate both this prediction and the expression for the LDOS shown above by computing the LDOS at the center of the cavity, the point of peak $|\vec{E}|$, in two ways. First, we compute the LDOS directly from the power radiated by a dipole, Fourier-transforming the result of a pulse using the `dft_ldos` command. Second, we compute the cavity mode and its lifetime $Q$ using Harminv and then compute the LDOS by the Purcell formula shown above. The latter technique is much more efficient for high Q (small $w$), since one must run the simulation for a very long time to directly accumulate the Fourier transform of a slowly-decaying mode. The two calculations, we will demonstrate, agree to within discretization error, verifying the LDOS analysis above, and $Q/V$ is asymptotically linear on a semilog scale versus $1/w$ as predicted.

A lossless localized mode yields a $\delta$-function spike in the LDOS, whereas a <i>lossy</i>, arising from either small absorption or radiation, localized mode &mdash; a resonant cavity mode &mdash; leads to a Lorentzian peak. The large enhancement in the LDOS at the resonant peak is known as a [Purcell effect](https://en.wikipedia.org/wiki/Purcell_effect), named after Purcell's proposal for enhancing spontaneous emission of an atom in a cavity. This is analogous to a microwave antenna resonating in a metal box. In this case, the resonant mode's contribution to the LDOS at $\omega^{(n)}$ can be shown to be:

$$\operatorname{resonant\ LDOS} \approx \frac{2}{\pi\omega^{(n)}} \frac{Q^{(n)}}{V^{(n)}}$$

where $Q^{(n)}=\omega^{(n)}/2\gamma^{(n)}$ is the dimensionless [quality factor](https://en.wikipedia.org/wiki/Q_factor) and $V^{(n)}$ is the modal volume. This represents another way to compute the LDOS. In this tutorial, we will verify this expression by comparing it to the earlier one.

The simulation script is [examples/metal-cavity-ldos.py](https://github.com/NanoComp/meep/blob/master/python/examples/metal-cavity-ldos.py). The notebook is [examples/metal-cavity-ldos.ipynb](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/metal-cavity-ldos.ipynb).

```py
import math
import meep as mp
import numpy as np
import matplotlib.pyplot as plt

def metal_cavity(w):
    resolution = 50
    sxy = 2
    dpml = 1
    sxy = sxy+2*dpml
    cell = mp.Vector3(sxy,sxy)

    pml_layers = [mp.PML(dpml)]
    a = 1
    t = 0.1
    geometry = [mp.Block(mp.Vector3(a+2*t,a+2*t,mp.inf), material=mp.metal),
                mp.Block(mp.Vector3(a,a,mp.inf), material=mp.air)]

    geometry.append(mp.Block(center=mp.Vector3(a/2),
                             size=mp.Vector3(2*t,w,mp.inf),
                             material=mp.air))

    fcen = math.sqrt(0.5)/a
    df = 0.2
    sources = [mp.Source(src=mp.GaussianSource(fcen,fwidth=df),
                         component=mp.Ez,
                         center=mp.Vector3())]

    symmetries = [mp.Mirror(mp.Y)]

    sim = mp.Simulation(cell_size=cell,
                        geometry=geometry,
                        boundary_layers=pml_layers,
                        sources=sources,
                        symmetries=symmetries,
                        resolution=resolution)

    h = mp.Harminv(mp.Ez, mp.Vector3(), fcen, df)
    sim.run(mp.after_sources(h), until_after_sources=500)

    m = h.modes[0]
    f = m.freq
    Q = m.Q
    Vmode = 0.25*a*a
    ldos_1 = Q / Vmode / (2 * math.pi * f * math.pi * 0.5)

    sim.reset_meep()

    T = 2*Q*(1/f)
    sim.run(mp.dft_ldos(f,0,1), until_after_sources=T)
    ldos_2 = sim.ldos_data[0]

    return ldos_1, ldos_2
```

We need to run for a sufficiently long time to ensure that the Fourier-transformed fields have converged. A suitable time interval is, due to the Fourier Uncertainty Principle, just one period of the decay which we can determine using the $Q$ we calculated previously. The smaller the notch size becomes and the higher the corresponding $Q$ of the mode, the longer the simulation has to run. This is why the former calculation is much more efficient for slowly-decaying modes.

We run several simulations spanning a number of different notch sizes and plot the result in the following figure which shows good agreement between the two methods.

```py
ws = np.arange(0.2,0.5,0.1)
ldos_1 = np.zeros(len(ws))
ldos_2 = np.zeros(len(ws))

for j in range(len(ws)):
    ldos_1[j], ldos_2[j] = metal_cavity(ws[j])
    print("ldos:, {}, {}".format(ldos_1[j],ldos_2[2]))
```

```py
plt.figure(dpi=150)
plt.semilogy(1/ws,ldos_1,'bo-',label="2Q/(πωV)")
plt.semilogy(1/ws,ldos_2,'rs-',label="LDOS")
plt.xlabel('a/w')
plt.ylabel('2Q/(πωW) or LDOS')
plt.show()
```

![](../images/Metalcavity_ldos.png#center)


Extraction Efficiency of a Light-Emitting Diode (LED)
-----------------------------------------------------

Another application of the LDOS feature involves computing the extraction efficiency of a dipole emitter. The extraction efficiency is used to determine the external [quantum efficiency](https://en.wikipedia.org/wiki/Quantum_efficiency) (EQE) of light-emitting diodes (LEDs). The spontaneous-emission rate, computed in the previous examples, is proportional to the internal quantum efficiency (IQE). The EQE is the product of the IQE and the extraction efficiency.

The extraction efficiency is defined as the ratio of the power extracted from the device to the total power generated by the dipole emitter. Computing the extracted power is straightforward: use a set of DFT monitors to capture all the outgoing flux. By Poynting's theorem, as long as all the outgoing flux is measured, the size and orientation of the DFT monitors is arbitrary. To compute the total power emitted by the dipole, one approach would be to enclose the point dipole with a six-sided box (in 3D) of DFT field monitors. The challenge, however, with this approach is that whenever the dipole is adjacent to an object (e.g., a surface or scatterer) then the DFT box must be made small enough to enclose *only* the dipole embedded within the homogeneous source medium. Because of this requirement, for cases in which the dipole is just a few nanometers from a nearby object (particularly if it is a lossy metal involving surface plasmon polaritons), the grid resolution may need to be significantly increased to obtain accurate results. Unfortunately, this can make 3D computations quite expensive.

A more efficient approach to computing the total emitted power from the dipole is to use the LDOS. In particular, the terms in the numerator from the LDOS expression above are exactly equivalent to the emitted power: $$-\operatorname{Re}[\hat{E}_{\ell}(\vec{x}_0,\omega)\hat{p}(\omega)^*].$$ The Fourier-transformed electric field $\hat{E}_{\ell}(\vec{x}_0,\omega)$ and current source $\hat{p}(\omega)$ can be obtained using the `ldos_Fdata` and `ldos_Jdata` members of the `Simulation` class object, respectively. The LDOS approach has the advantage that it involves just a single point monitor (colocated with the dipole source) and thus does not require the expense of collecting DFT fields for the flux through an enclosing box (and in some circumstances resolving the Poynting flux may also require higher resolution).

To demonstrate this feature of the LDOS, we will compute the extraction efficiency of an LED-like structure consisting of a point dipole embedded within a dielectric layer ($n$=2.4 and thickness 0.7 $\lambda/n$) surrounded by vacuum and positioned above an infinite ground plane of a lossless metal. We will compute the extraction efficiency (at a wavelength $\lambda$ = 1.0 μm for a dipole with polarization parallel to the ground plane) as a function of the height of the dipole source above the ground plane using two different coordinate systems — 3D Cartesian and cylindrical — and demonstrate that the results are equivalent (apart from discretization error). The key difference is that simulation using cylindrical (2D) coordinates is significantly faster and more memory efficient than 3D.

The simulation setup is shown in the figures below for 3D Cartesian (cross section in $xz$) and cylindrical coordinates. (Note that this single-dipole calculation differs from the somewhat related flux calculation in [Tutorials/Custom Source/Stochastic Dipole Emission in Light Emitting Diodes](Custom_Source.md#stochastic-dipole-emission-in-light-emitting-diodes) which involved modeling a *collection* of dipoles.) In this example, the point-dipole source is positioned at $r=0$ which involves a single simulation. Nonaxisymmetric dipoles positioned at $r>0$, however, are more challenging because they involve multiple simulations. For a demonstration, see [Cylindrical Coordinates/Nonaxisymmetric Dipole Sources](Cylindrical_Coordinates.md#nonaxisymmetric-dipole-sources).

Note: because of a [bug](https://github.com/NanoComp/meep/issues/2704) for an $E_r$ point source at $r = 0$ and $m = \pm 1$ simulation, it is necessary to slightly offset the source to $r = 1.5\Delta r$. This incurs a small error which decreases linearly with resolution.

![](../images/dipole_extraction_eff_3D.png#center)

![](../images/dipole_extraction_eff_cyl.png#center)

The total emitted power obtained from the LDOS terms of the formula above must be multiplied by $\Delta V$, the volume of the voxel. In cylindrical coordinates, $\Delta V = \Delta r \times \Delta z \times 2 \pi r$. Meep implements an $r = 0$ source at $r = 0.5 \Delta r$, corresponding to the smallest-$r$ $E_r$ Yee grid point. This means that for a source at $r = 0$, $\Delta V = \pi / resolution^3$ since $\Delta r = \Delta z = 1 / resolution$. In 3D, $\Delta V = \Delta x \times \Delta y \times \Delta z = 1 / resolution^3$ for every voxel in the cell.

As shown in the figure below, the results from the two coordinate systems have good agreement.

The simulation script is [examples/extraction_eff_ldos.py](https://github.com/NanoComp/meep/blob/master/python/examples/extraction_eff_ldos.py).

```py
import matplotlib.pyplot as plt
import meep as mp
import numpy as np


resolution = 80  # pixels/μm
dpml = 0.5  # thickness of PML
dair = 1.0  # thickness of air padding
L = 6.0  # length of non-PML region
n = 2.4  # refractive index of surrounding medium
wvl = 1.0  # wavelength (in vacuum)

fcen = 1 / wvl  # center frequency of source/monitor

# runtime termination criteria
tol = 1e-8


def extraction_eff_cyl(dmat: float, h: float) -> float:
    """Computes the extraction efficiency in cylindrical coordinates.

    Args:
      dmat: thickness of dielectric layer.
      h: height of dipole above ground plane as fraction of dmat.

    Returns:
      The extraction efficiency of the dipole within the dielecric layer.
    """
    sr = L + dpml
    sz = dmat + dair + dpml
    cell_size = mp.Vector3(sr, 0, sz)

    boundary_layers = [
        mp.PML(dpml, direction=mp.R),
        mp.PML(dpml, direction=mp.Z, side=mp.High),
    ]

    src_cmpt = mp.Er

    # Because (1) Er is not defined at r=0 on the Yee grid, and (2) there
    # seems to be a bug in the interpolation of an Er point source at r=0,
    # the source is placed at r=~Δr (just outside the first voxel).
    # This incurs a small error which decreases linearly with resolution.
    # Ref: https://github.com/NanoComp/meep/issues/2704
    src_pt = mp.Vector3(1.5 / resolution, 0, -0.5 * sz + h * dmat)

    sources = [
        mp.Source(
            src=mp.GaussianSource(fcen, fwidth=0.1 * fcen),
            component=src_cmpt,
            center=src_pt,
        )
    ]

    geometry = [
        mp.Block(
            material=mp.Medium(index=n),
            center=mp.Vector3(0, 0, -0.5 * sz + 0.5 * dmat),
            size=mp.Vector3(mp.inf, mp.inf, dmat),
        )
    ]

    sim = mp.Simulation(
        resolution=resolution,
        cell_size=cell_size,
        dimensions=mp.CYLINDRICAL,
        m=-1,
        boundary_layers=boundary_layers,
        sources=sources,
        geometry=geometry,
    )

    flux_air = sim.add_flux(
        fcen,
        0,
        1,
        mp.FluxRegion(
            center=mp.Vector3(0.5 * L, 0, 0.5 * sz - dpml),
            size=mp.Vector3(L, 0, 0),
        ),
        mp.FluxRegion(
            center=mp.Vector3(L, 0, 0.5 * sz - dpml - 0.5 * dair),
            size=mp.Vector3(0, 0, dair),
        ),
    )

    sim.run(
        mp.dft_ldos(fcen, 0, 1),
        until_after_sources=mp.stop_when_fields_decayed(20, src_cmpt, src_pt, tol),
    )

    out_flux = mp.get_fluxes(flux_air)[0]
    if src_pt.x == 0:
        dV = np.pi / (resolution**3)
    else:
        dV = 2 * np.pi * src_pt.x / (resolution**2)
    total_flux = -np.real(sim.ldos_Fdata[0] * np.conj(sim.ldos_Jdata[0])) * dV
    ext_eff = out_flux / total_flux
    print(f"extraction efficiency (cyl):, " f"{dmat:.4f}, {h:.4f}, {ext_eff:.6f}")

    return ext_eff


def extraction_eff_3D(dmat: float, h: float) -> float:
    """Computes the extraction efficiency in 3D Cartesian coordinates.

    Args:
      dmat: thickness of dielectric layer.
      h: height of dipole above ground plane as fraction of dmat.

    Returns:
      The extraction efficiency of the dipole within the dielecric layer.
    """
    sxy = L + 2 * dpml
    sz = dmat + dair + dpml
    cell_size = mp.Vector3(sxy, sxy, sz)

    symmetries = [
        mp.Mirror(direction=mp.X, phase=-1),
        mp.Mirror(direction=mp.Y),
    ]

    boundary_layers = [
        mp.PML(dpml, direction=mp.X),
        mp.PML(dpml, direction=mp.Y),
        mp.PML(dpml, direction=mp.Z, side=mp.High),
    ]

    src_cmpt = mp.Ex
    src_pt = mp.Vector3(0, 0, -0.5 * sz + h * dmat)
    sources = [
        mp.Source(
            src=mp.GaussianSource(fcen, fwidth=0.1 * fcen),
            component=src_cmpt,
            center=src_pt,
        )
    ]

    geometry = [
        mp.Block(
            material=mp.Medium(index=n),
            center=mp.Vector3(0, 0, -0.5 * sz + 0.5 * dmat),
            size=mp.Vector3(mp.inf, mp.inf, dmat),
        )
    ]

    sim = mp.Simulation(
        resolution=resolution,
        cell_size=cell_size,
        boundary_layers=boundary_layers,
        sources=sources,
        geometry=geometry,
        symmetries=symmetries,
    )

    flux_air = sim.add_flux(
        fcen,
        0,
        1,
        mp.FluxRegion(
            center=mp.Vector3(0, 0, 0.5 * sz - dpml),
            size=mp.Vector3(L, L, 0),
        ),
        mp.FluxRegion(
            center=mp.Vector3(0.5 * L, 0, 0.5 * sz - dpml - 0.5 * dair),
            size=mp.Vector3(0, L, dair),
        ),
        mp.FluxRegion(
            center=mp.Vector3(-0.5 * L, 0, 0.5 * sz - dpml - 0.5 * dair),
            size=mp.Vector3(0, L, dair),
            weight=-1.0,
        ),
        mp.FluxRegion(
            center=mp.Vector3(0, 0.5 * L, 0.5 * sz - dpml - 0.5 * dair),
            size=mp.Vector3(L, 0, dair),
        ),
        mp.FluxRegion(
            center=mp.Vector3(0, -0.5 * L, 0.5 * sz - dpml - 0.5 * dair),
            size=mp.Vector3(L, 0, dair),
            weight=-1.0,
        ),
    )

    sim.run(
        mp.dft_ldos(fcen, 0, 1),
        until_after_sources=mp.stop_when_fields_decayed(20, src_cmpt, src_pt, tol),
    )

    out_flux = mp.get_fluxes(flux_air)[0]
    dV = 1 / (resolution**3)
    total_flux = -np.real(sim.ldos_Fdata[0] * np.conj(sim.ldos_Jdata[0])) * dV
    ext_eff = out_flux / total_flux
    print(f"extraction efficiency (3D):, {dmat:.4f}, {h:.4f}, {ext_eff:.6f}")

    return ext_eff


if __name__ == "__main__":
    layer_thickness = 0.7 * wvl / n
    dipole_height = np.linspace(0.1, 0.9, 21)

    exteff_cyl = np.zeros(len(dipole_height))
    exteff_3D = np.zeros(len(dipole_height))
    for j in range(len(dipole_height)):
        exteff_cyl[j] = extraction_eff_cyl(layer_thickness, dipole_height[j])
        exteff_3D[j] = extraction_eff_3D(layer_thickness, dipole_height[j])

    plt.plot(dipole_height, exteff_cyl, "bo-", label="cylindrical")
    plt.plot(dipole_height, exteff_3D, "ro-", label="3D Cartesian")
    plt.xlabel("height of dipole above ground plane (fraction of layer thickness)")
    plt.ylabel("extraction efficiency")
    plt.legend()

    if mp.am_master():
        plt.savefig("extraction_eff_vs_dipole_height.png", dpi=150, bbox_inches="tight")
```


![](../images/extraction_eff_vs_dipole_height.png#center)
