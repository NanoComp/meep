---
# 2d Cell with Out-of-Plane Wavevector
---

A 2d cell with a `k_point` that has a *non-zero* $z$ component $\beta$ (e.g., planewaves incident on a 2d structure from an out-of-plane direction, propagating modes of fiber waveguides with 2d cladding cross section, etc.) could be accomplished by a "1-pixel-thick" 3d simulation with complex fields (this is done by the `kz_2d="3d"` option). However, Meep will by default model the $e^{i \beta z}$ dependence using a modified 2d simulation with complex fields (this is done by the `kz_2d="complex"` option), which improves performance with no loss in accuracy (as demonstrated below). As a further optimization, Meep can also model this problem with "real" fields via `kz_2d="real/imag"`, but the fields must be interpreted in a special way.

Mathematically, an $e^{i \beta z}$ dependence ($k_z = \beta$) of the fields can be treated by including an $i\beta\hat{z} \times {}$ cross-product in the curls of Maxwell's equations, which couples the $E_z$ and $H_z$ polarizations.   Since this term is complex ($\sim i \beta$), the electromagnetic fields are complex.

However, an additional trick is possible.  Since the $i\beta\hat{z} \times {}$ complex term only couples $E_z$ and $H_z$ polarizations to one another, then we can choose the $H_z$ polarization ($E_x, E_y, H_z$ fields) to be purely real while the $E_z$ polarization ($H_x, H_y, E_z$ fields) is purely *imaginary*.   The `kz_2d="real/imag"` option does precisely this, but stores the purely imaginary $E_z$ polarization as the "real" parts of the fields, both internally and in the output.  So, if you use the `kz_2d="real/imag"` and you output both $E_x$ and $E_z$, the output will be real numbers, but you should multiply $E_z$ by $i$ to get the actual fields.    This requires some care, which is why this option is not the default.  The good news is that calculations of flux, energy, forces and similar real quantities are insensitive to this implicit $i$ factor (the $i$ cancels), so the built-in calculations do the right thing with no modification.

As a demonstration of this feature, consider the example of computing the reflectance of a planar air/dielectric interface for a planewave incident from the out-of-plane direction. This is a slightly modified version of [Tutorial/Basics/Angular Reflectance Spectrum of a Planar Interface](Python_Tutorials/Basics.md#angular-reflectance-spectrum-of-a-planar-interface). The interface normal is along the $x$ axis. An $E_y$-polarized planewave is incident at an angle of 19.4° counter clockwise around the $y$-axis with 0° along the $+x$ direction (i.e., the plane of incidence is $xz$ which corresponds to the $\mathcal{S}$-polarization). The 2d cell is in the $xy$ plane (with a $y$ dimension of a single pixel). The reflectance at a wavelength of 1 μm is computed using separate simulations for `kz_2d` of `"real/imag"`, `"complex"`, and `"3d"`. The three results are compared with the analytic Fresnel equation.

The simulation script is [examples/refl-angular-kz2d.py](https://github.com/NanoComp/meep/blob/master/python/examples/refl-angular-kz2d.py). The Scheme version is [examples/refl-angular-kz2d.ctl](https://github.com/NanoComp/meep/blob/master/scheme/examples/refl-angular-kz2d.ctl).

```py
import meep as mp
import math

def refl_planar(theta, kz_2d):
    resolution = 100

    dpml = 1.0
    sx = 10
    sx = 10 + 2*dpml
    cell_size = mp.Vector3(sx)
    pml_layers = [mp.PML(dpml)]

    fcen = 1.0

    # plane of incidence is XZ
    k = mp.Vector3(z=math.sin(theta)).scale(fcen)

    sources = [mp.Source(mp.GaussianSource(fcen,fwidth=0.2*fcen),
                         component=mp.Ey,
                         center=mp.Vector3(-0.5*sx+dpml))]

    sim = mp.Simulation(cell_size=cell_size,
                        boundary_layers=pml_layers,
                        sources=sources,
                        k_point=k,
                        kz_2d=kz_2d,
                        resolution=resolution)

    refl_fr = mp.FluxRegion(center=mp.Vector3(-0.25*sx))
    refl = sim.add_flux(fcen, 0, 1, refl_fr)

    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ey, mp.Vector3(-0.5*sx+dpml), 1e-9))

    input_flux = mp.get_fluxes(refl)
    input_data = sim.get_flux_data(refl)
    sim.reset_meep()

    # add a block with n=3.5 for the air-dielectric interface
    geometry = [mp.Block(size=mp.Vector3(0.5*sx,mp.inf,mp.inf),
                         center=mp.Vector3(0.25*sx),
                         material=mp.Medium(index=3.5))]

    sim = mp.Simulation(cell_size=cell_size,
                        geometry=geometry,
                        boundary_layers=pml_layers,
                        sources=sources,
                        k_point=k,
                        kz_2d=kz_2d,
                        resolution=resolution)

    refl = sim.add_flux(fcen, 0, 1, refl_fr)
    sim.load_minus_flux_data(refl, input_data)

    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ey, mp.Vector3(-0.5*sx+dpml), 1e-9))

    refl_flux = mp.get_fluxes(refl)
    freqs = mp.get_flux_freqs(refl)

    Rmeep = -refl_flux[0]/input_flux[0]
    return Rmeep


# rotation angle of source: CCW around Y axis, 0 degrees along +X axis
theta_r = math.radians(19.4)

Rmeep_real_imag = refl_planar(theta_r,"real/imag")
Rmeep_complex = refl_planar(theta_r,"complex")
Rmeep_3d = refl_planar(theta_r,"3d")

n1=1
n2=3.5

# compute angle of refracted planewave in medium n2
# for incident planewave in medium n1 at angle theta_in
theta_out = lambda theta_in: math.asin(n1*math.sin(theta_in)/n2)

# compute Fresnel reflectance for S-polarization in medium n2
# for incident planewave in medium n1 at angle theta_in
Rfresnel = lambda theta_in: math.fabs((n2*math.cos(theta_out(theta_in))-n1*math.cos(theta_in))/(n2*math.cos(theta_out(theta_in))+n1*math.cos(theta_in)))**2

print("refl:, {} (real/imag), {} (complex), {} (3d), {} (analytic)".format(Rmeep_real_imag,Rmeep_complex,Rmeep_3d,Rfresnel(theta_r)))
```

The Meep results are identical to within six decimal digits and agree well with the analytic theory with an error of less than 0.7%.

```
refl:, 0.3272338236967464 (real/imag), 0.3272338244372344 (complex), 0.3272330216564413 (3d), 0.3293821216165117 (analytic)
```