import matplotlib.pyplot as plt
import numpy as np

import meep as mp

r = 0.7  # radius of cylinder
h = 2.3  # height of cylinder

wvl_min = 2 * np.pi * r / 10
wvl_max = 2 * np.pi * r / 2

frq_min = 1 / wvl_max
frq_max = 1 / wvl_min
frq_cen = 0.5 * (frq_min + frq_max)
dfrq = frq_max - frq_min
nfrq = 100

## at least 8 pixels per smallest wavelength, i.e. np.floor(8/wvl_min)
resolution = 25

dpml = 0.5 * wvl_max
dair = 1.0 * wvl_max

pml_layers = [mp.PML(thickness=dpml)]

sr = r + dair + dpml
sz = dpml + dair + h + dair + dpml
cell_size = mp.Vector3(sr, 0, sz)

sources = [
    mp.Source(
        mp.GaussianSource(frq_cen, fwidth=dfrq, is_integrated=True),
        component=mp.Er,
        center=mp.Vector3(0.5 * sr, 0, -0.5 * sz + dpml),
        size=mp.Vector3(sr),
    ),
    mp.Source(
        mp.GaussianSource(frq_cen, fwidth=dfrq, is_integrated=True),
        component=mp.Ep,
        center=mp.Vector3(0.5 * sr, 0, -0.5 * sz + dpml),
        size=mp.Vector3(sr),
        amplitude=-1j,
    ),
]

sim = mp.Simulation(
    cell_size=cell_size,
    boundary_layers=pml_layers,
    resolution=resolution,
    sources=sources,
    dimensions=mp.CYLINDRICAL,
    m=-1,
)

box_z1 = sim.add_flux(
    frq_cen,
    dfrq,
    nfrq,
    mp.FluxRegion(center=mp.Vector3(0.5 * r, 0, -0.5 * h), size=mp.Vector3(r)),
)
box_z2 = sim.add_flux(
    frq_cen,
    dfrq,
    nfrq,
    mp.FluxRegion(center=mp.Vector3(0.5 * r, 0, +0.5 * h), size=mp.Vector3(r)),
)
box_r = sim.add_flux(
    frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(r), size=mp.Vector3(z=h))
)

sim.run(until_after_sources=10)

freqs = mp.get_flux_freqs(box_z1)
box_z1_data = sim.get_flux_data(box_z1)
box_z2_data = sim.get_flux_data(box_z2)
box_r_data = sim.get_flux_data(box_r)

box_z1_flux0 = mp.get_fluxes(box_z1)

sim.reset_meep()

n_cyl = 2.0
geometry = [
    mp.Block(
        material=mp.Medium(index=n_cyl),
        center=mp.Vector3(0.5 * r),
        size=mp.Vector3(r, 0, h),
    )
]

sim = mp.Simulation(
    cell_size=cell_size,
    geometry=geometry,
    boundary_layers=pml_layers,
    resolution=resolution,
    sources=sources,
    dimensions=mp.CYLINDRICAL,
    m=-1,
)

box_z1 = sim.add_flux(
    frq_cen,
    dfrq,
    nfrq,
    mp.FluxRegion(center=mp.Vector3(0.5 * r, 0, -0.5 * h), size=mp.Vector3(r)),
)
box_z2 = sim.add_flux(
    frq_cen,
    dfrq,
    nfrq,
    mp.FluxRegion(center=mp.Vector3(0.5 * r, 0, +0.5 * h), size=mp.Vector3(r)),
)
box_r = sim.add_flux(
    frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(r), size=mp.Vector3(z=h))
)

sim.load_minus_flux_data(box_z1, box_z1_data)
sim.load_minus_flux_data(box_z2, box_z2_data)
sim.load_minus_flux_data(box_r, box_r_data)

sim.run(until_after_sources=100)

box_z1_flux = mp.get_fluxes(box_z1)
box_z2_flux = mp.get_fluxes(box_z2)
box_r_flux = mp.get_fluxes(box_r)

scatt_flux = np.asarray(box_z1_flux) - np.asarray(box_z2_flux) - np.asarray(box_r_flux)
intensity = np.asarray(box_z1_flux0) / (np.pi * r**2)
scatt_cross_section = np.divide(-scatt_flux, intensity)

if mp.am_master():
    plt.figure(dpi=150)
    plt.loglog(2 * np.pi * r * np.asarray(freqs), scatt_cross_section, "bo-")
    plt.grid(True, which="both", ls="-")
    plt.xlabel("(cylinder circumference)/wavelength, 2πr/λ")
    plt.ylabel("scattering cross section, σ")
    plt.title("Scattering Cross Section of a Lossless Dielectric Cylinder")
    plt.tight_layout()
    plt.savefig("cylinder_cross_section.png")
