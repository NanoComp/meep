import math

import matplotlib.pyplot as plt
import numpy as np

import meep as mp

resolution = 25  # pixels/μm

dpml = 1.0  # PML thickness
dsub = 2.0  # substrate thickness
dpad = 2.0  # padding betweeen zone plate and PML
zh = 0.5  # zone-plate height
zN = 25  # number of zones (odd zones: π phase shift, even zones: none)
focal_length = 200  # focal length of zone plate
spot_length = 100  # far-field line length
ff_res = 10  # far-field resolution

pml_layers = [mp.PML(thickness=dpml)]

wvl_cen = 0.5
frq_cen = 1 / wvl_cen
dfrq = 0.2 * frq_cen

## radii of zones
## ref: eq. 7 of http://zoneplate.lbl.gov/theory
r = [
    math.sqrt(n * wvl_cen * (focal_length + n * wvl_cen / 4)) for n in range(1, zN + 1)
]

sr = r[-1] + dpad + dpml
sz = dpml + dsub + zh + dpad + dpml
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

glass = mp.Medium(index=1.5)

geometry = [
    mp.Block(
        material=glass,
        size=mp.Vector3(sr, 0, dpml + dsub),
        center=mp.Vector3(0.5 * sr, 0, -0.5 * sz + 0.5 * (dpml + dsub)),
    )
]

geometry.extend(
    mp.Block(
        material=glass if n % 2 == 0 else mp.vacuum,
        size=mp.Vector3(r[n], 0, zh),
        center=mp.Vector3(0.5 * r[n], 0, -0.5 * sz + dpml + dsub + 0.5 * zh),
    )
    for n in range(zN - 1, -1, -1)
)

sim = mp.Simulation(
    cell_size=cell_size,
    boundary_layers=pml_layers,
    resolution=resolution,
    sources=sources,
    geometry=geometry,
    dimensions=mp.CYLINDRICAL,
    m=-1,
)

## near-field monitor
n2f_obj = sim.add_near2far(
    frq_cen,
    0,
    1,
    mp.Near2FarRegion(
        center=mp.Vector3(0.5 * (sr - dpml), 0, 0.5 * sz - dpml),
        size=mp.Vector3(sr - dpml),
    ),
    mp.Near2FarRegion(
        center=mp.Vector3(sr - dpml, 0, 0.5 * sz - dpml - 0.5 * (dsub + zh + dpad)),
        size=mp.Vector3(z=dsub + zh + dpad),
    ),
)

sim.plot2D()
if mp.am_master():
    plt.savefig("zone_plate_epsilon.png", bbox_inches="tight", dpi=150)

sim.run(until_after_sources=100)

ff_r = sim.get_farfields(
    n2f_obj,
    ff_res,
    center=mp.Vector3(
        0.5 * (sr - dpml), 0, -0.5 * sz + dpml + dsub + zh + focal_length
    ),
    size=mp.Vector3(sr - dpml),
)

ff_z = sim.get_farfields(
    n2f_obj,
    ff_res,
    center=mp.Vector3(z=-0.5 * sz + dpml + dsub + zh + focal_length),
    size=mp.Vector3(z=spot_length),
)

E2_r = (
    np.absolute(ff_r["Ex"]) ** 2
    + np.absolute(ff_r["Ey"]) ** 2
    + np.absolute(ff_r["Ez"]) ** 2
)
E2_z = (
    np.absolute(ff_z["Ex"]) ** 2
    + np.absolute(ff_z["Ey"]) ** 2
    + np.absolute(ff_z["Ez"]) ** 2
)

if mp.am_master():
    plt.figure(dpi=200)
    plt.subplot(1, 2, 1)
    plt.semilogy(np.linspace(0, sr - dpml, len(E2_r)), E2_r, "bo-")
    plt.xlim(-2, 20)
    plt.xticks(list(np.arange(0, 25, 5)))
    plt.grid(True, axis="y", which="both", ls="-")
    plt.xlabel(r"$r$ coordinate (μm)")
    plt.ylabel(r"energy density of far fields, |E|$^2$")
    plt.subplot(1, 2, 2)
    plt.semilogy(
        np.linspace(
            focal_length - 0.5 * spot_length,
            focal_length + 0.5 * spot_length,
            len(E2_z),
        ),
        E2_z,
        "bo-",
    )
    plt.grid(True, axis="y", which="both", ls="-")
    plt.xlabel(r"$z$ coordinate (μm)")
    plt.ylabel(r"energy density of far fields, |E|$^2$")
    plt.suptitle(f"binary-phase zone plate with focal length $z$ = {focal_length} μm")

    plt.tight_layout()
    plt.savefig("zone_plate_farfields.png")
