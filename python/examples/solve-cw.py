"""
Verifies that the relative error in the fields of a resonant mode
of a 2d ring resonator is monotonically decreasing with decreasing
tolerance of the CW solver. Also visualizes the fields of the resonant
mode in the time and frequency domains.
"""

import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
import meep as mp

resolution = 20  # pixels/μm
n = 3.4  # refractive index of ring
w = 1  # width of ring
r = 1  # inner radius of ring
pad = 4  # padding between outer ring and PML
dpml = 2  # PML thickness

sxy = 2 * (r + w + pad + dpml)
cell_size = mp.Vector3(sxy, sxy)

pml_layers = [mp.PML(dpml)]

nonpml_vol = mp.Volume(
    center=mp.Vector3(),
    size=mp.Vector3(sxy - 2 * dpml, sxy - 2 * dpml),
)

geometry = [
    mp.Cylinder(radius=r + w, material=mp.Medium(index=n)),
    mp.Cylinder(radius=r),
]

fcen = 0.118  # frequency of resonant mode

src = [
    mp.Source(
        mp.ContinuousSource(fcen),
        component=mp.Ez,
        center=mp.Vector3(r + 0.1),
    ),
    mp.Source(
        mp.ContinuousSource(fcen),
        component=mp.Ez,
        center=mp.Vector3(-(r + 0.1)),
        amplitude=-1,
    ),
]

symmetries = [
    mp.Mirror(mp.X, phase=-1),
    mp.Mirror(mp.Y, phase=+1),
]

sim = mp.Simulation(
    resolution=resolution,
    cell_size=cell_size,
    geometry=geometry,
    sources=src,
    force_complex_fields=True,
    symmetries=symmetries,
    boundary_layers=pml_layers,
)

# CW solver convergence properties
maxiters = 10000
L = 10
num_tols = 5
tols = np.logspace(-8, -8.0 - num_tols + 1, num_tols)

ez_dat = np.zeros(
    (
        int(nonpml_vol.size.x * resolution) + 2,
        int(nonpml_vol.size.y * resolution) + 2,
        num_tols,
    ),
    dtype=np.complex_,
)

for i in range(num_tols):
    sim.init_sim()
    sim.solve_cw(tols[i], maxiters, L)
    ez_dat[:, :, i] = sim.get_array(vol=nonpml_vol, component=mp.Ez)

err_dat = np.zeros(num_tols - 1)
for i in range(num_tols - 1):
    err_dat[i] = np.linalg.norm(ez_dat[:, :, i] - ez_dat[:, :, -1]) / np.linalg.norm(
        ez_dat[:, :, -1]
    )
    print(f"err:, {tols[i]}, {err_dat[i]}")

plt.figure(dpi=150)
plt.loglog(tols[: num_tols - 1], err_dat, "bo-")
plt.xlabel("frequency-domain solver tolerance")
plt.ylabel("relative error in fields of resonant mode")
plt.title("2d ring resonator")
plt.savefig("ring_err.png", dpi=150, bbox_inches="tight")

eps_data = sim.get_array(vol=nonpml_vol, component=mp.Dielectric)
ez_data = np.real(ez_dat[:, :, num_tols - 1])

plt.figure()
plt.imshow(
    eps_data.transpose(),
    interpolation="spline36",
    cmap="binary",
)
plt.imshow(
    ez_data.transpose(),
    interpolation="spline36",
    cmap="RdBu",
    alpha=0.9,
)
plt.title("time-domain fields ($E_z$)")
plt.axis("off")
plt.savefig("ring_ez.png", dpi=150, bbox_inches="tight")

if np.all(np.diff(err_dat) < 0):
    print(
        "PASSED solve_cw test: error in the fields is "
        "decreasing with increasing resolution."
    )
else:
    print(
        "FAILED solve_cw test: error in the fields is "
        "NOT decreasing with increasing resolution."
    )

sim.reset_meep()

df = 0.08  # frequency width of pulsed source
src = [
    mp.Source(
        mp.GaussianSource(fcen, fwidth=df),
        component=mp.Ez,
        center=mp.Vector3(r + 0.1),
    ),
    mp.Source(
        mp.GaussianSource(fcen, fwidth=df),
        component=mp.Ez,
        center=mp.Vector3(-(r + 0.1)),
        amplitude=-1,
    ),
]

sim = mp.Simulation(
    resolution=resolution,
    cell_size=mp.Vector3(sxy, sxy),
    geometry=geometry,
    sources=src,
    symmetries=symmetries,
    boundary_layers=pml_layers,
)

dft_obj = sim.add_dft_fields([mp.Ez], fcen, 0, 1, where=nonpml_vol)

sim.run(
    until_after_sources=mp.stop_when_fields_decayed(
        50,
        mp.Ez,
        mp.Vector3(r + 0.1523),
        1e-8,
    )
)

eps_data = sim.get_array(vol=nonpml_vol, component=mp.Dielectric)
ez_data = np.real(sim.get_dft_array(dft_obj, mp.Ez, 0))

plt.figure()
plt.imshow(
    eps_data.transpose(),
    interpolation="spline36",
    cmap="binary",
)
plt.imshow(
    ez_data.transpose(),
    interpolation="spline36",
    cmap="RdBu",
    alpha=0.9,
)
plt.title("DFT fields ($E_z$)")
plt.axis("off")
plt.savefig("ring_ez_dft.png", dpi=150, bbox_inches="tight")
