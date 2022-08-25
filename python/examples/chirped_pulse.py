## linear-chirped pulse planewave with higher frequencies at the front (down-chirp)
import numpy as np

import meep as mp

resolution = 40

dpml = 2
pml_layers = [mp.PML(thickness=dpml, direction=mp.X)]

sx = 40
sy = 6
cell_size = mp.Vector3(sx + 2 * dpml, sy)

v0 = 1.0  # pulse center frequency
a = 0.2  # Gaussian envelope half-width
b = -0.5  # linear chirp rate (positive: up-chirp, negative: down-chirp)
t0 = 15  # peak time

chirp = lambda t: np.exp(1j * 2 * np.pi * v0 * (t - t0)) * np.exp(
    -a * (t - t0) ** 2 + 1j * b * (t - t0) ** 2
)

sources = [
    mp.Source(
        src=mp.CustomSource(src_func=chirp),
        center=mp.Vector3(-0.5 * sx),
        size=mp.Vector3(y=sy),
        component=mp.Ez,
    )
]

sim = mp.Simulation(
    cell_size=cell_size,
    boundary_layers=pml_layers,
    resolution=resolution,
    k_point=mp.Vector3(),
    sources=sources,
    symmetries=[mp.Mirror(mp.Y)],
)

sim.run(
    mp.in_volume(
        mp.Volume(center=mp.Vector3(), size=mp.Vector3(sx, sy)),
        mp.at_every(2.7, mp.output_efield_z),
    ),
    until=t0 + 50,
)
