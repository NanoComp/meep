import matplotlib.pyplot as plt
import numpy as np

import meep as mp

resolution = 50  # pixels/μm

cell_size = mp.Vector3(14, 10, 0)

pml_layers = [mp.PML(thickness=2, direction=mp.X)]

# rotation angle (in degrees) of planewave, counter clockwise (CCW) around z-axis
rot_angle = np.radians(0)

fsrc = 1.0  # frequency of planewave (wavelength = 1/fsrc)

n = 1.5  # refractive index of homogeneous material
default_material = mp.Medium(index=n)

k_point = mp.Vector3(fsrc * n).rotate(mp.Vector3(z=1), rot_angle)

sources = [
    mp.EigenModeSource(
        src=mp.ContinuousSource(fsrc),
        center=mp.Vector3(),
        size=mp.Vector3(y=10),
        direction=mp.AUTOMATIC if rot_angle == 0 else mp.NO_DIRECTION,
        eig_kpoint=k_point,
        eig_band=1,
        eig_parity=mp.EVEN_Y + mp.ODD_Z if rot_angle == 0 else mp.ODD_Z,
        eig_match_freq=True,
    )
]

sim = mp.Simulation(
    cell_size=cell_size,
    resolution=resolution,
    boundary_layers=pml_layers,
    sources=sources,
    k_point=k_point,
    default_material=default_material,
    symmetries=[mp.Mirror(mp.Y)] if rot_angle == 0 else [],
)

sim.run(until=100)

nonpml_vol = mp.Volume(center=mp.Vector3(), size=mp.Vector3(10, 10, 0))

sim.plot2D(fields=mp.Ez, output_plane=nonpml_vol)

if mp.am_master():
    plt.axis("off")
    plt.savefig("pw.png", bbox_inches="tight")
