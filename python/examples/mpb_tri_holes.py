import math

import meep as mp
from meep import mpb

# 2d system: triangular lattice of air holes in dielectric
# This structure has a complete band gap (i.e. a gap in both TE and TM
# simultaneously) for a hole radius of 0.45a and a dielectric constant of
# 12.   (See, e.g., the book "Photonic Crystals" by Joannopoulos et al.)

# first, define the lattice vectors and k-points for a triangular lattice:

geometry_lattice = mp.Lattice(
    size=mp.Vector3(1, 1),
    basis1=mp.Vector3(math.sqrt(3) / 2, 0.5),
    basis2=mp.Vector3(math.sqrt(3) / 2, -0.5),
)

kz = 0  # use non-zero kz to consider vertical propagation

k_points = [
    mp.Vector3(z=kz),  # Gamma
    mp.Vector3(0, 0.5, kz),  # M
    mp.Vector3(1 / -3, 1 / 3, kz),  # K
    mp.Vector3(z=kz),  # Gamma
]

k_interp = 4
k_points = mp.interpolate(k_interp, k_points)

# Now, define the geometry, etcetera:

eps = 12  # the dielectric constant of the background
r = 0.45  # the hole radius

default_material = mp.Medium(epsilon=eps)
geometry = [mp.Cylinder(r, material=mp.air)]

resolution = 32
num_bands = 8

ms = mpb.ModeSolver(
    geometry_lattice=geometry_lattice,
    geometry=geometry,
    k_points=k_points,
    default_material=default_material,
    resolution=resolution,
    num_bands=num_bands,
)


def main():
    if kz == 0:
        ms.run_te()
        ms.run_tm()
    else:
        ms.run()  # if kz != 0 there are no purely te and tm bands


if __name__ == "__main__":
    main()
