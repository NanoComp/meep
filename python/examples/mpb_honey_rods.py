import math

import meep as mp
from meep import mpb

# A honeycomb lattice of dielectric rods in air.  (This structure has
# a complete (overlapping TE/TM) band gap.)  A honeycomb lattice is really
# just a triangular lattice with two rods per unit cell, so we just
# take the lattice, k_points, etcetera from mpb_tri_rods.py.

r = 0.14  # the rod radius
eps = 12  # the rod dielectric constant

# triangular lattice:
geometry_lattice = mp.Lattice(
    size=mp.Vector3(1, 1),
    basis1=mp.Vector3(math.sqrt(3) / 2, 0.5),
    basis2=mp.Vector3(math.sqrt(3) / 2, -0.5),
)

# Two rods per unit cell, at the correct positions to form a honeycomb
# lattice, and arranged to have inversion symmetry:
geometry = [
    mp.Cylinder(
        r,
        center=mp.Vector3(1 / 6, 1 / 6),
        height=mp.inf,
        material=mp.Medium(epsilon=eps),
    ),
    mp.Cylinder(
        r,
        center=mp.Vector3(1 / -6, 1 / -6),
        height=mp.inf,
        material=mp.Medium(epsilon=eps),
    ),
]

# The k_points list, for the Brillouin zone of a triangular lattice:
k_points = [
    mp.Vector3(),  # Gamma
    mp.Vector3(y=0.5),  # M
    mp.Vector3(1 / -3, 1 / 3),  # K
    mp.Vector3(),  # Gamma
]

k_interp = 4  # number of k_points to interpolate
k_points = mp.interpolate(k_interp, k_points)

resolution = 32
num_bands = 8

ms = mpb.ModeSolver(
    geometry_lattice=geometry_lattice,
    geometry=geometry,
    k_points=k_points,
    resolution=resolution,
    num_bands=num_bands,
)


def main():
    ms.run_tm()
    ms.run_te()


# Since there is a complete gap, we could instead see it just by using:
# run()
# The gap is between bands 12 and 13 in this case.  (Note that there is
# a false gap between bands 2 and 3, which disappears as you increase the
# k_point resolution.)

if __name__ == "__main__":
    main()
