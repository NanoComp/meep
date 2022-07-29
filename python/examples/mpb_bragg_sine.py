import math

import meep as mp
from meep import mpb

# Compute the band structure for a Bragg mirror consisting of a
# sinusoidally-varying dielectric index.

# The index will vary sinusoidally between index-min and index-max:
index_min = 1
index_max = 3


# Define a function of position p (in the lattice basis) that returns
# the material at that position.  In this case, we use the function:
#        index-min + 0.5 * (index-max - index-min)
#                        * (1 + cos(2*pi*x))
# This is periodic, and also has inversion symmetry.
def eps_func(p):
    return mp.Medium(
        index=index_min
        + 0.5 * (index_max - index_min) * (1 + math.cos(2 * math.pi * p.x))
    )


geometry_lattice = mp.Lattice(size=mp.Vector3(1))  # 1d cell

# We'll just make it the default material, so that it goes everywhere.
default_material = eps_func

k_points = mp.interpolate(9, [mp.Vector3(), mp.Vector3(x=0.5)])

resolution = 32
num_bands = 8

ms = mpb.ModeSolver(
    num_bands=num_bands,
    k_points=k_points,
    geometry_lattice=geometry_lattice,
    resolution=resolution,
    default_material=default_material,
)


def main():
    # the TM and TE bands are degenerate, so we only need TM:
    ms.run_tm()


if __name__ == "__main__":
    main()
