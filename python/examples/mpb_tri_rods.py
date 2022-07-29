import math

import meep as mp
from meep import mpb

# A triangular lattice of dielectric rods in air.  (This structure has
# a band_gap for TM fields.)  This file is used in the "Data Analysis
# Tutorial" section of the MPB manual.

num_bands = 8

geometry_lattice = mp.Lattice(
    size=mp.Vector3(1, 1),
    basis1=mp.Vector3(math.sqrt(3) / 2, 0.5),
    basis2=mp.Vector3(math.sqrt(3) / 2, -0.5),
)

geometry = [mp.Cylinder(0.2, material=mp.Medium(epsilon=12))]

k_points = [
    mp.Vector3(),  # Gamma
    mp.Vector3(y=0.5),  # M
    mp.Vector3(1 / -3, 1 / 3),  # K
    mp.Vector3(),  # Gamma
]

k_points = mp.interpolate(4, k_points)

resolution = 32

ms = mpb.ModeSolver(
    geometry=geometry,
    geometry_lattice=geometry_lattice,
    k_points=k_points,
    resolution=resolution,
    num_bands=num_bands,
)


def main():
    ms.run_tm(
        mpb.output_at_kpoint(
            mp.Vector3(1 / -3, 1 / 3), mpb.fix_efield_phase, mpb.output_efield_z
        )
    )
    ms.run_te()


if __name__ == "__main__":
    main()
