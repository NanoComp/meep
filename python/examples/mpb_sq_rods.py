import time

import meep as mp
from meep import mpb

# Compute band structure for a square lattice of dielectric rods
# in air.

# Define various parameters with define_param so that they are
# settable from the command_line (with mpb <param>=<value>):
r = 0.2  # radius of the rods
eps = 11.56  # dielectric constant
k_interp = 4  # number of k points to interpolate

GaAs = mp.Medium(epsilon=eps)

geometry_lattice = mp.Lattice(size=mp.Vector3(1, 1))  # 2d cell

geometry = [mp.Cylinder(r, material=GaAs)]

Gamma = mp.Vector3()
X = mp.Vector3(0.5, 0)
M = mp.Vector3(0.5, 0.5)
k_points = mp.interpolate(k_interp, [Gamma, X, M, Gamma])

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
    # Compute the TE and TM bands and report the total elapsed time:
    t0 = time.time()
    ms.run_te()
    ms.run_tm()
    print(f"total time for both TE and TM bands: {time.time() - t0:.2f} seconds")

    ms.display_eigensolver_stats()


if __name__ == "__main__":
    main()
