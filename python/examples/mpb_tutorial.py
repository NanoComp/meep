import math

from scipy.optimize import minimize_scalar, ridder

import meep as mp
from meep import mpb


def print_heading(h):
    stars = "*" * 10
    print("{0} {1} {0}".format(stars, h))


# Our First Band Structure

print_heading("Square lattice of rods in air")

num_bands = 8
k_points = [
    mp.Vector3(),  # Gamma
    mp.Vector3(0.5),  # X
    mp.Vector3(0.5, 0.5),  # M
    mp.Vector3(),
]  # Gamma

k_points = mp.interpolate(4, k_points)
geometry = [mp.Cylinder(0.2, material=mp.Medium(epsilon=12))]
geometry_lattice = mp.Lattice(size=mp.Vector3(1, 1))
resolution = 32

ms = mpb.ModeSolver(
    num_bands=num_bands,
    k_points=k_points,
    geometry=geometry,
    geometry_lattice=geometry_lattice,
    resolution=resolution,
)

print_heading("Square lattice of rods: TE bands")
ms.run_te()

print_heading("Square lattice of rods: TM bands")
ms.run_tm()

print_heading("Square lattice of rods: TM, w/efield")
ms.run_tm(mpb.output_efield_z)

print_heading("Square lattice of rods: TE, w/hfield & dpwr")
ms.run_te(mpb.output_at_kpoint(mp.Vector3(0.5), mpb.output_hfield_z, mpb.output_dpwr))

# Bands of a Triangular Lattice

print_heading("Triangular lattice of rods in air")

ms.geometry_lattice = mp.Lattice(
    size=mp.Vector3(1, 1),
    basis1=mp.Vector3(math.sqrt(3) / 2, 0.5),
    basis2=mp.Vector3(math.sqrt(3) / 2, -0.5),
)

ms.k_points = [
    mp.Vector3(),  # Gamma
    mp.Vector3(y=0.5),  # M
    mp.Vector3(-1 / 3, 1 / 3),  # K
    mp.Vector3(),
]  # Gamma

ms.k_points = mp.interpolate(4, k_points)
ms.run_tm()

# Maximizing the First TM Gap

print_heading("Maximizing the first TM gap")


def first_tm_gap(r):
    ms.geometry = [mp.Cylinder(r, material=mp.Medium(epsilon=12))]
    ms.run_tm()
    return -1 * ms.retrieve_gap(1)


ms.num_bands = 2
ms.mesh_size = 7

result = minimize_scalar(
    first_tm_gap, method="bounded", bounds=[0.1, 0.5], options={"xatol": 0.1}
)
print(f"radius at maximum: {result.x}")
print(f"gap size at maximum: {result.fun * -1}")

ms.mesh_size = 3  # Reset to default value of 3

# A Complete 2D Gap with an Anisotropic Dielectric

print_heading("Anisotropic complete 2d gap")

ms.geometry = [mp.Cylinder(0.3, material=mp.Medium(epsilon_diag=mp.Vector3(1, 1, 12)))]

ms.default_material = mp.Medium(epsilon_diag=mp.Vector3(12, 12, 1))
ms.num_bands = 8
ms.run()  # just use run, instead of run_te or run_tm, to find the complete gap

# Finding a Point-defect State

print_heading("5x5 point defect")

ms.geometry_lattice = mp.Lattice(size=mp.Vector3(5, 5))
ms.geometry = [mp.Cylinder(0.2, material=mp.Medium(epsilon=12))]

ms.geometry = mp.geometric_objects_lattice_duplicates(ms.geometry_lattice, ms.geometry)
ms.geometry.append(mp.Cylinder(0.2, material=mp.air))

ms.resolution = 16
ms.k_points = [mp.Vector3(0.5, 0.5)]

ms.num_bands = 50
ms.run_tm()

mpb.output_efield_z(ms, 25)

ms.get_dfield(25)  # compute the D field for band 25
ms.compute_field_energy()  # compute the energy density from D
c = mp.Cylinder(1.0, material=mp.air)
print(f"energy in cylinder: {ms.compute_energy_in_objects([c])}")

print_heading("5x5 point defect, targeted solver")

ms.num_bands = 1  # only need to compute a single band, now!
ms.target_freq = (0.2812 + 0.4174) / 2
ms.tolerance = 1e-8
ms.run_tm()

# Tuning the Point-defect Mode

print_heading("Tuning the 5x5 point defect")

old_geometry = ms.geometry  # save the 5x5 grid with a missing rod


def rootfun(eps):
    # add the cylinder of epsilon = eps to the old geometry:
    ms.geometry = old_geometry + [mp.Cylinder(0.2, material=mp.Medium(epsilon=eps))]
    ms.run_tm()  # solve for the mode (using the targeted solver)
    print(f"epsilon = {eps} gives freq. =  {ms.get_freqs()[0]}")
    return ms.get_freqs()[0] - 0.314159  # return 1st band freq. - 0.314159


rooteps = ridder(rootfun, 1, 12)
print(f"root (value of epsilon) is at: {rooteps}")

rootval = rootfun(rooteps)
print(f"root function at {rooteps} = {rootval}")
