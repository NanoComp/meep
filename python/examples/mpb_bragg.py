import meep as mp
from meep import mpb

# Compute the bands at the X point for a quarter-wave stack Bragg
# mirror (this is the point that defines the band gap edges).

# the high and low indices:
n_lo = 1.0
n_hi = 3.0

w_hi = n_lo / (n_hi + n_lo)  # a quarter_wave stack

geometry_lattice = mp.Lattice(size=mp.Vector3(1))  # 1d cell
default_material = mp.Medium(index=n_lo)
geometry = mp.Cylinder(
    material=mp.Medium(index=n_hi),
    center=mp.Vector3(),
    axis=mp.Vector3(1),
    radius=mp.inf,
    height=w_hi,
)

kx = 0.5
k_points = [mp.Vector3(kx)]

resolution = 32
num_bands = 8

ms = mpb.ModeSolver(
    num_bands=num_bands,
    k_points=k_points,
    geometry_lattice=geometry_lattice,
    geometry=[geometry],
    resolution=resolution,
    default_material=default_material,
)


def main():
    ms.run_tm(
        mpb.output_hfield_y
    )  # note that TM and TE bands are degenerate, so we only need TM


if __name__ == "__main__":
    main()
