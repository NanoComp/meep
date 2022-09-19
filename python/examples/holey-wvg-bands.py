# Meep Tutorial: Hz-polarized transmission and reflection through a cavity
# formed by a periodic sequence of holes in a dielectric waveguide,
# with a defect formed by a larger spacing between one pair of holes.

# This structure is based on one analyzed in:
#    S. Fan, J. N. Winn, A. Devenyi, J. C. Chen, R. D. Meade, and
#    J. D. Joannopoulos, "Guided and defect modes in periodic dielectric
#    waveguides," J. Opt. Soc. Am. B, 12 (7), 1267-1272 (1995).


import meep as mp


def main():
    # Some parameters to describe the geometry:
    eps = 13  # dielectric constant of waveguide
    w = 1.2  # width of waveguide
    r = 0.36  # radius of holes

    # The cell dimensions
    sy = 12  # size of cell in y direction (perpendicular to wvg.)
    dpml = 1  # PML thickness (y direction only!)

    cell = mp.Vector3(1, sy)

    b = mp.Block(size=mp.Vector3(mp.inf, w, mp.inf), material=mp.Medium(epsilon=eps))
    c = mp.Cylinder(radius=r)

    fcen = 0.25  # pulse center frequency
    df = 1.5  # pulse freq. width: large df = short impulse

    s = mp.Source(
        src=mp.GaussianSource(fcen, fwidth=df),
        component=mp.Hz,
        center=mp.Vector3(0.1234),
    )

    sym = mp.Mirror(direction=mp.Y, phase=-1)

    sim = mp.Simulation(
        cell_size=cell,
        geometry=[b, c],
        sources=[s],
        symmetries=[sym],
        boundary_layers=[mp.PML(dpml, direction=mp.Y)],
        resolution=20,
    )

    kx = False  # if true, do run at specified kx and get fields
    if kx:
        sim.k_point = mp.Vector3(kx)

        sim.run(
            mp.at_beginning(mp.output_epsilon),
            mp.after_sources(mp.Harminv(mp.Hz, mp.Vector3(0.1234), fcen, df)),
            until_after_sources=300,
        )

        sim.run(mp.at_every(1 / fcen / 20, mp.output_hfield_z), until=1 / fcen)

    else:
        k_interp = 19  # # k-points to interpolate, otherwise

        sim.run_k_points(300, mp.interpolate(k_interp, [mp.Vector3(), mp.Vector3(0.5)]))


if __name__ == "__main__":
    main()
