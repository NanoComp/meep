import argparse

import meep as mp

gdsII_file = "coupler.gds"
CELL_LAYER = 0
PORT1_LAYER = 1
PORT2_LAYER = 2
PORT3_LAYER = 3
PORT4_LAYER = 4
SOURCE_LAYER = 5
UPPER_BRANCH_LAYER = 31
LOWER_BRANCH_LAYER = 32
default_d = 0.3

t_oxide = 1.0
t_Si = 0.22
t_air = 0.78

dpml = 1
cell_thickness = dpml + t_oxide + t_Si + t_air + dpml

oxide = mp.Medium(epsilon=2.25)
silicon = mp.Medium(epsilon=12)

fcen = 1 / 1.55
df = 0.2 * fcen


def main(args):
    cell_zmax = 0.5 * cell_thickness if args.three_d else 0
    cell_zmin = -0.5 * cell_thickness if args.three_d else 0
    si_zmax = 0.5 * t_Si if args.three_d else 10
    si_zmin = -0.5 * t_Si if args.three_d else -10

    # read cell size, volumes for source region and flux monitors,
    # and coupler geometry from GDSII file
    upper_branch = mp.get_GDSII_prisms(
        silicon, gdsII_file, UPPER_BRANCH_LAYER, si_zmin, si_zmax
    )
    lower_branch = mp.get_GDSII_prisms(
        silicon, gdsII_file, LOWER_BRANCH_LAYER, si_zmin, si_zmax
    )

    cell = mp.GDSII_vol(gdsII_file, CELL_LAYER, cell_zmin, cell_zmax)
    p1 = mp.GDSII_vol(gdsII_file, PORT1_LAYER, si_zmin, si_zmax)
    p2 = mp.GDSII_vol(gdsII_file, PORT2_LAYER, si_zmin, si_zmax)
    p3 = mp.GDSII_vol(gdsII_file, PORT3_LAYER, si_zmin, si_zmax)
    p4 = mp.GDSII_vol(gdsII_file, PORT4_LAYER, si_zmin, si_zmax)
    src_vol = mp.GDSII_vol(gdsII_file, SOURCE_LAYER, si_zmin, si_zmax)

    # displace upper and lower branches of coupler (as well as source and flux regions)
    if args.d != default_d:
        delta_y = 0.5 * (args.d - default_d)
        delta = mp.Vector3(y=delta_y)
        p1.center += delta
        p2.center -= delta
        p3.center += delta
        p4.center -= delta
        src_vol.center += delta
        cell.size += 2 * delta
        for np in range(len(lower_branch)):
            lower_branch[np].center -= delta
            for nv in range(len(lower_branch[np].vertices)):
                lower_branch[np].vertices[nv] -= delta
        for np in range(len(upper_branch)):
            upper_branch[np].center += delta
            for nv in range(len(upper_branch[np].vertices)):
                upper_branch[np].vertices[nv] += delta

    geometry = upper_branch + lower_branch

    if args.three_d:
        oxide_center = mp.Vector3(z=-0.5 * t_oxide)
        oxide_size = mp.Vector3(cell.size.x, cell.size.y, t_oxide)
        oxide_layer = [mp.Block(material=oxide, center=oxide_center, size=oxide_size)]
        geometry = geometry + oxide_layer

    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            volume=src_vol,
            eig_parity=mp.NO_PARITY if args.three_d else mp.EVEN_Y + mp.ODD_Z,
        )
    ]

    sim = mp.Simulation(
        resolution=args.res,
        cell_size=cell.size,
        boundary_layers=[mp.PML(dpml)],
        sources=sources,
        geometry=geometry,
    )

    mode1 = sim.add_mode_monitor(fcen, 0, 1, mp.ModeRegion(volume=p1))
    mode2 = sim.add_mode_monitor(fcen, 0, 1, mp.ModeRegion(volume=p2))
    mode3 = sim.add_mode_monitor(fcen, 0, 1, mp.ModeRegion(volume=p3))
    mode4 = sim.add_mode_monitor(fcen, 0, 1, mp.ModeRegion(volume=p4))

    sim.run(until_after_sources=100)

    # S parameters
    p1_coeff = sim.get_eigenmode_coefficients(
        mode1, [1], eig_parity=mp.NO_PARITY if args.three_d else mp.EVEN_Y + mp.ODD_Z
    ).alpha[0, 0, 0]
    p2_coeff = sim.get_eigenmode_coefficients(
        mode2, [1], eig_parity=mp.NO_PARITY if args.three_d else mp.EVEN_Y + mp.ODD_Z
    ).alpha[0, 0, 1]
    p3_coeff = sim.get_eigenmode_coefficients(
        mode3, [1], eig_parity=mp.NO_PARITY if args.three_d else mp.EVEN_Y + mp.ODD_Z
    ).alpha[0, 0, 0]
    p4_coeff = sim.get_eigenmode_coefficients(
        mode4, [1], eig_parity=mp.NO_PARITY if args.three_d else mp.EVEN_Y + mp.ODD_Z
    ).alpha[0, 0, 0]

    # transmittance
    p2_trans = abs(p2_coeff) ** 2 / abs(p1_coeff) ** 2
    p3_trans = abs(p3_coeff) ** 2 / abs(p1_coeff) ** 2
    p4_trans = abs(p4_coeff) ** 2 / abs(p1_coeff) ** 2

    print(
        "trans:, {:.2f}, {:.6f}, {:.6f}, {:.6f}".format(
            args.d, p2_trans, p3_trans, p4_trans
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-res", type=int, default=50, help="resolution (default: 50 pixels/um)"
    )
    parser.add_argument(
        "-d", type=float, default=0.1, help="branch separation (default: 0.1 um)"
    )
    parser.add_argument(
        "--three_d",
        action="store_true",
        default=False,
        help="3d calculation? (default: False)",
    )
    args = parser.parse_args()
    main(args)
