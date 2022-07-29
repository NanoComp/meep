import argparse

from meep.materials import Al

import meep as mp


def main(args):

    resolution = 40
    cell_size = mp.Vector3(z=10)

    boundary_layers = [
        mp.PML(1, direction=mp.Z) if args.pml else mp.Absorber(1, direction=mp.Z)
    ]

    sources = [
        mp.Source(
            src=mp.GaussianSource(1 / 0.803, fwidth=0.1),
            center=mp.Vector3(),
            component=mp.Ex,
        )
    ]

    def print_stuff(sim):
        p = sim.get_field_point(mp.Ex, mp.Vector3())
        print(f"ex:, {sim.meep_time()}, {p.real}")

    sim = mp.Simulation(
        cell_size=cell_size,
        resolution=resolution,
        dimensions=1,
        default_material=Al,
        boundary_layers=boundary_layers,
        sources=sources,
    )

    sim.run(
        mp.at_every(10, print_stuff),
        until_after_sources=mp.stop_when_fields_decayed(50, mp.Ex, mp.Vector3(), 1e-6),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-pml", action="store_true", default=False, help="Use PML as boundary layer"
    )
    args = parser.parse_args()
    main(args)
