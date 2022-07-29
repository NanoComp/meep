import meep as mp


def main():

    c = mp.Cylinder(radius=3, material=mp.Medium(index=3.5))
    e = mp.Ellipsoid(size=mp.Vector3(1, 2, mp.inf))

    src_cmpt = mp.Hz
    sources = mp.Source(
        src=mp.GaussianSource(1, fwidth=0.1), component=src_cmpt, center=mp.Vector3()
    )

    if src_cmpt == mp.Ez:
        symmetries = [mp.Mirror(mp.X), mp.Mirror(mp.Y)]

    if src_cmpt == mp.Hz:
        symmetries = [mp.Mirror(mp.X, -1), mp.Mirror(mp.Y, -1)]

    sim = mp.Simulation(
        cell_size=mp.Vector3(10, 10),
        geometry=[c, e],
        boundary_layers=[mp.PML(1.0)],
        sources=[sources],
        symmetries=symmetries,
        resolution=100,
    )

    def print_stuff(sim_obj):
        v = mp.Vector3(4.13, 3.75, 0)
        p = sim.get_field_point(src_cmpt, v)
        print(f"t, Ez: {sim.round_time()} {p.real}+{p.imag}i")

    sim.run(
        mp.at_beginning(mp.output_epsilon),
        mp.at_every(0.25, print_stuff),
        mp.at_end(print_stuff),
        mp.at_end(mp.output_efield_z),
        until=23,
    )

    print(f"stopped at meep time = {sim.round_time()}")


if __name__ == "__main__":
    main()
