import unittest

import meep as mp


def dummy_eps(vec):
    return 1.0


class TestCylEllipsoid(unittest.TestCase):

    ref_Ez = -8.29555720049629e-5
    ref_Hz = -4.5623185899766e-5

    @classmethod
    def setUpClass(cls):
        cls.temp_dir = mp.make_output_directory()

    @classmethod
    def tearDownClass(cls):
        mp.delete_directory(cls.temp_dir)

    def init(self):

        c = mp.Cylinder(radius=3, material=mp.Medium(index=3.5))
        e = mp.Ellipsoid(size=mp.Vector3(1, 2, mp.inf))

        sources = mp.Source(
            src=mp.GaussianSource(1, fwidth=0.1),
            component=self.src_cmpt,
            center=mp.Vector3(),
        )

        if self.src_cmpt == mp.Ez:
            symmetries = [mp.Mirror(mp.X), mp.Mirror(mp.Y)]

        if self.src_cmpt == mp.Hz:
            symmetries = [mp.Mirror(mp.X, -1), mp.Mirror(mp.Y, -1)]

        self.sim = mp.Simulation(
            cell_size=mp.Vector3(10, 10),
            geometry=[c, e],
            boundary_layers=[mp.PML(1.0)],
            sources=[sources],
            symmetries=symmetries,
            resolution=100,
        )

        self.sim.use_output_directory(self.temp_dir)

        def print_stuff(sim_obj):
            v = mp.Vector3(4.13, 3.75, 0)
            p = self.sim.get_field_point(self.src_cmpt, v)
            print(f"t, Ez: {self.sim.round_time()} {p.real}+{p.imag}i")

        self.print_stuff = print_stuff

    def run_simulation(self):

        self.sim.run(
            mp.at_beginning(mp.output_epsilon),
            mp.at_every(0.25, self.print_stuff),
            mp.at_end(self.print_stuff),
            mp.at_end(mp.output_efield_z),
            until=23,
        )

        ref_out_field = self.ref_Ez if self.src_cmpt == mp.Ez else self.ref_Hz
        out_field = self.sim.fields.get_field(self.src_cmpt, mp.vec(4.13, 3.75)).real
        diff = abs(out_field - ref_out_field)

        self.assertTrue(abs(diff) <= 0.05 * abs(ref_out_field), "Field output differs")

    def test_ez_field(self):

        self.src_cmpt = mp.Ez
        self.init()
        self.run_simulation()

    def test_hz_field(self):

        self.src_cmpt = mp.Hz
        self.init()
        self.run_simulation()


if __name__ == "__main__":
    unittest.main()
