import unittest

import meep as mp


class TestForce(unittest.TestCase):
    def setUp(self):
        self.resolution = 20
        self.cell = mp.Vector3(10, 10)
        self.pml_layers = mp.PML(1.0)

        self.fcen = 1.0
        self.df = 1.0
        self.sources = mp.Source(
            src=mp.GaussianSource(frequency=self.fcen, fwidth=self.df),
            center=mp.Vector3(),
            component=mp.Ez,
        )

        self.fr = mp.ForceRegion(
            center=mp.Vector3(0, 1.27), direction=mp.Y, size=mp.Vector3(4.38, 0)
        )

        self.fr_box = [
            mp.ForceRegion(
                direction=mp.X,
                center=mp.Vector3(1.0, 0, 0),
                size=mp.Vector3(0, 2.0, 0),
                weight=1.0,
            ),
            mp.ForceRegion(
                direction=mp.X,
                center=mp.Vector3(-1.0, 0, 0),
                size=mp.Vector3(0, 2.0, 0),
                weight=-1.0,
            ),
            mp.ForceRegion(
                direction=mp.Y,
                center=mp.Vector3(0, 1.0, 0),
                size=mp.Vector3(2.0, 0, 0),
                weight=1.0,
            ),
            mp.ForceRegion(
                direction=mp.Y,
                center=mp.Vector3(0, -1.0, 0),
                size=mp.Vector3(2.0, 0, 0),
                weight=-1.0,
            ),
        ]

    def test_force(self):
        """Test force calculations using a golden master."""

        sim = mp.Simulation(
            resolution=self.resolution,
            cell_size=self.cell,
            boundary_layers=[self.pml_layers],
            sources=[self.sources],
        )

        force = sim.add_force(self.fcen, 0, 1, self.fr, decimation_factor=1)
        force_decimated = sim.add_force(self.fcen, 0, 1, self.fr, decimation_factor=10)

        sim.run(
            until_after_sources=mp.stop_when_fields_decayed(
                50, mp.Ez, mp.Vector3(), 1e-6
            )
        )

        sim.display_forces(force)
        f = mp.get_forces(force)
        f_decimated = mp.get_forces(force_decimated)

        self.assertAlmostEqual(f[0], -0.11039089113393187)

        places = 6 if mp.is_single_precision() else 7
        self.assertAlmostEqual(f[0], f_decimated[0], places=places)

        # Test store and load of force as numpy array
        fdata = sim.get_force_data(force)
        sim.load_force_data(force, fdata)
        f_fdata = mp.get_forces(force)
        self.assertAlmostEqual(f_fdata[0], f[0], places=7)

    def test_force_with_symmetry(self):
        """Test that force calculations produce correct results with mirror symmetry."""

        # Run without mirror symmetry.
        sim_nosym = mp.Simulation(
            resolution=self.resolution,
            cell_size=self.cell,
            boundary_layers=[self.pml_layers],
            sources=[self.sources],
        )

        force_nosym = sim_nosym.add_force(self.fcen, 0, 1, *self.fr_box)
        sim_nosym.run(
            until_after_sources=mp.stop_when_fields_decayed(
                50, mp.Ez, mp.Vector3(), 1e-6
            )
        )
        f_nosym = mp.get_forces(force_nosym)

        # Run with mirror symmetry.
        sim_sym = mp.Simulation(
            resolution=self.resolution,
            cell_size=self.cell,
            boundary_layers=[self.pml_layers],
            sources=[self.sources],
            symmetries=[mp.Mirror(mp.X), mp.Mirror(mp.Y)],
        )

        force_sym = sim_sym.add_force(self.fcen, 0, 1, *self.fr_box)
        sim_sym.run(
            until_after_sources=mp.stop_when_fields_decayed(
                50, mp.Ez, mp.Vector3(), 1e-6
            )
        )
        f_sym = mp.get_forces(force_sym)

        self.assertAlmostEqual(f_nosym[0], f_sym[0], places=7)


if __name__ == "__main__":
    unittest.main()
