import unittest

import meep as mp


class TestForce(unittest.TestCase):
    def setUp(self):

        resolution = 20
        cell = mp.Vector3(10, 10)
        pml_layers = mp.PML(1.0)
        fcen = 1.0
        df = 1.0
        sources = mp.Source(
            src=mp.GaussianSource(fcen, fwidth=df), center=mp.Vector3(), component=mp.Ez
        )

        self.sim = mp.Simulation(
            resolution=resolution,
            cell_size=cell,
            boundary_layers=[pml_layers],
            sources=[sources],
        )

        fr = mp.ForceRegion(mp.Vector3(y=1.27), direction=mp.Y, size=mp.Vector3(4.38))
        self.myforce = self.sim.add_force(fcen, 0, 1, fr, decimation_factor=1)
        self.myforce_decimated = self.sim.add_force(
            fcen, 0, 1, fr, decimation_factor=10
        )

    def test_force(self):

        self.sim.run(
            until_after_sources=mp.stop_when_fields_decayed(
                50, mp.Ez, mp.Vector3(), 1e-6
            )
        )

        # Test store and load of force as numpy array
        fdata = self.sim.get_force_data(self.myforce)
        self.sim.load_force_data(self.myforce, fdata)

        self.sim.display_forces(self.myforce)
        f = mp.get_forces(self.myforce)

        self.assertAlmostEqual(f[0], -0.11039089113393187)

        places = 6 if mp.is_single_precision() else 7
        self.assertAlmostEqual(
            f[0], mp.get_forces(self.myforce_decimated)[0], places=places
        )


if __name__ == "__main__":
    unittest.main()
