import math
import unittest

import meep as mp


class TestLDOS(unittest.TestCase):

    def setUp(self):
        resolution = 20

        cell = mp.Vector3(10, 10, 0)

        pml_layers = mp.PML(1.0)

        self.fcen = 1.0
        df = 1.0

        sources = mp.Source(src=mp.GaussianSource(self.fcen, fwidth=df), center=mp.Vector3(),
                            component=mp.Ez)

        symmetries = [mp.Mirror(mp.X), mp.Mirror(mp.Y)]

        self.sim = mp.Simulation(resolution=resolution,
                                 cell_size=cell,
                                 boundary_layers=[pml_layers],
                                 sources=[sources],
                                 symmetries=symmetries)

    def test_ldos(self):
        self.sim.run(
            mp.dft_ldos(self.fcen, 0, 1),
            until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(), 1e-6)
        )

        self.assertAlmostEqual(self.sim.ldos_data[0], 1.011459560620368)

    def test_ldos_user_object(self):
        ldos = mp.DftLdos(self.fcen, 0, 1)
        self.sim.run(
            mp.dft_ldos(ldos=ldos),
            until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(), 1e-6)
        )

        self.assertAlmostEqual(self.sim.ldos_data[0], 1.011459560620368)
        freqs = mp.get_ldos_freqs(ldos)
        self.assertEqual(ldos.omega_min, freqs[0] * 2 * math.pi)
        self.assertEqual(ldos.Nomega, 1)
        self.assertEqual(ldos.domega, 0)

    def test_invalid_dft_ldos(self):
        with self.assertRaises(ValueError):
            self.sim.run(mp.dft_ldos(mp.DftLdos(self.fcen, 0, 1)), until=200)


if __name__ == '__main__':
    unittest.main()
