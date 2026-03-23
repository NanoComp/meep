import unittest

import meep as mp


class TestHoleyWvgBands(unittest.TestCase):
    def setUp(self):
        cell = mp.Vector3(1, 12)
        b = mp.Block(
            size=mp.Vector3(mp.inf, 1.2, mp.inf), material=mp.Medium(epsilon=13)
        )
        c = mp.Cylinder(0.36)

        self.fcen = 0.25
        self.df = 1.5

        s = mp.Source(
            src=mp.GaussianSource(self.fcen, fwidth=self.df),
            component=mp.Hz,
            center=mp.Vector3(0.1234),
        )

        sym = mp.Mirror(direction=mp.Y, phase=-1)

        self.sim = mp.Simulation(
            cell_size=cell,
            geometry=[b, c],
            sources=[s],
            symmetries=[sym],
            boundary_layers=[mp.PML(1, direction=mp.Y)],
            resolution=20,
        )

    def test_run_k_points(self):
        all_freqs = self.sim.run_k_points(
            5, mp.interpolate(19, [mp.Vector3(), mp.Vector3(0.5)])
        )

        expected = [
            (0.19449071258888953, 0.0014994682279157307),
            (0.19805436668057624, -0.001207945176186543),
            (0.19302761095749146, 0.0010875391742976345),
            (0.1936140281065424, 0.0005279000201557691),
        ]

        self.assertTrue(any(all_freqs))
        for (r, i), f in zip(expected, all_freqs[17:21][0]):
            self.assertAlmostEqual(r, f.real)
            self.assertAlmostEqual(i, f.imag)

    def test_fields_at_kx(self):
        self.sim.k_point = mp.Vector3(3.5)
        h = mp.Harminv(mp.Hz, mp.Vector3(0.1234), self.fcen, self.df)
        self.sim.run(mp.after_sources(h), until_after_sources=300)

        expected = [
            (0.20007047969935934, 4.5993232095101406e-08),
            (0.3055995596872102, -5.4109280044667933e-08),
            (0.43969786779929054, 2.1860186449772234e-06),
            (0.4586461900867213, 5.635200431911015e-07),
            (0.5006545446251286, -0.0014323990588542283),
            (0.7415264941625447, -4.111723240258199e-05),
            (0.7626824561443201, -0.006619859775259674),
            (0.8244428928771226, -0.0004949457488574832),
            (0.8269995628775046, -0.001308667409393739),
            (0.9505882746507738, -0.005323624919387383),
            (0.9726976761260451, -0.003003111331548999),
            (0.986914116912195, -0.00417972044328907),
        ]

        self.assertTrue(h.modes)
        places = 4 if mp.is_single_precision() else 7
        for (r, i), m in zip(expected, h.modes):
            self.assertAlmostEqual(m.freq, r, places=places)
            self.assertAlmostEqual(m.decay, i, places=places)


if __name__ == "__main__":
    unittest.main()
