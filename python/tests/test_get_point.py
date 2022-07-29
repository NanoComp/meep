import math
import unittest

import numpy as np

import meep as mp


class TestGetPoint(unittest.TestCase):
    def test_get_point(self):
        sxy = 6  # cell size
        dpml = 1  # thickness of PML

        def sinusoid(p):
            r = (p.x**2 + p.y**2) ** 0.5
            return mp.Medium(index=1.0 + math.sin(2 * math.pi * r) ** 2)

        geometry = [
            mp.Block(center=mp.Vector3(), size=mp.Vector3(sxy, sxy), material=sinusoid)
        ]

        src = [
            mp.Source(
                mp.GaussianSource(1.0, fwidth=0.1), component=mp.Ez, center=mp.Vector3()
            )
        ]

        sim = mp.Simulation(
            cell_size=mp.Vector3(sxy, sxy),
            geometry=geometry,
            sources=src,
            k_point=mp.Vector3(),
            resolution=20,
            symmetries=[mp.Mirror(mp.X), mp.Mirror(mp.Y)],
            boundary_layers=[mp.PML(dpml)],
        )

        sim.run(until_after_sources=100)

        ## reference values for Ez and epsilon from serial run
        ez_ref = [
            -0.0002065983,
            -0.0001954795,
            -0.0000453570,
            0.0000311267,
            -0.0000121473,
            -0.0000410032,
            -0.0000341301,
            -0.0000275021,
            -0.0000397990,
            -0.0000351730,
            0.0000079602,
            0.0000227437,
            -0.0001092821,
            -0.0002202751,
            -0.0001408186,
            0.0006325076,
            0.0024890489,
            0.0027476069,
            0.0014815873,
            0.0004714913,
            -0.0004332029,
            -0.0007101315,
            -0.0003818581,
            -0.0000748507,
            0.0001408819,
            0.0001119776,
            0.0000395008,
            0.0000078844,
            -0.0000010431,
        ]

        eps_ref = [
            1.6458346134,
            1.2752837068,
            1.0974010956,
            1.0398089537,
            1.0465784716,
            1.0779924737,
            1.1059439286,
            1.1135579291,
            1.0971979186,
            1.0653178566,
            1.0391657283,
            1.0513779677,
            1.1466009312,
            1.3882154483,
            1.8496939317,
            2.5617731415,
            3.3788212533,
            3.9019494270,
            3.6743431894,
            2.7285622651,
            1.6635165033,
            1.0891237010,
            1.1485969863,
            1.9498398061,
            3.3100416367,
            3.9038800599,
            2.8471862395,
            1.4742605488,
            1.0370162714,
        ]

        x = np.linspace(-0.865692, 2.692867, 29)
        places = 5 if mp.is_single_precision() else 10
        for j in range(x.size):
            self.assertAlmostEqual(
                np.real(sim.get_field_point(mp.Ez, mp.Vector3(x[j], -0.394862))),
                ez_ref[j],
                places=places,
            )
            self.assertAlmostEqual(
                sim.get_epsilon_point(mp.Vector3(x[j], 2.967158)),
                eps_ref[j],
                places=places,
            )


if __name__ == "__main__":
    unittest.main()
