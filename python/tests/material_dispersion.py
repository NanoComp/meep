from __future__ import division

import unittest
import numpy as np
import meep as mp


class TestMaterialDispersion(unittest.TestCase):

    def setUp(self):
        susceptibilities = [
            mp.LorentzianSusceptibility(frequency=1.1, gamma=1e-5, sigma=0.5),
            mp.LorentzianSusceptibility(frequency=0.5, gamma=0.1, sigma=2e-5)
        ]

        default_material = mp.Medium(epsilon=2.25, E_susceptibilities=susceptibilities)

        fcen = 1.0
        df = 2.0

        sources = mp.Source(
            mp.GaussianSource(fcen, fwidth=df),
            component=mp.Ez,
            center=mp.Vector3()
        )

        kmin = 0.3
        kmax = 2.2
        k_interp = 5

        self.kpts = mp.interpolate(k_interp, [mp.Vector3(kmin), mp.Vector3(kmax)])

        self.sim = mp.Simulation(
            cell_size=mp.Vector3(),
            geometry=[],
            sources=[sources],
            default_material=default_material,
            resolution=20
        )

    def test_material_dispersion(self):
        all_freqs = self.sim.run(self.kpts, k_points=200)

        expected = [
            (0.18039256903351222, -9.43254394797079e-8),
            (1.2210179998635848, -4.995198119184291e-6),
            (0.3671711444546296, -9.116945419273966e-7),
            (1.2317965266011561, -4.900764116872976e-6),
            (0.5450978259018661, -4.937192953188026e-6),
            (1.2535585975158368, -4.683483806835463e-6),
            (0.7058362678935716, -1.311598513657167e-6),
            (1.2932545600392702, -4.2322347271138075e-6),
            (0.8372527022546067, -1.7723690643203654e-6),
            (1.3619816335728656, -3.4369744685704915e-6),
            (0.9296922571569821, -2.668853795179047e-6),
            (1.4685401195056051, -2.446704263810018e-6),
            (0.9866379289119844, -3.4474204597212896e-6),
            (1.6087133818517, -1.627744600303355e-6),
        ]

        res = [(f.real, f.imag) for fs in all_freqs[:10] for f in fs]

        np.testing.assert_allclose(expected, res)


if __name__ == '__main__':
    unittest.main()
