from __future__ import division

import unittest
import numpy as np
import meep as mp


susceptibilities = [
    mp.LorentzianSusceptibility(frequency=1.1, gamma=1e-5, sigma=0.5),
    mp.LorentzianSusceptibility(frequency=0.5, gamma=0.1, sigma=2e-5)
]


class MatFunc(object):

    def __call__(self, p):
        return mp.Medium(epsilon=2.25, E_susceptibilities=susceptibilities)


class TestMaterialDispersion(unittest.TestCase):

    def setUp(self):

        def mat_func(p):
            return mp.Medium(epsilon=2.25, E_susceptibilities=susceptibilities)

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
            material_function=mat_func,
            default_material=mp.air,
            resolution=20
        )

        self.expected = [
            0.1999342026399106,
            0.41053963810375294,
            0.6202409070451909,
            0.8285737385146619,
            1.0350739448523063,
            1.2392775309110078,
            1.4407208712852109,
        ]

    def test_material_dispersion_with_user_material_function(self):
        all_freqs = self.sim.run_k_points(200, self.kpts)
        res = [f.real for fs in all_freqs for f in fs]
        np.testing.assert_allclose(self.expected, res)

    def test_material_dispersion_with_user_callable_object(self):
        self.sim.material_function = MatFunc()
        all_freqs = self.sim.run_k_points(200, self.kpts)
        res = [f.real for fs in all_freqs for f in fs]
        np.testing.assert_allclose(self.expected, res)


if __name__ == '__main__':
    unittest.main()
