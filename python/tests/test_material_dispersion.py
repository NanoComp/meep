import unittest

import numpy as np

import meep as mp


class TestMaterialDispersion(unittest.TestCase):
    def test_material_dispersion_with_user_material(self):
        susceptibilities = [
            mp.LorentzianSusceptibility(frequency=1.1, gamma=1e-5, sigma=0.5),
            mp.LorentzianSusceptibility(frequency=0.5, gamma=0.1, sigma=2e-5),
        ]

        def mat_func(p):
            return mp.Medium(epsilon=2.25, E_susceptibilities=susceptibilities)

        fcen = 1.0
        df = 2.0

        sources = mp.Source(
            mp.GaussianSource(fcen, fwidth=df), component=mp.Ez, center=mp.Vector3()
        )

        kmin = 0.3
        kmax = 2.2
        k_interp = 5

        kpts = mp.interpolate(k_interp, [mp.Vector3(kmin), mp.Vector3(kmax)])

        self.sim = mp.Simulation(
            cell_size=mp.Vector3(),
            geometry=[],
            sources=[sources],
            material_function=mat_func,
            default_material=mp.air,
            resolution=20,
        )

        all_freqs = self.sim.run_k_points(200, kpts)
        res = [f.real for fs in all_freqs for f in fs]

        expected = [
            0.1999342026399106,
            0.41053963810375294,
            0.6202409070451909,
            0.8285737385146619,
            1.0350739448523063,
            1.2392775309110078,
            1.4407208712852109,
        ]

        np.testing.assert_allclose(expected, res)


if __name__ == "__main__":
    unittest.main()
