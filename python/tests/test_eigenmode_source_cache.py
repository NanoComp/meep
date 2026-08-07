"""Tests for EigenModeSource(eig_use_cache=True): re-placing an eigenmode
source from its cache (e.g. after reset_meep) must produce the same fields
as re-running the MPB eigensolver."""

import unittest

import numpy as np

import meep as mp


class TestEigenmodeSourceCache(unittest.TestCase):
    fcen = 0.15
    df = 0.1
    parity = mp.ODD_Z + mp.EVEN_Y

    def make_sim(self):
        return mp.Simulation(
            cell_size=mp.Vector3(12, 6),
            geometry=[
                mp.Block(
                    size=mp.Vector3(mp.inf, 1, mp.inf),
                    center=mp.Vector3(),
                    material=mp.Medium(epsilon=12),
                )
            ],
            boundary_layers=[mp.PML(1.0)],
            resolution=20,
        )

    def make_src(self, use_cache):
        return mp.EigenModeSource(
            mp.GaussianSource(self.fcen, fwidth=self.df),
            center=mp.Vector3(-4, 0),
            size=mp.Vector3(0, 4),
            eig_band=1,
            eig_parity=self.parity,
            eig_match_freq=True,
            eig_use_cache=use_cache,
        )

    def mode_coeff(self, sim):
        flux = sim.add_mode_monitor(
            self.fcen,
            0,
            1,
            mp.ModeRegion(center=mp.Vector3(3, 0), size=mp.Vector3(0, 4)),
        )
        sim.run(until_after_sources=mp.stop_when_dft_decayed(1e-8))
        res = sim.get_eigenmode_coefficients(flux, [1], eig_parity=self.parity)
        return res.alpha[0, 0, 0]

    def test_cache_matches_solver(self):
        sim = self.make_sim()
        sim.change_sources([self.make_src(use_cache=False)])
        alpha_ref = self.mode_coeff(sim)

        src = self.make_src(use_cache=True)
        sim = self.make_sim()
        sim.change_sources([src])
        alpha_cold = self.mode_coeff(sim)
        self.assertIsNotNone(src._eigenmode_data)

        # the per-iteration path of an optimization loop: reset + re-place
        sim.reset_meep()
        sim.change_sources([src])
        alpha_warm = self.mode_coeff(sim)

        self.assertAlmostEqual(
            abs(alpha_cold - alpha_ref) / abs(alpha_ref), 0, places=8
        )
        self.assertAlmostEqual(
            abs(alpha_warm - alpha_ref) / abs(alpha_ref), 0, places=8
        )

    def test_invalidate_cache(self):
        src = self.make_src(use_cache=True)
        sim = self.make_sim()
        sim.change_sources([src])
        sim.init_sim()
        first = src._eigenmode_data
        self.assertIsNotNone(first)

        src.invalidate_eigenmode_cache()
        self.assertIsNone(src._eigenmode_data)

        sim.reset_meep()
        sim.change_sources([src])
        sim.init_sim()
        self.assertIsNotNone(src._eigenmode_data)


if __name__ == "__main__":
    unittest.main()
