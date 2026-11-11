"""Tests the warnings for materials whose ε, μ, or σ tensors are not
positive-(semi)definite (issue #3231)."""

import contextlib
import io
import unittest

import meep as mp


class TestMaterialPositiveDefinite(unittest.TestCase):
    def init_sim_stderr(self, geometry=None, default_material=mp.air, **kwargs):
        """Initializes a small 2d simulation and returns what it printed to stderr."""
        sim = mp.Simulation(
            cell_size=mp.Vector3(2, 2),
            resolution=10,
            geometry=geometry or [],
            default_material=default_material,
            eps_averaging=False,
            **kwargs,
        )
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            sim.init_sim()
        return stderr.getvalue()

    def block(self, material, x):
        return mp.Block(
            center=mp.Vector3(x), size=mp.Vector3(0.5, 0.5), material=material
        )

    def test_physical_materials(self):
        drude = mp.Medium(
            epsilon=1,
            E_susceptibilities=[
                mp.DrudeSusceptibility(frequency=1, gamma=0.1, sigma=2)
            ],
        )
        uniaxial = mp.Medium(
            epsilon_diag=mp.Vector3(2, 3, 4),
            epsilon_offdiag=mp.Vector3(0.5, 0, 0),
            D_conductivity_diag=mp.Vector3(0.1, 0, 0),
            E_susceptibilities=[
                mp.LorentzianSusceptibility(
                    frequency=1, gamma=0.1, sigma_diag=mp.Vector3(1, 0, 0)
                )
            ],
        )
        geometry = [
            self.block(mp.Medium(index=3.5), -0.6),
            self.block(drude, 0),
            self.block(uniaxial, 0.6),
            self.block(mp.metal, 0),
        ]
        stderr = self.init_sim_stderr(geometry=geometry)
        self.assertNotIn("positive", stderr)

    def test_epsilon(self):
        # eigenvalues of epsilon are 2 ± 3, so it has a negative eigenvalue
        # even though all of its diagonal entries are positive
        bad = mp.Medium(
            epsilon_diag=mp.Vector3(2, 2, 2), epsilon_offdiag=mp.Vector3(3, 0, 0)
        )
        stderr = self.init_sim_stderr(geometry=[self.block(bad, 0)])
        if mp.am_master():
            self.assertIn("material epsilon tensor", stderr)
            self.assertIn("is not positive-definite", stderr)

    def test_mu_default_material(self):
        stderr = self.init_sim_stderr(default_material=mp.Medium(mu=-2))
        if mp.am_master():
            self.assertIn("material mu tensor", stderr)

    def test_susceptibility_sigma(self):
        # every 2x2 principal minor of sigma is positive but det(sigma) < 0
        sus = mp.LorentzianSusceptibility(
            frequency=1,
            gamma=0.1,
            sigma_diag=mp.Vector3(1, 1, 1),
            sigma_offdiag=mp.Vector3(0.9, 0.9, -0.9),
        )
        bad = mp.Medium(epsilon=2, E_susceptibilities=[sus])
        stderr = self.init_sim_stderr(extra_materials=[bad])
        if mp.am_master():
            self.assertIn("material E_susceptibility sigma tensor", stderr)
            self.assertIn("is not positive-semidefinite", stderr)

    def test_conductivity(self):
        bad = mp.Medium(epsilon=2, D_conductivity=-0.1)
        stderr = self.init_sim_stderr(geometry=[self.block(bad, 0)])
        if mp.am_master():
            self.assertIn("material D_conductivity tensor", stderr)

    def test_warn_once(self):
        bad = mp.Medium(epsilon=-2)
        geometry = [self.block(bad, -0.6), self.block(bad, 0.6)]
        stderr = self.init_sim_stderr(geometry=geometry, extra_materials=[bad])
        if mp.am_master():
            self.assertEqual(stderr.count("material epsilon tensor"), 1)


if __name__ == "__main__":
    unittest.main()
