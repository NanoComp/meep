import cmath
import math
import unittest

from utils import ApproxComparisonTestCase

import meep as mp


class TestPwSource(ApproxComparisonTestCase):
    def setUp(self):
        s = 11
        dpml = 1

        sxy = s + 2 * dpml
        cell = mp.Vector3(sxy, sxy, 0)

        pml_layers = [mp.PML(dpml)]
        resolution = 10

        def pw_amp(k, x0):
            def _pw_amp(x):
                return cmath.exp(1j * k.dot(x + x0))

            return _pw_amp

        fcen = 0.8
        df = 0.02
        kdir = mp.Vector3(1, 1)
        self.k = kdir.unit().scale(2 * math.pi * fcen)

        sources = [
            mp.Source(
                mp.ContinuousSource(fcen, fwidth=df),
                component=mp.Ez,
                center=mp.Vector3(-0.5 * s, 0),
                size=mp.Vector3(0, s),
                amp_func=pw_amp(self.k, mp.Vector3(x=-0.5 * s)),
            ),
            mp.Source(
                mp.ContinuousSource(fcen, fwidth=df),
                component=mp.Ez,
                center=mp.Vector3(0, -0.5 * s),
                size=mp.Vector3(s, 0),
                amp_func=pw_amp(self.k, mp.Vector3(y=-0.5 * s)),
            ),
        ]

        self.sim = mp.Simulation(
            cell_size=cell,
            sources=sources,
            boundary_layers=pml_layers,
            resolution=resolution,
        )
        self.sim.use_output_directory(self.temp_dir)
        self.s = s

    @classmethod
    def setUpClass(cls):
        cls.temp_dir = mp.make_output_directory()

    @classmethod
    def tearDownClass(cls):
        mp.delete_directory(cls.temp_dir)

    def test_pw_source(self):
        self.sim.run(mp.at_end(mp.output_efield_z), until=400)

        v1 = mp.Vector3(0.5 * self.s, 0)
        v2 = mp.Vector3(0.5 * self.s, 0.5 * self.s)

        pt1 = self.sim.get_field_point(mp.Ez, v1)
        pt2 = self.sim.get_field_point(mp.Ez, v2)

        tol = 1e-4 if mp.is_single_precision() else 1e-9
        self.assertClose(pt1 / pt2, 27.557668029008262, epsilon=tol)

        self.assertAlmostEqual(
            cmath.exp(1j * self.k.dot(v1 - v2)),
            0.7654030066070924 - 0.6435512702783076j,
        )


if __name__ == "__main__":
    unittest.main()
