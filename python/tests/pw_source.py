from __future__ import division

import cmath
import math
import unittest

import meep as mp


class TestPwSource(unittest.TestCase):

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
        k = kdir.unit().scale(2 * math.pi * fcen)

        sources = [
            mp.Source(
                mp.ContinuousSource(fcen, fwidth=df),
                component=mp.Ez,
                center=mp.Vector3(-0.5 * s, 0),
                size=mp.Vector3(0, s),
                amp_func=pw_amp(k, mp.Vector3(x=-0.5 * s))
            ),
            mp.Source(
                mp.ContinuousSource(fcen, fwidth=df),
                component=mp.Ez,
                center=mp.Vector3(0, -0.5 * s),
                size=mp.Vector3(s, 0),
                amp_func=pw_amp(k, mp.Vector3(y=-0.5 * s))
            )
        ]

        self.sim = mp.Simulation(
            cell_size=cell,
            sources=sources,
            boundary_layers=pml_layers,
            resolution=resolution
        )

    def test_pw_source(self):
        self.sim.run(mp.at_end(mp.output_efield_z), until=400)

if __name__ == '__main__':
    unittest.main()
