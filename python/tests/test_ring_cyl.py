import unittest

from utils import ApproxComparisonTestCase

import meep as mp


class TestRingCyl(ApproxComparisonTestCase):
    def setUp(self):
        n = 3.4
        w = 1
        self.r = 1
        pad = 4
        dpml = 2
        sr = self.r + w + pad + dpml
        dimensions = mp.CYLINDRICAL
        cell = mp.Vector3(sr, 0, 0)
        m = 3

        geometry = [
            mp.Block(
                center=mp.Vector3(self.r + (w / 2)),
                size=mp.Vector3(w, mp.inf, mp.inf),
                material=mp.Medium(index=n),
            )
        ]

        pml_layers = [mp.PML(dpml)]
        resolution = 10

        self.fcen = 0.15
        self.df = 0.1
        sources = [
            mp.Source(
                src=mp.GaussianSource(self.fcen, fwidth=self.df),
                component=mp.Ez,
                center=mp.Vector3(self.r + 0.1),
            )
        ]

        self.sim = mp.Simulation(
            cell_size=cell,
            geometry=geometry,
            boundary_layers=pml_layers,
            resolution=resolution,
            sources=sources,
            dimensions=dimensions,
            m=m,
            split_chunks_evenly=False,
        )

    def test_ring_cyl(self):
        expected = [
            0.11835455441250553,
            -6.907792691629741e-4,
            85.66741917133473,
            0.025701906263451237,
            -0.024027038833537524,
            -0.009126302124459489,
        ]

        h = mp.Harminv(mp.Ez, mp.Vector3(self.r + 0.1), self.fcen, self.df)
        self.sim.run(mp.after_sources(h), until_after_sources=200)

        m = h.modes[0]
        res = [m.freq, m.decay, m.Q, abs(m.amp), m.amp.real, m.amp.imag]

        tol = 1e-6 if mp.is_single_precision() else 1e-7
        self.assertClose(expected, res, epsilon=tol)


if __name__ == "__main__":
    unittest.main()
