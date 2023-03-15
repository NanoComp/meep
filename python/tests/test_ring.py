# Python port of meep/examples/ring.ctl
# Calculating 2d ring-resonator modes, from the Meep tutorial.
import unittest

import meep as mp
import numpy as np
from scipy.signal import find_peaks


class TestRing(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp_dir = mp.make_output_directory()

    @classmethod
    def tearDownClass(cls):
        mp.delete_directory(cls.temp_dir)

    def init(self):
        n = 3.4
        w = 1
        r = 1
        pad = 4
        dpml = 2
        sxy = 2 * (r + w + pad + dpml)

        dielectric = mp.Medium(epsilon=n**2)
        air = mp.Medium()

        c1 = mp.Cylinder(r + w, material=dielectric)
        c2 = mp.Cylinder(r, material=air)

        fcen = 0.15
        df = 0.1

        src = mp.Source(mp.GaussianSource(fcen, fwidth=df), mp.Ez, mp.Vector3(r + 0.1))

        self.sim = mp.Simulation(
            cell_size=mp.Vector3(sxy, sxy),
            geometry=[c1, c2],
            sources=[src],
            resolution=10,
            symmetries=[mp.Mirror(mp.Y)],
            boundary_layers=[mp.PML(dpml)],
        )

        self.sim.use_output_directory(self.temp_dir)
        self.h = mp.Harminv(mp.Ez, mp.Vector3(r + 0.1), fcen, df)
        self.p = mp.Pade(mp.Ez, mp.Vector3(r + 0.1), sample_rate=4)

    def test_harminv(self):
        self.init()

        self.sim.run(
            mp.at_beginning(mp.output_epsilon),
            mp.after_sources(self.h),
            until_after_sources=300,
        )

        m1 = self.h.modes[0]

        self.assertAlmostEqual(m1.freq, 0.118101315147, places=4)
        self.assertAlmostEqual(m1.decay, -0.000731513241623, places=4)
        self.assertAlmostEqual(abs(m1.amp), 0.00341267634436, places=4)
        self.assertAlmostEqual(m1.amp.real, -0.00304951667301, places=4)
        self.assertAlmostEqual(m1.amp.imag, -0.00153192946717, places=3)

        v = mp.Vector3(1, 1)
        fp = self.sim.get_field_point(mp.Ez, v)
        ep = self.sim.get_epsilon_point(v)

        places = 5 if mp.is_single_precision() else 7
        self.assertAlmostEqual(ep, 11.559999999999999, places=places)
        self.assertAlmostEqual(fp, -0.08185972142450348, places=places)

    def test_pade(self):
        self.init()

        self.sim.run(
            self.p,
            mp.after_sources(self.h),
            until_after_sources=300,
        )

        freqs = np.linspace(
            self.h.fcen - self.h.df / 2, self.h.fcen + self.h.df / 2, 1000
        )

        freq_domain = [self.p.freq_response(freq) for freq in freqs]

        idx = find_peaks(
            np.abs(freq_domain) ** 2 / max(np.abs(freq_domain) ** 2), prominence=1e-4
        )[0]

        self.assertAlmostEqual(freqs[idx[0]], self.h.modes[0].freq, places=3)


if __name__ == "__main__":
    unittest.main()
