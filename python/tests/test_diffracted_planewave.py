import cmath
import math
import unittest

import numpy as np

import meep as mp

# Computes the mode coefficient of the transmitted orders of
# a binary grating given an incident planewave and verifies
# that the results are the same when using either a band number
# or `DiffractedPlanewave` object in `get_eigenmode_coefficients`.


class TestDiffractedPlanewave(unittest.TestCase):
    @classmethod
    def setUp(cls):
        cls.resolution = 50  # pixels/um
        cls.dpml = 1.0  # PML thickness
        cls.dsub = 3.0  # substrate thickness
        cls.dpad = 3.0  # length of padding between grating and PML

        cls.wvl = 0.5  # center wavelength
        cls.fcen = 1 / cls.wvl  # center frequency

        cls.ng = 1.5
        cls.glass = mp.Medium(index=cls.ng)

        cls.pml_layers = [mp.PML(thickness=cls.dpml, direction=mp.X)]

    def run_binary_grating_diffraction(self, gp, gh, gdc, theta):
        sx = self.dpml + self.dsub + gh + self.dpad + self.dpml
        sy = gp
        cell_size = mp.Vector3(sx, sy, 0)

        # rotation angle of incident planewave
        # counter clockwise (CCW) about Z axis, 0 degrees along +X axis
        theta_in = math.radians(theta)

        # k (in source medium) with correct length (plane of incidence: XY)
        k = mp.Vector3(self.fcen * self.ng).rotate(mp.Vector3(z=1), theta_in)

        eig_parity = mp.ODD_Z
        if theta == 0:
            k = mp.Vector3()
            eig_parity += mp.EVEN_Y
            symmetries = [mp.Mirror(direction=mp.Y)]
        else:
            symmetries = []

        def pw_amp(k, x0):
            def _pw_amp(x):
                return cmath.exp(1j * 2 * math.pi * k.dot(x + x0))

            return _pw_amp

        src_pt = mp.Vector3(-0.5 * sx + self.dpml, 0, 0)
        sources = [
            mp.Source(
                mp.GaussianSource(self.fcen, fwidth=0.1 * self.fcen),
                component=mp.Ez,
                center=src_pt,
                size=mp.Vector3(0, sy, 0),
                amp_func=pw_amp(k, src_pt),
            )
        ]

        geometry = [
            mp.Block(
                material=self.glass,
                size=mp.Vector3(self.dpml + self.dsub, mp.inf, mp.inf),
                center=mp.Vector3(-0.5 * sx + 0.5 * (self.dpml + self.dsub), 0, 0),
            ),
            mp.Block(
                material=self.glass,
                size=mp.Vector3(gh, gdc * gp, mp.inf),
                center=mp.Vector3(-0.5 * sx + self.dpml + self.dsub + 0.5 * gh, 0, 0),
            ),
        ]

        sim = mp.Simulation(
            resolution=self.resolution,
            cell_size=cell_size,
            boundary_layers=self.pml_layers,
            geometry=geometry,
            k_point=k,
            sources=sources,
            symmetries=symmetries,
        )

        tran_pt = mp.Vector3(0.5 * sx - self.dpml, 0, 0)
        tran_flux = sim.add_mode_monitor(
            self.fcen, 0, 1, mp.FluxRegion(center=tran_pt, size=mp.Vector3(0, sy, 0))
        )

        sim.run(
            until_after_sources=mp.stop_when_fields_decayed(20, mp.Ez, src_pt, 1e-6)
        )

        m_plus = int(np.floor((self.fcen - k.y) * gp))
        m_minus = int(np.ceil((-self.fcen - k.y) * gp))

        if theta == 0:
            orders = range(m_plus)
        else:
            # ordering of the modes computed by MPB is according to *decreasing*
            # values of kx (i.e., closest to propagation direction of 0° or +x)
            ms = range(m_minus, m_plus + 1)
            kx = lambda m: np.power(self.fcen, 2) - np.power(k.y + m / gp, 2)
            kxs = [kx(m) for m in ms]
            ids = np.flip(np.argsort(kxs))
            orders = [ms[d] for d in ids]

        for band, order in enumerate(orders):
            res = sim.get_eigenmode_coefficients(
                tran_flux, [band + 1], eig_parity=eig_parity
            )
            tran_ref = abs(res.alpha[0, 0, 0]) ** 2
            if theta == 0:
                tran_ref = 0.5 * tran_ref
            vg_ref = res.vgrp[0]
            kdom_ref = res.kdom[0]

            res = sim.get_eigenmode_coefficients(
                tran_flux,
                mp.DiffractedPlanewave((0, order, 0), mp.Vector3(0, 1, 0), 1, 0),
            )
            tran_dp = abs(res.alpha[0, 0, 0]) ** 2
            if (theta == 0) and (order == 0):
                tran_dp = 0.5 * tran_dp
            vg_dp = res.vgrp[0]
            kdom_dp = res.kdom[0]

            err = abs(tran_ref - tran_dp) / tran_ref
            print(
                "tran:, {:2d} (band), {:2d} (order), "
                "{:10.8f} (band num.), {:10.8f} (diff. pw.), "
                "{:10.8f} (error)".format(band, order, tran_ref, tran_dp, err)
            )

            self.assertAlmostEqual(vg_ref, vg_dp, places=4)
            self.assertAlmostEqual(tran_ref, tran_dp, places=4)
            if theta == 0:
                self.assertAlmostEqual(abs(kdom_ref.x), kdom_dp.x, places=5)
                self.assertAlmostEqual(abs(kdom_ref.y), kdom_dp.y, places=5)
                self.assertAlmostEqual(abs(kdom_ref.z), kdom_dp.z, places=5)
            else:
                self.assertAlmostEqual(kdom_ref.x, kdom_dp.x, places=5)
                self.assertAlmostEqual(kdom_ref.y, kdom_dp.y, places=5)
                self.assertAlmostEqual(kdom_ref.z, kdom_dp.z, places=5)

        print("PASSED.")

    def test_diffracted_planewave(self):
        self.run_binary_grating_diffraction(2.6, 0.4, 0.6, 0)
        self.run_binary_grating_diffraction(2.6, 0.4, 0.6, 13.4)

        # self.run_binary_grating_diffraction(10.0,0.5,0.5,0)
        # self.run_binary_grating_diffraction(10.0,0.5,0.5,10.7)


if __name__ == "__main__":
    unittest.main()
