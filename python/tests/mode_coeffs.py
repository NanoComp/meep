from __future__ import division

import unittest
import numpy as np
import meep as mp


class TestModeCoeffs(unittest.TestCase):

    def run_mode_coeffs(self, mode_num, kpoint_func):

        resolution = 15

        w = 1   # width of waveguide
        L = 10  # length of waveguide

        Si = mp.Medium(epsilon=12.0)

        dair = 3.0
        dpml = 3.0

        sx = dpml + L + dpml
        sy = dpml + dair + w + dair + dpml
        cell_size = mp.Vector3(sx, sy, 0)

        prism_x = sx + 1
        prism_y = w / 2
        vertices = [mp.Vector3(-prism_x, prism_y),
                    mp.Vector3(prism_x, prism_y),
                    mp.Vector3(prism_x, -prism_y),
                    mp.Vector3(-prism_x, -prism_y)]

        geometry = [mp.Prism(vertices, height=mp.inf, material=Si)]

        boundary_layers = [mp.PML(dpml)]

        # mode frequency
        fcen = 0.20  # > 0.5/sqrt(11) to have at least 2 modes

        sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen, fwidth=0.5*fcen),
                                      eig_band=mode_num,
                                      size=mp.Vector3(0,sy-2*dpml,0),
                                      center=mp.Vector3(-0.5*sx+dpml,0,0),
                                      eig_match_freq=True,
                                      eig_resolution=32) ]

        sim = mp.Simulation(resolution=resolution,
                            cell_size=cell_size,
                            boundary_layers=boundary_layers,
                            geometry=geometry,
                            sources=sources,
                            symmetries=[mp.Mirror(mp.Y, phase=1 if mode_num % 2 == 1 else -1)])

        xm = 0.5*sx - dpml  # x-coordinate of monitor
        mflux = sim.add_mode_monitor(fcen, 0, 1, mp.ModeRegion(center=mp.Vector3(xm,0), size=mp.Vector3(0,sy-2*dpml)))
        mode_flux = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(xm,0), size=mp.Vector3(0,sy-2*dpml)))

        # sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(-0.5*sx+dpml,0), 1e-10))
        sim.run(until_after_sources=100)

        modes_to_check = [1, 2]  # indices of modes for which to compute expansion coefficients
        res = sim.get_eigenmode_coefficients(mflux, modes_to_check, kpoint_func=kpoint_func)

        self.assertTrue(res.kpoints[0].close(mp.Vector3(0.604301, 0, 0)))
        self.assertTrue(res.kpoints[1].close(mp.Vector3(0.494353, 0, 0), tol=1e-2))
        self.assertTrue(res.kdom[0].close(mp.Vector3(0.604301, 0, 0)))
        self.assertTrue(res.kdom[1].close(mp.Vector3(0.494353, 0, 0), tol=1e-2))

        mode_power = mp.get_fluxes(mode_flux)[0]

        TestPassed = True
        TOLERANCE = 5.0e-3
        c0 = res.alpha[mode_num - 1, 0, 0] # coefficient of forward-traveling wave for mode #mode_num
        for nm in range(1, len(modes_to_check)+1):
            if nm != mode_num:
                cfrel = np.abs(res.alpha[nm - 1, 0, 0]) / np.abs(c0)
                cbrel = np.abs(res.alpha[nm - 1, 0, 1]) / np.abs(c0)
                if cfrel > TOLERANCE or cbrel > TOLERANCE:
                    TestPassed = False

        self.sim = sim

        # test 1: coefficient of excited mode >> coeffs of all other modes
        self.assertTrue(TestPassed, msg="cfrel: {}, cbrel: {}".format(cfrel, cbrel))
        # test 2: |mode coeff|^2 = power
        self.assertAlmostEqual(mode_power / abs(c0**2), 1.0, places=1)

    def test_modes(self):
        self.run_mode_coeffs(1, None)
        self.run_mode_coeffs(2, None)

        # Test mp.get_eigenmode and EigenmodeData
        vol = mp.Volume(center=mp.Vector3(5), size=mp.Vector3(y=7))
        emdata = self.sim.get_eigenmode(0.2, mp.X, vol, 2, mp.Vector3())
        self.assertEqual(emdata.omega, 0.2)
        self.assertEqual(emdata.band_num, 2)

        eval_point = mp.Vector3(0.7, -0.2, 0.3)
        ex_at_eval_point = emdata.amplitude(eval_point, mp.Ex)
        hz_at_eval_point = emdata.amplitude(eval_point, mp.Hz)
        self.assertAlmostEqual(ex_at_eval_point, 2.4827383502368017e-05 + 0.16095103621219753j, places=2)
        self.assertAlmostEqual(hz_at_eval_point, 1.1381993864995965-3.7330585912764835e-06j, places=2)

    def test_kpoint_func(self):

        def kpoint_func(freq, mode):
            return mp.Vector3()

        self.run_mode_coeffs(1, kpoint_func)


if __name__ == '__main__':
    unittest.main()
