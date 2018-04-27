from __future__ import division

import argparse
import math
import unittest
import numpy as np
import meep as mp

class ModeCoeffs(unittest.TestCase):

    def run_mode_coeffs(self, mode_num):

        resolution = 15

        w = 1   # width of waveguide
        L = 10  # length of waveguide

        Si = mp.Medium(epsilon=12.0)

        dair = 3.0
        dpml = 3.0

        sx = dpml + L + dpml
        sy = dpml + dair + w + dair + dpml
        cell_size = mp.Vector3(sx, sy, 0)

        geometry = [ mp.Block(material=Si, center=mp.Vector3(), size=mp.Vector3(mp.inf,w,mp.inf)) ]

        boundary_layers = [ mp.PML(dpml) ]

        # mode frequency
        fcen = 0.20 # > 0.5/sqrt(11) to have at least 2 modes

        sources = [ mp.EigenModeSource(src=mp.GaussianSource(fcen, fwidth=0.5*fcen),
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
                            )

        xm = 0.5*sx - dpml  # x-coordinate of monitor
        mflux = sim.add_eigenmode(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(xm,0), size=mp.Vector3(0,sy-2*dpml)))
        mode_flux = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(xm,0), size=mp.Vector3(0,sy-2*dpml)))

        #sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(-0.5*sx+dpml,0), 1e-10))
        sim.run(until_after_sources=100)

        modes_to_check = [1,2]  # indices of modes for which to compute expansion coefficients
        alpha = sim.get_eigenmode_coefficients(mflux, modes_to_check)

        mode_power = mp.get_fluxes(mode_flux)[0]

        TestPassed=True
        TOLERANCE=5.0e-3
        c0 = alpha[2*(mode_num-1) + 0] # coefficient of forward-traveling wave for mode #mode_num
        for nm in range(1,len(modes_to_check)+1):
            if nm != mode_num:
                cfrel = np.abs(alpha[2*(nm-1)+0]) / np.abs(c0)
                cbrel = np.abs(alpha[2*(nm-1)+1]) / np.abs(c0)
                if cfrel > TOLERANCE or cbrel > TOLERANCE:
                    TestPassed=False
        self.assertTrue(TestPassed) # test 1: coefficient of excited mode >> coeffs of all other modes

        self.assertAlmostEqual(mode_power / abs(c0**2), 1.0, places=1) # test 2: |mode coeff|^2 = power

    def test_modes(self):
        self.run_mode_coeffs(1)
        self.run_mode_coeffs(2)


if __name__ == '__main__':
    unittest.main()
