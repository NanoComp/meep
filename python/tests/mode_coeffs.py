from __future__ import division

import unittest
import meep as mp
import numpy as np
import argparse

class ModeCoeffs(unittest.TestCase):

    def test_mode_coeffs(self):
        ##################################################
        # instantiate structure
        ##################################################
        resolution = 20
        fcen       = 0.15
        df         = 0.1
        mode_num   = 1
        plot       = False

        eps  = 12
        lx   = 1
        ly   = 3
        pad  = 4
        dpml = 2
        sx = (lx + pad + dpml)
        sy = (ly + pad + dpml)
        nfreq = 1

        x0        = mp.Vector3(0, 0)
        blk_size  = mp.Vector3(lx, ly)
        cell_size = mp.Vector3(sx, sy)

        blk = mp.Block(blk_size, center=x0, material=mp.Medium(epsilon=eps))

        src = mp.EigenModeSource(src=mp.GaussianSource(fcen, fwidth=df),
                                size=cell_size, center=x0,
                                eig_band=mode_num, direction=mp.Z)

        sim = mp.Simulation(cell_size=cell_size, geometry=[blk],
                            sources=[src], resolution=resolution,
                            boundary_layers=[mp.PML(dpml)])

        ##################################################
        # timestep, then compute mode-expansion coefficients
        ##################################################
        fr  = mp.FluxRegion(center=x0, size=cell_size, direction=mp.Z)
        flx = sim.add_flux(fcen, df, nfreq, fr)
        sim.run(until_after_sources=200)
        modes_to_check = [1,2,3,4]
        coeffs = sim.get_eigenmode_coefficients(flx, modes_to_check)

        ##################################################
        # confirm that the eigenmode coefficient for the
        # mode we excited is at least 3 orders of magnitude
        # larger than the coefficients of all other modes
        ##################################################
        TestPassed=True
        TOLERANCE=1.0e-3
        c0 = coeffs[2*(mode_num-1) + 0] # coefficient of forward-traveling wave for mode #mode_num
        for nm in range(1,len(modes_to_check)+1):
            if nm != mode_num:
                rx = np.abs(coeffs[2*(nm-1)+0]) / np.abs(c0)
                ix = np.abs(coeffs[2*(nm-1)+0]) / np.abs(c0)
                if rx > TOLERANCE or ix > TOLERANCE:
                    TestPassed=False;

        self.assertTrue(TestPassed)

        # generate plots of the structure and the eigenmode fields, useful for interactive use
        if plot:
            print(np.abs(coeffs))
            import matplotlib as mpl
            mpl.use('Qt4Agg')
            import matplotlib.pyplot as plt
            plt.ion()
            plt.figure()
            eps_slice = sim.get_array(sim.fields.total_volume(), component=mp.Dielectric)
            plt.subplot(1,3,1); plt.imshow(eps_slice.transpose())
            plt.colorbar();
            ex = sim.get_dft_array(flx, mp.Ex, 0)
            plt.subplot(1,3,2)
            plt.imshow(np.abs(ex.transpose()))
            plt.colorbar()
            ey = sim.get_dft_array(flx, mp.Ey, 0)
            plt.subplot(1,3,3)
            plt.imshow(np.abs(ey.transpose()))
            plt.colorbar();
            plt.show(block=True)

        #symmetries=[mp.Mirror(mp.Y)],
        #        self.h = mp.Harminv(mp.Ez, mp.Vector3(r + 0.1), fcen, df)

if __name__ == '__main__':
    unittest.main()
