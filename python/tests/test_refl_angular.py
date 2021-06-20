from __future__ import division

import math
import unittest
import numpy as np
import meep as mp


class TestReflAngular(unittest.TestCase):

    def test_refl_angular(self):
        resolution = 100

        dpml = 1.0
        sz = 10
        sz = sz + 2 * dpml
        cell_size = mp.Vector3(z=sz)
        pml_layers = [mp.PML(dpml)]

        wvl_min = 0.4
        wvl_max = 0.8
        fmin = 1 / wvl_max
        fmax = 1 / wvl_min
        fcen = 0.5 * (fmin + fmax)
        df = fmax - fmin
        nfreq = 50

        theta_r = math.radians(0)

        k = mp.Vector3(math.sin(theta_r), 0, math.cos(theta_r)).scale(fcen)

        dimensions = 1

        sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ex,
                             center=mp.Vector3(z=-0.5 * sz + dpml))]

        sim = mp.Simulation(cell_size=cell_size,
                            boundary_layers=pml_layers,
                            sources=sources,
                            k_point=k,
                            dimensions=dimensions,
                            resolution=resolution)

        refl_fr = mp.FluxRegion(center=mp.Vector3(z=-0.25 * sz))
        refl = sim.add_flux(fcen, df, nfreq, refl_fr)

        sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ex, mp.Vector3(z=-0.5 * sz + dpml), 1e-9))

        empty_data = sim.get_flux_data(refl)
        sim.reset_meep()

        geometry = [mp.Block(mp.Vector3(mp.inf, mp.inf, 0.5 * sz), center=mp.Vector3(z=0.25 * sz),
                             material=mp.Medium(index=3.5))]

        sim = mp.Simulation(cell_size=cell_size,
                            geometry=geometry,
                            boundary_layers=pml_layers,
                            sources=sources,
                            k_point=k,
                            dimensions=dimensions,
                            resolution=resolution)

        refl = sim.add_flux(fcen, df, nfreq, refl_fr)
        sim.load_minus_flux_data(refl, empty_data)

        sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ex, mp.Vector3(z=-0.5 * sz + dpml), 1e-9))

        refl_flux = mp.get_fluxes(refl)
        freqs = mp.get_flux_freqs(refl)

        expected = [
            (1.25, -1.123696883299492e-6),
            (1.2755102040816326, -2.5749667658387866e-6),
            (1.3010204081632653, -5.70480204599006e-6),
            (1.3265306122448979, -1.2220464827582253e-5),
            (1.3520408163265305, -2.531247480206961e-5),
            (1.3775510204081631, -5.069850309492639e-5),
            (1.4030612244897958, -9.819256552437341e-5),
            (1.4285714285714284, -1.8390448659017395e-4),
            (1.454081632653061, -3.330762066794769e-4),
            (1.4795918367346936, -5.833650417163753e-4),
            (1.5051020408163263, -9.8807834237052e-4),
            (1.5306122448979589, -0.001618472171445976),
            (1.5561224489795915, -0.0025638388059825985),
            (1.5816326530612241, -0.003927863989816029),
            (1.6071428571428568, -0.005819831283556752),
            (1.6326530612244894, -0.008339881000982728),
            (1.658163265306122, -0.011558769654206626),
            (1.6836734693877546, -0.015494308354153143),
            (1.7091836734693873, -0.02008850084337135),
            (1.73469387755102, -0.025190871516857616),
            (1.7602040816326525, -0.030553756123198477),
            (1.7857142857142851, -0.03584404966066722),
            (1.8112244897959178, -0.040672967700428275),
            (1.8367346938775504, -0.04464118393086191),
            (1.862244897959183, -0.047392712128477496),
            (1.8877551020408156, -0.048667403362887635),
            (1.9132653061224483, -0.048341494285878264),
            (1.938775510204081, -0.04644739000778679),
            (1.9642857142857135, -0.043168390293742316),
            (1.9897959183673462, -0.0388094755730579),
            (2.0153061224489788, -0.03375052221907117),
            (2.0408163265306114, -0.02839209067703472),
            (2.066326530612244, -0.023104245646230648),
            (2.0918367346938767, -0.01818725699718267),
            (2.1173469387755093, -0.013849270759480073),
            (2.142857142857142, -0.010201733597436358),
            (2.1683673469387745, -0.007269616609175294),
            (2.193877551020407, -0.005011210495189995),
            (2.21938775510204, -0.0033417192031464896),
            (2.2448979591836724, -0.0021557351734376256),
            (2.270408163265305, -0.0013453062176115673),
            (2.2959183673469377, -8.121742663131631e-4),
            (2.3214285714285703, -4.7433135191915683e-4),
            (2.346938775510203, -2.6799188013374266e-4),
            (2.3724489795918355, -1.464781343401766e-4),
            (2.397959183673468, -7.745339273024636e-5),
            (2.423469387755101, -3.9621374769542025e-5),
            (2.4489795918367334, -1.9608458558430508e-5),
            (2.474489795918366, -9.38818477949983e-6),
            (2.4999999999999987, -4.3484671364929225e-6),
        ]

        np.testing.assert_allclose(expected, list(zip(freqs, refl_flux)), rtol=1e-6)


if __name__ == '__main__':
    unittest.main()
