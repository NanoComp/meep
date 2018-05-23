from __future__ import division

import unittest
import numpy as np
import meep as mp
from utils import compare_arrays


class TestBendFlux(unittest.TestCase):

    def init(self, no_bend=False):
        sx = 16
        sy = 32
        cell = mp.Vector3(sx, sy, 0)
        pad = 4
        w = 1
        wvg_ycen = -0.5 * (sy - w - (2 * pad))
        wvg_xcen = 0.5 * (sx - w - (2 * pad))

        if no_bend:
            no_bend_vertices = [mp.Vector3(-0.5 * sx, wvg_ycen - 0.5 * w),
                                mp.Vector3(+0.5 * sx, wvg_ycen - 0.5 * w),
                                mp.Vector3(+0.5 * sx, wvg_ycen + 0.5 * w),
                                mp.Vector3(-0.5 * sx, wvg_ycen + 0.5 * w)]

            geometry = [mp.Prism(no_bend_vertices, 0, material=mp.Medium(epsilon=12))]
        else:
            bend_vertices = [mp.Vector3(-0.5 * sx, wvg_ycen - 0.5 * w),
                             mp.Vector3(wvg_xcen + 0.5 * w, wvg_ycen - 0.5 * w),
                             mp.Vector3(wvg_xcen + 0.5 * w, 0.5 * sy),
                             mp.Vector3(wvg_xcen - 0.5 * w, 0.5 * sy),
                             mp.Vector3(wvg_xcen - 0.5 * w, wvg_ycen + 0.5 * w),
                             mp.Vector3(-0.5 * sx, wvg_ycen + 0.5 * w)]

            geometry = [mp.Prism(bend_vertices, 0, material=mp.Medium(epsilon=12))]

        fcen = 0.15
        df = 0.1
        sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ez,
                             center=mp.Vector3(1 + (-0.5 * sx), wvg_ycen), size=mp.Vector3(0, w))]

        pml_layers = [mp.PML(1.0)]
        resolution = 10
        nfreq = 100

        self.sim = mp.Simulation(cell_size=cell,
                                 boundary_layers=pml_layers,
                                 geometry=geometry,
                                 sources=sources,
                                 resolution=resolution)

        if no_bend:
            fr = mp.FluxRegion(center=mp.Vector3((sx / 2) - 1.5, wvg_ycen), size=mp.Vector3(0, w * 2))
        else:
            fr = mp.FluxRegion(center=mp.Vector3(wvg_xcen, (sy / 2) - 1.5), size=mp.Vector3(w * 2, 0))

        self.trans = self.sim.add_flux(fcen, df, nfreq, fr)
        refl_fr = mp.FluxRegion(center=mp.Vector3((-0.5 * sx) + 1.5, wvg_ycen),
                                size=mp.Vector3(0, w * 2))

        self.refl = self.sim.add_flux(fcen, df, nfreq, refl_fr)

        if no_bend:
            self.pt = mp.Vector3((sx / 2) - 1.5, wvg_ycen)
        else:
            self.pt = mp.Vector3(wvg_xcen, (sy / 2) - 1.5)

    def test_bend_flux(self):
        # Normalization run
        self.init(no_bend=True)
        self.sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, self.pt, 1e-3))
        # Save flux data for use in real run below
        fdata = self.sim.get_flux_data(self.refl)

        expected = [
            (0.09999999999999999, 3.759029112132568e-5, 3.808240613412702e-5),
            (0.10101010101010101, 5.718360627275768e-5, 5.787858557476235e-5),
            (0.10202020202020202, 8.626178149785403e-5, 8.723095735323989e-5),
            (0.10303030303030304, 1.2904302255339786e-4, 1.303734865376113e-4),
            (0.10404040404040406, 1.91440977959859e-4, 1.9323229749290062e-4),
            (0.10505050505050507, 2.816640015532028e-4, 2.8401853836185545e-4),
            (0.10606060606060609, 4.1099104864621264e-4, 4.139936984202774e-4),
            (0.10707070707070711, 5.947620615041214e-4, 5.984438388719244e-4),
            (0.10808080808080812, 8.536232917474e-4, 8.579019620596675e-4),
            (0.10909090909090914, 0.0012150684695463483, 0.0012196587415521445),
            (0.11010101010101016, 0.0017153140696826942, 0.0017196000673953989),
            (0.11111111111111117, 0.0024015400405837696, 0.0024044040886885415),
            (0.11212121212121219, 0.0033345151119150893, 0.0033341192817500924),
            (0.11313131313131321, 0.004591607098594056, 0.0045851272360724615),
            (0.11414141414141422, 0.006270152592938964, 0.006253468410682328),
            (0.11515151515151524, 0.008491126571916682, 0.008458475808916224),
            (0.11616161616161626, 0.01140301079799375, 0.011346620857916318),
            (0.11717171717171727, 0.015185711272838182, 0.0150954288696048),
            (0.11818181818181829, 0.020054321141823886, 0.01991726961986203),
            (0.11919191919191931, 0.02626246907430508, 0.026062772997176274),
        ]

        res = list(zip(mp.get_flux_freqs(self.trans), mp.get_fluxes(self.trans), mp.get_fluxes(self.refl)))

        compare_arrays(self, np.array(expected), np.array(res[:20]))

        # Real run
        self.sim = None
        self.init()
        # Load flux data obtained from normalization run
        self.sim.load_minus_flux_data(self.refl, fdata)
        self.sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, self.pt, 1e-3))

        expected = [
            (0.09999999999999999, 2.0170954669155e-5, -7.145323406243266e-6),
            (0.10101010101010101, 3.070001955299509e-5, -1.0941385713179106e-5),
            (0.10202020202020202, 4.626587530507024e-5, -1.6558923559642594e-5),
            (0.10303030303030304, 6.892295567480947e-5, -2.4836579153858833e-5),
            (0.10404040404040406, 1.0141364425635142e-4, -3.699982489006669e-5),
            (0.10505050505050507, 1.474365085880599e-4, -5.473578454006719e-5),
            (0.10606060606060609, 2.1207288913368974e-4, -8.028391176449571e-5),
            (0.10707070707070711, 3.0239431640933715e-4, -1.1664290104872828e-4),
            (0.10808080808080812, 4.2819888830440315e-4, -1.679368899871326e-4),
            (0.10909090909090914, 6.027275359145561e-4, -2.398486110177475e-4),
            (0.11010101010101016, 8.431482588216411e-4, -3.39975596024738e-4),
            (0.11111111111111117, 0.0011706479532335638, -4.78110538203311e-4),
            (0.11212121212121219, 0.0016101926665890486, -6.66666253433329e-4),
            (0.11313131313131321, 0.0021903704321079953, -9.214760453657109e-4),
            (0.11414141414141422, 0.002944040018564378, -0.0012628973860187313),
            (0.11515151515151524, 0.003910493690885223, -0.0017168242626265833),
            (0.11616161616161626, 0.005139265971645378, -0.0023153058291112685),
            (0.11717171717171727, 0.006694611064215856, -0.0030970004778165016),
            (0.11818181818181829, 0.008658482837267306, -0.004108081960927432),
            (0.11919191919191931, 0.0111293966280746, -0.005403865420469877),
        ]
        res = list(zip(mp.get_flux_freqs(self.trans), mp.get_fluxes(self.trans), mp.get_fluxes(self.refl)))

        compare_arrays(self, np.array(expected), np.array(res[:20]))


if __name__ == '__main__':
    unittest.main()
