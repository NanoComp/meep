import os
import unittest

import numpy as np
from utils import ApproxComparisonTestCase

import meep as mp


class TestBendFlux(ApproxComparisonTestCase):
    def init(self, no_bend=False, gdsii=False):
        sx = 16
        sy = 32
        cell = mp.Vector3(sx, sy, 0)
        pad = 4
        w = 1
        wvg_ycen = -0.5 * (sy - w - (2 * pad))
        wvg_xcen = 0.5 * (sx - w - (2 * pad))
        height = mp.inf
        data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "data"))
        gdsii_file = os.path.join(data_dir, "bend-flux.gds")

        if no_bend:
            if gdsii:
                geometry = mp.get_GDSII_prisms(
                    mp.Medium(epsilon=12), gdsii_file, 1, 0, height
                )
            else:
                no_bend_vertices = [
                    mp.Vector3(-0.5 * sx - 5, wvg_ycen - 0.5 * w),
                    mp.Vector3(+0.5 * sx + 5, wvg_ycen - 0.5 * w),
                    mp.Vector3(+0.5 * sx + 5, wvg_ycen + 0.5 * w),
                    mp.Vector3(-0.5 * sx - 5, wvg_ycen + 0.5 * w),
                ]

                geometry = [
                    mp.Prism(no_bend_vertices, height, material=mp.Medium(epsilon=12))
                ]
        elif gdsii:
            geometry = mp.get_GDSII_prisms(
                mp.Medium(epsilon=12), gdsii_file, 2, 0, height
            )
        else:
            bend_vertices = [
                mp.Vector3(-0.5 * sx, wvg_ycen - 0.5 * w),
                mp.Vector3(wvg_xcen + 0.5 * w, wvg_ycen - 0.5 * w),
                mp.Vector3(wvg_xcen + 0.5 * w, 0.5 * sy),
                mp.Vector3(wvg_xcen - 0.5 * w, 0.5 * sy),
                mp.Vector3(wvg_xcen - 0.5 * w, wvg_ycen + 0.5 * w),
                mp.Vector3(-0.5 * sx, wvg_ycen + 0.5 * w),
            ]

            geometry = [mp.Prism(bend_vertices, height, material=mp.Medium(epsilon=12))]

        fcen = 0.15
        df = 0.1
        sources = [
            mp.Source(
                mp.GaussianSource(fcen, fwidth=df),
                component=mp.Ez,
                center=mp.Vector3(1 + (-0.5 * sx), wvg_ycen),
                size=mp.Vector3(0, w),
            )
        ]

        pml_layers = [mp.PML(1.0)]
        resolution = 10
        nfreq = 100

        self.sim = mp.Simulation(
            cell_size=cell,
            boundary_layers=pml_layers,
            geometry=geometry,
            sources=sources,
            resolution=resolution,
        )

        if no_bend:
            fr = mp.FluxRegion(
                center=mp.Vector3((sx / 2) - 1.5, wvg_ycen), size=mp.Vector3(0, w * 2)
            )
        else:
            fr = mp.FluxRegion(
                center=mp.Vector3(wvg_xcen, (sy / 2) - 1.5), size=mp.Vector3(w * 2, 0)
            )

        self.trans = self.sim.add_flux(fcen, df, nfreq, fr, decimation_factor=1)
        self.trans_decimated = self.sim.add_flux(
            fcen, df, nfreq, fr, decimation_factor=5
        )

        refl_fr = mp.FluxRegion(
            center=mp.Vector3((-0.5 * sx) + 1.5, wvg_ycen), size=mp.Vector3(0, w * 2)
        )
        self.refl = self.sim.add_flux(
            np.linspace(fcen - 0.5 * df, fcen + 0.5 * df, nfreq),
            refl_fr,
            decimation_factor=1,
        )
        self.refl_decimated = self.sim.add_flux(
            np.linspace(fcen - 0.5 * df, fcen + 0.5 * df, nfreq),
            refl_fr,
            decimation_factor=10,
        )

        if no_bend:
            self.pt = mp.Vector3((sx / 2) - 1.5, wvg_ycen)
        else:
            self.pt = mp.Vector3(wvg_xcen, (sy / 2) - 1.5)

    def run_bend_flux(self, from_gdsii_file):
        # Normalization run
        self.init(no_bend=True, gdsii=from_gdsii_file)
        self.sim.run(until_after_sources=mp.stop_when_energy_decayed(100, 1e-3))
        # Save flux data for use in real run below
        fdata = self.sim.get_flux_data(self.refl)
        fdata_decimated = self.sim.get_flux_data(self.refl_decimated)

        expected = [
            (0.1, 3.65231563251e-05, 3.68932495077e-05),
            (0.10101010101, 5.55606718876e-05, 5.6065539588e-05),
            (0.10202020202, 8.38211697478e-05, 8.44909864736e-05),
            (0.10303030303, 0.000125411162229, 0.000126268639045),
            (0.10404040404, 0.000186089117531, 0.000187135303398),
            (0.105050505051, 0.000273848867869, 0.000275039134667),
            (0.106060606061, 0.000399674037745, 0.000400880269423),
            (0.107070707071, 0.00057849953593, 0.000579454087881),
            (0.108080808081, 0.000830418432986, 0.000830635406881),
            (0.109090909091, 0.00118217282661, 0.00118084271347),
            (0.110101010101, 0.00166896468348, 0.00166481944189),
            (0.111111111111, 0.00233661613864, 0.00232776318321),
            (0.112121212121, 0.00324409729096, 0.00322782257917),
            (0.113131313131, 0.00446642217385, 0.00443896468822),
            (0.114141414141, 0.0060978895019, 0.0060541922825),
            (0.115151515152, 0.00825561352398, 0.00818906047274),
            (0.116161616162, 0.0110832518495, 0.010985404883),
            (0.117171717172, 0.0147547920552, 0.0146151488236),
            (0.118181818182, 0.0194782085272, 0.0192840042241),
            (0.119191919192, 0.0254987474079, 0.0252348211592),
        ]

        res = list(
            zip(
                mp.get_flux_freqs(self.trans),
                mp.get_fluxes(self.trans),
                mp.get_fluxes(self.refl),
            )
        )

        res_decimated = list(
            zip(
                mp.get_flux_freqs(self.trans_decimated),
                mp.get_fluxes(self.trans_decimated),
                mp.get_fluxes(self.refl_decimated),
            )
        )

        tol = 1e-6 if mp.is_single_precision() else 1e-8
        self.assertClose(np.array(expected), np.array(res[:20]), epsilon=tol)
        self.assertClose(np.array(expected), np.array(res_decimated[:20]), epsilon=tol)

        # Real run
        self.sim = None
        self.init(gdsii=from_gdsii_file)
        # Load flux data obtained from normalization run
        self.sim.load_minus_flux_data(self.refl, fdata)
        self.sim.load_minus_flux_data(self.refl_decimated, fdata_decimated)
        self.sim.run(until_after_sources=mp.stop_when_energy_decayed(100, 1e-3))

        expected = [
            (0.09999999999999999, 1.8392235204829767e-5, -7.259467687598002e-6),
            (0.10101010101010101, 2.7629932558236724e-5, -1.1107162110079347e-5),
            (0.10202020202020202, 4.1001228946782745e-5, -1.687561915798036e-5),
            (0.10303030303030304, 6.018966076122556e-5, -2.5425779493709066e-5),
            (0.10404040404040406, 8.758554071933231e-5, -3.794958119189475e-5),
            (0.10505050505050507, 1.2656696778129198e-4, -5.612512808928115e-5),
            (0.10606060606060609, 1.817948859871414e-4, -8.232188174309142e-5),
            (0.10707070707070711, 2.594514094902856e-4, -1.1981531280672989e-4),
            (0.10808080808080812, 3.6736164837695035e-4, -1.7300125173897737e-4),
            (0.10909090909090914, 5.150131339048232e-4, -2.476730940385436e-4),
            (0.11010101010101016, 7.136181099374187e-4, -3.5145561406042276e-4),
            (0.11111111111111117, 9.76491765781944e-4, -4.944142331545938e-4),
            (0.11212121212121219, 0.001320033637882244, -6.897357105189368e-4),
            (0.11313131313131321, 0.0017653940714397098, -9.543556354451615e-4),
            (0.11414141414141422, 0.0023404727796352857, -0.0013095604571818236),
            (0.11515151515151524, 0.0030813962415392098, -0.00178176942635486),
            (0.11616161616161626, 0.00403238648982478, -0.0024036650652026112),
            (0.11717171717171727, 0.005243320443599316, -0.003215529845495731),
            (0.11818181818181829, 0.0067654019326068, -0.004266367104375331),
            (0.11919191919191931, 0.008646855439680507, -0.005614491919262783),
        ]

        res = list(
            zip(
                mp.get_flux_freqs(self.trans),
                mp.get_fluxes(self.trans),
                mp.get_fluxes(self.refl),
            )
        )

        res_decimated = list(
            zip(
                mp.get_flux_freqs(self.trans_decimated),
                mp.get_fluxes(self.trans_decimated),
                mp.get_fluxes(self.refl_decimated),
            )
        )

        tol = 1e-3
        self.assertClose(np.array(expected), np.array(res[:20]), epsilon=tol)
        self.assertClose(np.array(expected), np.array(res_decimated[:20]), epsilon=tol)

    def test_bend_flux(self):
        self.run_bend_flux(False)
        if mp.with_libGDSII():
            self.run_bend_flux(True)


if __name__ == "__main__":
    unittest.main()
