import unittest

from utils import ApproxComparisonTestCase

import meep as mp


class TestHoleyWvgCavity(ApproxComparisonTestCase):
    def setUp(self):
        eps = 13
        self.w = 1.2
        r = 0.36
        d = 1.4
        N = 3
        sy = 6
        pad = 2
        self.dpml = 1
        self.sx = (2 * (pad + self.dpml + N)) + d - 1
        self.fcen = 0.25
        self.df = 0.2
        self.nfreq = 500

        cell = mp.Vector3(self.sx, sy, 0)

        blk = mp.Block(
            size=mp.Vector3(mp.inf, self.w, mp.inf), material=mp.Medium(epsilon=eps)
        )

        geometry = [blk]

        geometry.extend(mp.Cylinder(r, center=mp.Vector3(d / 2 + i)) for i in range(3))
        for i in range(3):
            geometry.append(mp.Cylinder(r, center=mp.Vector3(d / -2 - i)))

        self.sim = mp.Simulation(
            cell_size=cell,
            geometry=geometry,
            sources=[],
            boundary_layers=[mp.PML(self.dpml)],
            resolution=20,
        )

    @classmethod
    def setUpClass(cls):
        cls.temp_dir = mp.make_output_directory()

    @classmethod
    def tearDownClass(cls):
        mp.delete_directory(cls.temp_dir)

    def test_resonant_modes(self):
        self.sim.sources = [
            mp.Source(mp.GaussianSource(self.fcen, fwidth=self.df), mp.Hz, mp.Vector3())
        ]

        self.sim.symmetries = [mp.Mirror(mp.Y, phase=-1), mp.Mirror(mp.X, phase=-1)]

        self.sim.use_output_directory(self.temp_dir)
        h = mp.Harminv(mp.Hz, mp.Vector3(), self.fcen, self.df)
        self.sim.run(
            mp.at_beginning(mp.output_epsilon),
            mp.after_sources(h),
            until_after_sources=400,
        )

        expected = [
            0.23445415346009466,
            -3.147812367338531e-4,
            372.40808234438254,
            5.8121430334347135,
            -3.763107485715599,
            -4.429450156854109,
        ]

        m = h.modes[0]
        res = [m.freq, m.decay, m.Q, abs(m.amp), m.amp.real, m.amp.imag]

        tol = 1e-6 if mp.is_single_precision() else 1e-8
        self.assertClose(expected, res, epsilon=tol)

    def test_transmission_spectrum(self):
        expected = [
            (0.15, 7.218492264696595e-6),
            (0.1504008016032064, 6.445696315927592e-6),
            (0.1508016032064128, 5.140949243632777e-6),
            (0.15120240480961922, 3.6159747936427164e-6),
            (0.15160320641282563, 2.263940553705969e-6),
            (0.15200400801603203, 1.4757165844336744e-6),
            (0.15240480961923844, 1.5491803919142815e-6),
            (0.15280561122244485, 2.612053246626972e-6),
            (0.15320641282565126, 4.577504371188737e-6),
            (0.15360721442885766, 7.1459089162998185e-6),
            (0.15400801603206407, 9.856622013418823e-6),
            (0.15440881763527048, 1.2182309227954296e-5),
            (0.1548096192384769, 1.3647726444709649e-5),
            (0.1552104208416833, 1.3947420613633674e-5),
            (0.1556112224448897, 1.303466755716231e-5),
            (0.1560120240480961, 1.115807915037775e-5),
            (0.15641282565130252, 8.832335196969796e-6),
            (0.15681362725450892, 6.743645773127985e-6),
            (0.15721442885771533, 5.605913756087576e-6),
            (0.15761523046092174, 5.996668564026961e-6),
            (0.15801603206412815, 8.209400611614078e-6),
            (0.15841683366733456, 1.2158641936828497e-5),
            (0.15881763527054096, 1.73653230513453e-5),
            (0.15921843687374737, 2.303382576477893e-5),
            (0.15961923847695378, 2.821180350795834e-5),
            (0.1600200400801602, 3.200359292911769e-5),
            (0.1604208416833666, 3.3792624373001934e-5),
            (0.160821643286573, 3.342171394788991e-5),
            (0.1612224448897794, 3.1284866146526904e-5),
            (0.16162324649298582, 2.830022088581398e-5),
            (0.16202404809619222, 2.5758413657344014e-5),
            (0.16242484969939863, 2.506899997971769e-5),
            (0.16282565130260504, 2.7453508915303887e-5),
            (0.16322645290581145, 3.365089813497114e-5),
            (0.16362725450901786, 4.370486834112e-5),
            (0.16402805611222426, 5.689050715055283e-5),
            (0.16442885771543067, 7.181133157470506e-5),
            (0.16482965931863708, 8.666168027415369e-5),
            (0.16523046092184349, 9.961094123261317e-5),
            (0.1656312625250499, 1.0923388232657953e-4),
            (0.1660320641282563, 1.1489334204708105e-4),
            (0.1664328657314627, 1.1698318060032011e-4),
            (0.16683366733466912, 1.169621456132733e-4),
            (0.16723446893787552, 1.1714995241571987e-4),
            (0.16763527054108193, 1.2030783847222252e-4),
            (0.16803607214428834, 1.2907652919660887e-4),
        ]

        self.sim.sources = [
            mp.Source(
                mp.GaussianSource(self.fcen, fwidth=self.df),
                mp.Ey,
                mp.Vector3(self.dpml + (-0.5 * self.sx)),
                size=mp.Vector3(0, self.w),
            )
        ]

        self.sim.symmetries = [mp.Mirror(mp.Y, phase=-1)]

        freg = mp.FluxRegion(
            center=mp.Vector3((0.5 * self.sx) - self.dpml - 0.5),
            size=mp.Vector3(0, 2 * self.w),
        )

        trans = self.sim.add_flux(
            self.fcen, self.df, self.nfreq, freg, decimation_factor=1
        )

        self.sim.run(
            until_after_sources=mp.stop_when_fields_decayed(
                50, mp.Ey, mp.Vector3((0.5 * self.sx) - self.dpml - 0.5, 0), 1e-1
            )
        )

        res = zip(mp.get_flux_freqs(trans), mp.get_fluxes(trans))

        tol = 1e-8 if mp.is_single_precision() else 1e-10
        for e, r in zip(expected, res):
            self.assertClose(e, r, epsilon=tol)


if __name__ == "__main__":
    unittest.main()
