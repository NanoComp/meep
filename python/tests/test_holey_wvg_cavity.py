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
            0.23460397880408404,
            -0.000313881191161031,
            373.7146178404248,
            5.840690675863821,
            -3.8828446356969186,
            -4.363162282813034,
        ]

        m = h.modes[0]
        res = [m.freq, m.decay, m.Q, abs(m.amp), m.amp.real, m.amp.imag]

        tol = 1e-6 if mp.is_single_precision() else 1e-8
        self.assertClose(expected, res, epsilon=tol)

    def test_transmission_spectrum(self):
        expected = [
            (0.15, 7.641150802552664e-06),
            (0.1504008016032064, 7.4099162291116155e-06),
            (0.1508016032064128, 6.470424836107634e-06),
            (0.15120240480961925, 5.036688729610434e-06),
            (0.15160320641282565, 3.4596969804911187e-06),
            (0.15200400801603206, 2.1579966191987915e-06),
            (0.15240480961923847, 1.5263243157861142e-06),
            (0.15280561122244488, 1.841647499043871e-06),
            (0.15320641282565128, 3.1883721049726715e-06),
            (0.15360721442885772, 5.422004251871202e-06),
            (0.15400801603206413, 8.183334469188715e-06),
            (0.15440881763527053, 1.0964489123629207e-05),
            (0.15480961923847694, 1.3216147259725806e-05),
            (0.15521042084168335, 1.4474552106315964e-05),
            (0.15561122244488979, 1.4480288057536148e-05),
            (0.1560120240480962, 1.326008843753029e-05),
            (0.1564128256513026, 1.1148917390095156e-05),
            (0.156813627254509, 8.741462998544e-06),
            (0.15721442885771542, 6.777759641756439e-06),
            (0.15761523046092185, 5.98361992265642e-06),
            (0.15801603206412826, 6.899190056394783e-06),
            (0.15841683366733467, 9.734963866380887e-06),
            (0.15881763527054107, 1.429195266661998e-05),
            (0.15921843687374748, 1.9971198114791834e-05),
            (0.1596192384769539, 2.5879277159802237e-05),
            (0.16002004008016033, 3.101454849773109e-05),
            (0.16042084168336673, 3.449843135700863e-05),
            (0.16082164328657314, 3.580191484366967e-05),
            (0.16122244488977955, 3.4913677664286064e-05),
            (0.16162324649298596, 3.240453702696152e-05),
            (0.16202404809619236, 2.9362669845099298e-05),
            (0.1624248496993988, 2.720166132620488e-05),
            (0.1628256513026052, 2.737326994739111e-05),
            (0.16322645290581161, 3.104213382715732e-05),
            (0.16362725450901802, 3.8794169323131676e-05),
            (0.16402805611222443, 5.044972561862373e-05),
            (0.16442885771543087, 6.50353631708701e-05),
            (0.16482965931863727, 8.093681722135633e-05),
            (0.16523046092184368, 9.621616842969988e-05),
            (0.1656312625250501, 0.00010903679275505979),
            (0.1660320641282565, 0.00011810937312947968),
            (0.16643286573146293, 0.0001230589288665698),
            (0.16683366733466934, 0.00012462108367484838),
            (0.16723446893787575, 0.00012460575950881206),
            (0.16763527054108215, 0.0001256133166935415),
            (0.16803607214428856, 0.00013054275124547288),
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
