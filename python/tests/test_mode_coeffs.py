import unittest

import numpy as np

import meep as mp


class TestModeCoeffs(unittest.TestCase):
    def run_mode_coeffs(self, mode_num, kpoint_func, nf=1, resolution=15):

        w = 1  # width of waveguide
        L = 10  # length of waveguide

        Si = mp.Medium(epsilon=12.0)

        dair = 3.0
        dpml = 3.0

        sx = dpml + L + dpml
        sy = dpml + dair + w + dair + dpml
        cell_size = mp.Vector3(sx, sy, 0)

        prism_x = sx + 1
        prism_y = w / 2
        vertices = [
            mp.Vector3(-prism_x, prism_y),
            mp.Vector3(prism_x, prism_y),
            mp.Vector3(prism_x, -prism_y),
            mp.Vector3(-prism_x, -prism_y),
        ]

        geometry = [mp.Prism(vertices, height=mp.inf, material=Si)]

        boundary_layers = [mp.PML(dpml)]

        # mode frequency
        fcen = 0.20  # > 0.5/sqrt(11) to have at least 2 modes
        df = 0.5 * fcen

        source = mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            eig_band=mode_num,
            size=mp.Vector3(0, sy - 2 * dpml, 0),
            center=mp.Vector3(-0.5 * sx + dpml, 0, 0),
        )

        symmetries = [mp.Mirror(mp.Y, phase=1 if mode_num % 2 == 1 else -1)]

        sim = mp.Simulation(
            resolution=resolution,
            cell_size=cell_size,
            boundary_layers=boundary_layers,
            geometry=geometry,
            sources=[source],
            symmetries=symmetries,
        )

        xm = 0.5 * sx - dpml  # x-coordinate of monitor
        mflux = sim.add_mode_monitor(
            fcen,
            df,
            nf,
            mp.ModeRegion(center=mp.Vector3(xm, 0), size=mp.Vector3(0, sy - 2 * dpml)),
        )
        mode_flux = sim.add_flux(
            fcen,
            df,
            nf,
            mp.FluxRegion(center=mp.Vector3(xm, 0), size=mp.Vector3(0, sy - 2 * dpml)),
        )

        # sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(-0.5*sx+dpml,0), 1e-10))
        sim.run(until_after_sources=100)

        ##################################################
        # If the number of analysis frequencies is >1, we
        # are testing the unit-power normalization
        # of the eigenmode source: we observe the total
        # power flux through the mode_flux monitor (which
        # equals the total power emitted by the source as
        # there is no scattering in this ideal waveguide)
        # and check that it agrees with the prediction
        # of the eig_power() class method in EigenmodeSource.
        ##################################################
        if nf > 1:
            power_observed = mp.get_fluxes(mode_flux)
            freqs = mp.get_flux_freqs(mode_flux)
            power_expected = [source.eig_power(f) for f in freqs]
            return freqs, power_expected, power_observed

        modes_to_check = [
            1,
            2,
        ]  # indices of modes for which to compute expansion coefficients
        res = sim.get_eigenmode_coefficients(
            mflux, modes_to_check, kpoint_func=kpoint_func
        )

        self.assertTrue(res.kpoints[0].close(mp.Vector3(0.604301, 0, 0)))
        self.assertTrue(res.kpoints[1].close(mp.Vector3(0.494353, 0, 0), tol=1e-2))
        self.assertTrue(res.kdom[0].close(mp.Vector3(0.604301, 0, 0)))
        self.assertTrue(res.kdom[1].close(mp.Vector3(0.494353, 0, 0), tol=1e-2))
        self.assertAlmostEqual(res.cscale[0], 0.50000977, places=5)
        self.assertAlmostEqual(res.cscale[1], 0.50096888, places=5)
        mode_power = mp.get_fluxes(mode_flux)[0]

        TestPassed = True
        TOLERANCE = 5.0e-3
        c0 = res.alpha[
            mode_num - 1, 0, 0
        ]  # coefficient of forward-traveling wave for mode #mode_num
        for nm in range(1, len(modes_to_check) + 1):
            if nm != mode_num:
                cfrel = np.abs(res.alpha[nm - 1, 0, 0]) / np.abs(c0)
                cbrel = np.abs(res.alpha[nm - 1, 0, 1]) / np.abs(c0)
                if cfrel > TOLERANCE or cbrel > TOLERANCE:
                    TestPassed = False

        self.sim = sim

        # test 1: coefficient of excited mode >> coeffs of all other modes
        self.assertTrue(TestPassed, msg=f"cfrel: {cfrel}, cbrel: {cbrel}")
        # test 2: |mode coeff|^2 = power
        self.assertAlmostEqual(mode_power / abs(c0**2), 1.0, places=1)

        return res

    def test_modes(self):
        self.run_mode_coeffs(1, None)
        res = self.run_mode_coeffs(2, None)

        # Test mp.get_eigenmode and EigenmodeData
        vol = mp.Volume(center=mp.Vector3(5), size=mp.Vector3(y=7))
        emdata = self.sim.get_eigenmode(0.2, mp.X, vol, 2, mp.Vector3())
        self.assertEqual(emdata.freq, 0.2)
        self.assertEqual(emdata.band_num, 2)
        self.assertTrue(emdata.kdom.close(res.kdom[1]))

        eval_point = mp.Vector3(0.7, -0.2, 0.3)
        ex_at_eval_point = emdata.amplitude(eval_point, mp.Ex)
        hz_at_eval_point = emdata.amplitude(eval_point, mp.Hz)

        places = 5 if mp.is_single_precision() else 7
        self.assertAlmostEqual(
            ex_at_eval_point, 0.4887779638178009 + 0.484240145324284j, places=places
        )
        self.assertAlmostEqual(
            hz_at_eval_point, 3.4249236584603495 - 3.455974863884166j, places=places
        )

    def test_kpoint_func(self):
        def kpoint_func(freq, mode):
            return mp.Vector3()

        self.run_mode_coeffs(1, kpoint_func)

    def test_eigensource_normalization(self):
        f, p_exp, p_obs = self.run_mode_coeffs(1, None, nf=51, resolution=15)
        # self.assertAlmostEqual(max(p_exp),max(p_obs),places=1)
        max_exp, max_obs = max(p_exp), max(p_obs)
        self.assertLess(abs(max_exp - max_obs), 0.5 * max(abs(max_exp), abs(max_obs)))

    def test_reciprocity_kpoint(self):
        resolution = 40

        sx = 7.0
        sy = 5.0
        cell_size = mp.Vector3(sx, sy)

        dpml = 1.0
        pml_layers = [mp.PML(thickness=dpml)]

        w = 1.0
        geometry = [
            mp.Block(
                center=mp.Vector3(),
                size=mp.Vector3(mp.inf, w, mp.inf),
                material=mp.Medium(epsilon=12),
            )
        ]

        fsrc = 0.15
        sources = [
            mp.EigenModeSource(
                src=mp.GaussianSource(fsrc, fwidth=0.2 * fsrc),
                center=mp.Vector3(x=-0.5 * sx + dpml),
                size=mp.Vector3(y=sy),
                eig_parity=mp.EVEN_Y + mp.ODD_Z,
            )
        ]

        symmetries = [mp.Mirror(mp.Y)]

        sim = mp.Simulation(
            cell_size=cell_size,
            resolution=resolution,
            boundary_layers=pml_layers,
            sources=sources,
            geometry=geometry,
            symmetries=symmetries,
        )

        tran = sim.add_mode_monitor(
            fsrc,
            0,
            1,
            mp.ModeRegion(center=mp.Vector3(x=0.5 * sx - dpml), size=mp.Vector3(y=sy)),
            yee_grid=False,
        )

        sim.run(until_after_sources=50)

        res_fwd = sim.get_eigenmode_coefficients(
            tran,
            [1],
            eig_parity=mp.EVEN_Y + mp.ODD_Z,
            direction=mp.NO_DIRECTION,
            kpoint_func=lambda f, n: mp.Vector3(+1, 0, 0),
        )

        res_bwd = sim.get_eigenmode_coefficients(
            tran,
            [1],
            eig_parity=mp.EVEN_Y + mp.ODD_Z,
            direction=mp.NO_DIRECTION,
            kpoint_func=lambda f, n: mp.Vector3(-1, 0, 0),
        )

        print(f"S11:, {res_fwd.alpha[0,0,1]}, {res_bwd.alpha[0,0,0]}")
        print(f"S21:, {res_fwd.alpha[0,0,0]}, {res_bwd.alpha[0,0,1]}")

        # |S11|^2
        self.assertAlmostEqual(
            abs(res_fwd.alpha[0, 0, 1]) ** 2, abs(res_bwd.alpha[0, 0, 0]) ** 2, places=4
        )

        # |S21|^2
        self.assertAlmostEqual(
            abs(res_fwd.alpha[0, 0, 0]) ** 2 / abs(res_bwd.alpha[0, 0, 1]) ** 2,
            1.00,
            places=2,
        )

    def test_reciprocity_monitor(self):
        resolution = 25

        sx = 7.0
        sy = 5.0
        cell_size = mp.Vector3(sx, sy)

        dpml = 1.0
        pml_layers = [mp.PML(thickness=dpml)]

        w = 1.0
        geometry = [
            mp.Block(
                center=mp.Vector3(),
                size=mp.Vector3(mp.inf, w, mp.inf),
                material=mp.Medium(epsilon=12),
            )
        ]

        fsrc = 0.15

        # source is at the left edge of the waveguide
        sources = [
            mp.EigenModeSource(
                src=mp.GaussianSource(fsrc, fwidth=0.2 * fsrc),
                center=mp.Vector3(x=-0.5 * sx + dpml),
                size=mp.Vector3(y=sy),
                eig_parity=mp.EVEN_Y + mp.ODD_Z,
            )
        ]

        symmetries = [mp.Mirror(mp.Y)]

        sim = mp.Simulation(
            cell_size=cell_size,
            resolution=resolution,
            boundary_layers=pml_layers,
            sources=sources,
            geometry=geometry,
            symmetries=symmetries,
        )

        # monitor is at the right edge of the waveguide
        tran = sim.add_mode_monitor(
            fsrc,
            0,
            1,
            mp.ModeRegion(center=mp.Vector3(x=0.5 * sx - dpml), size=mp.Vector3(y=sy)),
            yee_grid=False,
        )

        sim.run(until_after_sources=50)

        res_fwd = sim.get_eigenmode_coefficients(
            tran, [1], eig_parity=mp.EVEN_Y + mp.ODD_Z
        )

        print(f"S11:, {res_fwd.alpha[0,0,1]}")
        print(f"S21:, {res_fwd.alpha[0,0,0]}")

        sim.reset_meep()

        # source is at the right edge of the waveguide
        sources = [
            mp.EigenModeSource(
                src=mp.GaussianSource(fsrc, fwidth=0.2 * fsrc),
                center=mp.Vector3(x=0.5 * sx - dpml),
                size=mp.Vector3(y=sy),
                direction=mp.NO_DIRECTION,
                eig_kpoint=mp.Vector3(-1, 0, 0),
                eig_parity=mp.EVEN_Y + mp.ODD_Z,
            )
        ]

        sim = mp.Simulation(
            cell_size=cell_size,
            resolution=resolution,
            boundary_layers=pml_layers,
            sources=sources,
            geometry=geometry,
            symmetries=symmetries,
        )

        # monitor is at the left edge of the waveguide
        tran = sim.add_mode_monitor(
            fsrc,
            0,
            1,
            mp.ModeRegion(center=mp.Vector3(x=-0.5 * sx + dpml), size=mp.Vector3(y=sy)),
            yee_grid=False,
        )

        sim.run(until_after_sources=50)

        res_bwd = sim.get_eigenmode_coefficients(
            tran, [1], eig_parity=mp.EVEN_Y + mp.ODD_Z
        )

        print(f"S12:, {res_bwd.alpha[0,0,1]}")
        print(f"S22:, {res_bwd.alpha[0,0,0]}")

        # |S21|^2 = |S12|^2
        self.assertAlmostEqual(
            abs(res_fwd.alpha[0, 0, 0]) ** 2 / abs(res_bwd.alpha[0, 0, 1]) ** 2,
            1.00,
            places=2,
        )

        # |S11|^2 = |S22|^2
        self.assertAlmostEqual(
            abs(res_fwd.alpha[0, 0, 1]) ** 2, abs(res_bwd.alpha[0, 0, 0]) ** 2, places=2
        )


if __name__ == "__main__":
    unittest.main()
