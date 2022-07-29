import cmath
import math
import unittest

import numpy as np
import parameterized

import meep as mp


class TestEigCoeffs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.resolution = 30  # pixels/μm

        cls.dpml = 1.0  # PML thickness
        cls.dsub = 1.0  # substrate thickness
        cls.dpad = 1.0  # padding thickness between grating and PML
        cls.gp = 6.0  # grating period
        cls.gh = 0.5  # grating height
        cls.gdc = 0.5  # grating duty cycle

        cls.sx = cls.dpml + cls.dsub + cls.gh + cls.dpad + cls.dpml
        cls.sy = cls.gp

        cls.cell_size = mp.Vector3(cls.sx, cls.sy, 0)

        # replace anisotropic PML with isotropic Absorber to
        # attenuate parallel-directed fields of oblique source
        cls.abs_layers = [mp.Absorber(thickness=cls.dpml, direction=mp.X)]

        wvl = 0.5  # center wavelength
        cls.fcen = 1 / wvl  # center frequency
        cls.df = 0.05 * cls.fcen  # frequency width

        cls.ng = 1.5
        cls.glass = mp.Medium(index=cls.ng)

        cls.geometry = [
            mp.Block(
                material=cls.glass,
                size=mp.Vector3(cls.dpml + cls.dsub, mp.inf, mp.inf),
                center=mp.Vector3(-0.5 * cls.sx + 0.5 * (cls.dpml + cls.dsub), 0, 0),
            ),
            mp.Block(
                material=cls.glass,
                size=mp.Vector3(cls.gh, cls.gdc * cls.gp, mp.inf),
                center=mp.Vector3(
                    -0.5 * cls.sx + cls.dpml + cls.dsub + 0.5 * cls.gh, 0, 0
                ),
            ),
        ]

    @parameterized.parameterized.expand([(0,), (10.7,)])
    def test_binary_grating_oblique(self, theta):
        # rotation angle of incident planewave
        # counterclockwise (CCW) about Z axis, 0 degrees along +X axis
        theta_in = math.radians(theta)

        # k (in source medium) with correct length (plane of incidence: XY)
        k = mp.Vector3(self.fcen * self.ng).rotate(mp.Vector3(0, 0, 1), theta_in)

        symmetries = []
        eig_parity = mp.ODD_Z
        if theta_in == 0:
            symmetries = [mp.Mirror(mp.Y)]
            eig_parity += mp.EVEN_Y

        def pw_amp(k, x0):
            def _pw_amp(x):
                return cmath.exp(1j * 2 * math.pi * k.dot(x + x0))

            return _pw_amp

        src_pt = mp.Vector3(-0.5 * self.sx + self.dpml, 0, 0)
        sources = [
            mp.Source(
                mp.GaussianSource(self.fcen, fwidth=self.df),
                component=mp.Ez,  # S polarization
                center=src_pt,
                size=mp.Vector3(0, self.sy, 0),
                amp_func=pw_amp(k, src_pt),
            )
        ]

        sim = mp.Simulation(
            resolution=self.resolution,
            cell_size=self.cell_size,
            boundary_layers=self.abs_layers,
            k_point=k,
            default_material=self.glass,
            sources=sources,
            symmetries=symmetries,
        )

        refl_pt = mp.Vector3(-0.5 * self.sx + self.dpml + 0.5 * self.dsub, 0, 0)
        refl_flux = sim.add_flux(
            self.fcen,
            0,
            1,
            mp.FluxRegion(center=refl_pt, size=mp.Vector3(0, self.sy, 0)),
        )

        sim.run(until_after_sources=mp.stop_when_dft_decayed())

        input_flux = mp.get_fluxes(refl_flux)
        input_flux_data = sim.get_flux_data(refl_flux)

        sim.reset_meep()

        sim = mp.Simulation(
            resolution=self.resolution,
            cell_size=self.cell_size,
            boundary_layers=self.abs_layers,
            geometry=self.geometry,
            k_point=k,
            sources=sources,
            symmetries=symmetries,
        )

        refl_flux = sim.add_flux(
            self.fcen,
            0,
            1,
            mp.FluxRegion(center=refl_pt, size=mp.Vector3(0, self.sy, 0)),
        )

        sim.load_minus_flux_data(refl_flux, input_flux_data)

        tran_pt = mp.Vector3(0.5 * self.sx - self.dpml - 0.5 * self.dpad, 0, 0)
        tran_flux = sim.add_flux(
            self.fcen,
            0,
            1,
            mp.FluxRegion(center=tran_pt, size=mp.Vector3(0, self.sy, 0)),
        )

        sim.run(until_after_sources=mp.stop_when_dft_decayed())

        # number of reflected orders
        nm_r = np.floor((self.fcen * self.ng - k.y) * self.gp) - np.ceil(
            (-self.fcen * self.ng - k.y) * self.gp
        )
        if theta_in == 0:
            nm_r = nm_r / 2  # since eig_parity removes degeneracy in y-direction
        nm_r = int(nm_r)

        res = sim.get_eigenmode_coefficients(
            refl_flux, range(1, nm_r + 1), eig_parity=eig_parity
        )
        r_coeffs = res.alpha

        Rsum = 0
        for nm in range(nm_r):
            Rsum += abs(r_coeffs[nm, 0, 1]) ** 2 / input_flux[0]

        # number of transmitted orders
        nm_t = np.floor((self.fcen - k.y) * self.gp) - np.ceil(
            (-self.fcen - k.y) * self.gp
        )
        if theta_in == 0:
            nm_t = nm_t / 2  # since eig_parity removes degeneracy in y-direction
        nm_t = int(nm_t)

        res = sim.get_eigenmode_coefficients(
            tran_flux, range(1, nm_t + 1), eig_parity=eig_parity
        )
        t_coeffs = res.alpha

        Tsum = 0
        for nm in range(nm_t):
            Tsum += abs(t_coeffs[nm, 0, 0]) ** 2 / input_flux[0]

        r_flux = mp.get_fluxes(refl_flux)
        t_flux = mp.get_fluxes(tran_flux)
        Rflux = -r_flux[0] / input_flux[0]
        Tflux = t_flux[0] / input_flux[0]

        print(f"refl:, {Rsum}, {Rflux}")
        print(f"tran:, {Tsum}, {Tflux}")
        print(f"sum:,  {Rsum + Tsum}, {Rflux + Tflux}")

        self.assertAlmostEqual(Rsum, Rflux, places=2)
        self.assertAlmostEqual(Tsum, Tflux, places=2)
        self.assertAlmostEqual(Rsum + Tsum, 1.00, places=2)

    @parameterized.parameterized.expand(
        [(13.2, "real/imag"), (17.7, "complex"), (21.2, "3d")]
    )
    def test_binary_grating_special_kz(self, theta, kz_2d):
        # rotation angle of incident planewave
        # counterclockwise (CCW) about Y axis, 0 degrees along +X axis
        theta_in = math.radians(theta)

        # k (in source medium) with correct length (plane of incidence: XZ)
        k = mp.Vector3(self.fcen * self.ng).rotate(mp.Vector3(0, 1, 0), theta_in)

        symmetries = [mp.Mirror(mp.Y)]

        def pw_amp(k, x0):
            def _pw_amp(x):
                return cmath.exp(1j * 2 * math.pi * k.dot(x + x0))

            return _pw_amp

        src_pt = mp.Vector3(-0.5 * self.sx + self.dpml, 0, 0)
        sources = [
            mp.Source(
                mp.GaussianSource(self.fcen, fwidth=self.df),
                component=mp.Ez,
                center=src_pt,
                size=mp.Vector3(0, self.sy, 0),
                amp_func=pw_amp(k, src_pt),
            )
        ]

        sim = mp.Simulation(
            resolution=self.resolution,
            cell_size=self.cell_size,
            boundary_layers=self.abs_layers,
            k_point=k,
            default_material=self.glass,
            sources=sources,
            symmetries=symmetries,
            kz_2d=kz_2d,
        )

        refl_pt = mp.Vector3(-0.5 * self.sx + self.dpml + 0.5 * self.dsub, 0, 0)
        refl_flux = sim.add_mode_monitor(
            self.fcen,
            0,
            1,
            mp.FluxRegion(center=refl_pt, size=mp.Vector3(0, self.sy, 0)),
        )

        sim.run(until_after_sources=mp.stop_when_dft_decayed())

        input_flux = mp.get_fluxes(refl_flux)
        input_flux_data = sim.get_flux_data(refl_flux)

        sim.reset_meep()

        sim = mp.Simulation(
            resolution=self.resolution,
            cell_size=self.cell_size,
            boundary_layers=self.abs_layers,
            geometry=self.geometry,
            k_point=k,
            sources=sources,
            symmetries=symmetries,
            kz_2d=kz_2d,
        )

        refl_flux = sim.add_mode_monitor(
            self.fcen,
            0,
            1,
            mp.FluxRegion(center=refl_pt, size=mp.Vector3(0, self.sy, 0)),
        )

        sim.load_minus_flux_data(refl_flux, input_flux_data)

        tran_pt = mp.Vector3(0.5 * self.sx - self.dpml - 0.5 * self.dpad, 0, 0)
        tran_flux = sim.add_mode_monitor(
            self.fcen,
            0,
            1,
            mp.FluxRegion(center=tran_pt, size=mp.Vector3(0, self.sy, 0)),
        )

        sim.run(until_after_sources=mp.stop_when_dft_decayed())

        # number of reflected orders
        nm_r = np.ceil(
            (np.sqrt((self.fcen * self.ng) ** 2 - k.z**2) - k.y) * self.gp
        ) - np.floor((-np.sqrt((self.fcen * self.ng) ** 2 - k.z**2) - k.y) * self.gp)
        nm_r = int(nm_r / 2)

        Rsum = 0
        for nm in range(nm_r):
            for S_pol in [False, True]:
                res = sim.get_eigenmode_coefficients(
                    refl_flux,
                    mp.DiffractedPlanewave(
                        [0, nm, 0],
                        mp.Vector3(1, 0, 0),
                        1 if S_pol else 0,
                        0 if S_pol else 1,
                    ),
                )
                r_coeffs = res.alpha
                Rmode = abs(r_coeffs[0, 0, 1]) ** 2 / input_flux[0]
                print(
                    "refl-order:, {}, {}, {}".format("s" if S_pol else "p", nm, Rmode)
                )
                Rsum += Rmode if nm == 0 else 2 * Rmode

        # number of transmitted orders
        nm_t = np.ceil((np.sqrt(self.fcen**2 - k.z**2) - k.y) * self.gp) - np.floor(
            (-np.sqrt(self.fcen**2 - k.z**2) - k.y) * self.gp
        )
        nm_t = int(nm_t / 2)

        Tsum = 0
        for nm in range(nm_t):
            for S_pol in [False, True]:
                res = sim.get_eigenmode_coefficients(
                    tran_flux,
                    mp.DiffractedPlanewave(
                        [0, nm, 0],
                        mp.Vector3(1, 0, 0),
                        1 if S_pol else 0,
                        0 if S_pol else 1,
                    ),
                )
                t_coeffs = res.alpha
                Tmode = abs(t_coeffs[0, 0, 0]) ** 2 / input_flux[0]
                print(
                    "tran-order:, {}, {}, {}".format("s" if S_pol else "p", nm, Tmode)
                )
                Tsum += Tmode if nm == 0 else 2 * Tmode

        r_flux = mp.get_fluxes(refl_flux)
        t_flux = mp.get_fluxes(tran_flux)
        Rflux = -r_flux[0] / input_flux[0]
        Tflux = t_flux[0] / input_flux[0]

        print(f"refl:, {Rsum}, {Rflux}")
        print(f"tran:, {Tsum}, {Tflux}")
        print(f"sum:,  {Rsum + Tsum}, {Rflux + Tflux}")

        self.assertAlmostEqual(Rsum, Rflux, places=2)
        self.assertAlmostEqual(Tsum, Tflux, places=2)
        self.assertAlmostEqual(Rsum + Tsum, 1.00, places=2)


if __name__ == "__main__":
    unittest.main()
