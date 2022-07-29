import cmath
import math
import unittest
from time import time

import parameterized

import meep as mp


class TestSpecialKz(unittest.TestCase):
    def refl_planar(self, theta, kz_2d):
        resolution = 100  # pixels/um

        dpml = 1.0
        sx = 3.0 + 2 * dpml
        sy = 1 / resolution
        cell_size = mp.Vector3(sx, sy)
        pml_layers = [mp.PML(dpml, direction=mp.X)]

        fcen = 1.0

        # plane of incidence is XZ
        k_point = mp.Vector3(1, 0, 0).rotate(mp.Vector3(0, 1, 0), theta).scale(fcen)

        sources = [
            mp.Source(
                mp.GaussianSource(fcen, fwidth=0.2 * fcen),
                component=mp.Ez,  # P-polarization
                center=mp.Vector3(-0.5 * sx + dpml),
                size=mp.Vector3(y=sy),
            )
        ]

        sim = mp.Simulation(
            cell_size=cell_size,
            boundary_layers=pml_layers,
            sources=sources,
            k_point=k_point,
            kz_2d=kz_2d,
            resolution=resolution,
        )

        refl_fr = mp.FluxRegion(center=mp.Vector3(-0.25 * sx), size=mp.Vector3(y=sy))
        refl = sim.add_flux(fcen, 0, 1, refl_fr)

        sim.run(
            until_after_sources=mp.stop_when_fields_decayed(
                50, mp.Ez, mp.Vector3(), 1e-9
            )
        )

        empty_flux = mp.get_fluxes(refl)
        empty_data = sim.get_flux_data(refl)
        sim.reset_meep()

        geometry = [
            mp.Block(
                material=mp.Medium(index=3.5),
                size=mp.Vector3(0.5 * sx, mp.inf, mp.inf),
                center=mp.Vector3(0.25 * sx),
            )
        ]

        sim = mp.Simulation(
            cell_size=cell_size,
            boundary_layers=pml_layers,
            geometry=geometry,
            sources=sources,
            k_point=k_point,
            kz_2d=kz_2d,
            resolution=resolution,
        )

        refl = sim.add_flux(fcen, 0, 1, refl_fr)
        sim.load_minus_flux_data(refl, empty_data)

        sim.run(
            until_after_sources=mp.stop_when_fields_decayed(
                50, mp.Ez, mp.Vector3(), 1e-9
            )
        )

        refl_flux = mp.get_fluxes(refl)

        return -refl_flux[0] / empty_flux[0]

    def test_special_kz(self):
        n1 = 1
        n2 = 3.5

        # compute angle of refracted planewave in medium n2
        # for incident planewave in medium n1 at angle theta_in
        theta_out = lambda theta_in: math.asin(n1 * math.sin(theta_in) / n2)

        # compute Fresnel reflectance for P-polarization in medium n2
        # for incident planewave in medium n1 at angle theta_in
        Rfresnel = (
            lambda theta_in: math.fabs(
                (n1 * math.cos(theta_out(theta_in)) - n2 * math.cos(theta_in))
                / (n1 * math.cos(theta_out(theta_in)) + n2 * math.cos(theta_in))
            )
            ** 2
        )

        # angle of incident planewave; clockwise (CW) about Y axis, 0 degrees along +X axis
        theta = math.radians(23)

        start = time()
        Rmeep_complex = self.refl_planar(theta, "complex")
        t_complex = time() - start

        start = time()
        Rmeep_real_imag = self.refl_planar(theta, "real/imag")
        t_real_imag = time() - start

        Rfres = Rfresnel(theta)

        self.assertAlmostEqual(Rmeep_complex, Rfres, places=2)
        self.assertAlmostEqual(Rmeep_real_imag, Rfres, places=2)

        # the real/imag algorithm should be faster, but on CI machines performance is too variable
        # for this to reliably hold
        # self.assertLess(t_real_imag,t_complex)

    @parameterized.parameterized.expand([("complex",), ("real/imag",)])
    def test_eigsrc_kz(self, kz_2d):
        resolution = 30  # pixels/um

        cell_size = mp.Vector3(14, 14)

        pml_layers = [mp.PML(thickness=2)]

        geometry = [
            mp.Block(
                center=mp.Vector3(),
                size=mp.Vector3(mp.inf, 1, mp.inf),
                material=mp.Medium(epsilon=12),
            )
        ]

        fsrc = 0.3  # frequency of eigenmode or constant-amplitude source
        bnum = 1  # band number of eigenmode
        kz = 0.2  # fixed out-of-plane wavevector component

        sources = [
            mp.EigenModeSource(
                src=mp.GaussianSource(fsrc, fwidth=0.2 * fsrc),
                center=mp.Vector3(),
                size=mp.Vector3(y=14),
                eig_band=bnum,
                eig_parity=mp.EVEN_Y,
                eig_match_freq=True,
            )
        ]

        sim = mp.Simulation(
            cell_size=cell_size,
            resolution=resolution,
            boundary_layers=pml_layers,
            sources=sources,
            geometry=geometry,
            symmetries=[mp.Mirror(mp.Y)],
            k_point=mp.Vector3(z=kz),
            kz_2d=kz_2d,
        )

        tran = sim.add_flux(
            fsrc, 0, 1, mp.FluxRegion(center=mp.Vector3(x=5), size=mp.Vector3(y=14))
        )

        sim.run(until_after_sources=50)

        res = sim.get_eigenmode_coefficients(tran, [1, 2], eig_parity=mp.EVEN_Y)

        total_flux = mp.get_fluxes(tran)[0]
        mode1_flux = abs(res.alpha[0, 0, 0]) ** 2
        mode2_flux = abs(res.alpha[1, 0, 0]) ** 2

        mode1_frac = 0.99
        self.assertGreater(mode1_flux / total_flux, mode1_frac)
        self.assertLess(mode2_flux / total_flux, 1 - mode1_frac)

        d = 3.5
        ez1 = sim.get_field_point(mp.Ez, mp.Vector3(2.3, -5.7, 4.8))
        ez2 = sim.get_field_point(mp.Ez, mp.Vector3(2.3, -5.7, 4.8 + d))
        ratio_ez = ez2 / ez1
        phase_diff = cmath.exp(1j * 2 * cmath.pi * kz * d)
        self.assertAlmostEqual(ratio_ez.real, phase_diff.real, places=10)
        self.assertAlmostEqual(ratio_ez.imag, phase_diff.imag, places=10)


if __name__ == "__main__":
    unittest.main()
