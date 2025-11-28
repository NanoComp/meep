"""Verfies boundary conditions for a 1D cell using a 3D simulation."""


import math
from typing import Tuple
import unittest

import meep as mp
import numpy as np


class TestBoundaries1D(unittest.TestCase):
    @classmethod
    def setUp(cls):
        cls.resolution_um = 100
        cls.pml_um = 1.0
        cls.bulk_um = 10.0
        cls.n_background = 1.2
        cls.antenna_height_um = 1.25
        cls.wavelength_um = 0.65
        cls.field_decay_threshold = 1e-6
        cls.field_decay_period = 25
        cls.frequency = 1 / cls.wavelength_um

    def planewave_above_ground_plane(self, kx: float, ky: float, kz: float) -> float:
        """
        Returns the radiated flux from a planewave source above a ground plane.

        Args:
            kx, ky, kz: the wavevector components of the planewave.

        Returns:
            The radiated flux from a linearly polarized planewave.
        """
        size_z_um = self.bulk_um + self.pml_um
        pml_layers = [mp.PML(thickness=self.pml_um, direction=mp.Z, side=mp.High)]
        src_pt = mp.Vector3(0, 0, -0.5 * size_z_um + self.antenna_height_um)

        cell_size = mp.Vector3(0, 0, size_z_um)

        src_cmpt = mp.Ey
        sources = [
            mp.Source(
                src=mp.GaussianSource(self.frequency, fwidth=0.1 * self.frequency),
                component=src_cmpt,
                center=src_pt,
                size=mp.Vector3(),
            )
        ]

        sim = mp.Simulation(
            resolution=self.resolution_um,
            default_material=mp.Medium(index=self.n_background),
            cell_size=cell_size,
            sources=sources,
            boundary_layers=pml_layers,
            k_point=mp.Vector3(kx, ky, kz),
        )

        sim.set_boundary(mp.Low, mp.Z, mp.Metallic)
        sim.set_boundary(mp.High, mp.Z, mp.Metallic)

        mon_pt = mp.Vector3(0, 0, 0.5 * size_z_um - self.pml_um)
        dft_flux_z = sim.add_flux(
            self.frequency,
            0,
            1,
            mp.FluxRegion(center=mon_pt, size=mp.Vector3(), direction=mp.Z),
        )

        sim.run(
            until_after_sources=mp.stop_when_fields_decayed(
                self.field_decay_period, src_cmpt, mon_pt, self.field_decay_threshold
            )
        )

        flux_z = mp.get_fluxes(dft_flux_z)[0]

        return flux_z

    def planewave_wavevector(self, polar_rad: float) -> Tuple[float, float, float]:
        """Returns the wavevector of the outgoing planewave in free space.

        Args:
            polar_rad: polar angle of the wavevector for φ = 0 (xz plane). 0 is
                the +z direction.

        Returns:
            The wavevector of the outgoing planewave.
        """
        kx = self.n_background * self.frequency * np.sin(polar_rad)
        ky = 0
        kz = self.n_background * self.frequency * np.cos(polar_rad)

        return kx, ky, kz

    def test_boundaries_1D(self):
        # Find the smallest angle at which the radiation pattern is a maximum using
        # the formula for the radiation pattern. This angle is used to compute the
        # radial flux which is then used for normalization of the radiation pattern.
        k_free_space = 2 * math.pi / (self.wavelength_um / self.n_background)
        polar_rad_max = math.acos(
            9 * math.pi / (2 * k_free_space * self.antenna_height_um)
        )
        kx, ky, kz = self.planewave_wavevector(polar_rad_max)
        flux_z = self.planewave_above_ground_plane(kx, ky, kz)
        radial_flux_meep_max = math.cos(polar_rad_max) * flux_z

        for polar_rad in np.deg2rad([0, 10.62, 26.7, 66.2]):
            kx, ky, kz = self.planewave_wavevector(polar_rad)

            # Skip wavevectors which are close to the light cone
            # due to poor absorption by PML (i.e. glancing-angle waves).
            if kx > (0.95 * self.n_background * self.frequency):
                continue

            flux_z = self.planewave_above_ground_plane(kx, ky, kz)
            radial_flux_meep = (math.cos(polar_rad) * flux_z) / radial_flux_meep_max

            # The radiation pattern of a two-element antenna array is equivalent
            # to the radiation pattern of a single antenna multiplied by its
            # array factor as derived in Section 6.2 "Two-Element Array" of
            # Antenna Theory: Analysis and Design, Fourth Edition (2016) by C.A.
            # Balanis.
            radial_flux_analytic = (
                math.sin(k_free_space * self.antenna_height_um * math.cos(polar_rad))
                ** 2
            )

            delta = abs(radial_flux_meep - radial_flux_analytic)
            print(
                f"radial_flux:, {np.rad2deg(polar_rad):.2f}°, "
                f"{radial_flux_meep:.6f} (meep), {radial_flux_analytic:.6f} "
                f"(theory), {delta:.6f}"
            )

            self.assertAlmostEqual(radial_flux_meep, radial_flux_analytic, delta=0.01)


if __name__ == "__main__":
    unittest.main()
