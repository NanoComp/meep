import math
from typing import List, Tuple
import unittest

import numpy as np
import parameterized
from utils import ApproxComparisonTestCase

import meep as mp


class TestReflectanceAngular(ApproxComparisonTestCase):
    @classmethod
    def setUpClass(cls):
        cls.resolution = 200  # pixels/μm

        cls.n1 = 1.4  # refractive index of medium 1
        cls.n2 = 3.5  # refractive index of medium 2

        cls.t_pml = 1.0
        cls.length_z = 7.0
        cls.size_z = cls.length_z + 2 * cls.t_pml

        cls.wavelength_min = 0.4
        cls.wavelength_max = 0.8
        cls.frequency_min = 1 / cls.wavelength_max
        cls.frequency_max = 1 / cls.wavelength_min
        cls.frequency_center = 0.5 * (cls.frequency_min + cls.frequency_max)
        cls.frequency_width = cls.frequency_max - cls.frequency_min
        cls.num_freq = 11

    def reflectance_angular(
        self, theta_deg: float, use_bfast: bool
    ) -> Tuple[List, List, np.ndarray]:
        """Computes properties of the incident and reflected planewave.

        Args:
          theta_deg: angle of incident planewave.
          use_bfast: whether to use the same angle for the incident planewave
            for all frequencies. If False, the incident angle is frequency
            dependent.

        Returns:
          A 3-tuple comprising the frequencies of the incident planewave,
          angles of the incident planewave, and the reflectance.
        """
        theta_rad = math.radians(theta_deg)

        if use_bfast:
            bfast_scaled_k = (self.n1 * np.sin(theta_rad), 0, 0)

            Courant = (1 - bfast_scaled_k[0]) / 3**0.5

            k = mp.Vector3()
        else:
            bfast_scaled_k = (0, 0, 0)

            Courant = 0.5

            # Wavevector in source medium with refractive index n1.
            # Plane of incidence is XZ. Rotation is counter clockwise about
            # Y axis. A rotation angle of zero is the +Z axis.
            k = (
                mp.Vector3(0, 0, 1)
                .rotate(mp.Vector3(0, 1, 0), theta_rad)
                .scale(self.n1 * self.frequency_min)
            )

        dimensions = 1 if theta_deg == 0 else 3
        cell_size = mp.Vector3(z=self.size_z)
        pml_layers = [mp.PML(self.t_pml)]

        # P polarization.
        source_component = mp.Ex

        sources = [
            mp.Source(
                mp.GaussianSource(self.frequency_center, fwidth=self.frequency_width),
                component=source_component,
                center=mp.Vector3(z=-0.5 * self.size_z + self.t_pml),
            )
        ]

        sim = mp.Simulation(
            resolution=self.resolution,
            cell_size=cell_size,
            dimensions=dimensions,
            default_material=mp.Medium(index=self.n1),
            sources=sources,
            boundary_layers=pml_layers,
            k_point=k,
            bfast_scaled_k=bfast_scaled_k,
            Courant=Courant,
        )

        monitor_point = -0.5 * self.size_z + self.t_pml + 0.25 * self.length_z
        monitor_region = mp.FluxRegion(center=mp.Vector3(z=monitor_point))
        flux_monitor = sim.add_flux(
            self.frequency_center, self.frequency_width, self.num_freq, monitor_region
        )

        termination_criteria = mp.stop_when_fields_decayed(
            50, source_component, mp.Vector3(z=monitor_point), 1e-6
        )
        sim.run(until_after_sources=termination_criteria)

        empty_data = sim.get_flux_data(flux_monitor)
        empty_flux = mp.get_fluxes(flux_monitor)

        sim.reset_meep()

        geometry = [
            mp.Block(
                size=mp.Vector3(mp.inf, mp.inf, 0.5 * self.size_z),
                center=mp.Vector3(z=0.25 * self.size_z),
                material=mp.Medium(index=self.n2),
            )
        ]

        sim = mp.Simulation(
            resolution=self.resolution,
            cell_size=cell_size,
            dimensions=dimensions,
            default_material=mp.Medium(index=self.n1),
            sources=sources,
            boundary_layers=pml_layers,
            k_point=k,
            bfast_scaled_k=bfast_scaled_k,
            Courant=Courant,
            geometry=geometry,
        )

        flux_monitor = sim.add_flux(
            self.frequency_center, self.frequency_width, self.num_freq, monitor_region
        )
        sim.load_minus_flux_data(flux_monitor, empty_data)

        sim.run(until_after_sources=termination_criteria)

        flux_monitor_flux = mp.get_fluxes(flux_monitor)
        freqs = mp.get_flux_freqs(flux_monitor)

        reflectance = -np.array(flux_monitor_flux) / np.array(empty_flux)

        if use_bfast:
            theta_in_rad = [theta_rad] * self.num_freq
        else:
            # Returns the angle of the incident planewave in medium n1 based
            # on its frequency given a fixed wavevector component in X.
            theta_in_rad = [
                math.asin(k.x / (self.n1 * freqs[i])) for i in range(self.num_freq)
            ]

        return freqs, theta_in_rad, reflectance

    @parameterized.parameterized.expand([(0, False), (20.6, False), (35.7, True)])
    def test_reflectance_angular(self, theta_deg: float, use_bfast: bool):
        (
            frequency_meep,
            theta_in_rad_meep,
            reflectance_meep,
        ) = self.reflectance_angular(theta_deg, use_bfast)

        # Returns angle of refracted planewave in medium n2 given
        # an incident planewave in medium n1 at angle theta_in_rad.
        theta_out = lambda theta_in_rad: math.asin(
            self.n1 * math.sin(theta_in_rad) / self.n2
        )

        # Returns Fresnel reflectance for P polarization in medium n2
        # for an incident planewave in medium n1 at angle theta_in_rad.
        reflectance_fresnel = lambda theta_in_rad: (
            math.fabs(
                (
                    self.n1 * math.cos(theta_out(theta_in_rad))
                    - self.n2 * math.cos(theta_in_rad)
                )
                / (
                    self.n1 * math.cos(theta_out(theta_in_rad))
                    + self.n2 * math.cos(theta_in_rad)
                )
            )
            ** 2
        )

        reflectance_analytic = np.empty((self.num_freq,))
        print(
            "refl:, wavelength (μm), incident angle (°), reflectance (Meep), "
            "reflectance (analytic), error"
        )
        for i in range(self.num_freq):
            reflectance_analytic[i] = reflectance_fresnel(theta_in_rad_meep[i])
            err = (
                abs(reflectance_meep[i] - reflectance_analytic[i])
                / reflectance_analytic[i]
            )
            print(
                "refl:, {:4.2f}, {:4.2f}, {:8.6f}, {:8.6f}, {:6.4f}".format(
                    1 / frequency_meep[i],
                    math.degrees(theta_in_rad_meep[i]),
                    reflectance_meep[i],
                    reflectance_analytic[i],
                    err,
                )
            )

        tol = 0.03
        self.assertClose(reflectance_meep, reflectance_analytic, epsilon=tol)


if __name__ == "__main__":
    unittest.main()
