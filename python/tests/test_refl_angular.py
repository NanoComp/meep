import math
import unittest

import numpy as np
import parameterized
from utils import ApproxComparisonTestCase

import meep as mp


class TestReflAngular(ApproxComparisonTestCase):
    @classmethod
    def setUpClass(cls):
        cls.resolution = 400  # pixels/μm

        cls.n1 = 1.4  # refractive index of medium 1
        cls.n2 = 3.5  # refractive index of medium 2

        cls.dpml = 1.0
        cls.dz = 7.0
        cls.sz = cls.dz + 2 * cls.dpml

        cls.wvl_min = 0.4
        cls.wvl_max = 0.8
        cls.fmin = 1 / cls.wvl_max
        cls.fmax = 1 / cls.wvl_min
        cls.fcen = 0.5 * (cls.fmin + cls.fmax)
        cls.df = cls.fmax - cls.fmin
        cls.nfreq = 11

    def refl_angular(self, theta):
        theta_r = math.radians(theta)

        # wavevector (in source medium); plane of incidence is XZ
        k = (
            mp.Vector3(0, 0, 1)
            .rotate(mp.Vector3(0, 1, 0), theta_r)
            .scale(self.n1 * self.fmin)
        )

        dimensions = 1 if theta == 0 else 3
        cell_size = mp.Vector3(z=self.sz)
        pml_layers = [mp.PML(self.dpml)]

        sources = [
            mp.Source(
                mp.GaussianSource(self.fcen, fwidth=self.df),
                component=mp.Ex,  # P polarization
                center=mp.Vector3(z=-0.5 * self.sz + self.dpml),
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
        )

        mon_pt = -0.5 * self.sz + self.dpml + 0.25 * self.dz
        refl_fr = mp.FluxRegion(center=mp.Vector3(z=mon_pt))
        refl = sim.add_flux(self.fcen, self.df, self.nfreq, refl_fr)

        termination_cond = mp.stop_when_fields_decayed(
            50, mp.Ex, mp.Vector3(z=mon_pt), 1e-9
        )
        sim.run(until_after_sources=termination_cond)

        empty_data = sim.get_flux_data(refl)
        empty_flux = mp.get_fluxes(refl)
        sim.reset_meep()

        geometry = [
            mp.Block(
                size=mp.Vector3(mp.inf, mp.inf, 0.5 * self.sz),
                center=mp.Vector3(z=0.25 * self.sz),
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
            geometry=geometry,
        )

        refl = sim.add_flux(self.fcen, self.df, self.nfreq, refl_fr)
        sim.load_minus_flux_data(refl, empty_data)

        sim.run(until_after_sources=termination_cond)

        refl_flux = mp.get_fluxes(refl)
        freqs = mp.get_flux_freqs(refl)

        Rs = -np.array(refl_flux) / np.array(empty_flux)

        thetas = [math.asin(k.x / (self.n1 * freqs[i])) for i in range(self.nfreq)]
        return freqs, thetas, Rs

    @parameterized.parameterized.expand([(0,), (20.6,)])
    def test_refl_angular(self, theta):
        fmeep, tmeep, Rmeep = self.refl_angular(theta)

        # angle of refracted planewave in medium n2 for an
        # incident planewave in medium n1 at angle theta_in
        theta_out = lambda theta_in: math.asin(self.n1 * math.sin(theta_in) / self.n2)

        # Fresnel reflectance for P polarization in medium n2 for
        # an incident planewave in medium n1 at angle theta_in
        Rfresnel = lambda theta_in: (
            math.fabs(
                (self.n1 * math.cos(theta_out(theta_in)) - self.n2 * math.cos(theta_in))
                / (
                    self.n1 * math.cos(theta_out(theta_in))
                    + self.n2 * math.cos(theta_in)
                )
            )
            ** 2
        )

        Ranalytic = np.empty((self.nfreq,))
        print(
            "refl:, wavelength (μm), incident angle (°), reflectance (Meep), reflectance (analytic), error"
        )
        for i in range(self.nfreq):
            Ranalytic[i] = Rfresnel(tmeep[i])
            err = abs(Rmeep[i] - Ranalytic[i]) / Ranalytic[i]
            print(
                "refl:, {:4.2f}, {:4.2f}, {:8.6f}, {:8.6f}, {:6.4f}".format(
                    1 / fmeep[i], math.degrees(tmeep[i]), Rmeep[i], Ranalytic[i], err
                )
            )

        tol = 0.005 if mp.is_single_precision() else 0.004
        self.assertClose(Rmeep, Ranalytic, epsilon=tol)


if __name__ == "__main__":
    unittest.main()
