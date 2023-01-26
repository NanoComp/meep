import unittest
import meep as mp
import numpy as np


class TestPMLCylindrical(unittest.TestCase):
    @classmethod
    def setUp(cls):
        cls.resolution = 25  # pixels/um
        cls.s = 5.0
        cls.dpml_r = 1.0
        cls.dpml_z = 1.0
        cls.cell_size = mp.Vector3(
            cls.s + cls.dpml_r,
            0,
            cls.s + 2 * cls.dpml_z,
        )
        cls.fcen = 1.0

    def test_pml_cyl_coords(self):
        """Verifies that the z-PML in cylindrical coordinates properly
        attenuates fields at r=0.
        """

        pml_layers = [
            mp.PML(self.dpml_r, direction=mp.R),
            mp.PML(self.dpml_z, direction=mp.Z),
        ]

        sources = [
            mp.Source(
                src=mp.GaussianSource(self.fcen, fwidth=0.1 * self.fcen),
                center=mp.Vector3(),
                component=mp.Er,
            ),
        ]

        sim = mp.Simulation(
            resolution=self.resolution,
            cell_size=self.cell_size,
            dimensions=mp.CYLINDRICAL,
            m=-1,
            sources=sources,
            boundary_layers=pml_layers,
        )

        flux_plus_z = sim.add_flux(
            self.fcen,
            0,
            1,
            mp.FluxRegion(
                center=mp.Vector3(
                    0.5 * self.s,
                    0,
                    0.5 * self.s,
                ),
                size=mp.Vector3(self.s, 0, 0),
            ),
        )

        flux_plus_r = sim.add_flux(
            self.fcen,
            0,
            1,
            mp.FluxRegion(
                center=mp.Vector3(self.s, 0, 0),
                size=mp.Vector3(0, 0, self.s),
            ),
        )

        flux_minus_z = sim.add_flux(
            self.fcen,
            0,
            1,
            mp.FluxRegion(
                center=mp.Vector3(
                    0.5 * self.s,
                    0,
                    -0.5 * self.s,
                ),
                size=mp.Vector3(self.s, 0, 0),
                weight=-1.0,
            ),
        )

        sim.run(until_after_sources=50.94)

        prev_flux_plusz = mp.get_fluxes(flux_plus_z)[0]
        prev_flux_plusr = mp.get_fluxes(flux_plus_r)[0]
        prev_flux_minusz = mp.get_fluxes(flux_minus_z)[0]

        for t in [142.15, 214.64, 365.32]:
            sim.run(until_after_sources=t)

            flux_plusz = mp.get_fluxes(flux_plus_z)[0]
            flux_plusr = mp.get_fluxes(flux_plus_r)[0]
            flux_minusz = mp.get_fluxes(flux_minus_z)[0]
            flux_tot = flux_plusz + flux_plusr + flux_minusz

            print(
                f"flux:, {sim.meep_time()}, {flux_plusz:.6f}, "
                f"{flux_plusr:.6f}, {flux_minusz:.6f}, {flux_tot:.6f}"
            )

            places = 6 if mp.is_single_precision() else 9
            self.assertAlmostEqual(flux_plusz, prev_flux_plusz, places=places)
            self.assertAlmostEqual(flux_plusr, prev_flux_plusr, places=places)
            self.assertAlmostEqual(flux_minusz, prev_flux_minusz, places=places)

            prev_flux_plusz = flux_plusz
            prev_flux_plusr = flux_plusr
            prev_flux_minusz = flux_minusz


if __name__ == "__main__":
    unittest.main()
