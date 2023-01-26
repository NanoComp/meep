import unittest
import meep as mp
import numpy as np


class TestPMLCylindrical(unittest.TestCase):
    def test_pml_cyl_coords(self):
        """Verifies that the z-PML in cylindrical coordinates properly
        attenuates fields at r=0.
        """

        resolution = 25  # pixels/um
        s = 5.0
        dpml_r = 1.0
        dpml_z = 1.0
        cell_size = mp.Vector3(s + dpml_r, 0, s + 2 * dpml_z)

        pml_layers = [
            mp.PML(dpml_r, direction=mp.R),
            mp.PML(dpml_z, direction=mp.Z),
        ]

        fcen = 1.0

        sources = [
            mp.Source(
                src=mp.GaussianSource(fcen, fwidth=0.1 * fcen),
                center=mp.Vector3(),
                component=mp.Er,
            ),
        ]

        sim = mp.Simulation(
            resolution=resolution,
            cell_size=cell_size,
            dimensions=mp.CYLINDRICAL,
            m=-1,
            accurate_fields_near_cylorigin=False,
            sources=sources,
            boundary_layers=pml_layers,
        )

        flux_plus_z = sim.add_flux(
            fcen,
            0,
            1,
            mp.FluxRegion(
                center=mp.Vector3(0.5 * s, 0, 0.5 * s), size=mp.Vector3(s, 0, 0)
            ),
        )

        flux_plus_r = sim.add_flux(
            fcen,
            0,
            1,
            mp.FluxRegion(center=mp.Vector3(s, 0, 0), size=mp.Vector3(0, 0, s)),
        )

        flux_minus_z = sim.add_flux(
            fcen,
            0,
            1,
            mp.FluxRegion(
                center=mp.Vector3(0.5 * s, 0, -0.5 * s),
                size=mp.Vector3(s, 0, 0),
                weight=-1.0,
            ),
        )

        sim.run(until_after_sources=50)

        prev_flux_plusz = mp.get_fluxes(flux_plus_z)[0]
        prev_flux_plusr = mp.get_fluxes(flux_plus_r)[0]
        prev_flux_minusz = mp.get_fluxes(flux_minus_z)[0]

        for t in [100, 200, 400]:
            sim.run(until_after_sources=t)

            flux_plusz = mp.get_fluxes(flux_plus_z)[0]
            flux_plusr = mp.get_fluxes(flux_plus_r)[0]
            flux_minusz = mp.get_fluxes(flux_minus_z)[0]
            flux_tot = flux_plusz + flux_plusr + flux_minusz

            print(
                f"flux:, {sim.meep_time()}, {flux_plusz:.6f}, "
                f"{flux_plusr:.6f}, {flux_minusz:.6f}, {flux_tot:.6f}"
            )

            self.assertAlmostEqual(flux_plusz, prev_flux_plusz, places=9)
            self.assertAlmostEqual(flux_plusr, prev_flux_plusr, places=9)
            self.assertAlmostEqual(flux_minusz, prev_flux_minusz, places=9)

            prev_flux_plusz = flux_plusz
            prev_flux_plusr = flux_plusr
            prev_flux_minusz = flux_minusz


if __name__ == "__main__":
    unittest.main()
