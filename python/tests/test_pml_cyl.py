import unittest
import parameterized
import numpy as np
import meep as mp


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

    @parameterized.parameterized.expand(
        [
            (-1.0, 0.0, False),
            (2.0, 0.14, False),
            (3.0, 0.17, True),
        ]
    )
    def test_pml_cyl(
        self, m: float, rpos: float, accurate_fields_near_cylorigin: bool = False
    ):
        """Verifies that the z-PML in cylindrical coordinates properly
        attenuates fields at r=0.

        Args:
           m: exp(imϕ) angular dependence of the fields.
           rpos: position of source along R direction.
           accurate_fields_near_cylorigin: whether to compute more accurate
              fields near the origin r=0.
        """

        pml_layers = [
            mp.PML(self.dpml_r, direction=mp.R),
            mp.PML(self.dpml_z, direction=mp.Z),
        ]

        sources = [
            mp.Source(
                src=mp.GaussianSource(self.fcen, fwidth=0.1 * self.fcen),
                center=mp.Vector3(rpos, 0, 0),
                component=mp.Er,
            ),
        ]

        sim = mp.Simulation(
            resolution=self.resolution,
            cell_size=self.cell_size,
            dimensions=mp.CYLINDRICAL,
            m=m,
            sources=sources,
            boundary_layers=pml_layers,
            accurate_fields_near_cylorigin=accurate_fields_near_cylorigin,
        )

        if accurate_fields_near_cylorigin and abs(m) > 1:
            sim.Courant = 1 / (abs(m) + 0.6)

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

        prev_flux = [
            mp.get_fluxes(flux_plus_z)[0],
            mp.get_fluxes(flux_plus_r)[0],
            mp.get_fluxes(flux_minus_z)[0],
        ]

        for t in [142.15, 214.64, 365.32]:
            sim.run(until_after_sources=t)

            cur_flux = [
                mp.get_fluxes(flux_plus_z)[0],
                mp.get_fluxes(flux_plus_r)[0],
                mp.get_fluxes(flux_minus_z)[0],
            ]
            cur_flux_str = ", ".join(f"{c:.6f}" for c in cur_flux)
            flux_tot = sum(cur_flux)

            print(f"flux:, {sim.meep_time()}, {cur_flux_str}, {flux_tot:.6f}")

            places = 6 if mp.is_single_precision() else 8
            for i in range(len(cur_flux)):
                self.assertAlmostEqual(
                    prev_flux[i],
                    cur_flux[i],
                    places=places,
                )
                prev_flux[i] = cur_flux[i]


if __name__ == "__main__":
    unittest.main()
