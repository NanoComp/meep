import sys
import unittest

import meep as mp


class TestWvgSrc(unittest.TestCase):
    def setUp(self):

        cell = mp.Vector3(16, 8)

        geometry = [
            mp.Block(
                center=mp.Vector3(),
                size=mp.Vector3(mp.inf, 1, mp.inf),
                material=mp.Medium(epsilon=12),
            ),
            mp.Block(
                center=mp.Vector3(y=0.3),
                size=mp.Vector3(mp.inf, 0.1, mp.inf),
                material=mp.Medium(),
            ),
        ]

        sources = [
            mp.EigenModeSource(
                src=mp.ContinuousSource(0.15),
                size=mp.Vector3(y=6),
                center=mp.Vector3(x=-5),
                component=mp.Dielectric,
                eig_parity=mp.ODD_Z,
            )
        ]

        pml_layers = [mp.PML(1.0)]

        self.sim = mp.Simulation(
            cell_size=cell,
            geometry=geometry,
            sources=sources,
            boundary_layers=pml_layers,
            force_complex_fields=True,
            resolution=10,
        )

    def test_wvg_src(self):
        self.sim.run(until=200)

        flux1 = self.sim.flux_in_box(
            mp.X, mp.Volume(center=mp.Vector3(-6.0), size=mp.Vector3(1.8, 6))
        )
        flux2 = self.sim.flux_in_box(
            mp.X, mp.Volume(center=mp.Vector3(6.0), size=mp.Vector3(1.8, 6))
        )

        self.assertAlmostEqual(flux1, -1.775216564842667e-03)
        places = 5 if mp.is_single_precision() else 7
        self.assertAlmostEqual(flux2, 7.215785537102116, places=places)


if __name__ == "__main__":
    unittest.main()
