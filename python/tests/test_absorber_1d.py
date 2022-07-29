import unittest

from meep.materials import Al

import meep as mp


class TestAbsorber(unittest.TestCase):
    def setUp(self):

        resolution = 40
        cell_size = mp.Vector3(z=10)

        absorber_layers = [mp.Absorber(1, direction=mp.Z)]

        sources = [
            mp.Source(
                src=mp.GaussianSource(1 / 0.803, fwidth=0.1),
                center=mp.Vector3(),
                component=mp.Ex,
            )
        ]

        self.sim = mp.Simulation(
            cell_size=cell_size,
            resolution=resolution,
            dimensions=1,
            default_material=Al,
            boundary_layers=absorber_layers,
            sources=sources,
        )

    def test_absorber(self):

        self.sim.run(
            until_after_sources=mp.stop_when_fields_decayed(
                50, mp.Ex, mp.Vector3(), 1e-6
            )
        )

        f = self.sim.get_field_point(mp.Ex, mp.Vector3())
        self.assertAlmostEqual(f.real, 3.218846961494622e-13, places=6)

    def test_absorber_2d(self):
        source = mp.Source(
            src=mp.GaussianSource(frequency=0.1, fwidth=0.1),
            component=mp.Hz,
            center=mp.Vector3(),
        )

        sim = mp.Simulation(
            cell_size=mp.Vector3(20, 20, 0),
            resolution=10,
            sources=[source],
            boundary_layers=[mp.Absorber(5)],
        )

        sim.run(until_after_sources=1000)
        v = mp.Vector3(4.13, 3.75, 0)
        p = sim.get_field_point(mp.Hz, v)

        self.assertAlmostEqual(-4.058476603571745e-11, p.real)


if __name__ == "__main__":
    unittest.main()
