from __future__ import division

import unittest
import meep as mp
from meep import mpb


class TestMPBWrappers(unittest.TestCase):

    def setUp(self):
        self.num_bands = 8
        self.k_points = [mp.Vector3(),
                         mp.Vector3(0.5),
                         mp.Vector3(0.5, 0.5),
                         mp.Vector3()]

        self.k_points = mp.interpolate(4, self.k_points)
        self.geometry = [mp.Cylinder(0.2, material=mp.Medium(epsilon=12))]
        self.geometry_lattice = mpb.Lattice(size=mp.Vector3(1, 1))
        self.resolution = 32

    # def get_mode_solver(self):
    #     return mpb.ModeSolver(
    #         num_bands=self.num_bands,
    #         k_points=self.k_points,
    #         geometry=self.geometry,
    #         geometry_lattice=self.geometry_lattice,
    #         resolution=self.resolution
    #     )

    def test_mode_solver_constructor(self):
        mpb.mode_solver(self.num_bands, 0, self.resolution, self.geometry_lattice,
                        1.0e-7, mp.Medium(), self.geometry)

    def test_mode_solver_init(self):
        ms = mpb.mode_solver(self.num_bands, 0, self.resolution, self.geometry_lattice,
                             1.0e-7, mp.Medium(), self.geometry)
        ms.init(0, True)


class TestModeSolver(unittest.TestCase):

    def test_list_split(self):
        k_points = [
            mp.Vector3(),
            mp.Vector3(0.5),
            mp.Vector3(0.5, 0.5),
            mp.Vector3()
        ]

        k_points = mp.interpolate(4, k_points)

        k_split = mp.list_split(k_points, 1, 0)

        expected = [
            (0, [mp.Vector3(),
                 mp.Vector3(0.10000000000000003),
                 mp.Vector3(0.20000000000000004),
                 mp.Vector3(0.30000000000000004),
                 mp.Vector3(0.4),
                 mp.Vector3(0.5),
                 mp.Vector3(0.5, 0.10000000000000003),
                 mp.Vector3(0.5, 0.20000000000000004),
                 mp.Vector3(0.5, 0.30000000000000004),
                 mp.Vector3(0.5, 0.4),
                 mp.Vector3(0.5, 0.5),
                 mp.Vector3(0.4, 0.4),
                 mp.Vector3(0.30000000000000004, 0.30000000000000004),
                 mp.Vector3(0.2, 0.2),
                 mp.Vector3(0.1, 0.1),
                 mp.Vector3(0.0, 0.0)]),
        ]

        indx = k_split[0][0]
        split_list = k_split[0][1]
        self.assertEqual(indx, 0)
        for res, exp in zip(split_list, expected[0][1]):
            self.assertEqual(res, exp)


if __name__ == '__main__':
    unittest.main()
