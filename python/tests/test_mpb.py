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


# class TestLattice(unittest.TestCase):

#     def test_lattice(self):
#         pass


# class TestMatrix(unittest.TestCase):

#     def test_matrix(self):
#         pass


if __name__ == '__main__':
    unittest.main()
