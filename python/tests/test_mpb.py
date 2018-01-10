from __future__ import division

import unittest
import meep as mp
from meep import mpb


class TestMPB(unittest.TestCase):

    def test_initialization(self):
        num_bands = 8
        k_points = [mp.Vector3(),
                    mp.Vector3(0.5),
                    mp.Vector3(0.5, 0.5),
                    mp.Vector3()]

        k_points = mp.interpolate(4, k_points)
        geometry = [mp.Cylinder(0.2, material=mp.Medium(epsilon=12))]
        geometry_lattice = mpb.Lattice(size=mp.Vector3(1, 1))
        resolution = 32

        ms = mpb.ModeSolver(num_bands=num_bands,
                            k_points=k_points,
                            geometry=geometry,
                            geometry_lattice=geometry_lattice,
                            resolution=resolution)

        ms.mode_solver.init(0, True)

if __name__ == '__main__':
    unittest.main()
