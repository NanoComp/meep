import os
import sys
import unittest

# import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import meep as mp
import geom as gm


class TestTypemaps(unittest.TestCase):

    def test_point_in_objectp(self):
        s = gm.Sphere(center=gm.Vector3(0, 0, 0), radius=2.0)

        # Accepts Vector3
        self.assertTrue(mp.point_in_objectp(gm.Vector3(0, 0, 0), s))
        self.assertFalse(mp.point_in_objectp(gm.Vector3(10, 10, 10), s))

        # # Accepts tuple
        # self.assertTrue(mp.point_in_objectp((0, 0, 0), s))
        # self.assertFalse(mp.point_in_objectp((10, 10, 10), s))

        # # Accepts list
        # self.assertTrue(mp.point_in_objectp([0, 0, 0], s))
        # self.assertFalse(mp.point_in_objectp([10, 10, 10], s))

        # # Accepts numpy array
        # self.assertTrue(mp.point_in_objectp(np.array([0, 0, 0]), s))
        # self.assertFalse(mp.point_in_objectp(np.array([10, 10, 10]), s))


if __name__ == '__main__':
    unittest.main()
