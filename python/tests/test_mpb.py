from __future__ import division

import unittest

import meep as mp
from mpb import MaterialGrid

unit_vector = mp.Vector3(1, 1, 1)


class TestMaterialGrid(unittest.TestCase):

    def test_material_grid_kind_constraints(self):
        with self.assertRaises(ValueError) as ctx:
            MaterialGrid(0, 1, unit_vector, material_grid_kind=-1)
            self.assertIn("Got -1", ctx.exception)

        with self.assertRaises(ValueError) as ctx:
            MaterialGrid(0, 1, unit_vector, material_grid_kind=10)
            self.assertIn("Got 10", ctx.exception)

        mg = MaterialGrid(0, 1, unit_vector)

        with self.assertRaises(ValueError) as ctx:
            mg.material_grid_kind = -1
            self.assertIn("Got -1", ctx.exception)

        with self.assertRaises(ValueError) as ctx:
            mg.material_grid_kind = 10
            self.assertIn("Got 10", ctx.exception)

        mg.material_grid_kind = 1
        self.assertEqual(mg.material_grid_kind, 1)

    def test_size_constraints(self):
        with self.assertRaises(ValueError) as ctx:
            MaterialGrid(0, 1, mp.Vector3(-1, 1, 1))
            self.assertIn("Got -1", ctx.exception)

        with self.assertRaises(ValueError) as ctx:
            MaterialGrid(0, 1, mp.Vector3(1, 1, 1.5))
            self.assertIn("Got 1.5", ctx.exception)

        MaterialGrid(0, 1, mp.Vector3(1, 1, 1))

    def test_matgrid_init_constraints(self):

        def wrong(x, y, z, oops):
            return 1.0

        with self.assertRaises(ValueError) as ctx:
            MaterialGrid(0, 1, unit_vector, matgrid_init=wrong)
            self.assertIn("Got 4", ctx.exception)

if __name__ == '__main__':
    unittest.main()
