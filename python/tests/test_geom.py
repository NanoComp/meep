import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import meep as mp
import geom as gm


class TestSphere(unittest.TestCase):

    def test_kwargs_passed_to_parent(self):
        s = gm.Sphere()
        self.assertEqual(s.material.data, None)
        self.assertEqual(s.center, None)
        self.assertEqual(s.radius, None)

        s = gm.Sphere(radius=1.0)
        self.assertEqual(s.material.data, None)
        self.assertEqual(s.center, None)
        self.assertEqual(s.radius, 1.0)

        s = gm.Sphere(center=(1, 1, 1))
        self.assertEqual(s.material.data, None)
        self.assertEqual(s.center, (1, 1, 1))
        self.assertEqual(s.radius, None)

    def test_invalid_kwarg_raises_exception(self):
        with self.assertRaises(TypeError):
            gm.Sphere(invalid='This is not allowed')
        with self.assertRaises(TypeError):
            gm.Sphere(radius=1.0, oops='Nope')

    def test_non_neg_radius_constructor(self):
        gm.Sphere(radius=0.0)
        gm.Sphere(radius=1.0)

        with self.assertRaises(ValueError) as ctx:
            gm.Sphere(radius=-1)
            self.assertIn("Got -1", ctx.exception)

    def test_non_neg_radius_setter(self):
        s = gm.Sphere(radius=0.0)
        s.radius = 1.0

        with self.assertRaises(ValueError) as ctx:
            s.radius = -1.0
            self.assertIn("Got -1.0", ctx.exception)

    def test_contains_point(self):
        s = gm.Sphere(center=gm.Vector3(0, 0, 0), radius=2.0)
        point = gm.Vector3(1, 1, 1)
        self.assertTrue(point in s)
        self.assertIn(point, s)
        self.assertFalse(gm.Vector3(10, 10, 10) in s)


class TestCylinder(unittest.TestCase):

    def test_non_neg_height_constructor(self):
        gm.Cylinder(height=0.0)
        gm.Cylinder(height=1.0)

        with self.assertRaises(ValueError) as ctx:
            gm.Cylinder(height=-1)
            self.assertIn("Got -1", ctx.exception)

    def test_non_neg_height_setter(self):
        s = gm.Cylinder(height=0.0)
        s.height = 1.0

        with self.assertRaises(ValueError) as ctx:
            s.height = -1.0
            self.assertIn("Got -1.0", ctx.exception)

    def test_contains_point(self):
        c = gm.Cylinder(center=gm.Vector3(0, 0, 0), radius=2.0, height=4.0, axis=gm.Vector3(0, 0, 1))

        self.assertIn(gm.Vector3(0, 0, 0), c)
        self.assertIn(gm.Vector3(2, 0, 0), c)
        self.assertIn(gm.Vector3(2, 0, 2), c)

        self.assertNotIn(gm.Vector3(2.0001, 0, 0), c)
        self.assertNotIn(gm.Vector3(10, 10, 10), c)


class TestWedge(unittest.TestCase):

    def test_default_properties(self):
        import math
        w = gm.Wedge(center=gm.Vector3(0, 0, 0), radius=2.0, height=4.0, axis=gm.Vector3(0, 0, 1))
        self.assertEqual(w.wedge_angle, 8 * math.atan(1))

    def test_contains_point(self):
        w = gm.Wedge(center=gm.Vector3(0, 0, 0), radius=2.0, height=4.0, axis=gm.Vector3(0, 0, 1))
        self.assertIn(gm.Vector3(2.0, 0, 0), w)


class TestCone(unittest.TestCase):

    def test_contains_point(self):
        c = gm.Cone(center=gm.Vector3(0, 0, 0), radius=2.0, height=3.0, axis=gm.Vector3(0, 0, 1))
        self.assertIn(gm.Vector3(0, 0, 1), c)


class TestBlock(unittest.TestCase):

    def test_contains_point(self):
        b = gm.Block(size=gm.Vector3(1, 1, 1), center=gm.Vector3(0, 0, 0))
        self.assertIn(gm.Vector3(0, 0, 0), b)


class TestEllipsoid(unittest.TestCase):

    def test_contains_point(self):
        e = gm.Ellipsoid(size=gm.Vector3(1, 1, 1), center=gm.Vector3(0, 0, 0))
        self.assertIn(gm.Vector3(0, 0, 0), e)

    # TODO(chogan): Allow python to read this member after it's computed in C
    def test_inverse_semi_axes(self):
        pass


class TestCompoundGeometricObject(unittest.TestCase):

    def test_convert_pylist_to_cgo_list(self):
        s = gm.Sphere(center=gm.Vector3(0, 0, 0), radius=2.0)
        c = gm.Cylinder(center=gm.Vector3(0, 0, 0), radius=2.0, height=4.0, axis=gm.Vector3(0, 0, 1))
        b = gm.Block(size=gm.Vector3(1, 1, 1), center=gm.Vector3(0, 0, 0))

        cgo = gm.CompoundGeometricObject(center=gm.Vector3(0, 0, 0), component_objects=[s, c, b])

        sxy = 2.0 * (1 + 1 + 4 + 2)
        resolution = 10.0
        gv = mp.voltwo(sxy, sxy, resolution)

        def dummy_eps(v):
            return 1.0

        dpml = 2.0
        sym = mp.mirror(mp.Y, gv)

        the_structure = mp.structure(gv, dummy_eps, mp.pml(dpml), sym)

        mp.set_materials_from_geometry(the_structure, cgo)


if __name__ == '__main__':
    unittest.main()
