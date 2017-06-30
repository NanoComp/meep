import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import geom as gm


def zeros():
    return gm.Vector3(0, 0, 0)


def ones():
    return gm.Vector3(1, 1, 1)


class TestSphere(unittest.TestCase):

    def test_kwargs_passed_to_parent(self):
        s = gm.Sphere()
        self.assertEqual(s.material.epsilon_diag, ones())
        self.assertEqual(s.center, zeros())
        self.assertEqual(s.radius, 1)

        s = gm.Sphere(radius=2.0)
        self.assertEqual(s.material.epsilon_diag, ones())
        self.assertEqual(s.center, zeros())
        self.assertEqual(s.radius, 2.0)

        s = gm.Sphere(center=ones())
        self.assertEqual(s.material.epsilon_diag, ones())
        self.assertEqual(s.center, ones())
        self.assertEqual(s.radius, 1)

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
        s = gm.Sphere(center=zeros(), radius=2.0)
        point = ones()
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
        c = gm.Cylinder(center=zeros(), radius=2.0, height=4.0)

        self.assertIn(zeros(), c)
        self.assertIn(gm.Vector3(2, 0, 0), c)
        self.assertIn(gm.Vector3(2, 0, 2), c)

        self.assertNotIn(gm.Vector3(2.0001, 0, 0), c)
        self.assertNotIn(gm.Vector3(10, 10, 10), c)

    def test_missing_required_arg_throws(self):
        c = gm.Cylinder(radius=2.0, height=4.0, center=None)

        with self.assertRaises(ValueError) as ctx:
            self.assertIn(zeros(), c)
            self.assertIn("Vector3 is not initialized", ctx.exception)


class TestWedge(unittest.TestCase):

    def test_default_properties(self):
        import math
        w = gm.Wedge(center=zeros(), radius=2.0, height=4.0, axis=gm.Vector3(0, 0, 1))
        self.assertEqual(w.wedge_angle, 8 * math.atan(1))

    def test_contains_point(self):
        w = gm.Wedge(center=zeros(), radius=2.0, height=4.0, axis=gm.Vector3(0, 0, 1))
        self.assertIn(gm.Vector3(2.0, 0, 0), w)


class TestCone(unittest.TestCase):

    def test_contains_point(self):
        c = gm.Cone(center=zeros(), radius=2.0, height=3.0, axis=gm.Vector3(0, 0, 1))
        self.assertIn(gm.Vector3(0, 0, 1), c)


class TestBlock(unittest.TestCase):

    def test_contains_point(self):
        b = gm.Block(size=ones(), center=zeros())
        self.assertIn(zeros(), b)


class TestEllipsoid(unittest.TestCase):

    def test_contains_point(self):
        e = gm.Ellipsoid(size=ones(), center=zeros())
        self.assertIn(zeros(), e)

    # TODO(chogan): Allow python to read this member after it's computed in C
    def test_inverse_semi_axes(self):
        pass


# TODO(chogan): Need to call a method that takes a CGO to test its typemap
# class TestCompoundGeometricObject(unittest.TestCase):

#     def test_convert_pylist_to_cgo_list(self):
#         s = gm.Sphere(center=zeros(), radius=2.0)
#         c = gm.Cylinder(center=zeros(), radius=2.0, height=4.0, axis=gm.Vector3(0, 0, 1))
#         b = gm.Block(size=ones(), center=zeros())

#         cgo = gm.CompoundGeometricObject(center=zeros(), component_objects=[s, c, b])


if __name__ == '__main__':
    unittest.main()
