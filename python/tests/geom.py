import unittest
from math import pi
import numpy as np
import meep as mp
import meep.geom as gm


def zeros():
    return gm.Vector3(0, 0, 0)


def ones():
    return gm.Vector3(1, 1, 1)


class TestGeom(unittest.TestCase):

    def test_geometric_object_duplicates_x(self):
        rad = 1
        s = mp.Sphere(rad)
        res = mp.geometric_object_duplicates(mp.Vector3(x=1), 1, 5, s)

        expected = [
            mp.Sphere(rad, center=mp.Vector3(x=1)),
            mp.Sphere(rad, center=mp.Vector3(x=2)),
            mp.Sphere(rad, center=mp.Vector3(x=3)),
            mp.Sphere(rad, center=mp.Vector3(x=4)),
            mp.Sphere(rad, center=mp.Vector3(x=5))
        ]

        for r, e in zip(res, expected):
            self.assertEqual(r.center, e.center)

    def test_geometric_object_duplicates_xyz(self):
        rad = 1
        s = mp.Sphere(rad)
        res = mp.geometric_object_duplicates(mp.Vector3(1, 1, 1), 1, 5, s)

        expected = [
            mp.Sphere(rad, center=mp.Vector3(1, 1, 1)),
            mp.Sphere(rad, center=mp.Vector3(2, 2, 2)),
            mp.Sphere(rad, center=mp.Vector3(3, 3, 3)),
            mp.Sphere(rad, center=mp.Vector3(4, 4, 4)),
            mp.Sphere(rad, center=mp.Vector3(5, 5, 5))
        ]

        for r, e in zip(res, expected):
            self.assertEqual(r.center, e.center)

    def test_geometric_object_duplicates_multiple_objs(self):
        rad = 1
        s = mp.Sphere(rad)
        c = mp.Cylinder(rad)

        res = mp.geometric_object_duplicates(mp.Vector3(1, 1, 1), 1, 5, s, c)

        expected = [
            mp.Sphere(rad, center=mp.Vector3(1, 1, 1)),
            mp.Sphere(rad, center=mp.Vector3(2, 2, 2)),
            mp.Sphere(rad, center=mp.Vector3(3, 3, 3)),
            mp.Sphere(rad, center=mp.Vector3(4, 4, 4)),
            mp.Sphere(rad, center=mp.Vector3(5, 5, 5)),
            mp.Cylinder(rad, center=mp.Vector3(1, 1, 1)),
            mp.Cylinder(rad, center=mp.Vector3(2, 2, 2)),
            mp.Cylinder(rad, center=mp.Vector3(3, 3, 3)),
            mp.Cylinder(rad, center=mp.Vector3(4, 4, 4)),
            mp.Cylinder(rad, center=mp.Vector3(5, 5, 5))
        ]

        for r, e in zip(res, expected):
            self.assertEqual(r.center, e.center)


class TestSphere(unittest.TestCase):

    def test_kwargs_passed_to_parent(self):
        s = gm.Sphere(1.0)
        self.assertEqual(s.material.epsilon_diag, ones())
        self.assertEqual(s.center, zeros())
        self.assertEqual(s.radius, 1)

        s = gm.Sphere(radius=2.0)
        self.assertEqual(s.material.epsilon_diag, ones())
        self.assertEqual(s.center, zeros())
        self.assertEqual(s.radius, 2.0)

        s = gm.Sphere(1.0, center=ones())
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
        self.assertTrue(mp.is_point_in_periodic_object(mp.Vector3(), s))
        self.assertIn(point, s)
        self.assertFalse(gm.Vector3(10, 10, 10) in s)

    def test_shift(self):
        s = gm.Sphere(center=zeros(), radius=2.0)
        self.assertEqual(s.center, gm.Vector3())

        s.shift(gm.Vector3(10))
        self.assertEqual(s.center, gm.Vector3(10, 0, 0))

        s = gm.Sphere(center=gm.Vector3(10, 10), radius=2.0)
        s.shift(gm.Vector3(-10, -10))
        self.assertEqual(s.center, gm.Vector3())

    def test_info(self):
        # Sanity test to ensure that display_geometric_object_info is callable
        s = gm.Sphere(2)
        s.info()
        s.info(2)


class TestCylinder(unittest.TestCase):

    def test_non_neg_height_constructor(self):
        gm.Cylinder(radius=1.0, height=0.0)
        gm.Cylinder(radius=1.0, height=1.0)

        with self.assertRaises(ValueError) as ctx:
            gm.Cylinder(radius=1.0, height=-1)
            self.assertIn("Got -1", ctx.exception)

    def test_non_neg_height_setter(self):
        s = gm.Cylinder(radius=1.0, height=0.0)
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

    def test_wrong_type_exception(self):
        """Test for Issue 180"""
        with self.assertRaises(TypeError) as ctx:
            c = gm.Cylinder(radius=mp.Vector3())


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


class TestMedium(unittest.TestCase):

    def test_D_conductivity(self):
        m = gm.Medium(D_conductivity=2)
        self.assertEqual(m.D_conductivity_diag.x, 2)
        self.assertEqual(m.D_conductivity_diag.y, 2)
        self.assertEqual(m.D_conductivity_diag.z, 2)

    def test_B_conductivity(self):
        m = gm.Medium(B_conductivity=2)
        self.assertEqual(m.B_conductivity_diag.x, 2)
        self.assertEqual(m.B_conductivity_diag.y, 2)
        self.assertEqual(m.B_conductivity_diag.z, 2)

    def test_E_chi2(self):
        m = gm.Medium(E_chi2=2)
        self.assertEqual(m.E_chi2_diag.x, 2)
        self.assertEqual(m.E_chi2_diag.y, 2)
        self.assertEqual(m.E_chi2_diag.z, 2)

    def test_E_chi3(self):
        m = gm.Medium(E_chi3=2)
        self.assertEqual(m.E_chi3_diag.x, 2)
        self.assertEqual(m.E_chi3_diag.y, 2)
        self.assertEqual(m.E_chi3_diag.z, 2)

    def test_H_chi2(self):
        m = gm.Medium(H_chi2=2)
        self.assertEqual(m.H_chi2_diag.x, 2)
        self.assertEqual(m.H_chi2_diag.y, 2)
        self.assertEqual(m.H_chi2_diag.z, 2)

    def test_H_chi3(self):
        m = gm.Medium(H_chi3=2)
        self.assertEqual(m.H_chi3_diag.x, 2)
        self.assertEqual(m.H_chi3_diag.y, 2)
        self.assertEqual(m.H_chi3_diag.z, 2)


class TestVector3(unittest.TestCase):

    def test_use_as_numpy_array(self):
        v = gm.Vector3(10, 10, 10)
        res = np.add(v, np.array([10, 10, 10]))

        self.assertTrue(type(res) is np.ndarray)
        np.testing.assert_array_equal(np.array([20, 20, 20]), res)

    def test_cross(self):
        v1 = mp.Vector3(x=1)
        v2 = mp.Vector3(z=1)
        self.assertEqual(v1.cross(v2), mp.Vector3(y=-1))

        v1 = mp.Vector3(1, 1)
        v2 = mp.Vector3(0, 1, 1)
        self.assertEqual(v1.cross(v2), mp.Vector3(1, -1, 1))

    def test_cdot(self):
        complex_vec1 = mp.Vector3(complex(1, 1), complex(1, 1), complex(1, 1))
        complex_vec2 = mp.Vector3(complex(2, 2), complex(2, 2), complex(2, 2))

        self.assertEqual(complex_vec1.cdot(complex_vec2), 12 + 0j)

    def test_rotate(self):
        axis = mp.Vector3(z=1)
        v = mp.Vector3(x=1)
        res = v.rotate(axis, pi)
        self.assertTrue(res.close(mp.Vector3(x=-1)))

    def test_close(self):
        v1 = mp.Vector3(1e-7)
        v2 = mp.Vector3(1e-8)
        self.assertTrue(v1.close(v2))

        v1 = mp.Vector3(1e-6)
        v2 = mp.Vector3(1e-7)
        self.assertFalse(v1.close(v2))

        v1 = mp.Vector3(1e-10)
        v2 = mp.Vector3(1e-11)
        self.assertTrue(v1.close(v2, tol=1e-10))

    def test_mul_operator(self):

        self.assertEqual(mp.Vector3(2, 2, 2) * 0.5, mp.Vector3(1, 1, 1))
        self.assertEqual(mp.Vector3(1, 1, 1) * mp.Vector3(1, 1, 1), 3)
        self.assertEqual(0.5 * mp.Vector3(2, 2, 2), mp.Vector3(1, 1, 1))


class TestLattice(unittest.TestCase):

    def test_lattice(self):
        # TODO
        pass


class TestMatrix(unittest.TestCase):

    def test_matrix(self):
        # TODO
        pass


if __name__ == '__main__':
    unittest.main()
