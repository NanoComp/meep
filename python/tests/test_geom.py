import math
import unittest
import warnings

import meep.geom as gm
import numpy as np

import meep as mp


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
            mp.Sphere(rad, center=mp.Vector3(x=5)),
            mp.Sphere(rad, center=mp.Vector3(x=4)),
            mp.Sphere(rad, center=mp.Vector3(x=3)),
            mp.Sphere(rad, center=mp.Vector3(x=2)),
            mp.Sphere(rad, center=mp.Vector3(x=1)),
        ]

        for r, e in zip(res, expected):
            self.assertEqual(r.center, e.center)

    def test_geometric_object_duplicates_xyz(self):
        rad = 1
        s = mp.Sphere(rad)
        res = mp.geometric_object_duplicates(mp.Vector3(1, 1, 1), 1, 5, s)

        expected = [
            mp.Sphere(rad, center=mp.Vector3(5, 5, 5)),
            mp.Sphere(rad, center=mp.Vector3(4, 4, 4)),
            mp.Sphere(rad, center=mp.Vector3(3, 3, 3)),
            mp.Sphere(rad, center=mp.Vector3(2, 2, 2)),
            mp.Sphere(rad, center=mp.Vector3(1, 1, 1)),
        ]

        for r, e in zip(res, expected):
            self.assertEqual(r.center, e.center)

    def test_geometric_objects_duplicates(self):
        rad = 1
        s = mp.Sphere(rad)
        c = mp.Cylinder(rad)

        res = mp.geometric_objects_duplicates(mp.Vector3(1, 1, 1), 1, 5, [s, c])

        expected = [
            mp.Sphere(rad, center=mp.Vector3(5, 5, 5)),
            mp.Sphere(rad, center=mp.Vector3(4, 4, 4)),
            mp.Sphere(rad, center=mp.Vector3(3, 3, 3)),
            mp.Sphere(rad, center=mp.Vector3(2, 2, 2)),
            mp.Sphere(rad, center=mp.Vector3(1, 1, 1)),
            mp.Cylinder(rad, center=mp.Vector3(5, 5, 5)),
            mp.Cylinder(rad, center=mp.Vector3(4, 4, 4)),
            mp.Cylinder(rad, center=mp.Vector3(3, 3, 3)),
            mp.Cylinder(rad, center=mp.Vector3(2, 2, 2)),
            mp.Cylinder(rad, center=mp.Vector3(1, 1, 1)),
        ]
        for r, e in zip(res, expected):
            self.assertEqual(r.center, e.center)

    def test_geometric_objects_lattice_duplicates(self):
        geometry_lattice = mp.Lattice(
            size=mp.Vector3(1, 7),
            basis1=mp.Vector3(math.sqrt(3) / 2, 0.5),
            basis2=mp.Vector3(math.sqrt(3) / 2, -0.5),
        )
        eps = 12
        r = 0.2

        geometry = [mp.Cylinder(r, material=mp.Medium(epsilon=eps))]
        geometry = mp.geometric_objects_lattice_duplicates(geometry_lattice, geometry)

        med = mp.Medium(epsilon=12)

        expected = [
            mp.Cylinder(0.2, material=med, center=mp.Vector3(y=3.0)),
            mp.Cylinder(0.2, material=med, center=mp.Vector3(y=2.0)),
            mp.Cylinder(0.2, material=med, center=mp.Vector3(y=1.0)),
            mp.Cylinder(0.2, material=med, center=mp.Vector3(y=0.0)),
            mp.Cylinder(0.2, material=med, center=mp.Vector3(y=-1.0)),
            mp.Cylinder(0.2, material=med, center=mp.Vector3(y=-2.0)),
            mp.Cylinder(0.2, material=med, center=mp.Vector3(y=-3.0)),
        ]

        for exp, res in zip(expected, geometry):
            self.assertEqual(exp.center, res.center)
            self.assertEqual(exp.radius, res.radius)

    def test_cartesian_to_lattice(self):
        lattice = mp.Lattice(
            size=mp.Vector3(1, 7),
            basis1=mp.Vector3(math.sqrt(3) / 2, 0.5),
            basis2=mp.Vector3(math.sqrt(3) / 2, -0.5),
        )
        res = mp.cartesian_to_lattice(lattice.basis * mp.Vector3(1), lattice)
        self.assertEqual(res, mp.Vector3(1))


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
            gm.Sphere(invalid="This is not allowed")
        with self.assertRaises(TypeError):
            gm.Sphere(radius=1.0, oops="Nope")

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

        shifted = s.shift(gm.Vector3(10))
        self.assertEqual(shifted.center, gm.Vector3(10, 0, 0))
        self.assertEqual(shifted.radius, 2.0)

        s = gm.Sphere(center=gm.Vector3(10, 10), radius=2.0)
        shifted = s.shift(gm.Vector3(-10, -10))
        self.assertEqual(shifted.center, gm.Vector3())

        s = gm.Sphere(center=zeros(), radius=2.0)
        s += mp.Vector3(5, 5)
        self.assertEqual(s.center, mp.Vector3(5, 5))

        s = gm.Sphere(center=zeros(), radius=2.0)
        new_sphere = s + mp.Vector3(5, 5)
        self.assertEqual(new_sphere.center, mp.Vector3(5, 5))
        self.assertEqual(s.center, zeros())

        s = gm.Sphere(center=zeros(), radius=2.0)
        new_sphere = mp.Vector3(5, 5) + s
        self.assertEqual(new_sphere.center, mp.Vector3(5, 5))
        self.assertEqual(s.center, zeros())

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

    def test_wrong_type_exception(self):
        """Test for Issue 180"""
        with self.assertRaises(TypeError):
            gm.Cylinder(radius=mp.Vector3())


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


class TestPrism(unittest.TestCase):
    def test_contains_point(self):
        vertices = [
            gm.Vector3(-1, 1),
            gm.Vector3(1, 1),
            gm.Vector3(1, -1),
            gm.Vector3(-1, -1),
        ]
        p = gm.Prism(vertices, height=1)
        self.assertIn(zeros(), p)
        self.assertNotIn(gm.Vector3(2, 2), p)


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

    def test_check_material_frequencies(self):
        mat = mp.Medium(valid_freq_range=mp.FreqRange(min=10, max=20))
        invalid_sources = [
            [mp.Source(mp.GaussianSource(5, fwidth=1), mp.Ez, mp.Vector3())],
            [mp.Source(mp.ContinuousSource(10, fwidth=1), mp.Ez, mp.Vector3())],
            [mp.Source(mp.GaussianSource(10, width=1), mp.Ez, mp.Vector3())],
            [mp.Source(mp.GaussianSource(20, width=1), mp.Ez, mp.Vector3())],
        ]

        cell_size = mp.Vector3(5, 5)
        resolution = 5

        def check_warnings(sim, should_warn=True):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                sim.run(until=5)

                if should_warn:
                    self.assertEqual(len(w), 1)
                    self.assertIn("material", str(w[-1].message))
                else:
                    self.assertEqual(len(w), 0)

        geom = [mp.Sphere(0.2, material=mat)]

        for s in invalid_sources:
            # Check for invalid extra_materials
            sim = mp.Simulation(
                cell_size=cell_size,
                resolution=resolution,
                sources=s,
                extra_materials=[mat],
            )
            check_warnings(sim)

            # Check for invalid geometry materials
            sim = mp.Simulation(
                cell_size=cell_size, resolution=resolution, sources=s, geometry=geom
            )
            check_warnings(sim)

        valid_sources = [
            [mp.Source(mp.GaussianSource(15, fwidth=1), mp.Ez, mp.Vector3())],
            [mp.Source(mp.ContinuousSource(15, width=5), mp.Ez, mp.Vector3())],
        ]

        for s in valid_sources:
            sim = mp.Simulation(
                cell_size=cell_size,
                resolution=resolution,
                sources=s,
                extra_materials=[mat],
            )
            check_warnings(sim, False)

        # Check DFT frequencies

        # Invalid extra_materials
        sim = mp.Simulation(
            cell_size=cell_size,
            resolution=resolution,
            sources=valid_sources[0],
            extra_materials=[mat],
        )
        fregion = mp.FluxRegion(
            center=mp.Vector3(0, 1), size=mp.Vector3(2, 2), direction=mp.X
        )
        sim.add_flux(18, 6, 2, fregion, decimation_factor=1)
        check_warnings(sim)

        # Invalid geometry material
        sim = mp.Simulation(
            cell_size=cell_size,
            resolution=resolution,
            sources=valid_sources[0],
            geometry=geom,
        )
        sim.add_flux(18, 6, 2, fregion, decimation_factor=1)
        check_warnings(sim)

    def test_transform(self):

        e_sus = [
            mp.LorentzianSusceptibility(
                sigma_diag=mp.Vector3(1, 2, 3), sigma_offdiag=mp.Vector3(12, 13, 14)
            ),
            mp.DrudeSusceptibility(
                sigma_diag=mp.Vector3(1, 2, 3), sigma_offdiag=mp.Vector3(12, 13, 14)
            ),
        ]

        mat = mp.Medium(
            epsilon_diag=mp.Vector3(1, 2, 3),
            epsilon_offdiag=mp.Vector3(12, 13, 14),
            E_susceptibilities=e_sus,
        )

        rot_angle = math.radians(23.9)
        rot_matrix = mp.Matrix(
            mp.Vector3(math.cos(rot_angle), math.sin(rot_angle), 0),
            mp.Vector3(-math.sin(rot_angle), math.cos(rot_angle), 0),
            mp.Vector3(0, 0, 1),
        )
        mat.transform(rot_matrix)

        expected_diag = mp.Vector3(-7.72552, 10.72552, 3)
        expected_offdiag = mp.Vector3(7.69024, 6.21332, 18.06640)

        self.assertTrue(mat.epsilon_diag.close(expected_diag, tol=4))
        self.assertTrue(mat.epsilon_offdiag.close(expected_offdiag, tol=4))
        self.assertEqual(mat.mu_diag, mp.Vector3(1, 1, 1))
        self.assertEqual(mat.mu_offdiag, mp.Vector3())
        self.assertEqual(len(mat.E_susceptibilities), 2)
        self.assertTrue(
            mat.E_susceptibilities[0].sigma_diag.close(expected_diag, tol=4)
        )
        self.assertTrue(
            mat.E_susceptibilities[0].sigma_offdiag.close(expected_offdiag, tol=4)
        )
        self.assertTrue(
            mat.E_susceptibilities[1].sigma_diag.close(expected_diag, tol=4)
        )
        self.assertTrue(
            mat.E_susceptibilities[1].sigma_offdiag.close(expected_offdiag, tol=4)
        )


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
        res = v.rotate(axis, math.pi)
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

    def test_rotate_lattice(self):
        axis = mp.Vector3(1)
        v = mp.Vector3(2, 2, 2)
        lattice = mp.Lattice(size=mp.Vector3(1, 1))
        res = v.rotate_lattice(axis, 3, lattice)
        self.assertTrue(
            res.close(mp.Vector3(2.0, -2.262225009320625, -1.6977449770811563))
        )

    def test_rotate_reciprocal(self):
        axis = mp.Vector3(1)
        v = mp.Vector3(2, 2, 2)
        lattice = mp.Lattice(size=mp.Vector3(1, 1))
        res = v.rotate_reciprocal(axis, 3, lattice)
        self.assertTrue(
            res.close(mp.Vector3(2.0, -2.262225009320625, -1.6977449770811563))
        )

    def test_complex_norm(self):
        # issue #722
        v = mp.Vector3(1, 1j, 0)
        self.assertAlmostEqual(v.norm(), math.sqrt(2))


class TestLattice(unittest.TestCase):
    def test_basis(self):
        lattice = mp.Lattice(
            size=mp.Vector3(1, 7),
            basis1=mp.Vector3(math.sqrt(3) / 2, 0.5),
            basis2=mp.Vector3(math.sqrt(3) / 2, -0.5),
        )
        b = lattice.basis
        exp = mp.Matrix(
            mp.Vector3(0.8660254037844388, 0.5000000000000001),
            mp.Vector3(0.8660254037844388, -0.5000000000000001),
            mp.Vector3(z=1.0),
        )

        for e, r in zip([exp.c1, exp.c2, exp.c3], [b.c1, b.c2, b.c3]):
            self.assertTrue(e.close(r))


class TestMatrix(unittest.TestCase):

    identity = mp.Matrix(mp.Vector3(1), mp.Vector3(y=1), mp.Vector3(z=1))

    def matrix_eq(self, exp, res):
        for e, r in zip([exp.c1, exp.c2, exp.c3], [res.c1, res.c2, res.c3]):
            self.assertEqual(e, r)

    def matrix_close(self, exp, res):
        for e, r in zip([exp.c1, exp.c2, exp.c3], [res.c1, res.c2, res.c3]):
            self.assertTrue(e.close(r))

    def test_indexing(self):
        self.assertEqual(self.identity[0][0], 1)
        self.assertEqual(self.identity[1][1], 1)
        self.assertEqual(self.identity[2][2], 1)
        self.assertEqual(self.identity[0][1], 0)

    def test_row(self):
        self.assertEqual(self.identity.row(0), self.identity.c1)
        self.assertEqual(self.identity.row(1), self.identity.c2)
        self.assertEqual(self.identity.row(2), self.identity.c3)

    def test_mm_mult(self):
        m1 = mp.Matrix(mp.Vector3(1, 2, 3), mp.Vector3(4, 5, 6), mp.Vector3(7, 8, 9))
        m2 = mp.Matrix(mp.Vector3(9, 8, 7), mp.Vector3(6, 5, 4), mp.Vector3(3, 2, 1))
        res = m1 * m2
        exp = mp.Matrix(
            mp.Vector3(90.0, 114.0, 138.0),
            mp.Vector3(54.0, 69.0, 84.0),
            mp.Vector3(18.0, 24.0, 30.0),
        )

        self.matrix_eq(exp, res)

    def test_add(self):
        result = self.identity + self.identity
        self.assertEqual(result.row(0), mp.Vector3(x=2))
        self.assertEqual(result.row(1), mp.Vector3(y=2))
        self.assertEqual(result.row(2), mp.Vector3(z=2))

    def test_sub(self):
        ones_matrix = mp.Matrix(ones(), ones(), ones())
        result = ones_matrix - ones_matrix
        self.assertEqual(result.row(0), zeros())
        self.assertEqual(result.row(1), zeros())
        self.assertEqual(result.row(2), zeros())

    def test_mv_mult(self):
        lattice = mp.Lattice(
            size=mp.Vector3(1, 7),
            basis1=mp.Vector3(math.sqrt(3) / 2, 0.5),
            basis2=mp.Vector3(math.sqrt(3) / 2, -0.5),
        )
        res = lattice.basis * mp.Vector3(1)
        exp = mp.Vector3(0.8660254037844388, 0.5000000000000001)
        self.assertTrue(res.close(exp))

    def test_scale(self):
        m = mp.Matrix(
            mp.Vector3(90.0, 114.0, 138.0),
            mp.Vector3(54.0, 69.0, 84.0),
            mp.Vector3(18.0, 24.0, 30.0),
        )
        res = m.scale(0.5)
        exp = mp.Matrix(
            mp.Vector3(45.0, 57.0, 69.0),
            mp.Vector3(27.0, 34.5, 42.0),
            mp.Vector3(9.0, 12.0, 15.0),
        )
        self.matrix_eq(exp, res)

        self.matrix_eq(exp, m * 0.5)
        self.matrix_eq(exp, 0.5 * m)

    def test_determinant(self):
        m = mp.Matrix(mp.Vector3(2), mp.Vector3(y=2), mp.Vector3(z=2))

        m1 = mp.Matrix(mp.Vector3(1, 2, 3), mp.Vector3(4, 5, 6), mp.Vector3(7, 8, 9))

        self.assertEqual(8, m.determinant())
        self.assertEqual(0, m1.determinant())

    def test_transpose(self):
        m = mp.Matrix(mp.Vector3(1, 2, 3), mp.Vector3(4, 5, 6), mp.Vector3(7, 8, 9))
        exp = mp.Matrix(mp.Vector3(1, 4, 7), mp.Vector3(2, 5, 8), mp.Vector3(3, 6, 9))

        self.matrix_eq(exp, m.transpose())

    def test_inverse(self):
        self.matrix_eq(self.identity, self.identity.inverse())

        lattice = mp.Lattice(
            size=mp.Vector3(1, 7),
            basis1=mp.Vector3(math.sqrt(3) / 2, 0.5),
            basis2=mp.Vector3(math.sqrt(3) / 2, -0.5),
        )

        res = lattice.basis.inverse()
        exp = mp.Matrix(
            mp.Vector3(0.5773502691896256, 0.5773502691896256, -0.0),
            mp.Vector3(0.9999999999999998, -0.9999999999999998, -0.0),
            mp.Vector3(-0.0, -0.0, 1.0),
        )

        self.matrix_close(exp, res)

    def test_get_rotation_matrix(self):
        result = mp.get_rotation_matrix(ones(), 5)
        self.assertTrue(
            result.c1.close(
                mp.Vector3(0.5224414569754843, -0.3148559165969717, 0.7924144596214877)
            )
        )
        self.assertTrue(
            result.c2.close(
                mp.Vector3(0.7924144596214877, 0.5224414569754843, -0.3148559165969717)
            )
        )
        self.assertTrue(
            result.c3.close(
                mp.Vector3(-0.3148559165969717, 0.7924144596214877, 0.5224414569754843)
            )
        )

    def test_conj(self):
        m = mp.Matrix(mp.Vector3(x=1 + 1j), mp.Vector3(y=1 + 1j), mp.Vector3(z=1 + 1j))
        result = m.conj()
        self.assertEqual(result.c1, mp.Vector3(x=1 - 1j))
        self.assertEqual(result.c2, mp.Vector3(y=1 - 1j))
        self.assertEqual(result.c3, mp.Vector3(z=1 - 1j))

    def test_adjoint(self):
        m = mp.Matrix(mp.Vector3(1 + 1j), mp.Vector3(1 + 1j), mp.Vector3(1 + 1j))
        getH_result = m.getH()
        H_result = m.H
        self.assertEqual(getH_result.c1, mp.Vector3(1 - 1j, 1 - 1j, 1 - 1j))
        self.assertEqual(getH_result.c2, mp.Vector3())
        self.assertEqual(getH_result.c3, mp.Vector3())
        np.testing.assert_allclose(getH_result, H_result)

    def test_to_numpy_array(self):
        m = mp.Matrix(mp.Vector3(1 + 1j), mp.Vector3(1 + 1j), mp.Vector3(1 + 1j))
        adjoint = m.H
        m_arr = np.array(m)
        np_adjoint = m_arr.conj().transpose()
        np.testing.assert_allclose(adjoint, np_adjoint)


if __name__ == "__main__":
    unittest.main()
