import unittest

import meep as mp
from meep.geom import Cylinder, Vector3
from meep.source import EigenModeSource, ContinuousSource, GaussianSource


class TestEigenModeSource(unittest.TestCase):

    def test_eig_lattice_defaults(self):
        src = ContinuousSource(5.0)
        center = Vector3()

        default_lattice = EigenModeSource(src, center)
        self.assertEqual(default_lattice.eig_lattice_size, Vector3())
        self.assertEqual(default_lattice.eig_lattice_center, Vector3())

        elc = Vector3(1, 1, 1)
        els = Vector3(1, 1, 1)
        custom_lattice = EigenModeSource(src, center, eig_lattice_center=elc, eig_lattice_size=els)
        self.assertEqual(custom_lattice.eig_lattice_size, els)
        self.assertEqual(custom_lattice.eig_lattice_center, elc)


class TestSourceTypemaps(unittest.TestCase):

    expected_msg = "Expected a meep.source.SourceTime or a meep.src_time\n"

    def setUp(self):

        def dummy_eps(v):
            return 1.0

        gv = mp.voltwo(16, 16, 10)
        gv.center_origin()
        sym = mp.mirror(mp.Y, gv)
        the_structure = mp.structure(gv, dummy_eps, mp.pml(2), sym)
        objects = []
        objects.append(Cylinder(1))
        mp.set_materials_from_geometry(the_structure, objects)
        self.f = mp.fields(the_structure)
        self.v = mp.volume(mp.vec(1.1, 0.0), mp.vec(0.0, 0.0))

    def test_typemap_swig(self):
        src = mp.gaussian_src_time(0.15, 0.1)
        self.f.add_volume_source(mp.Ez, src, self.v)

    def test_typemap_py(self):
        src = GaussianSource(0.15, 0.1)
        self.f.add_volume_source(mp.Ez, src, self.v)

    def test_typemap_swig_raises(self):
        src = mp.gaussian_src_time(0.15, 0.1)
        self.assertTrue(src.is_equal(src))

        with self.assertRaises(TypeError) as error:
            src.is_equal(mp.vec())
            self.assertEqual(error.exception.message, self.expected_msg)

    def test_typemap_py_raises(self):
        src = GaussianSource(0.15, 0.1)
        self.assertTrue(src.swigobj.is_equal(src))

        with self.assertRaises(TypeError) as error:
            src.swigobj.is_equal(Vector3())
            self.assertEqual(error.exception.message, self.expected_msg)


if __name__ == '__main__':
    unittest.main()
    mp.all_wait()
