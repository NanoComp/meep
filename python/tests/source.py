from __future__ import division

import math
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


class TestSourceTime(unittest.TestCase):

    def test_source_wavelength(self):
        g_src = GaussianSource(wavelength=10)
        c_src = ContinuousSource(wavelength=10)

        self.assertAlmostEqual(1. / 10., g_src.frequency)
        self.assertAlmostEqual(1. / 10., c_src.frequency)

    def test_source_frequency(self):
        g_src = GaussianSource(10)
        c_src = ContinuousSource(10)

        self.assertEqual(10, g_src.frequency)
        self.assertEqual(10, c_src.frequency)

        with self.assertRaises(ValueError):
            GaussianSource()

        with self.assertRaises(ValueError):
            ContinuousSource()


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

    def test_custom_source(self):
        n = 3.4
        w = 1
        r = 1
        pad = 4
        dpml = 2
        sxy = 2 * (r + w + pad + dpml)

        cell = mp.Vector3(sxy, sxy)

        geometry = [
            mp.Cylinder(r + w, material=mp.Medium(index=n)),
            mp.Cylinder(r, material=mp.air)
        ]

        boundary_layers = [mp.PML(dpml)]
        resolution = 10
        fcen = 0.15
        df = 0.1

        # Bump function
        def my_src_func(t):
            if t > 0 and t < 2:
                return math.exp(-1 / (1 - ((t - 1)**2)))
            return 0j

        sources = [mp.Source(src=mp.CustomSource(src_func=my_src_func, end_time=100),
                             component=mp.Ez, center=mp.Vector3(r + 0.1))]

        symmetries = [mp.Mirror(mp.Y)]

        sim = mp.Simulation(cell_size=cell,
                            resolution=resolution,
                            geometry=geometry,
                            boundary_layers=boundary_layers,
                            sources=sources,
                            symmetries=symmetries)

        h = mp.Harminv(mp.Ez, mp.Vector3(r + 0.1), fcen, df)
        sim.run(mp.after_sources(h), until_after_sources=200)
        fp = sim.get_field_point(mp.Ez, mp.Vector3(1))

        self.assertAlmostEqual(fp, -0.021997617628500023 + 0j)

if __name__ == '__main__':
    unittest.main()
