import os
import shutil
import sys
import unittest
import warnings
import h5py
import numpy as np
import meep as mp


try:
    unicode
except NameError:
    unicode = str


examples_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', 'examples'))
sys.path.insert(0, examples_dir)


class TestSimulation(unittest.TestCase):

    fname = 'simulation-ez-000200.00.h5'

    def test_interpolate_numbers(self):

        expected = [
            1.0, 1.0909090909090908, 1.1818181818181819, 1.2727272727272727, 1.3636363636363635, 1.4545454545454546, 1.5454545454545454, 1.6363636363636365, 1.7272727272727273, 1.8181818181818181, 1.9090909090909092,
            2.0, 2.090909090909091, 2.1818181818181817, 2.272727272727273, 2.3636363636363638, 2.4545454545454546, 2.5454545454545454, 2.6363636363636362, 2.727272727272727, 2.8181818181818183, 2.909090909090909,
            3.0, 3.090909090909091, 3.1818181818181817, 3.272727272727273, 3.3636363636363638, 3.4545454545454546, 3.5454545454545454, 3.6363636363636362, 3.727272727272727, 3.8181818181818183, 3.909090909090909,
            4.0, 4.090909090909091, 4.181818181818182, 4.2727272727272725, 4.363636363636363, 4.454545454545454, 4.545454545454546, 4.636363636363637, 4.7272727272727275, 4.818181818181818, 4.909090909090909,
            5.0, 5.090909090909091, 5.181818181818182, 5.2727272727272725, 5.363636363636363, 5.454545454545454, 5.545454545454546, 5.636363636363637, 5.7272727272727275, 5.818181818181818, 5.909090909090909,
            6.0, 6.090909090909091, 6.181818181818182, 6.2727272727272725, 6.363636363636363, 6.454545454545454, 6.545454545454546, 6.636363636363637, 6.7272727272727275, 6.818181818181818, 6.909090909090909,
            7.0, 7.090909090909091, 7.181818181818182, 7.2727272727272725, 7.363636363636363, 7.454545454545454, 7.545454545454546, 7.636363636363637, 7.7272727272727275, 7.818181818181818, 7.909090909090909,
            8.0, 8.090909090909092, 8.181818181818182, 8.272727272727273, 8.363636363636363, 8.454545454545455, 8.545454545454545, 8.636363636363637, 8.727272727272727, 8.818181818181818, 8.909090909090908,
            9.0, 9.090909090909092, 9.181818181818182, 9.272727272727273, 9.363636363636363, 9.454545454545455, 9.545454545454545, 9.636363636363637, 9.727272727272727, 9.818181818181818, 9.909090909090908,
            10.0
        ]

        nums = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        result = mp.interpolate(10, nums)

        np.testing.assert_allclose(expected, result)

    def test_interpolate_vectors(self):

        expected = [
            mp.Vector3(),
            mp.Vector3(0.024999999999999842),
            mp.Vector3(0.04999999999999984),
            mp.Vector3(0.07499999999999984),
            mp.Vector3(0.09999999999999984),
            mp.Vector3(0.12499999999999983),
            mp.Vector3(0.14999999999999983),
            mp.Vector3(0.17499999999999982),
            mp.Vector3(0.19999999999999982),
            mp.Vector3(0.2249999999999998),
            mp.Vector3(0.2499999999999998),
            mp.Vector3(0.2749999999999998),
            mp.Vector3(0.2999999999999998),
            mp.Vector3(0.32499999999999984),
            mp.Vector3(0.34999999999999987),
            mp.Vector3(0.3749999999999999),
            mp.Vector3(0.3999999999999999),
            mp.Vector3(0.42499999999999993),
            mp.Vector3(0.44999999999999996),
            mp.Vector3(0.475),
            mp.Vector3(0.5)
        ]

        res = mp.interpolate(19, [mp.Vector3(), mp.Vector3(0.5)])

        np.testing.assert_allclose([v.x for v in expected], [v.x for v in res])
        np.testing.assert_allclose([v.y for v in expected], [v.y for v in res])
        np.testing.assert_allclose([v.z for v in expected], [v.z for v in res])

    def test_arith_sequence(self):

        expected = [
            0.15,
            0.15040080160320599,
            0.15080160320641198,
            0.15120240480961797,
            0.15160320641282396,
            0.15200400801602995,
            0.15240480961923594,
            0.15280561122244193,
            0.15320641282564793,
            0.15360721442885392
        ]

        res = np.linspace(0.15, 0.15 + 0.000400801603206 * 10, num=10, endpoint=False)

        self.assertEqual(len(expected), len(res))
        np.testing.assert_allclose(expected, res)

    def init_simple_simulation(self, **kwargs):
        resolution = 20

        cell = mp.Vector3(10, 10)

        pml_layers = mp.PML(1.0)

        fcen = 1.0
        df = 1.0

        sources = mp.Source(src=mp.GaussianSource(fcen, fwidth=df), center=mp.Vector3(),
                            component=mp.Ez)

        symmetries = [mp.Mirror(mp.X), mp.Mirror(mp.Y)]

        return mp.Simulation(resolution=resolution,
                             cell_size=cell,
                             boundary_layers=[pml_layers],
                             sources=[sources],
                             symmetries=symmetries,
                             **kwargs)

    @unittest.skipIf(not mp.with_mpi(), "MPI specific test")
    def test_mpi(self):
        self.assertGreater(mp.comm.Get_size(), 1)

    def test_use_output_directory_default(self):
        sim = self.init_simple_simulation()
        sim.use_output_directory()
        sim.run(mp.at_end(mp.output_efield_z), until=200)

        output_dir = 'simulation-out'
        self.assertTrue(os.path.exists(os.path.join(output_dir, self.fname)))

        mp.all_wait()
        if mp.am_master():
            shutil.rmtree(output_dir)

    def test_use_output_directory_custom(self):
        sim = self.init_simple_simulation()
        sim.use_output_directory('custom_dir')
        sim.run(mp.at_end(mp.output_efield_z), until=200)

        output_dir = 'custom_dir'
        self.assertTrue(os.path.exists(os.path.join(output_dir, self.fname)))

        mp.all_wait()
        if mp.am_master():
            shutil.rmtree(output_dir)

    def test_at_time(self):
        sim = self.init_simple_simulation()
        sim.run(mp.at_time(100, mp.output_efield_z), until=200)

        fname = 'simulation-ez-000100.00.h5'
        self.assertTrue(os.path.exists(fname))

        mp.all_wait()
        if mp.am_master():
            os.remove(fname)

    def test_after_sources_and_time(self):
        sim = self.init_simple_simulation()

        done = [False]

        def _done(sim, todo):
            done[0] = True

        sim.run(mp.after_sources_and_time(1, _done), until_after_sources=2)

        self.assertTrue(done[0])

    def test_with_prefix(self):
        sim = self.init_simple_simulation()
        sim.run(mp.with_prefix('test_prefix-', mp.at_end(mp.output_efield_z)), until=200)

        fname = 'test_prefix-simulation-ez-000200.00.h5'
        self.assertTrue(os.path.exists(fname))

        mp.all_wait()
        if mp.am_master():
            os.remove(fname)

    def test_extra_materials(self):
        sim = self.init_simple_simulation()
        sim.extra_materials = [mp.Medium(epsilon=5), mp.Medium(epsilon=10)]
        sim.run(mp.at_end(lambda sim: None), until=5)

    def test_require_dimensions(self):
        sim = self.init_simple_simulation()
        self.assertIsNone(sim.structure)
        self.assertEqual(sim.dimensions, 3)

        sim.require_dimensions()
        sim._init_structure(k=mp.Vector3())
        self.assertEqual(sim.structure.gv.dim, mp.D2)

    def test_infer_dimensions(self):
        sim = self.init_simple_simulation()
        self.assertEqual(sim.dimensions, 3)
        sim._init_structure()
        self.assertEqual(sim.dimensions, 2)

    def test_in_volume(self):
        sim = self.init_simple_simulation()
        sim.filename_prefix = 'test_in_volume'
        vol = mp.Volume(mp.Vector3(), size=mp.Vector3(x=2))
        sim.run(mp.at_end(mp.in_volume(vol, mp.output_efield_z)), until=200)

    def test_in_point(self):
        sim = self.init_simple_simulation(filename_prefix='test_in_point')
        fn = sim.filename_prefix + '-ez-000200.00.h5'
        pt = mp.Vector3()
        sim.run(mp.at_end(mp.in_point(pt, mp.output_efield_z)), until=200)
        self.assertTrue(os.path.exists(fn))

        mp.all_wait()
        if mp.am_master():
            os.remove(fn)

    def test_epsilon_input_file(self):
        sim = self.init_simple_simulation()
        eps_input_fname = 'cyl-ellipsoid-eps-ref.h5'
        eps_input_dir = os.path.join(os.path.abspath(os.path.realpath(os.path.dirname(__file__))),
                                     '..', '..', 'libmeepgeom')
        eps_input_path = os.path.join(eps_input_dir, eps_input_fname)
        sim.epsilon_input_file = eps_input_path

        sim.run(until=200)
        fp = sim.get_field_point(mp.Ez, mp.Vector3(x=1))

        self.assertAlmostEqual(fp, -0.002989654055823199 + 0j)

        # Test unicode file name for Python 2
        if sys.version_info[0] == 2:
            sim = self.init_simple_simulation(epsilon_input_file=unicode(eps_input_path))
            sim.run(until=200)
            fp = sim.get_field_point(mp.Ez, mp.Vector3(x=1))
            self.assertAlmostEqual(fp, -0.002989654055823199 + 0j)

    def test_set_materials(self):

        def change_geom(sim):
            t = sim.meep_time()
            fn = t * 0.02
            geom = [mp.Cylinder(radius=3, material=mp.Medium(index=3.5), center=mp.Vector3(fn, fn)),
                    mp.Ellipsoid(size=mp.Vector3(1, 2, mp.inf), center=mp.Vector3(fn, fn))]

            sim.set_materials(geometry=geom)

        c = mp.Cylinder(radius=3, material=mp.Medium(index=3.5))
        e = mp.Ellipsoid(size=mp.Vector3(1, 2, mp.inf))

        sources = mp.Source(src=mp.GaussianSource(1, fwidth=0.1), component=mp.Hz, center=mp.Vector3())
        symmetries = [mp.Mirror(mp.X, -1), mp.Mirror(mp.Y, -1)]

        sim = mp.Simulation(cell_size=mp.Vector3(10, 10),
                            geometry=[c, e],
                            boundary_layers=[mp.PML(1.0)],
                            sources=[sources],
                            symmetries=symmetries,
                            resolution=16)

        eps = {'arr1': None, 'arr2': None}

        def get_arr1(sim):
            eps['arr1'] = sim.get_array(mp.Volume(mp.Vector3(), mp.Vector3(10, 10)),
                                        component=mp.Dielectric)

        def get_arr2(sim):
            eps['arr2'] = sim.get_array(mp.Volume(mp.Vector3(), mp.Vector3(10, 10)),
                                        component=mp.Dielectric)

        sim.run(mp.at_time(50, get_arr1), mp.at_time(100, change_geom),
                mp.at_end(get_arr2), until=200)

        self.assertFalse(np.array_equal(eps['arr1'], eps['arr2']))

    def test_modal_volume_in_box(self):
        sim = self.init_simple_simulation()
        sim.run(until=200)
        vol = sim.fields.total_volume()
        self.assertAlmostEqual(sim.fields.modal_volume_in_box(vol),
                               sim.modal_volume_in_box())

        vol = mp.Volume(mp.Vector3(), size=mp.Vector3(1, 1, 1))
        self.assertAlmostEqual(sim.fields.modal_volume_in_box(vol.swigobj),
                               sim.modal_volume_in_box(vol))

    def test_in_box_volumes(self):
        sim = self.init_simple_simulation()
        sim.run(until=200)

        tv = sim.fields.total_volume()
        v = mp.Volume(mp.Vector3(), size=mp.Vector3(5, 5))

        sim.electric_energy_in_box(tv)
        sim.electric_energy_in_box(v)
        sim.flux_in_box(mp.X, tv)
        sim.flux_in_box(mp.X, v)
        sim.magnetic_energy_in_box(tv)
        sim.magnetic_energy_in_box(v)
        sim.field_energy_in_box(tv)
        sim.field_energy_in_box(v)

    def test_load_dump_structure(self):
        resolution = 10
        cell = mp.Vector3(10, 10)
        pml_layers = mp.PML(1.0)
        fcen = 1.0
        df = 1.0
        sources = mp.Source(src=mp.GaussianSource(fcen, fwidth=df), center=mp.Vector3(),
                            component=mp.Hz)
        geometry = mp.Cylinder(0.2, material=mp.Medium(index=3))

        sim = mp.Simulation(resolution=resolution,
                            cell_size=cell,
                            default_material=mp.Medium(index=1),
                            geometry=[geometry],
                            boundary_layers=[pml_layers],
                            sources=[sources])

        sim.run(until=200)
        ref_field = sim.get_field_point(mp.Hz, mp.Vector3(z=2))
        dump_fn = 'test_load_dump_structure.h5'
        sim.dump_structure(dump_fn)

        sim = mp.Simulation(resolution=resolution,
                            cell_size=cell,
                            default_material=mp.Medium(index=1),
                            geometry=[],
                            boundary_layers=[pml_layers],
                            sources=[sources],
                            load_structure=dump_fn)
        sim.run(until=200)
        field = sim.get_field_point(mp.Hz, mp.Vector3(z=2))

        self.assertAlmostEqual(ref_field, field)

        mp.all_wait()
        if mp.am_master():
            os.remove(dump_fn)

    def test_get_array_output(self):
        sim = self.init_simple_simulation()
        sim.symmetries = []
        sim.geometry = [mp.Cylinder(0.2, material=mp.Medium(index=3))]
        sim.filename_prefix = 'test_get_array_output'
        sim.run(until=20)

        mp.output_epsilon(sim)
        mp.output_efield_z(sim)
        mp.output_tot_pwr(sim)
        mp.output_efield(sim)

        eps_arr = sim.get_epsilon()
        efield_z_arr = sim.get_efield_z()
        energy_arr = sim.get_tot_pwr()
        efield_arr = sim.get_efield()

        fname_fmt = "test_get_array_output-{}-000020.00.h5"

        with h5py.File(fname_fmt.format('eps'), 'r') as f:
            eps = f['eps'].value

        with h5py.File(fname_fmt.format('ez'), 'r') as f:
            efield_z = f['ez'].value

        with h5py.File(fname_fmt.format('energy'), 'r') as f:
            energy = f['energy'].value

        with h5py.File(fname_fmt.format('e'), 'r') as f:
            ex = f['ex'].value
            ey = f['ey'].value
            ez = f['ez'].value
            efield = np.stack([ex, ey, ez], axis=-1)

        np.testing.assert_allclose(eps, eps_arr)
        np.testing.assert_allclose(efield_z, efield_z_arr)
        np.testing.assert_allclose(energy, energy_arr)
        np.testing.assert_allclose(efield, efield_arr)

    def test_synchronized_magnetic(self):
        # Issue 309
        cell = mp.Vector3(16, 8, 0)

        geometry = [mp.Block(mp.Vector3(1e20, 1, 1e20),
                             center=mp.Vector3(0, 0),
                             material=mp.Medium(epsilon=12))]

        sources = [mp.Source(mp.ContinuousSource(frequency=0.15),
                             component=mp.Ez,
                             center=mp.Vector3(-7, 0))]

        pml_layers = [mp.PML(1.0)]
        resolution = 10

        sim = mp.Simulation(
            cell_size=cell,
            boundary_layers=pml_layers,
            geometry=geometry,
            sources=sources,
            resolution=resolution
        )

        sim.run(mp.synchronized_magnetic(mp.output_bfield_y), until=10)

    def test_harminv_warnings(self):

        def check_warnings(sim, h, should_warn=True):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                sim.run(mp.after_sources(h), until_after_sources=5)

                if should_warn:
                    self.assertEqual(len(w), 1)
                    self.assertIn("Harminv", str(w[-1].message))
                else:
                    self.assertEqual(len(w), 0)

        sources = [mp.Source(src=mp.GaussianSource(1, fwidth=1), center=mp.Vector3(), component=mp.Ez)]
        sim = mp.Simulation(cell_size=mp.Vector3(10, 10), resolution=10, sources=sources)
        h = mp.Harminv(mp.Ez, mp.Vector3(), 1.4, 0.5)
        check_warnings(sim, h)

        sim = mp.Simulation(cell_size=mp.Vector3(10, 10), resolution=10, sources=sources)
        h = mp.Harminv(mp.Ez, mp.Vector3(), 0.5, 0.5)
        check_warnings(sim, h)

        sim = mp.Simulation(cell_size=mp.Vector3(10, 10), resolution=10, sources=sources)
        h = mp.Harminv(mp.Ez, mp.Vector3(), 1, 1)
        check_warnings(sim, h, should_warn=False)

    def test_vec_constructor(self):

        def assert_one(v):
            self.assertEqual(v.z(), 1)

        def assert_two(v):
            self.assertEqual(v.x(), 1)
            self.assertEqual(v.y(), 2)

        def assert_three(v):
            assert_two(v)
            self.assertEqual(v.z(), 3)

        def assert_raises(it, err):
            with self.assertRaises(err):
                mp.vec(it)

        v1 = mp.vec(1)
        assert_one(v1)
        v2 = mp.vec(1, 2)
        assert_two(v2)
        v3 = mp.vec(1, 2, 3)
        assert_three(v3)
        mp.vec()

        with self.assertRaises(TypeError):
            mp.vec(1, 2, 3, 4)

        def check_iterable(one, two, three, four):
            v1 = mp.vec(one)
            assert_one(v1)
            v2 = mp.vec(two)
            assert_two(v2)
            v3 = mp.vec(three)
            assert_three(v3)
            assert_raises(four, NotImplementedError)

        check_iterable([1], [1, 2], [1, 2, 3], [1, 2, 3, 4])
        check_iterable((1,), (1, 2), (1, 2, 3), (1, 2, 3, 4))
        check_iterable(np.array([1.]),
                       np.array([1., 2.]),
                       np.array([1., 2., 3.]),
                       np.array([1., 2., 3., 4.]))

        with self.assertRaises(TypeError):
            mp.vec([1, 2], 3)

        with self.assertRaises(TypeError):
            mp.vec(1, [2, 3])

    def test_epsilon_warning(self):

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            from materials_library import Si
            self.assertEqual(len(w), 0)

        from materials_library import Mo
        geom = [mp.Sphere(radius=0.2, material=Mo)]
        sim = self.init_simple_simulation(geometry=geom)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            sim.run(until=5)
            self.assertGreater(len(w), 0)
            self.assertIn("Epsilon", str(w[0].message))

        from materials_library import SiO2
        geom = [mp.Sphere(radius=0.2, material=SiO2)]
        sim = self.init_simple_simulation(geometry=geom)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            sim.run(until=5)
            self.assertEqual(len(w), 1)
            self.assertNotIn("Epsilon", str(w[0].message))


if __name__ == '__main__':
    unittest.main()
