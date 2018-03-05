from __future__ import division, print_function

import glob
import math
import os
import re
import sys
import unittest

import h5py
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.optimize import ridder
import meep as mp
from meep import mpb


class TestModeSolver(unittest.TestCase):

    data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))
    examples_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', 'examples'))
    sys.path.insert(0, examples_dir)

    def setUp(self):
        """Store the test name and register a function to clean up all the
        generated h5 files."""

        self.filename_prefix = self.id().split('.')[-1]
        print()
        print(self.filename_prefix)
        print('=' * 24)

        def rm_h5():
            mp.all_wait()
            if mp.am_master():
                for f in glob.glob("{}*.h5".format(self.filename_prefix)):
                    os.remove(f)

        self.addCleanup(rm_h5)

    def init_solver(self, geom=True):
        num_bands = 8
        k_points = [
            mp.Vector3(),
            mp.Vector3(0.5),
            mp.Vector3(0.5, 0.5),
            mp.Vector3()
        ]

        geometry = [mp.Cylinder(0.2, material=mp.Medium(epsilon=12))] if geom else []
        k_points = mp.interpolate(4, k_points)
        geometry_lattice = mp.Lattice(size=mp.Vector3(1, 1))
        resolution = 32

        return mpb.ModeSolver(
            num_bands=num_bands,
            k_points=k_points,
            geometry=geometry,
            geometry_lattice=geometry_lattice,
            resolution=resolution,
            filename_prefix=self.filename_prefix,
            deterministic=True,
            tolerance=1e-12
        )

    def test_resolution(self):
        ms = self.init_solver()
        self.assertEqual([32, 32, 32], ms.resolution)

        ms.resolution = mp.Vector3(16, 16, 32)
        self.assertEqual([16, 16, 32], ms.resolution)

        with self.assertRaises(TypeError):
            ms.resolution = [32, 32, 32]

    def test_list_split(self):
        k_points = [
            mp.Vector3(),
            mp.Vector3(0.5),
            mp.Vector3(0.5, 0.5),
            mp.Vector3()
        ]

        k_points = mp.interpolate(4, k_points)
        ms = mpb.ModeSolver()
        k_split = ms.list_split(k_points, 1, 0)

        expected_list = [
            mp.Vector3(),
            mp.Vector3(0.10000000000000003),
            mp.Vector3(0.20000000000000004),
            mp.Vector3(0.30000000000000004),
            mp.Vector3(0.4),
            mp.Vector3(0.5),
            mp.Vector3(0.5, 0.10000000000000003),
            mp.Vector3(0.5, 0.20000000000000004),
            mp.Vector3(0.5, 0.30000000000000004),
            mp.Vector3(0.5, 0.4),
            mp.Vector3(0.5, 0.5),
            mp.Vector3(0.4, 0.4),
            mp.Vector3(0.30000000000000004, 0.30000000000000004),
            mp.Vector3(0.2, 0.2),
            mp.Vector3(0.1, 0.1),
            mp.Vector3(0.0, 0.0),
        ]

        self.assertEqual(k_split[0], 0)

        for res, exp in zip(k_split[1], expected_list):
            self.assertTrue(res.close(exp))

    def test_first_brillouin_zone(self):
        ms = self.init_solver()

        res = []
        for k in ms.k_points:
            res.append(ms.first_brillouin_zone(k))

        expected = [
            mp.Vector3(0.0, 0.0, 0.0),
            mp.Vector3(0.10000000000000003, 0.0, 0.0),
            mp.Vector3(0.20000000000000004, 0.0, 0.0),
            mp.Vector3(0.30000000000000004, 0.0, 0.0),
            mp.Vector3(0.4, 0.0, 0.0),
            mp.Vector3(0.5, 0.0, 0.0),
            mp.Vector3(0.5, 0.10000000000000003, 0.0),
            mp.Vector3(0.5, 0.20000000000000004, 0.0),
            mp.Vector3(0.5, 0.30000000000000004, 0.0),
            mp.Vector3(0.5, 0.4, 0.0),
            mp.Vector3(0.5, 0.5, 0.0),
            mp.Vector3(0.4, 0.4, 0.0),
            mp.Vector3(0.30000000000000004, 0.30000000000000004, 0.0),
            mp.Vector3(0.2, 0.2, 0.0),
            mp.Vector3(0.1, 0.1, 0.0),
            mp.Vector3(0.0, 0.0, 0.0),
        ]

        for e, r in zip(expected, res):
            self.assertTrue(e.close(r))

    def check_band_range_data(self, expected_brd, result, places=3, tol=1e-3):
        for exp, res in zip(expected_brd, result):
            # Compare min freqs
            self.assertAlmostEqual(exp[0][0], res[0][0], places=places)
            # Compare min k
            msg = "expected {}, got {}"
            self.assertTrue(exp[0][1].close(res[0][1], tol=tol),
                            msg=msg.format(exp[0][1], res[0][1]))
            # Compare max freqs
            self.assertAlmostEqual(exp[1][0], res[1][0], places=places)
            # Compare max k
            self.assertTrue(exp[1][1].close(res[1][1], tol=tol),
                            msg=msg.format(exp[1][1], res[1][1]))

    def check_freqs(self, expected_freqs, result):
        for exp, res in zip(expected_freqs, result):
            for r, e in zip(res, exp):
                self.assertAlmostEqual(r, e, places=3)

    def check_gap_list(self, expected_gap_list, result):
        self.check_freqs(expected_gap_list, result)

    def check_fields_against_h5(self, ref_path, field):
        with h5py.File(ref_path, 'r') as ref:
            # Reshape the reference data into a component-wise 1d array like
            # [x1,y1,z1,x2,y2,z2,etc.] to match the layout of `field`
            ref_x = ref['x.r'].value + ref['x.i'].value * 1j
            ref_y = ref['y.r'].value + ref['y.i'].value * 1j
            ref_z = ref['z.r'].value + ref['z.i'].value * 1j

            ref_arr = np.zeros(field.shape[0], dtype=np.complex128)
            ref_arr[0::3] = ref_x.ravel()
            ref_arr[1::3] = ref_y.ravel()
            ref_arr[2::3] = ref_z.ravel()

            self.compare_arrays(ref_arr, field)

    def compare_arrays(self, exp, res, tol=1e-3):
        exp_1d = exp.ravel()
        res_1d = res.ravel()

        norm_exp = np.linalg.norm(exp_1d)
        norm_res = np.linalg.norm(res_1d)

        if norm_exp == 0:
            self.assertEqual(norm_res, 0)
        else:
            diff = np.linalg.norm(res_1d - exp_1d) / norm_exp
            self.assertLess(diff, tol)

    def compare_h5_files(self, ref_path, res_path, tol=1e-3):
        mp.all_wait()
        with h5py.File(ref_path) as ref:
            with h5py.File(res_path, 'r') as res:
                for k in ref.keys():
                    if k == 'description':
                        self.assertEqual(ref[k].value, res[k].value)
                    else:
                        self.compare_arrays(ref[k].value, res[k].value, tol=tol)

    def test_update_band_range_data(self):
        brd = []
        freqs = [0.0, 1.0000000001231053, 1.0000000001577114, 1.000000000183077,
                 1.0000000003647922, 1.4142135627385737, 1.4142135630373556,
                 1.4142135634172286]
        kpoint = mp.Vector3()

        expected = [
            ((0.0, mp.Vector3()), (0.0, mp.Vector3())),
            ((1.0000000001231053, mp.Vector3()), (1.0000000001231053, mp.Vector3())),
            ((1.0000000001577114, mp.Vector3()), (1.0000000001577114, mp.Vector3())),
            ((1.000000000183077, mp.Vector3()), (1.000000000183077, mp.Vector3())),
            ((1.0000000003647922, mp.Vector3()), (1.0000000003647922, mp.Vector3())),
            ((1.4142135627385737, mp.Vector3()), (1.4142135627385737, mp.Vector3())),
            ((1.4142135630373556, mp.Vector3()), (1.4142135630373556, mp.Vector3())),
            ((1.4142135634172286, mp.Vector3()), (1.4142135634172286, mp.Vector3())),
        ]

        ms = mpb.ModeSolver()
        res = ms.update_band_range_data(brd, freqs, kpoint)
        self.check_band_range_data(expected, res)

    def test_run_te_no_geometry(self):

        expected_freqs = [0.0, 1.0, 1.0000000000000004, 1.0000000000000013,
                          1.0000000000000016, 1.4142135623730958, 1.4142135623730965,
                          1.414213562373097]

        expected_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.7071067811865476, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.5000000000350678, mp.Vector3(0.5, 0.0, 0.0)),
             (1.0000000000464613, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.7071067811884221, mp.Vector3(0.5, 0.5, 0.0)),
             (1.1180339887555244, mp.Vector3(0.5, 0.0, 0.0))),
            ((0.7071067811901163, mp.Vector3(0.5, 0.5, 0.0)),
             (1.1313708499266775, mp.Vector3(0.2, 0.2, 0.0))),
            ((1.000000000001028, mp.Vector3(0.0, 0.0, 0.0)),
             (1.5811388300846416, mp.Vector3(0.5, 0.5, 0.0))),
            ((1.1180339887687103, mp.Vector3(0.5, 0.0, 0.0)),
             (1.5811388300858549, mp.Vector3(0.5, 0.5, 0.0))),
            ((1.2806248475163993, mp.Vector3(0.20000000000000004, 0.0, 0.0)),
             (1.5811388300892486, mp.Vector3(0.5, 0.5, 0.0))),
            ((1.4142135623752818, mp.Vector3(0.0, 0.0, 0.0)),
             (1.8027756376524453, mp.Vector3(0.5, 0.0, 0.0))),
        ]

        ms = self.init_solver(geom=False)
        ms.tolerance = 1e-7
        ms.run_te()

        self.check_band_range_data(expected_brd, ms.band_range_data)

        for e, r in zip(expected_freqs, ms.all_freqs[-1]):
            self.assertAlmostEqual(e, r, places=3)

        self.assertEqual(len(ms.gap_list), 0)

    def test_run_te(self):

        expected_freqs = [0.0, 0.5527092320101986, 0.7732265593069759, 0.773229883948054,
                          0.9229965195855876, 1.0001711176882833, 1.0001720032257042,
                          1.092820931747826]

        # Obtained by passing NULL to set_maxwell_dielectric for epsilon_mean_func
        # in mpb/mpb/medium.c.
        expected_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.49683586474489877, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.4415884497225449, mp.Vector3(0.5, 0.0, 0.0)),
             (0.5931405141160885, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.5931535863117832, mp.Vector3(0.5, 0.5, 0.0)),
             (0.7732265593069759, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.6791690130757013, mp.Vector3(0.5, 0.5, 0.0)),
             (0.80968915516771, mp.Vector3(0.30000000000000004, 0.30000000000000004, 0.0))),
            ((0.8241814443502151, mp.Vector3(0.5, 0.30000000000000004, 0.0)),
             (0.9229965195855876, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.8819770916660669, mp.Vector3(0.5, 0.5, 0.0)),
             (1.0291597050650205, mp.Vector3(0.5, 0.0, 0.0))),
            ((0.8819818134421844, mp.Vector3(0.5, 0.5, 0.0)),
             (1.086072932359415, mp.Vector3(0.5, 0.0, 0.0))),
            ((1.0878689635052692, mp.Vector3(0.5, 0.0, 0.0)),
             (1.1119173707556929, mp.Vector3(0.5, 0.5, 0.0))),
        ]

        expected_gap_list = [
            (0.0022038709776893727, 0.5931405141160885, 0.5931535863117832),
            (1.7739824912427062, 0.80968915516771, 0.8241814443502151),
            (0.1652326724344101, 1.086072932359415, 1.0878689635052692),
        ]

        ms = self.init_solver()
        ms.run_te()

        self.check_band_range_data(expected_brd, ms.band_range_data)
        for e, r in zip(expected_freqs, ms.all_freqs[-1]):
            self.assertAlmostEqual(e, r, places=3)

        self.check_gap_list(expected_gap_list, ms.gap_list)

    def test_run_tm(self):
        expected_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.28094795352537366, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.4171142493246634, mp.Vector3(0.5, 0.0, 0.0)),
             (0.5460267793370319, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.49619745276546123, mp.Vector3(0.5, 0.5, 0.0)),
             (0.5576576362977246, mp.Vector3(0.5, 0.0, 0.0))),
            ((0.5520955864542503, mp.Vector3(0.0, 0.0, 0.0)),
             (0.7133951516423513, mp.Vector3(0.5, 0.0, 0.0))),
            ((0.7413109657068678, mp.Vector3(0.5, 0.0, 0.0)),
             (0.8594741069571248, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.8295176150251525, mp.Vector3(0.30000000000000004, 0.30000000000000004, 0.0)),
             (0.8783155473463833, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.8625159053811312, mp.Vector3(0.5, 0.0, 0.0)),
             (0.9511074539064021, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.8793510958294801, mp.Vector3(0.5, 0.5, 0.0)),
             (1.0825923841450287, mp.Vector3(0.0, 0.0, 0.0))),
        ]

        ms = self.init_solver()
        ms.run_tm()

        self.check_band_range_data(expected_brd, ms.band_range_data)

    def _test_get_field(self, field):
        ms = self.init_solver()
        ms.run_te()
        get_field_func = getattr(ms, "get_{}field".format(field))
        fix_phase_func = getattr(mpb, "fix_{}field_phase".format(field))
        fix_phase_func(ms, ms.num_bands)
        fields = get_field_func(ms.num_bands)

        ref_fname = "tutorial-{}.k16.b08.te.h5".format(field)
        ref_path = os.path.join(self.data_dir, ref_fname)

        self.check_fields_against_h5(ref_path, fields)

    def test_get_bfield(self):
        self._test_get_field('b')

    def test_get_efield(self):
        self._test_get_field('e')

    def test_get_dfield(self):
        self._test_get_field('d')

    def test_get_hfield(self):
        self._test_get_field('h')

    def test_output_field_to_file(self):
        fname = 'tutorial-epsilon.h5'
        data_path = os.path.join(self.data_dir, fname)

        ms = self.init_solver()
        ms.run_te()

        res_path = self.filename_prefix + '-epsilon.h5'
        self.compare_h5_files(data_path, res_path)

    def test_compute_field_energy(self):
        ms = self.init_solver()
        ms.run_te()
        ms.get_dfield(8)
        res = ms.compute_field_energy()

        expected = [1.0000000000000002, 1.726755206037815e-5, 0.4999827324479414,
                    1.7267552060375955e-5, 0.4999827324479377, 0.0, 0.0]

        self.compare_arrays(np.array(expected), np.array(res))

    def test_output_efield_z(self):
        ms = self.init_solver()
        ms.run_tm()
        mpb.fix_efield_phase(ms, 8)
        mpb.output_efield_z(ms, 8)

        ref_fname = 'tutorial-e.k16.b08.z.tm.h5'
        ref_path = os.path.join(self.data_dir, ref_fname)
        res_path = re.sub('tutorial', ms.filename_prefix, ref_fname)
        self.compare_h5_files(ref_path, res_path)

    def test_triangular_lattice(self):

        expected_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.2746902258623623, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0))),
            ((0.44533108084715683, mp.Vector3(0.0, 0.5, 0.0)),
             (0.5605181423162835, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.4902389149027666, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0)),
             (0.5605607947797747, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.5932960873585144, mp.Vector3(0.0, 0.0, 0.0)),
             (0.7907195974443698, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0))),
            ((0.790832076332758, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0)),
             (0.8374511167537562, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.8375948528443267, mp.Vector3(0.0, 0.0, 0.0)),
             (0.867200926490345, mp.Vector3(-0.2, 0.39999999999999997, 0.0))),
            ((0.8691349955739203, mp.Vector3(-0.13333333333333336, 0.4333333333333333, 0.0)),
             (0.9941291022664892, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.8992499095547049, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0)),
             (1.098318352915696, mp.Vector3(0.0, 0.0, 0.0))),
        ]

        ms = self.init_solver()
        ms.geometry_lattice = mp.Lattice(
            size=mp.Vector3(1, 1),
            basis1=mp.Vector3(math.sqrt(3) / 2, 0.5),
            basis2=mp.Vector3(math.sqrt(3) / 2, -0.5)
        )

        k_points = [
            mp.Vector3(),
            mp.Vector3(y=0.5),
            mp.Vector3(-1 / 3, 1 / 3),
            mp.Vector3()
        ]

        ms.k_points = mp.interpolate(4, k_points)
        ms.run_tm()

        self.check_band_range_data(expected_brd, ms.band_range_data)

    def test_maximize_first_tm_gap(self):

        def first_tm_gap(r):
            ms.geometry = [mp.Cylinder(r, material=mp.Medium(epsilon=12))]
            ms.run_tm()
            return -1 * ms.retrieve_gap(1)

        ms = self.init_solver()
        ms.num_bands = 2
        ms.mesh_size = 7

        result = minimize_scalar(first_tm_gap, method='bounded', bounds=[0.1, 0.5], tol=0.1)
        expected = 39.10325687542367
        self.assertAlmostEqual(expected, result.fun * -1, places=2)

    def test_anisotropic_2d_gap(self):

        expected_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.2213165540404889, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.23068427462181276, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.21691192680454566, mp.Vector3(0.5, 0.0, 0.0)),
             (0.319020283148369, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.30110868065428525, mp.Vector3(0.5, 0.0, 0.0)),
             (0.3648353129125716, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.30701621910773247, mp.Vector3(0.5, 0.5, 0.0)),
             (0.3852513546698513, mp.Vector3(0.1, 0.1, 0.0))),
            ((0.341835260571013, mp.Vector3(0.5, 0.30000000000000004, 0.0)),
             (0.391421048600237, mp.Vector3(0.5, 0.10000000000000003, 0.0))),
            ((0.34982139739904294, mp.Vector3(0.5, 0.5, 0.0)),
             (0.4075642914057991, mp.Vector3(0.4, 0.0, 0.0))),
            ((0.3963465468598276, mp.Vector3(0.1, 0.1, 0.0)),
             (0.4772237204606825, mp.Vector3(0.5, 0.5, 0.0))),

        ]

        ms = self.init_solver()
        ms.geometry = [mp.Cylinder(0.3, material=mp.Medium(epsilon_diag=mp.Vector3(1, 1, 12)))]
        ms.default_material = mp.Medium(epsilon_diag=mp.Vector3(12, 12, 1))
        ms.num_bands = 8
        ms.run()

        self.check_band_range_data(expected_brd, ms.band_range_data)

    def test_point_defect_state(self):

        ms = self.init_solver()
        ms.geometry_lattice = mp.Lattice(size=mp.Vector3(5, 5))
        ms.geometry = [mp.Cylinder(0.2, material=mp.Medium(epsilon=12))]

        ms.geometry = mp.geometric_objects_lattice_duplicates(ms.geometry_lattice, ms.geometry)
        ms.geometry.append(mp.Cylinder(0.2, material=mp.air))

        ms.resolution = 16
        ms.k_points = [mp.Vector3(0.5, 0.5)]

        ms.num_bands = 50
        ms.run_tm()

        mpb.fix_efield_phase(ms, 25)
        mpb.output_efield_z(ms, 25)

        mpb.fix_dfield_phase(ms, 25)
        ms.get_dfield(25)
        ms.compute_field_energy()
        c = mp.Cylinder(1.0, material=mp.air)
        e = ms.compute_energy_in_objects([c])
        self.assertAlmostEqual(0.6227482574427817, e, places=3)

        ms.num_bands = 1
        ms.target_freq = (0.2812 + 0.4174) / 2
        ms.tolerance = 1e-8
        ms.run_tm()

        expected_brd = [
            ((0.37730041222979477, mp.Vector3(0.5, 0.5, 0.0)),
             (0.37730041222979477, mp.Vector3(0.5, 0.5, 0.0))),
        ]

        self.check_band_range_data(expected_brd, ms.band_range_data)

        old_geometry = ms.geometry  # save the 5x5 grid with a missing rod

        def rootfun(eps):
            ms.geometry = old_geometry + [mp.Cylinder(0.2, material=mp.Medium(epsilon=eps))]
            ms.run_tm()
            return ms.get_freqs()[0] - 0.314159

        rooteps = ridder(rootfun, 1, 12)
        rootval = rootfun(rooteps)

        self.assertAlmostEqual(5.288830111797463, rooteps, places=3)
        self.assertAlmostEqual(9.300716530269426e-9, rootval, places=3)

    def test_output_charge_density(self):
        ms = self.init_solver()
        ms.run_te()
        mpb.fix_efield_phase(ms, 8)
        mpb.output_charge_density(ms, 8)

        ref_fn = 'tutorial-C.k16.b08.te.h5'
        ref_path = os.path.join(self.data_dir, ref_fn)

        res_path = re.sub('tutorial', ms.filename_prefix, ref_fn)
        self.compare_h5_files(ref_path, res_path)

    def test_bragg_sine(self):
        from mpb_bragg_sine import ms
        ms.deterministic = True
        ms.tolerance = 1e-12
        ms.filename_prefix = self.filename_prefix
        ms.run_tm()

        expected_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.19477466366820298, mp.Vector3(0.5, 0.0, 0.0))),
            ((0.306403026230844, mp.Vector3(0.5, 0.0, 0.0)),
             (0.4687748525867193, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.5466257501317459, mp.Vector3(0.0, 0.0, 0.0)),
             (0.7316504426541637, mp.Vector3(0.5, 0.0, 0.0))),
            ((0.7842615905093812, mp.Vector3(0.5, 0.0, 0.0)),
             (0.9893486155437277, mp.Vector3(0.0, 0.0, 0.0))),
            ((1.0240548648147831, mp.Vector3(0.0, 0.0, 0.0)),
             (1.244098004202588, mp.Vector3(0.5, 0.0, 0.0))),
            ((1.266656686185507, mp.Vector3(0.5, 0.0, 0.0)),
             (1.4970379696966822, mp.Vector3(0.0, 0.0, 0.0))),
            ((1.5115800994652884, mp.Vector3(0.0, 0.0, 0.0)),
             (1.7488359039910502, mp.Vector3(0.5, 0.0, 0.0))),
            ((1.7581683208483643, mp.Vector3(0.5, 0.0, 0.0)),
             (1.9999072007119787, mp.Vector3(0.0, 0.0, 0.0))),
        ]

        self.check_band_range_data(expected_brd, ms.band_range_data)

    def test_bragg(self):
        from mpb_bragg import ms
        ms.deterministic = True
        ms.tolerance = 1e-12
        ms.filename_prefix = self.filename_prefix
        ms.run_tm()
        mpb.fix_hfield_phase(ms, 8)
        mpb.output_hfield_y(ms, 8)

        ref_fn = 'bragg-h.k01.b08.y.tm.h5'
        ref_path = os.path.join(self.data_dir, ref_fn)
        res_path = re.sub('bragg', self.filename_prefix, ref_fn)

        self.compare_h5_files(ref_path, res_path)

    def test_diamond(self):
        from mpb_diamond import ms
        ms.deterministic = True
        ms.filename_prefix = self.filename_prefix
        ms.tolerance = 1e-12

        ms.run(mpb.output_at_kpoint(mp.Vector3(0, 0.625, 0.375), mpb.fix_dfield_phase,
                                    mpb.output_dpwr))

        expected_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.39455107895905156, mp.Vector3(0.25, 0.75, 0.5))),
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.3967658014080592, mp.Vector3(0.29999999999999993, 0.75, 0.45000000000000007))),
            ((0.4423707668172989, mp.Vector3(0.0, 0.5, 0.0)),
             (0.5955899630254676, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.44355516512198145, mp.Vector3(0.0, 0.5, 0.0)),
             (0.5958191312898851, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.5030135895148902, mp.Vector3(0.0, 0.6, 0.4)),
             (0.5958386856926985, mp.Vector3(0.0, 0.0, 0.0))),


        ]

        self.check_band_range_data(expected_brd, ms.band_range_data)

        ref_fn = 'diamond-dpwr.k06.b05.h5'
        ref_path = os.path.join(self.data_dir, ref_fn)
        res_path = re.sub('diamond', self.filename_prefix, ref_fn)
        self.compare_h5_files(ref_path, res_path)

    def test_hole_slab(self):
        from mpb_hole_slab import ms
        ms.deterministic = True
        ms.filename_prefix = self.filename_prefix
        ms.k_points = [mp.Vector3(1 / -3, 1 / 3)]
        ms.tolerance = 1e-12

        ms.run_zeven()
        mpb.fix_hfield_phase(ms, 9)
        mpb.output_hfield_z(ms, 9)

        ref_fn = 'hole-slab-h.k01.b09.z.zeven.h5'
        ref_path = os.path.join(self.data_dir, ref_fn)
        res_path = re.sub('hole-slab', self.filename_prefix, ref_fn)
        ms.display_eigensolver_stats()
        self.compare_h5_files(ref_path, res_path)

    def test_honey_rods(self):
        from mpb_honey_rods import ms
        ms.deterministic = True
        ms.filename_prefix = self.filename_prefix
        ms.tolerance = 1e-12

        expected_tm_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.3351167660354989, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0))),
            ((0.3351850759916969, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0)),
             (0.42984811237816406, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.5751709345431462, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0)),
             (0.7255897672261712, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.6918627724774271, mp.Vector3(0.0, 0.5, 0.0)),
             (0.747622077830657, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0))),
            ((0.7443374497087805, mp.Vector3(-0.06666666666666667, 0.06666666666666667, 0.0)),
             (0.7793792212092525, mp.Vector3(0.0, 0.5, 0.0))),
            ((0.7852786984418492, mp.Vector3(0.0, 0.0, 0.0)),
             (0.8193652861712535, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0))),
            ((0.7856577771856611, mp.Vector3(0.0, 0.0, 0.0)),
             (0.9122560439014182, mp.Vector3(0.0, 0.5, 0.0))),
            ((1.0540350508135123, mp.Vector3(0.0, 0.5, 0.0)),
             (1.1492769389234725, mp.Vector3(0.0, 0.0, 0.0))),

        ]

        ms.run_tm()
        self.check_band_range_data(expected_tm_brd, ms.band_range_data)

        expected_te_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.5535093489972593, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0))),
            ((0.5203183590840945, mp.Vector3(0.0, 0.5, 0.0)),
             (0.7278447515454929, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.576335859651312, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0)),
             (0.7880878930618354, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.8161730293674944, mp.Vector3(0.0, 0.5, 0.0)),
             (0.9209611432140968, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.8385562359606971, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0)),
             (0.9220849425898038, mp.Vector3(0.0, 0.0, 0.0))),
            ((1.0168656683915511, mp.Vector3(0.0, 0.0, 0.0)),
             (1.1083536673418435, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0))),
            ((1.0184507253059425, mp.Vector3(0.0, 0.0, 0.0)),
             (1.159370227370719, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0))),
            ((1.1636719050364361, mp.Vector3(-0.2, 0.2, 0.0)),
             (1.2433411839870618, mp.Vector3(0.0, 0.0, 0.0))),
        ]

        ms.run_te()
        self.check_band_range_data(expected_te_brd, ms.band_range_data)

    def test_line_defect(self):
        from mpb_line_defect import ms, k_points
        ms.deterministic = True
        ms.filename_prefix = self.filename_prefix
        ms.tolerance = 1e-12

        ms.run_tm(mpb.output_at_kpoint(k_points[len(k_points) // 2]),
                  mpb.fix_efield_phase, mpb.output_efield_z)

        ref_fn = 'line-defect-e.k04.b12.z.tm.h5'
        ref_path = os.path.join(self.data_dir, ref_fn)
        res_path = re.sub('line-defect', self.filename_prefix, ref_fn)
        self.compare_h5_files(ref_path, res_path)

    def test_sq_rods(self):
        from mpb_sq_rods import ms
        ms.deterministic = True
        ms.filename_prefix = self.filename_prefix
        ms.tolerance = 1e-12

        ms.run_te()

        expected_te_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.5036058015219026, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.4446229134706744, mp.Vector3(0.5, 0.0, 0.0)),
             (0.5943440245519593, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.5943566394470321, mp.Vector3(0.5, 0.5, 0.0)),
             (0.7808428121911926, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.6793887413076383, mp.Vector3(0.5, 0.5, 0.0)),
             (0.8173893719542897, mp.Vector3(0.30000000000000004, 0.30000000000000004, 0.0))),
            ((0.83045738223392, mp.Vector3(0.30000000000000004, 0.30000000000000004, 0.0)),
             (0.9243716830585584, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.8957817684117546, mp.Vector3(0.5, 0.5, 0.0)),
             (1.0331104139200438, mp.Vector3(0.5, 0.0, 0.0))),
            ((0.8957868745330811, mp.Vector3(0.5, 0.5, 0.0)),
             (1.0958021492221048, mp.Vector3(0.5, 0.0, 0.0))),
            ((1.097416809585406, mp.Vector3(0.5, 0.0, 0.0)),
             (1.1280127648119964, mp.Vector3(0.5, 0.5, 0.0))),
        ]

        self.check_band_range_data(expected_te_brd, ms.band_range_data)

        ms.run_tm()

        expected_tm_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.285905779119655, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.42065733839975095, mp.Vector3(0.5, 0.0, 0.0)),
             (0.5503360754831277, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.5029830978365978, mp.Vector3(0.5, 0.5, 0.0)),
             (0.5671632878128118, mp.Vector3(0.5, 0.0, 0.0))),
            ((0.5613397939889757, mp.Vector3(0.0, 0.0, 0.0)),
             (0.7200918204563993, mp.Vector3(0.5, 0.0, 0.0))),
            ((0.7472029910480836, mp.Vector3(0.5, 0.0, 0.0)),
             (0.874359380500508, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.8404509697526803, mp.Vector3(0.30000000000000004, 0.30000000000000004, 0.0)),
             (0.8833173725822788, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.8770118718330763, mp.Vector3(0.5, 0.0, 0.0)),
             (0.9653253808981632, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.8929933495598104, mp.Vector3(0.5, 0.5, 0.0)),
             (1.089377682009333, mp.Vector3(0.0, 0.0, 0.0))),
        ]

        self.check_band_range_data(expected_tm_brd, ms.band_range_data)

    def test_strip(self):
        from mpb_strip import ms
        ms.deterministic = True
        ms.filename_prefix = self.filename_prefix
        ms.tolerance = 1e-12

        ms.run(mpb.display_zparities, mpb.display_yparities)

        y_parities = ms.mode_solver.compute_yparities()
        z_parities = ms.mode_solver.compute_zparities()

        expected_y_parities = [-0.9997979443175137, 1.0000061871738222, -1.000010781704281, -0.9997880312884855]
        expected_z_parities = [0.9992335747085693, -0.9955122771195629, -0.9970929846091117, -0.995110556144587]

        for e, r in zip(expected_y_parities, y_parities):
            self.assertAlmostEqual(e, r, places=3)

        for e, r in zip(expected_z_parities, z_parities):
            self.assertAlmostEqual(e, r, places=3)

        omega = 1 / 1.55

        kvals = ms.find_k(mp.NO_PARITY, omega, 1, ms.num_bands, mp.Vector3(1), 1e-3,
                          omega * 3.45, omega * 0.1, omega * 4, mpb.output_poynting_x,
                          mpb.display_yparities, mpb.display_group_velocities)

        expected_kvals = [
            1.066321795284513,
            1.0186792189943261,
            0.8398943502679427,
            0.7990426389486213
        ]

        for e, r in zip(expected_kvals, kvals):
            self.assertAlmostEqual(e, r, places=3)

        ref_fn = 'strip-flux.v.k01.b04.x.h5'
        ref_path = os.path.join(self.data_dir, ref_fn)
        res_path = re.sub('strip', self.filename_prefix, ref_fn)

        self.compare_h5_files(ref_path, res_path)

    def test_tri_holes(self):
        from mpb_tri_holes import ms
        ms.deterministic = True
        ms.filename_prefix = self.filename_prefix
        ms.tolerance = 1e-12

        ms.run_te()

        expected_te_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.2993049473117476, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0))),
            ((0.4924342823622065, mp.Vector3(0.0, 0.5, 0.0)),
             (0.6568362683499375, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.5269710506448809, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0)),
             (0.7156232200212518, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.6568031427446027, mp.Vector3(0.0, 0.5, 0.0)),
             (0.7578382217502109, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0))),
            ((0.7383774303752574, mp.Vector3(0.0, 0.0, 0.0)),
             (0.7988168792802597, mp.Vector3(0.0, 0.5, 0.0))),
            ((0.8259787164701536, mp.Vector3(0.0, 0.0, 0.0)),
             (0.9629215441012396, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0))),
            ((0.8271634538840886, mp.Vector3(0.0, 0.0, 0.0)),
             (0.9834563303529568, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0))),
            ((0.9984200611839882, mp.Vector3(-0.26666666666666666, 0.26666666666666666, 0.0)),
             (1.0411551252079034, mp.Vector3(0.0, 0.0, 0.0))),

        ]

        self.check_band_range_data(expected_te_brd, ms.band_range_data)

        ms.run_tm()

        expected_tm_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.28009156817399916, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0))),
            ((0.28015523913784407, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0)),
             (0.3985126081046686, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.4390817228448606, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0)),
             (0.49288810189980625, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.49336847788268695, mp.Vector3(0.0, 0.0, 0.0)),
             (0.5808701365268192, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0))),
            ((0.581035246804371, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0)),
             (0.6824860801372987, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.682531744671499, mp.Vector3(0.0, 0.0, 0.0)),
             (0.7011061593213783, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0))),
            ((0.6920145742134771, mp.Vector3(0.0, 0.0, 0.0)),
             (0.7841042622393081, mp.Vector3(0.0, 0.4, 0.0))),
            ((0.7980077872594108, mp.Vector3(0.0, 0.4, 0.0)),
             (0.8982239424823442, mp.Vector3(0.0, 0.0, 0.0))),

        ]

        self.check_band_range_data(expected_tm_brd, ms.band_range_data)

    def test_tri_rods(self):
        from mpb_tri_rods import ms
        ms.deterministic = True
        ms.tolerance = 1e-12
        ms.filename_prefix = self.filename_prefix

        ms.run_tm(mpb.output_at_kpoint(mp.Vector3(1 / -3, 1 / 3), mpb.fix_efield_phase,
                  mpb.output_efield_z))

        ref_fn = 'tri-rods-e.k11.b08.z.tm.h5'
        ref_path = os.path.join(self.data_dir, ref_fn)
        res_path = re.sub('tri-rods', self.filename_prefix, ref_fn)

        self.compare_h5_files(ref_path, res_path)

        ms.run_te()

        expected_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.49123581186757737, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0))),
            ((0.4730722390280754, mp.Vector3(0.0, 0.5, 0.0)),
             (0.5631059378714038, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0))),
            ((0.5631505198559038, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0)),
             (0.7939289395839766, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.7676614799039024, mp.Vector3(0.0, 0.5, 0.0)),
             (0.8214230044191525, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0))),
            ((0.865194814441525, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0)),
             (1.0334130018594276, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.8652307994936862, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0)),
             (1.0334230910419813, mp.Vector3(0.0, 0.0, 0.0))),
            ((1.021367669109368, mp.Vector3(0.0, 0.5, 0.0)),
             (1.115966990757518, mp.Vector3(0.0, 0.0, 0.0))),
            ((1.108662658537423, mp.Vector3(-0.26666666666666666, 0.26666666666666666, 0.0)),
             (1.1168107191255379, mp.Vector3(0.0, 0.0, 0.0))),
        ]

        self.check_band_range_data(expected_brd, ms.band_range_data)

    def test_subpixel_averaging(self):
        ms = self.init_solver()
        ms.tolerance = 1e-12
        ms.run_te()

        expected_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.49683586474489877, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.4415884497225449, mp.Vector3(0.5, 0.0, 0.0)),
             (0.5931405141160885, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.5931535863117832, mp.Vector3(0.5, 0.5, 0.0)),
             (0.7732265593069759, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.6791690130757013, mp.Vector3(0.5, 0.5, 0.0)),
             (0.80968915516771, mp.Vector3(0.30000000000000004, 0.30000000000000004, 0.0))),
            ((0.8241814443502151, mp.Vector3(0.5, 0.30000000000000004, 0.0)),
             (0.9229965195855876, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.8819770916660669, mp.Vector3(0.5, 0.5, 0.0)),
             (1.0291597050650205, mp.Vector3(0.5, 0.0, 0.0))),
            ((0.8819818134421844, mp.Vector3(0.5, 0.5, 0.0)),
             (1.086072932359415, mp.Vector3(0.5, 0.0, 0.0))),
            ((1.0878689635052692, mp.Vector3(0.5, 0.0, 0.0)),
             (1.1119173707556929, mp.Vector3(0.5, 0.5, 0.0))),
        ]

        expected_gap_list = [
            (0.0022038709776893727, 0.5931405141160885, 0.5931535863117832),
            (1.7739824912427062, 0.80968915516771, 0.8241814443502151),
            (0.1652326724344101, 1.086072932359415, 1.0878689635052692)
        ]

        ref_fn = 'subpixel_avg-epsilon.h5'
        ref_path = os.path.join(self.data_dir, ref_fn)
        res_path = re.sub('subpixel_avg', self.filename_prefix, ref_fn)

        self.compare_h5_files(ref_path, res_path)
        self.check_band_range_data(expected_brd, ms.band_range_data)
        self.check_gap_list(expected_gap_list, ms.gap_list)

    def test_run_te_with_mu_material(self):
        ms = self.init_solver(geom=False)
        ms.geometry = [mp.Cylinder(0.2, material=mp.Medium(mu=5))]

        expected_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.4165291233037574, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.47328232348733695, mp.Vector3(0.5, 0.0, 0.0)),
             (0.6699867281290507, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.6301802646818523, mp.Vector3(0.5, 0.5, 0.0)),
             (0.8037365323032135, mp.Vector3(0.5, 0.0, 0.0))),
            ((0.7017932556977557, mp.Vector3(0.5, 0.5, 0.0)),
             (0.8863448167711359, mp.Vector3(0.5, 0.10000000000000003, 0.0))),
            ((0.9047498485809726, mp.Vector3(0.5, 0.10000000000000003, 0.0)),
             (1.0557468193007016, mp.Vector3(0.5, 0.5, 0.0))),
            ((1.0077925606103986, mp.Vector3(0.2, 0.2, 0.0)),
             (1.1815403744341757, mp.Vector3(0.0, 0.0, 0.0))),
            ((1.122424251973878, mp.Vector3(0.20000000000000004, 0.0, 0.0)),
             (1.2351567679231688, mp.Vector3(0.30000000000000004, 0.30000000000000004, 0.0))),
            ((1.2059728636717586, mp.Vector3(0.0, 0.0, 0.0)),
             (1.3135062523646421, mp.Vector3(0.30000000000000004, 0.0, 0.0))),
        ]

        ms.run_te()
        self.check_band_range_data(expected_brd, ms.band_range_data)

if __name__ == '__main__':
    unittest.main()
