from __future__ import division, print_function

import glob
import math
import os
import re
import sys
import time
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

        self.start = time.time()

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

    def tearDown(self):
        end = time.time() - self.start
        print("{}: {:.2f}s".format(self.filename_prefix, end))

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

    def check_fields_against_h5(self, ref_path, field, suffix=''):
        with h5py.File(ref_path, 'r') as ref:
            # Reshape the reference data into a component-wise 1d array like
            # [x1,y1,z1,x2,y2,z2,etc.]
            ref_x = ref["x.r{}".format(suffix)].value + ref["x.i{}".format(suffix)].value * 1j
            ref_y = ref["y.r{}".format(suffix)].value + ref["y.i{}".format(suffix)].value * 1j
            ref_z = ref["z.r{}".format(suffix)].value + ref["z.i{}".format(suffix)].value * 1j

            ref_arr = np.zeros(np.prod(field.shape), dtype=np.complex128)
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

        pt = ms.get_epsilon_point(mp.Vector3(0.5, 0.5))
        self.assertEqual(pt, 1.0)

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
        field_pt = ms.get_field_point(mp.Vector3(0.5, 0.5))
        bloch_field_pt = ms.get_bloch_field_point(mp.Vector3(0.5, 0.5))
        eps_inv_tensor = ms.get_epsilon_inverse_tensor_point(mp.Vector3(0.5, 0.5))

        energy = ms.compute_field_energy()
        pt = ms.get_energy_point(mp.Vector3(0.5, 0.5))

        self.assertAlmostEqual(pt, 1.330368347216153e-9)

        expected_fp = mp.Vector3(
            2.5823356723958247e-5 + 6.713243287584132e-12j,
            -2.575955745071957e-5 - 6.696552990958943e-12j,
            0.0 - 0.0j
        )

        expected_bloch_fp = mp.Vector3(
            2.5823356723958247e-5 + 6.713243287584132e-12j,
            -2.575955745071957e-5 - 6.696552990958943e-12j,
            -0.0 - 0.0j
        )

        expected_eps_inv_tensor = mp.Matrix(
            mp.Vector3(1.0 + 0.0j, 0.0 - 0.0j, 0.0 - 0.0j),
            mp.Vector3(0.0 + 0.0j, 1.0 + 0.0j, 0.0 - 0.0j),
            mp.Vector3(0.0 + 0.0j, 0.0 + 0.0j, 1.0 + 0.0j)
        )

        self.assertTrue(expected_fp.close(field_pt))
        self.assertTrue(expected_bloch_fp.close(bloch_field_pt))
        self.assertEqual(expected_eps_inv_tensor.c1, eps_inv_tensor.c1)
        self.assertEqual(expected_eps_inv_tensor.c2, eps_inv_tensor.c2)
        self.assertEqual(expected_eps_inv_tensor.c3, eps_inv_tensor.c3)

        energy_in_dielectric = ms.compute_energy_in_dielectric(0, 1)

        expected_energy = [1.0000000000000002, 1.726755206037815e-5, 0.4999827324479414,
                           1.7267552060375955e-5, 0.4999827324479377, 0.0, 0.0]

        expected_energy_in_dielectric = 0.6990769686037558

        self.compare_arrays(np.array(expected_energy), np.array(energy))
        self.assertAlmostEqual(expected_energy_in_dielectric, energy_in_dielectric, places=3)

    def test_compute_group_velocity(self):
        ms = self.init_solver()
        ms.run_te()
        res1 = ms.compute_group_velocity_component(mp.Vector3(0.5, 0.5))
        res2 = ms.compute_one_group_velocity(8)
        res3 = ms.compute_one_group_velocity_component(mp.Vector3(0.5, 0.5), 8)

        expected1 = [
            0.0, 1.470762578355642e-4, -1.4378185933055663e-4, 1.1897503996483383e-4,
            -4.892687048681629e-4, 1.1240346140784176e-4, 1.5842474585356007e-4,
            4.496945573323881e-5,
        ]

        expected2 = mp.Vector3(3.180010979062989e-5, 3.179611968757397e-5)
        expected3 = 4.496932512216992e-5

        for e, r in zip(expected1, res1):
            self.assertAlmostEqual(e, r, places=4)

        self.assertTrue(expected2.close(res2, tol=3))
        self.assertAlmostEqual(expected3, res3, places=3)

    def test_output_efield_z(self):
        ms = self.init_solver()
        ms.run_tm()
        mpb.fix_efield_phase(ms, 8)
        mpb.output_efield_z(ms, 8)

        ref_fname = 'tutorial-e.k16.b08.z.tm.h5'
        ref_path = os.path.join(self.data_dir, ref_fname)
        res_path = re.sub('tutorial', ms.filename_prefix, ref_fname)
        self.compare_h5_files(ref_path, res_path)

    def test_output_dpwr_in_objects(self):
        ms = self.init_solver()
        ms.run_te(mpb.output_dpwr_in_objects(mpb.output_dfield, 0.85, ms.geometry))

        ref_fname1 = 'tutorial-d.k01.b02.te.h5'
        ref_fname2 = 'tutorial-d.k16.b02.te.h5'

        ref_path1 = os.path.join(self.data_dir, ref_fname1)
        ref_path2 = os.path.join(self.data_dir, ref_fname2)

        res_path1 = re.sub('tutorial', ms.filename_prefix, ref_fname1)
        res_path2 = re.sub('tutorial', ms.filename_prefix, ref_fname2)

        self.compare_h5_files(ref_path1, res_path1)
        self.compare_h5_files(ref_path2, res_path2)

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

        # Test MPBData
        with h5py.File(ref_path, 'r') as f:
            efield_re = f['z.r'].value
            efield_im = f['z.i'].value
            efield = np.vectorize(complex)(efield_re, efield_im)

        # rectangularize
        md = mpb.MPBData(ms.get_lattice(), rectify=True, resolution=32, periods=3, verbose=True)
        new_efield = md.convert(efield, ms.k_points[10])
        # check with ref file

        ref_fn = 'tri-rods-e.k11.b08.z.tm-r-m3-n32.h5'
        ref_path = os.path.join(self.data_dir, ref_fn)

        with h5py.File(ref_path, 'r') as f:
            expected_re = f['z.r-new'].value
            expected_im = f['z.i-new'].value
            expected = np.vectorize(complex)(expected_re, expected_im)
            self.compare_arrays(expected, new_efield)

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

        # Test MPBData
        eps = ms.get_epsilon()
        md = mpb.MPBData(ms.get_lattice(), rectify=True, resolution=32, periods=3, verbose=True)
        new_eps = md.convert(eps)
        ref_fn = 'tri-rods-epsilon-r-m3-n32.h5'
        ref_path = os.path.join(self.data_dir, ref_fn)

        with h5py.File(ref_path, 'r') as f:
            ref = f['data-new'].value
            self.compare_arrays(ref, new_eps, tol=1e-3)

    def test_subpixel_averaging(self):
        ms = self.init_solver()
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

    def test_output_tot_pwr(self):
        ms = self.init_solver()
        ms.run_te()
        mpb.output_tot_pwr(ms, 8)

        ref_fname = 'tutorial-tot.rpwr.k16.b08.te.h5'
        ref_path = os.path.join(self.data_dir, ref_fname)
        res_path = re.sub('tutorial', self.filename_prefix, ref_fname)

        self.compare_h5_files(ref_path, res_path)

    def test_get_eigenvectors(self):
        ms = self.init_solver()
        ms.run_te(mpb.fix_hfield_phase)

        def compare_eigenvectors(ref_fn, start, cols):
            with h5py.File(os.path.join(self.data_dir, ref_fn), 'r') as f:
                expected = f['rawdata'].value
                # Reshape the last dimension of 2 reals into one complex
                expected = np.vectorize(complex)(expected[..., 0], expected[..., 1])
                ev = ms.get_eigenvectors(start, cols)
                np.testing.assert_allclose(expected, ev, rtol=1e-3)

        # Get all columns
        compare_eigenvectors('tutorial-te-eigenvectors.h5', 1, 8)
        # Get last column
        compare_eigenvectors('tutorial-te-eigenvectors-8-1.h5', 8, 1)
        # Get columns 3,4, and 5
        compare_eigenvectors('tutorial-te-eigenvectors-3-3.h5', 3, 3)

    def test_set_eigenvectors(self):
        ms = self.init_solver()

        def set_H_to_zero_and_check(start, num_bands):
            ev = ms.get_eigenvectors(start, num_bands)
            self.assertNotEqual(np.count_nonzero(ev), 0)
            zeros = np.zeros(ev.shape, dtype=np.complex128)
            ms.set_eigenvectors(zeros, start)
            new_ev = ms.get_eigenvectors(start, num_bands)
            self.assertEqual(np.count_nonzero(new_ev), 0)

        ms.run_te()
        set_H_to_zero_and_check(8, 1)
        ms.run_te()
        set_H_to_zero_and_check(1, 8)
        ms.run_te()
        set_H_to_zero_and_check(3, 3)

    def test_load_and_save_eigenvectors(self):
        ms = self.init_solver()
        ms.run_te()
        fn = self.filename_prefix + '.h5'

        ev = ms.get_eigenvectors(8, 1)
        zeros = np.zeros(ev.shape, dtype=np.complex128)
        ms.set_eigenvectors(zeros, 8)
        ms.save_eigenvectors(fn)

        ms.run_te()
        ms.load_eigenvectors(fn)
        new_ev = ms.get_eigenvectors(8, 1)
        self.assertEqual(np.count_nonzero(new_ev), 0)

    def test_handle_cvector(self):
        from mpb_tri_rods import ms
        ms.deterministic = True
        ms.tolerance = 1e-12
        ms.filename_prefix = self.filename_prefix
        efields = []

        def get_efields(ms, band):
            efields.append(ms.get_efield(8, output=True))

        k = mp.Vector3(1 / -3, 1 / 3)
        ms.run_tm(mpb.output_at_kpoint(k, mpb.fix_efield_phase, get_efields))

        lat = ms.get_lattice()
        md = mpb.MPBData(lat, rectify=True, periods=3, resolution=32, verbose=True)
        result = md.convert(efields[-1], ms.k_points[10])

        ref_fn = 'converted-tri-rods-e.k11.b08.tm.h5'
        ref_path = os.path.join(self.data_dir, ref_fn)

        self.check_fields_against_h5(ref_path, result.ravel(), suffix='-new')

    def test_epsilon_input_file(self):
        ms = self.init_solver(geom=False)
        eps_fn = 'eps_input_file_test.h5'
        ms.epsilon_input_file = os.path.join(self.data_dir, eps_fn)

        ms.run_te()

        expected_freqs = np.array([
            0.0, 0.5543986136451342, 0.7613327775255415, 0.7613339178956054,
            0.8940893915924257, 0.998342969572652, 0.9983441882455961, 1.0747466061007138
        ])

        expected_gap_list = [
            (3.848610367089048e-5, 0.5781812856814899, 0.5781815082009817),
            (1.4651880980150234, 0.8051999839699242, 0.8170847453549156),
            (0.75255857475812, 1.0580309832489785, 1.0660233597945266),
        ]

        expected_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.4970977843772992, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.4402896410505961, mp.Vector3(0.5, 0.0, 0.0)),
             (0.5781812856814899, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.5781815082009817, mp.Vector3(0.5, 0.5, 0.0)),
             (0.761332777525562, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.6689126424359774, mp.Vector3(0.5, 0.5, 0.0)),
             (0.8051999839699242, mp.Vector3(0.30000000000000004, 0.30000000000000004, 0.0))),
            ((0.8170847453549156, mp.Vector3(0.30000000000000004, 0.30000000000000004, 0.0)),
             (0.8940893915924548, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.8826671164993868, mp.Vector3(0.30000000000000004, 0.30000000000000004, 0.0)),
             (1.0014926328155058, mp.Vector3(0.5, 0.0, 0.0))),
            ((0.8832199143682116, mp.Vector3(0.5, 0.5, 0.0)),
             (1.0580309832489785, mp.Vector3(0.5, 0.0, 0.0))),
            ((1.0660233597945266, mp.Vector3(0.2, 0.2, 0.0)),
             (1.087345829555555, mp.Vector3(0.5, 0.5, 0.0))),
        ]

        self.check_band_range_data(expected_brd, ms.band_range_data)
        self.compare_arrays(expected_freqs, ms.all_freqs[-1])

        self.check_gap_list(expected_gap_list, ms.gap_list)

        pt = ms.get_epsilon_point(mp.Vector3(0.5, 0.5))
        self.assertEqual(pt, 1.0)

    def test_hermitian_eps(self):
        ms = self.init_solver()
        ms.num_bands = 10
        ms.geometry_lattice = mp.Lattice(size=mp.Vector3(1, 1))
        ms.k_points = mp.interpolate(2, [mp.Vector3(), mp.Vector3(0.5),
                                         mp.Vector3(0.5, 0.5), mp.Vector3()])
        if mpb.with_hermitian_epsilon():
            mu_offdiag = mp.Vector3(0 + 12.4j, 0j, 0j)
        else:
            mu_offdiag = mp.Vector3(1.1 + 0j, 0j, 0j)

        mat = mp.Medium(epsilon=15, mu_diag=mp.Vector3(14, 14, 1), mu_offdiag=mu_offdiag)
        ms.geometry = [mp.Cylinder(0.11, material=mat)]

        ms.run_tm()

        expected_freqs_with_herm_eps = [
            (0.0, 0.4632939699892961, 0.5786056046494645, 0.6510415824942094, 0.895078332795855,
             0.9065629078750282, 0.9718615669186841, 1.0031527446201098, 1.0142458802909737, 1.0248230445372033),
            (0.12937769848279101, 0.4619527556321284, 0.5815596420443466, 0.6509684019890601,
             0.8515085512592131, 0.8991327667388375, 0.9483381392427291, 0.9868373331614826,
             1.0201294704318276, 1.1113539456722876),
            (0.240618523732184, 0.45558173972666255, 0.5943681967388351, 0.650900939169996,
             0.7483572681471874, 0.9008690943027022, 0.9516677787271349, 0.9850364677965752,
             1.0516821906747753, 1.1214943845916658),
            (0.2922690373273036, 0.4479821144386504, 0.613801526149148, 0.6509434790185421,
             0.6876982556529582, 0.9015193133195613, 0.9570671227792547, 0.9847711283698155,
             1.0828576982954996, 1.1106109238677777),
            (0.30027099141005154, 0.46956697245972573, 0.5973756776701357, 0.6743055424064167,
             0.6920274562199666, 0.8979504826615435, 0.9268408824618528, 0.9691495966240021,
             0.9983260085044127, 1.1065083471050117),
            (0.31666208802132145, 0.5137036663942733, 0.5882926042080675, 0.6899229118026092,
             0.730249913793744, 0.8458381750994521, 0.8595877992328249, 0.9137388415537298,
             0.9866008233455089, 1.137383764975547),
            (0.3251247277636798, 0.5292011796591106, 0.6018330031246529, 0.7028040151334913,
             0.7794097325510528, 0.7819161956650196, 0.8016335408886606, 0.910192351683647,
             0.98598162196522, 1.1535093885340242),
            (0.29904860785910004, 0.49821749875617755, 0.5818628214691952, 0.6702162529015839,
             0.792305404698029, 0.8010082951265327, 0.9011331789530838, 0.9347832593312477,
             0.9878915728570912, 1.1274287845362319),
            (0.17916974535365005, 0.46423758343427207, 0.5818626159151191, 0.6522567432851836,
             0.8578543711269236, 0.8654932847250656, 0.914866301019591, 0.9836433091978996,
             1.0332068416637614, 1.1280615056475125),
            (0.0, 0.46329396998948535, 0.5786056046501502, 0.6510415824943165, 0.8950783327952294,
             0.9065629078748172, 0.9718615669185836, 1.0031527446209287, 1.01424588029229, 1.0248230445379556)
        ]

        expected_freqs = [
            (0.0, 0.2652450752888948, 0.36211931924886037, 0.3621209970866579,
             0.5052984215604861, 0.5076357969986622, 0.5339446396841523, 0.6471305539348171,
             0.64718320186009, 0.6571600477142543),
            (0.12326602984561506, 0.27771339530892053, 0.3621959012731131, 0.36302565951244253,
             0.505387092405693, 0.5076195849020112, 0.5322316893298896, 0.6470796178474362,
             0.647223221628587, 0.654828206364591),
            (0.1943638388282225, 0.33137494705193943, 0.3623389900641213, 0.37550977433277527,
             0.5054646014393102, 0.5075913386897623, 0.5205986685505688, 0.6412030761801905,
             0.6471864573522761, 0.6473418100828358),
            (0.20839485247782383, 0.3504460140954443, 0.3624169907450078, 0.44431780234615237,
             0.4776254910634774, 0.5057427518243194, 0.5075800427433821, 0.6158600411915314,
             0.647002575504946, 0.6474104093525016),
            (0.2102544365791045, 0.35331213233768893, 0.3612597560359688, 0.45411907757436154,
             0.5002241723073143, 0.5058251418661104, 0.5083203737653703, 0.6162003560297682,
             0.647044345633276, 0.6473905401985292),
            (0.21363558887664455, 0.3572681960432815, 0.3594435751607656, 0.47402103865478207,
             0.5059401251635635, 0.5072182933392233, 0.5640174836064029, 0.6158978436810754,
             0.647079534025694, 0.6473213496606646),
            (0.2151912435586334, 0.358668029449057, 0.3586787257531444, 0.4841939607451779,
             0.5060858864792117, 0.5072753097601064, 0.615065597358213, 0.6150871735511474,
             0.6470586621215278, 0.6472657348303715),
            (0.2104321015504009, 0.3552493202332755, 0.359768447238033, 0.43918040591847685,
             0.5057723660215659, 0.5073426514067818, 0.541478117963846, 0.6308592295716045,
             0.6471525585808731, 0.6473326177488542),
            (0.16216053877476436, 0.2952804946342481, 0.36144739185221775, 0.36540213680816513,
             0.5054838640368768, 0.5075880729152275, 0.5318237561722323, 0.6470214642130181,
             0.6471660484741889, 0.6514980509005492),
            (0.0, 0.2652450753607585, 0.36211931929038943, 0.3621209970933658,
             0.5052984215620827, 0.5076357970004499, 0.5339446396882321, 0.647130554181123,
             0.6471832006793388, 0.6571643599704929)
        ]

        if mpb.with_hermitian_epsilon():
            self.check_freqs(expected_freqs_with_herm_eps, ms.all_freqs)
        else:
            self.check_freqs(expected_freqs, ms.all_freqs)


if __name__ == '__main__':
    unittest.main()
