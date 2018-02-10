from __future__ import division

import glob
import math
import os
import re
import unittest

import h5py
# TODO: Importing numpy loads MKL which breaks zdotc_
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.optimize import ridder
import meep as mp
from meep import mpb


class TestModeSolver(unittest.TestCase):

    data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))

    def setUp(self):
        """Store the test name and register a function to clean up all the
        generated h5 files."""

        self.filename_prefix = self.id().split('.')[-1]

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
            deterministic=True
        )

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

    def check_band_range_data(self, expected_brd, result):
        for exp, res in zip(expected_brd, result):
            # Compare min freqs
            self.assertAlmostEqual(exp[0][0], res[0][0])
            # Compare min k
            self.assertTrue(exp[0][1].close(res[0][1]))
            # Compare max freqs
            self.assertAlmostEqual(exp[1][0], res[1][0])
            # Compare max k
            self.assertTrue(exp[1][1].close(res[1][1]))

    def check_freqs(self, expected_freqs, result):
        for exp, res in zip(expected_freqs, result):
            for r, e in zip(res, exp):
                self.assertAlmostEqual(r, e)

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

            np.testing.assert_allclose(ref_arr, field)

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

        expected_freqs = [
            (0.0, 1.0000000000464613, 1.0000000000543054, 1.0000000000615705,
             1.000000000146055, 1.4142135624581003, 1.4142135625089993, 1.4142135625528056),
            (0.1, 0.9000000001550924, 1.0049875621904198, 1.0049875622885223,
             1.1000000000933774, 1.3453624048874995, 1.3453624052436801, 1.4866068751122161),
            (0.2, 0.8000000000116375, 1.0198039027271282, 1.019803902730517,
             1.200000000010754, 1.2806248474929454, 1.2806248475163993, 1.5620499352326278),
            (0.30000000000000004, 0.7000000000011868, 1.0440306508918364, 1.0440306508920387,
             1.2206555615737378, 1.2206555615763066, 1.3000000000014578, 1.6401219466969021),
            (0.4000000000000003, 0.6000000000004628, 1.0770329614271468, 1.0770329614271774,
             1.1661903789691512, 1.166190378969408, 1.4000000000003532, 1.7204650534095256),
            (0.5, 0.5000000000350678, 1.1180339887555244, 1.1180339887627508,
             1.1180339887649489, 1.1180339887687103, 1.5000000000111153, 1.8027756376524453),
            (0.509901951359278, 0.5099019514458716, 1.0295630141348482, 1.0295630142004515,
             1.2083045973874058, 1.2083045974222109, 1.5033296378612853, 1.7492855555308888),
            (0.5385164807134508, 0.5667749251868665, 0.9569710818749532, 0.9598001998613896,
             1.3060978520614315, 1.3097390466091283, 1.5200437573891414, 1.5961061937551981),
            (0.58309518948453, 0.5830951911457611, 0.8602325273332244, 0.8602325277393493,
             1.3928388281223145, 1.3928388365843103, 1.5297058564565378, 1.5297058995212813),
            (0.6403124237432845, 0.6403124247682638, 0.781024968061956, 0.7810249682383915,
             1.4866068748830337, 1.486606876659867, 1.5524174699814146, 1.5524174856540145),
            (0.7071067811865476, 0.7071067811880555, 0.7071067811884221, 0.7071067811901163,
             1.5811388300846416, 1.5811388300858549, 1.5811388300892486, 1.581138830114666),
            (0.5656854249492381, 0.7211102567829735, 0.7211102586521964, 0.8485281398532061,
             1.4560219792287434, 1.4560219922982203, 1.5231546215078169, 1.5231546281251418),
            (0.4242640687119286, 0.761577311243524, 0.7615773118865434, 0.9899494940070525,
             1.3341664068961834, 1.3341664087834073, 1.47648230615162, 1.4764823071570101),
            (0.282842712474619, 0.8246211251927931, 0.8246211252532635, 1.1313708499266775,
             1.2165525061076403, 1.216552506264088, 1.442220510231749, 1.4422205103997088),
            (0.14142135623730961, 0.9055385138206714, 0.9055385138262214, 1.1045361017236859,
             1.1045361017248032, 1.2727922061382657, 1.4212670403652352, 1.4212670403803151),
            (0.0, 1.0000000000003246, 1.0000000000004943, 1.0000000000005942,
             1.000000000001028, 1.4142135623732919, 1.414213562374368, 1.4142135623752818)
        ]

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
            ((1.2806248475163993, mp.Vector3(0.2000000000000004, 0.0, 0.0)),
             (1.5811388300892486, mp.Vector3(0.5, 0.5, 0.0))),
            ((1.4142135623752818, mp.Vector3(0.0, 0.0, 0.0)),
             (1.8027756376524453, mp.Vector3(0.5, 0.0, 0.0))),

        ]

        ms = self.init_solver(False)
        ms.run_te()

        self.check_band_range_data(expected_brd, ms.band_range_data)
        self.check_freqs(expected_freqs, ms.all_freqs)
        self.assertEqual(len(ms.gap_list), 0)

    def test_run_te(self):

        expected_freqs = [
            (0.0, 0.5543123040724638, 0.7744332061802593, 0.7744499175611301,
             0.9235813808951666, 1.0008467617941965, 1.000853657978345, 1.0937421750914536),
            (0.08975526912402747, 0.5527259748791986, 0.7625580252843803, 0.775927071773125,
             0.9083025904214346, 1.0025041485479258, 1.0040131495259865, 1.1151425798690515),
            (0.17874194256279652, 0.5464061721502865, 0.7290550665093838, 0.7798577251387895,
             0.8820233099342442, 1.0074514510489234, 1.0172937305690262, 1.1371240270581144),
            (0.265889317531334, 0.529456650246997, 0.6858506574662463, 0.7848207792703579,
             0.8628144586715519, 1.015371497818431, 1.0388248806041873, 1.0978603871063697),
            (0.3491098642318947, 0.491614128634062, 0.6534277127248703, 0.7889250441298363,
             0.8521868291383247, 1.0243890854263964, 1.0648782507133647, 1.0918380516599975),
            (0.41281930737977607, 0.4424096429060285, 0.6427858365432216, 0.7905166037050794,
             0.8487837791081655, 1.028981766498537, 1.0865991750573356, 1.088248151258527),
            (0.4237693545494147, 0.4466340166909042, 0.6389975571819423, 0.7946968579934286,
             0.8449963934371933, 0.9865104944996825, 1.050849906850091, 1.1042021258981072),
            (0.45480717337736154, 0.458384478782914, 0.6285854682590182, 0.8067270292869435,
             0.8274711584785993, 0.928705439321874, 1.005398394501796, 1.10918073349481),
            (0.4748865013922755, 0.50131006667719, 0.6141096017945592, 0.7818160269356271,
             0.8255441878500366, 0.8967928418027945, 0.9600200633481636, 1.111104880065585),
            (0.4909986510517719, 0.5561139188988492, 0.5998745448965, 0.7195409655301487,
             0.850456196583463, 0.8866021931139869, 0.9182966027540888, 1.1120385722459878),
            (0.49832784140942327, 0.5932903636375229, 0.5933300432825276, 0.6789328334547384,
             0.8786729653444213, 0.8843176495094911, 0.8843467579864347, 1.112337628001184),
            (0.4718692662905701, 0.5466449267509191, 0.6069228332255978, 0.7435877802945671,
             0.8404058252567083, 0.8845903173166171, 0.9379844985289354, 1.1116214392114732),
            (0.3726512334709978, 0.5374910186290283, 0.6448252155990313, 0.8107898648188872,
             0.8261977806810092, 0.8858710104130919, 1.010853588402447, 1.1076219716334672),
            (0.2521701793341113, 0.5443731521852448, 0.6985070869163165, 0.7905695248250529,
             0.8904434436928694, 0.8983876721337736, 1.0797388604006626, 1.090638038627879),
            (0.12687260757609428, 0.551484835775791, 0.7512227072577523, 0.778461951871455,
             0.9043827134262497, 0.9626336794092828, 1.0422273104938644, 1.101944971644869),
            (0.0, 0.5543122885042596, 0.7744332085266307, 0.7744499167032798,
             0.9235813819246859, 1.0008467515107635, 1.0008536928499856, 1.0937421550082655),
        ]

        # Obtained by passing NULL to set_maxwell_dielectric for epsilon_mean_func
        # in mpb/mpb/medium.c.
        expected_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.49832784140942327, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.4424096429060285, mp.Vector3(0.5, 0.0, 0.0)),
             (0.5932903636375229, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.5933300432825276, mp.Vector3(0.5, 0.5, 0.0)),
             (0.7744332085266307, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.6789328334547384, mp.Vector3(0.5, 0.5, 0.0)),
             (0.8107898648188872, mp.Vector3(0.30000000000000004, 0.30000000000000004, 0.0))),
            ((0.8255441878500366, mp.Vector3(0.5, 0.30000000000000004, 0.0)),
             (0.9235813819246859, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.8843176495094911, mp.Vector3(0.5, 0.5, 0.0)),
             (1.028981766498537, mp.Vector3(0.5, 0.0, 0.0))),
            ((0.8843467579864347, mp.Vector3(0.5, 0.5, 0.0)),
             (1.0865991750573356, mp.Vector3(0.5, 0.0, 0.0))),
            ((1.088248151258527, mp.Vector3(0.5, 0.0, 0.0)),
             (1.1371240270581144, mp.Vector3(0.20000000000000004, 0.0, 0.0))),
        ]

        expected_gap_list = [
            (0.006687841330442318, 0.5932903636375229, 0.5933300432825276),
            (1.8033387506768026, 0.8107898648188872, 0.8255441878500366),
            (0.15164063989582524, 1.0865991750573356, 1.088248151258527),
        ]

        ms = self.init_solver()
        ms.run_te()

        self.check_band_range_data(expected_brd, ms.band_range_data)
        self.check_freqs(expected_freqs, ms.all_freqs)
        self.check_gap_list(expected_gap_list, ms.gap_list)

    def test_run_tm(self):
        expected_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.2812417875105901, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.41769014812162525, mp.Vector3(0.5, 0.0, 0.0)),
             (0.5460171106569119, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.49700744460873947, mp.Vector3(0.5, 0.5, 0.0)),
             (0.5585358608818134, mp.Vector3(0.5, 0.0, 0.0))),
            ((0.5529886297870371, mp.Vector3(0.0, 0.0, 0.0)),
             (0.71347812930332, mp.Vector3(0.5, 0.0, 0.0))),
            ((0.741761112489577, mp.Vector3(0.5, 0.0, 0.0)),
             (0.8593597943953297, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.8295123038043076, mp.Vector3(0.30000000000000004, 0.30000000000000004, 0.0)),
             (0.878308388274097, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.8627300031476208, mp.Vector3(0.5, 0.0, 0.0)),
             (0.9522930815607895, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.8809282925679489, mp.Vector3(0.5, 0.5, 0.0)),
             (1.0840705439027458, mp.Vector3(0.0, 0.0, 0.0)))
        ]

        ms = self.init_solver()
        ms.run_tm()

        self.check_band_range_data(expected_brd, ms.band_range_data)

    def _test_get_field(self, field):
        ms = self.init_solver()
        ms.run_te()
        get_field_func = getattr(ms, "get_{}field".format(field))
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
        fname = 'mpb-tutorial-epsilon.h5'
        data_path = os.path.join(self.data_dir, fname)

        ms = self.init_solver(False)
        ms.run_te()

        with h5py.File('test_output_field_to_file-epsilon.h5', 'r') as f:
            with h5py.File(data_path) as ref:
                for k in ref.keys():
                    if k == 'description':
                        self.assertEqual(ref[k].value, f[k].value)
                    else:
                        np.testing.assert_array_equal(ref[k].value, f[k].value)

    def test_compute_field_energy(self):
        ms = self.init_solver()
        ms.run_te()
        ms.get_dfield(8)
        res = ms.compute_field_energy()

        expected = [1.0000000000000004, 1.725992366452563e-5, 0.49998273982991737,
                    1.725992370353374e-5, 0.49998274032271406, 0.0, 0.0]

        np.testing.assert_allclose(expected, res)

    def test_output_efield_z(self):
        ms = self.init_solver()
        ms.run_tm()
        mpb.output_efield_z(ms, 8)

        ref_fname = 'tutorial-e.k16.b08.z.tm.h5'
        ref_path = os.path.join(self.data_dir, ref_fname)

        with h5py.File(ref_path, 'r') as ref:
            with h5py.File(re.sub('tutorial', ms.filename_prefix, ref_fname)) as res:
                for k in ref.keys():
                    if k == 'description':
                        self.assertEqual(ref[k].value, res[k].value)
                    elif k == 'z.i' or k == 'z.r':
                        for e, r in zip(ref[k].value.ravel(), res[k].value.ravel()):
                            self.assertAlmostEqual(e, r, places=3)
                    else:
                        np.testing.assert_allclose(ref[k].value, res[k].value)


    def test_triangular_lattice(self):

        expected_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)),
             (0.27418366655053683, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0))),
            ((0.4443772345549731, mp.Vector3(0.0, 0.5, 0.0)),
             (0.5589242526869668, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.4890716313301616, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0)),
             (0.5590346475838462, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.5933240874910938, mp.Vector3(0.0, 0.0, 0.0)),
             (0.7895347583021823, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0))),
            ((0.7896048828436291, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0)),
             (0.8347553850744481, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.8348898399363115, mp.Vector3(0.0, 0.0, 0.0)),
             (0.8651618823872929, mp.Vector3(-0.2, 0.39999999999999997, 0.0))),
            ((0.8669299160898801, mp.Vector3(-0.13333333333333336, 0.4333333333333333, 0.0)),
             (0.9918057091113671, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.8973819286907876, mp.Vector3(-0.3333333333333333, 0.3333333333333333, 0.0)),
             (1.0937048661670459, mp.Vector3(0.0, 0.0, 0.0)))
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
        expected = 39.10065045011942
        self.assertAlmostEqual(expected, result.fun * -1)

    def test_anisotropic_2d_gap(self):

        expected_brd = [
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)), (0.22131658647819047, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.0, mp.Vector3(0.0, 0.0, 0.0)), (0.23068427740922454, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.2169119289078649, mp.Vector3(0.5, 0.0, 0.0)), (0.3190202969087597, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.30110868688276016, mp.Vector3(0.5, 0.0, 0.0)), (0.36483532251792833, mp.Vector3(0.0, 0.0, 0.0))),
            ((0.3070162288067315, mp.Vector3(0.5, 0.5, 0.0)), (0.3852513891468177, mp.Vector3(0.1, 0.1, 0.0))),
            ((0.341835282475095, mp.Vector3(0.5, 0.30000000000000004, 0.0)), (0.3914210645066325, mp.Vector3(0.5, 0.10000000000000003, 0.0))),
            ((0.3498214344265213, mp.Vector3(0.5, 0.5, 0.0)), (0.40756434315471013, mp.Vector3(0.4, 0.0, 0.0))),
            ((0.3963466078874464, mp.Vector3(0.1, 0.1, 0.0)), (0.477223780015547, mp.Vector3(0.5, 0.5, 0.0))),
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

        ms.geometry = mp.geometric_object_lattice_duplicates(ms.geometry_lattice, ms.geometry)
        ms.geometry.append(mp.Cylinder(0.2, material=mp.air))

        ms.resolution = 16
        ms.k_points = [mp.Vector3(0.5, 0.5)]

        ms.num_bands = 50
        ms.run_tm()

        mpb.output_efield_z(ms, 25)

        ms.get_dfield(25)
        ms.compute_field_energy()
        c = mp.Cylinder(1.0, material=mp.air)
        e = ms.compute_energy_in_objects([c])
        self.assertAlmostEqual(0.6202776947564257, e)

        ms.num_bands = 1
        ms.target_freq = (0.2812 + 0.4174) / 2
        ms.tolerance = 1e-8
        ms.run_tm()

        expected_brd = [
            ((0.37677410124156546, mp.Vector3(0.5, 0.5, 0.0)),
             (0.37677410124156546, mp.Vector3(0.5, 0.5, 0.0)))
        ]

        self.check_band_range_data(expected_brd, ms.band_range_data)

        old_geometry = ms.geometry  # save the 5x5 grid with a missing rod

        def rootfun(eps):
            ms.geometry = old_geometry + [mp.Cylinder(0.2, material=mp.Medium(epsilon=eps))]
            ms.run_tm()
            return ms.get_freqs()[0] - 0.314159

        rooteps = ridder(rootfun, 1, 12)
        rootval = rootfun(rooteps)

        self.assertAlmostEqual(5.221288723644521, rooteps, places=5)
        self.assertAlmostEqual(1.1874461136596182e-8, rootval)

if __name__ == '__main__':
    unittest.main()
