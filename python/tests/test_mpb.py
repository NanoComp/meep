from __future__ import division

import unittest
# import numpy as np
import meep as mp
from meep import mpb


# FIXME: Can only run one test becuase random initialization of fields requires
# setting the seed before each test, and there's currently no way to do that.

# class TestMPBWrappers(unittest.TestCase):

#     def setUp(self):
#         self.num_bands = 8
#         self.k_points = [mp.Vector3(),
#                          mp.Vector3(0.5),
#                          mp.Vector3(0.5, 0.5),
#                          mp.Vector3()]

#         self.k_points = mp.interpolate(4, self.k_points)
#         self.geometry = []  # [mp.Cylinder(0.2, material=mp.Medium(epsilon=12))]
#         self.geometry_lattice = mp.Lattice(size=mp.Vector3(1, 1))
#         self.resolution = 32

#     def test_mode_solver_constructor(self):
#         mpb.mode_solver(self.num_bands, 0, self.resolution, self.geometry_lattice,
#                         1.0e-7, mp.Medium(), self.geometry, True)


class TestModeSolver(unittest.TestCase):

    # def test_list_split(self):
    #     k_points = [
    #         mp.Vector3(),
    #         mp.Vector3(0.5),
    #         mp.Vector3(0.5, 0.5),
    #         mp.Vector3()
    #     ]

    #     k_points = mp.interpolate(4, k_points)

    #     k_split = mp.list_split(k_points, 1, 0)

    #     expected = [
    #         (0, [mp.Vector3(),
    #              mp.Vector3(0.10000000000000003),
    #              mp.Vector3(0.20000000000000004),
    #              mp.Vector3(0.30000000000000004),
    #              mp.Vector3(0.4),
    #              mp.Vector3(0.5),
    #              mp.Vector3(0.5, 0.10000000000000003),
    #              mp.Vector3(0.5, 0.20000000000000004),
    #              mp.Vector3(0.5, 0.30000000000000004),
    #              mp.Vector3(0.5, 0.4),
    #              mp.Vector3(0.5, 0.5),
    #              mp.Vector3(0.4, 0.4),
    #              mp.Vector3(0.30000000000000004, 0.30000000000000004),
    #              mp.Vector3(0.2, 0.2),
    #              mp.Vector3(0.1, 0.1),
    #              mp.Vector3(0.0, 0.0)]),
    #     ]

    #     indx = k_split[0][0]
    #     split_list = k_split[0][1]
    #     self.assertEqual(indx, 0)
    #     for res, exp in zip(split_list, expected[0][1]):
    #         self.assertEqual(res, exp)

    # def test_update_band_range_data(self):
    #     brd = []
    #     freqs = [0.0, 1.0000000001231053, 1.0000000001577114, 1.000000000183077,
    #              1.0000000003647922, 1.4142135627385737, 1.4142135630373556,
    #              1.4142135634172286]
    #     kpoint = mp.Vector3()

    #     expected = [
    #         ((0.0, mp.Vector3()), (0.0, mp.Vector3())),
    #         ((1.0000000001231053, mp.Vector3()), (1.0000000001231053, mp.Vector3())),
    #         ((1.0000000001577114, mp.Vector3()), (1.0000000001577114, mp.Vector3())),
    #         ((1.000000000183077, mp.Vector3()), (1.000000000183077, mp.Vector3())),
    #         ((1.0000000003647922, mp.Vector3()), (1.0000000003647922, mp.Vector3())),
    #         ((1.4142135627385737, mp.Vector3()), (1.4142135627385737, mp.Vector3())),
    #         ((1.4142135630373556, mp.Vector3()), (1.4142135630373556, mp.Vector3())),
    #         ((1.4142135634172286, mp.Vector3()), (1.4142135634172286, mp.Vector3())),
    #     ]

    #     ms = mpb.ModeSolver()
    #     res = ms.update_band_range_data(brd, freqs, kpoint)
    #     self.assertEqual(expected, res)

    def test_no_geometry(self):
        num_bands = 8
        k_points = [
            mp.Vector3(),
            mp.Vector3(0.5),
            mp.Vector3(0.5, 0.5),
            mp.Vector3()
        ]

        k_points = mp.interpolate(4, k_points)
        geometry_lattice = mp.Lattice(size=mp.Vector3(1, 1))
        resolution = 32

        ms = mpb.ModeSolver(
            num_bands=num_bands,
            k_points=k_points,
            geometry=[],
            geometry_lattice=geometry_lattice,
            resolution=resolution
        )

        ms.run_te()

        expected_freqs = [
            (0.0,
             1.0000000001231053,
             1.0000000001577114,
             1.000000000183077,
             1.0000000003647922,
             1.4142135627385737,
             1.4142135630373556,
             1.4142135634172286),
            (0.10000000000000007,
             0.9000000005231413,
             1.004987562417261,
             1.0049875626187106,
             1.1000000004086365,
             1.34536240526912,
             1.3453624069375645,
             1.4866068773856784),
            (0.20000000000000004,
             0.8000000000806242,
             1.0198039027587842,
             1.0198039027805834,
             1.2000000000589406,
             1.2806248475110769,
             1.280624847574805,
             1.562049935491582),
            (0.30000000000000004,
             0.7000000000058019,
             1.0440306508934871,
             1.0440306508948238,
             1.220655561574748,
             1.2206555615791594,
             1.300000000006018,
             1.6401219467453945),
            (0.4,
             0.6000000000045848,
             1.0770329614284637,
             1.0770329614291014,
             1.1661903789696544,
             1.1661903789710815,
             1.4000000000020343,
             1.7204650533727484),
            (0.5,
             0.5500015585576113,
             1.123621141563782,
             1.1322013728463447,
             1.136006887102713,
             1.1389977759520866,
             1.5140363419950278,
             1.6930083259493984),
            (0.5099019513592784,
             0.5099019517772004,
             1.0295630152934836,
             1.029563016461934,
             1.2083045989260803,
             1.2083045997195248,
             1.5033296404803116,
             1.503329646066483),
            (0.5385164807134505,
             0.5385164807531467,
             0.9433981132894825,
             0.9433981133212272,
             1.30000000010291,
             1.3000000001270966,
             1.5132745955004465,
             1.5132745964408887),
            (0.5830951894845302,
             0.583095189493567,
             0.860232526727502,
             0.8602325267368164,
             1.3928388277377648,
             1.392838827741362,
             1.5297058541754436,
             1.5297058543682367),
            (0.6403124237432852,
             0.6403124237442485,
             0.7810249675940715,
             0.7810249675954682,
             1.4866068747336605,
             1.4866068747342678,
             1.5524174696456052,
             1.5524174697169943),
            (0.7071067811865476,
             0.707106781186665,
             0.7071067811871665,
             0.7071067811874233,
             1.581138830084377,
             1.5811388300844045,
             1.581138830084745,
             1.5811388300879181),
            (0.5656854249492379,
             0.7211102550956993,
             0.7211102551151036,
             0.8485281374378021,
             1.456021977862241,
             1.456021977888115,
             1.5231546211779559,
             1.649242250170493),
            (0.4242640687119283,
             0.7615773106248578,
             0.7615773109059598,
             0.9899494938352668,
             1.3341664065094723,
             1.3341664066931676,
             1.4764823061004433,
             1.7262676489509081),
            (0.28284271247461856,
             0.8249001570384027,
             0.8271093536413782,
             1.132601027465129,
             1.2174674299504669,
             1.2192193433933378,
             1.4428744775576858,
             1.8060315065555141),
            (0.14142135623730953,
             0.9055385141690545,
             0.9055385172046134,
             1.1045361022245819,
             1.10453610404424,
             1.2727922109822825,
             1.4212670410466024,
             1.4212671028465198),
            (0.0,
             1.0000000000091356,
             1.0000000000136406,
             1.0000000000238949,
             1.0000000005045615,
             1.4142135623739474,
             1.4142135624052996,
             1.4142135625149197)
        ]

        expected_brd = [
            ((0.0, mp.Vector3()), (0.7071067811865476, mp.Vector3(0.5, 0.5, 0.0))),
            ((0.5099019517772004, mp.Vector3(0.5, 0.10000000000000003, 0.0)),
             (1.0000000001231053, mp.Vector3())),
            ((0.7071067811871665, mp.Vector3(0.5, 0.5, 0.0)),
             (1.123621141563782, mp.Vector3(0.5, 0.0, 0.0))),
            ((0.7071067811874233, mp.Vector3(0.5, 0.5, 0.0)),
             (1.132601027465129, mp.Vector3(0.2, 0.2, 0.0))),
            ((1.0000000003647922, mp.Vector3()),
             (1.581138830084377, mp.Vector3(0.5, 0.5, 0.0))),
            ((1.1389977759520866, mp.Vector3(0.5, 0.0, 0.0)),
             (1.5811388300844045, mp.Vector3(0.5, 0.5, 0.0))),
            ((1.280624847574805, mp.Vector3(0.20000000000000004, 0.0, 0.0)),
             (1.581138830084745, mp.Vector3(0.5, 0.5, 0.0))),
            ((1.4142135625149197, mp.Vector3()),
             (1.8060315065555141, mp.Vector3(0.2, 0.2, 0.0)))
        ]

        for exp, res in zip(expected_brd, ms.band_range_data):
            # Compare min freqs
            self.assertAlmostEqual(exp[0][0], res[0][0])
            # Compare min k
            self.assertTrue(exp[0][1].close(res[0][1]))
            # Compare max freqs
            self.assertAlmostEqual(exp[1][0], res[1][0])
            # Compare max k
            self.assertTrue(exp[1][1].close(res[1][1]))

        for res, exp in zip(ms.all_freqs, expected_freqs):
            for r, e in zip(res, exp):
                self.assertAlmostEqual(r, e)

        gaps = ms.output_gaps(ms.band_range_data)
        self.assertEqual(len(gaps), 0)

if __name__ == '__main__':
    unittest.main()
