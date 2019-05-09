
# medium_evaluations.py - Tests the evaluation of material permitivity profiles.
# Checks materials with lorentizian, drude, and non uniform diagonals.
# The extracted values are compared against actual datapoints pulled from
#   refractiveindex.info.
# TODO: 
#  * check materials with off diagonal components
#  * check magnetic profiles

from __future__ import division

import unittest
import meep as mp
import numpy as np


class TestMediumEvaluations(unittest.TestCase):

    def test_medium_evaluations(self):
        from meep.materials import Si, Au, LiNbO3

        # Check Silicon
        w0 = 1/1.357
        w1 = 1/11.04
        eps = Si.epsilon([w0,w1])
        self.assertAlmostEqual(np.real(np.sqrt(eps[0,0,0])), 3.4975134228247, places=6)
        self.assertAlmostEqual(np.real(np.sqrt(eps[1,0,0])), 3.4175012372164, places=6)

        # Check Gold
        w0 = 1/0.24797
        w1 = 1/6.1992
        eps = Au.epsilon([w0,w1])
        self.assertAlmostEqual(np.real(np.sqrt(eps[0,0,0])), 1.0941, places=4)
        self.assertAlmostEqual(np.real(np.sqrt(eps[1,0,0])), 5.0974, places=4)

        # Check Lithium Niobate
        w0 = 1/0.4
        w1 = 1/5.0
        eps = LiNbO3.epsilon([w0,w1])
        self.assertAlmostEqual(np.real(np.sqrt(eps[0,0,0])), 2.4392710888511, places=6)
        self.assertAlmostEqual(np.real(np.sqrt(eps[1,0,0])), 2.0508048794821, places=6)
        self.assertAlmostEqual(np.real(np.sqrt(eps[0,2,2])), 2.332119801394, places=6)
        self.assertAlmostEqual(np.real(np.sqrt(eps[1,2,2])), 2.0025312699385, places=6)

if __name__ == '__main__':
    unittest.main()
