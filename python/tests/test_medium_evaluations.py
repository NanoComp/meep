# medium_evaluations.py - Tests the evaluation of material permitivity profiles.
# Checks materials with lorentizian, drude, and non uniform diagonals.
# The extracted values are compared against actual datapoints pulled from
#   refractiveindex.info.
# TODO:
#  * check materials with off diagonal components
#  * check magnetic profiles
import unittest

import numpy as np

import meep as mp


class TestMediumEvaluations(unittest.TestCase):
    def test_medium_evaluations(self):
        from meep.materials import Ag, LiNbO3, Si, fused_quartz

        # Check that scalars work
        w0 = LiNbO3.valid_freq_range.min
        eps = LiNbO3.epsilon(w0)
        self.assertAlmostEqual(np.real(np.sqrt(eps[0, 0])), 2.0508, places=4)

        # Check numpy arrays
        try:
            w0 = Si.valid_freq_range.min
            w1 = Si.valid_freq_range.max
            eps = Si.epsilon(np.linspace(w0, w1, 100))
        except ExceptionType:
            self.fail("myFunc() raised ExceptionType unexpectedly!")

        # Check that regions outside of domain don't work
        self.assertRaises(ValueError, LiNbO3.epsilon, -1.0)
        self.assertRaises(ValueError, LiNbO3.epsilon, 10000.0)

        # Check complex vs non complex numbers
        self.assertTrue(np.iscomplex(Ag.epsilon(1.0)[0, 0]))
        self.assertFalse(np.iscomplex(fused_quartz.epsilon(1.0)[0, 0]))

        # Check Silicon
        w0 = Si.valid_freq_range.min
        w1 = Si.valid_freq_range.max
        eps = Si.epsilon([w0, w1])
        self.assertAlmostEqual(np.real(np.sqrt(eps[0, 0, 0])), 3.4175, places=4)
        self.assertAlmostEqual(np.real(np.sqrt(eps[1, 0, 0])), 3.4971, places=4)

        # Check Silver
        w0 = Ag.valid_freq_range.min
        w1 = Ag.valid_freq_range.max
        eps = Ag.epsilon([w0, w1])
        self.assertAlmostEqual(np.real(np.sqrt(eps[0, 0, 0])), 17.485, places=2)
        self.assertAlmostEqual(np.real(np.sqrt(eps[1, 0, 0])), 0.44265, places=4)

        # Check Lithium Niobate
        w0 = LiNbO3.valid_freq_range.min
        w1 = LiNbO3.valid_freq_range.max
        eps = LiNbO3.epsilon([w0, w1])
        self.assertAlmostEqual(np.real(np.sqrt(eps[0, 0, 0])), 2.0508, places=4)
        self.assertAlmostEqual(np.real(np.sqrt(eps[1, 0, 0])), 2.4393, places=4)
        self.assertAlmostEqual(np.real(np.sqrt(eps[0, 2, 2])), 2.0025, places=4)
        self.assertAlmostEqual(np.real(np.sqrt(eps[1, 2, 2])), 2.3321, places=4)


if __name__ == "__main__":
    unittest.main()
