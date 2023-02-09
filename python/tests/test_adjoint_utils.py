"""test_adjoint_utils.py

Check various components of the adjoint solver codebase, like
filters, which may not need explicit gradient computation
(i.e. forward and adjoint runs).

"""
import meep as mp

try:
    import meep.adjoint as mpa
except:
    import adjoint as mpa

import unittest
from enum import Enum

import numpy as np
import parameterized
from autograd import numpy as npa
from autograd import tensor_jacobian_product
from utils import ApproxComparisonTestCase

_TOL = 1e-6 if mp.is_single_precision() else 1e-14

## ensure reproducible results
rng = np.random.RandomState(9861548)


class TestAdjointUtils(ApproxComparisonTestCase):
    @parameterized.parameterized.expand(
        [
            ("1.0_1.0_20_conic", 1.0, 1.0, 20, 0.24, mpa.conic_filter),
            ("1.0_1.0_23_conic", 1.0, 1.0, 23, 0.24, mpa.conic_filter),
            ("0.887_1.56_conic", 0.887, 1.56, 20, 0.24, mpa.conic_filter),
            ("0.887_1.56_conic", 0.887, 1.56, 31, 0.24, mpa.conic_filter),
            ("0.887_1.56_gaussian", 0.887, 1.56, 20, 0.24, mpa.gaussian_filter),
            ("0.887_1.56_cylindrical", 0.887, 1.56, 20, 0.24, mpa.cylindrical_filter),
        ]
    )
    def test_filter_offset(self, test_name, Lx, Ly, resolution, radius, filter_func):
        """ensure that the filters are indeed zero-phase"""
        print("Testing ", test_name)
        Nx, Ny = int(round(resolution * Lx)), int(round(resolution * Ly))
        x = np.random.rand(Nx, Ny)
        x = x + np.fliplr(x)
        x = x + np.flipud(x)
        y = filter_func(x, radius, Lx, Ly, resolution)
        self.assertClose(y, np.fliplr(y), epsilon=_TOL)
        self.assertClose(y, np.flipud(y), epsilon=_TOL)


if __name__ == "__main__":
    unittest.main()
