import meep as mp
try:
    import meep.adjoint as mpa
except:
    import adjoint as mpa
import numpy as np
from autograd import numpy as npa
from autograd import tensor_jacobian_product
import unittest
from enum import Enum
from utils import ApproxComparisonTestCase

## ensure reproducible results
rng = np.random.RandomState(9861548)

class TestAdjointUtils(ApproxComparisonTestCase):
    def test_filter_offset(self):
        '''ensure that the filters are indeed zero-phase'''



if __name__ == '__main__':
    unittest.main()