from typing import Union
import unittest

import numpy as np


class ApproxComparisonTestCase(unittest.TestCase):
    """A mixin for adding correct scalar/vector comparison."""

    def assertClose(
        self,
        x: Union[float, np.ndarray],
        y: Union[float, np.ndarray],
        epsilon: float = 1e-2,
        msg: str = "",
    ):
        """Checks if two scalars or vectors satisfy ‖x-y‖ ≤ ε * max(‖x‖, ‖y‖).

        Args:
            x, y: two quantities to be compared (scalars or 1d arrays).
            epsilon: threshold value (maximum) of the relative error.
            msg: a string to display if the inequality is violated.
        """
        x = np.atleast_1d(x).ravel()
        y = np.atleast_1d(y).ravel()
        x_norm = np.linalg.norm(x, ord=np.inf)
        y_norm = np.linalg.norm(y, ord=np.inf)
        diff_norm = np.linalg.norm(x - y, ord=np.inf)
        self.assertLessEqual(diff_norm, epsilon * np.maximum(x_norm, y_norm), msg)
