from autograd import numpy as npa
from autograd import value_and_grad
from typing import Callable, List

from scipy.optimize import minimize


def unfilter_design(target: List[float], processing: Callable, maxiter: int = 100):
    """Given a processing function, uses optimization to compute x that minimizes
    the frobenius norm ||target-processing(x)||_F

    Args:
        target: 1D array, the target design weight after processing function
        processing: a differentiable function (e.g. filter and projection) that processes the input
            design variable and outputs the actual design weights for the structure
            For example, we normally have some mapping function for filtering and projection

            def mapping(x, eta, beta):
                filtered_field = mpa.conic_filter(x, ...)
                projected_field = mpa.tanh_projection(filtered_field, beta, eta)
                return projected_field.flatten()

            If eta=0.5, and the initial beta is 8, then we can pass the following processing
            function to find the desired initialization

            processing = lambda x: mapping(x, 0.5, 8)

        maxiter: maximum number of iterations for the optimization

    Returns:
        Optimized design variables x
    """

    def design_diff(x):
        return npa.sum((processing(x) - target) ** 2)

    # scipy's `jac=True` takes the objective and its gradient from a single
    # call, which also avoids evaluating `processing` twice per iteration.
    f = value_and_grad(design_diff)

    n = len(target)
    x = target
    ftol = 1e-5
    # L-BFGS-B is the gradient-based, box-constrained solver in scipy, and the
    # design weights are bounded to [0,1]. Its `ftol` is the relative decrease
    # in the objective, the same convergence criterion used previously.
    result = minimize(
        f,
        x,
        jac=True,
        method="L-BFGS-B",
        bounds=[(0.0, 1.0)] * n,
        options={"maxiter": maxiter, "ftol": ftol},
    )
    x[:] = result.x
    return x
