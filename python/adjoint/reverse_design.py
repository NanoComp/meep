from autograd import numpy as npa
from autograd import grad
from typing import Callable, List


def reverse_design(target: List[float], processing: Callable, maxiter: int = 50):
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

    def f(x, gradient):
        gradient[:] = grad(design_diff, 0)(x)
        return design_diff(x)
        
    import nlopt
    algorithm = nlopt.LD_CCSAQ
    n = len(target)
    x = target
    lb, ub = np.zeros((n,)), np.ones((n,))
    ftol = 1e-5
    solver = nlopt.opt(algorithm, n)
    solver.set_lower_bounds(lb)
    solver.set_upper_bounds(ub)
    solver.set_min_objective(f)
    solver.set_maxeval(maxiter)
    solver.set_ftol_rel(ftol)
    x[:] = solver.optimize(x)
    return x
