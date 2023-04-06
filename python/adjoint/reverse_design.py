import nlopt
from autograd import numpy as npa
from autograd import grad


def reverse_design(target, processing, maxiter=50):
    def design_diff(x):
        return npa.sum((processing(x) - target) ** 2)

    def f(x, grad):
        grad[:] = grad(design_diff, 0)(x)
        return design_diff(x)

    algorithm = nlopt.LD_MMA
    n = len(target)
    x = target
    lb, ub = np.zeros((n,)), np.ones((n,))
    ftol = 1e-5
    for iters in range(num_betas):
        solver = nlopt.opt(algorithm, n)
        solver.set_lower_bounds(lb)
        solver.set_upper_bounds(ub)
        solver.set_min_objective(f)
        solver.set_maxeval(maxiter)
        solver.set_ftol_rel(ftol)
        x[:] = solver.optimize(x)
    return x
