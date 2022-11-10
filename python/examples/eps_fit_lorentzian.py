# fit complex refractive index profile over broad bandwidth imported from
# a file to a sum of Lorentzian polarizability terms using gradient-based
# optimization via NLopt (nlopt.readthedocs.io)
import matplotlib
import nlopt
import numpy as np

matplotlib.use("agg")
import matplotlib.pyplot as plt

import meep as mp


def lorentzfunc(p, x):
    """Return the complex ε profile given a set of Lorentzian parameters p
    (σ_0, ω_0, γ_0, σ_1, ω_1, γ_1, ...) for a set of frequencies x.
    """
    N = len(p) // 3
    y = np.zeros(len(x))
    for n in range(N):
        A_n = p[3 * n + 0]
        x_n = p[3 * n + 1]
        g_n = p[3 * n + 2]
        y = y + A_n / (np.square(x_n) - np.square(x) - 1j * x * g_n)
    return y


def lorentzerr(p, x, y, grad):
    """Return the error (or residual or loss funnction) as the L2 norm
    of the difference of ε(p,x) and y over a set of frequencies x as
    well as the gradient of this error with respect to each Lorentzian
    polarizability parameter in p and saving the result in grad.
    """
    N = len(p) // 3
    yp = lorentzfunc(p, x)
    val = np.sum(np.square(abs(y - yp)))
    for n in range(N):
        A_n = p[3 * n + 0]
        x_n = p[3 * n + 1]
        g_n = p[3 * n + 2]
        d = 1 / (np.square(x_n) - np.square(x) - 1j * x * g_n)
        if grad.size > 0:
            grad[3 * n + 0] = 2 * np.real(np.dot(np.conj(yp - y), d))
            grad[3 * n + 1] = (
                -4 * x_n * A_n * np.real(np.dot(np.conj(yp - y), np.square(d)))
            )
            grad[3 * n + 2] = (
                -2 * A_n * np.imag(np.dot(np.conj(yp - y), x * np.square(d)))
            )
    return val


def lorentzfit(p0, x, y, alg=nlopt.LD_LBFGS, tol=1e-25, maxeval=10000):
    """Return the optimal Lorentzian polarizability parameters and error
    which minimize the error in ε(p0,x) relative to y for an initial
    set of Lorentzian polarizability parameters p0 over a set of
    frequencies x using the NLopt algorithm alg for a relative
    tolerance tol and a maximum number of iterations maxeval.
    """
    opt = nlopt.opt(alg, len(p0))
    opt.set_ftol_rel(tol)
    opt.set_maxeval(maxeval)
    opt.set_lower_bounds(np.zeros(len(p0)))
    opt.set_upper_bounds(float("inf") * np.ones(len(p0)))
    opt.set_min_objective(lambda p, grad: lorentzerr(p, x, y, grad))
    local_opt = nlopt.opt(nlopt.LD_LBFGS, len(p0))
    local_opt.set_ftol_rel(1e-10)
    local_opt.set_xtol_rel(1e-8)
    opt.set_local_optimizer(local_opt)
    popt = opt.optimize(p0)
    minf = opt.last_optimum_value()
    return popt, minf


# format of input data file is three comma-separated columns:
# wavelength (nm), real(n), complex(n)
mydata = np.genfromtxt("mymaterial.csv", delimiter=",")
n = mydata[:, 1] + 1j * mydata[:, 2]
eps_inf = 1.1  # should be >= 1.0 for stability and chosen such that np.amin(np.real(eps)) is ~1.0
eps = np.square(n) - eps_inf
wl = mydata[:, 0]
wl_min = 399  # minimum wavelength (units of nm)
wl_max = 701  # maximum wavelength (units of nm)
start_idx = np.where(wl > wl_min)
idx_start = start_idx[0][0]
end_idx = np.where(wl < wl_max)
idx_end = end_idx[0][-1] + 1
freqs = 1000 / wl  # units of 1/μm
freqs_reduced = freqs[idx_start:idx_end]
wl_reduced = wl[idx_start:idx_end]
eps_reduced = eps[idx_start:idx_end]

num_lorentzians = 2
num_repeat = (
    30  # number of times to repeat local optimization with random initial values
)
ps = np.zeros((num_repeat, 3 * num_lorentzians))
mins = np.zeros(num_repeat)
for m in range(num_repeat):
    # specify Lorentzian polarizability terms each consisting of three parameters (σ_0, ω_0, γ_0)
    # with random initial values
    # note: for the case of no absorption, set γ (every third parameter) to zero
    p_rand = [10 ** (np.random.random()) for _ in range(3 * num_lorentzians)]
    ps[m, :], mins[m] = lorentzfit(
        p_rand, freqs_reduced, eps_reduced, nlopt.LD_MMA, 1e-25, 50000
    )
    print(f"iteration{m}:, {ps[m,:]}, {mins[m]}")

# find the best performing set of parameters
idx_opt = np.where(np.min(mins) == mins)
print(f"optimal:, {ps[idx_opt]}, {mins[idx_opt]}")


# define a Meep `Medium` class object using the optimal fitting parameters
E_susceptibilities = []

for n in range(num_lorentzians):
    mymaterial_freq = ps[idx_opt][0][3 * n + 1]
    mymaterial_gamma = ps[idx_opt][0][3 * n + 2]

    if mymaterial_freq == 0:
        mymaterial_sigma = ps[idx_opt][0][3 * n + 0]
        E_susceptibilities.append(
            mp.DrudeSusceptibility(
                frequency=1.0, gamma=mymaterial_gamma, sigma=mymaterial_sigma
            )
        )
    else:
        mymaterial_sigma = ps[idx_opt][0][3 * n + 0] / mymaterial_freq**2
        E_susceptibilities.append(
            mp.LorentzianSusceptibility(
                frequency=mymaterial_freq,
                gamma=mymaterial_gamma,
                sigma=mymaterial_sigma,
            )
        )

mymaterial = mp.Medium(epsilon=eps_inf, E_susceptibilities=E_susceptibilities)


# plot the fit and the actual data for comparison

mymaterial_eps = [mymaterial.epsilon(f)[0][0] for f in freqs_reduced]

plt.subplot(1, 2, 1)
plt.plot(wl_reduced, np.real(eps_reduced) + eps_inf, "bo-", label="actual")
plt.plot(wl_reduced, np.real(mymaterial_eps), "ro-", label="fit")
plt.xlabel("wavelength (nm)")
plt.ylabel(r"real($\epsilon$)")
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(wl_reduced, np.imag(eps_reduced), "bo-", label="actual")
plt.plot(wl_reduced, np.imag(mymaterial_eps), "ro-", label="fit")
plt.xlabel("wavelength (nm)")
plt.ylabel(r"imag($\epsilon$)")
plt.legend()

plt.suptitle(
    "Comparison of Actual Material Data and Fit using Drude-Lorentzian Susceptibility"
)
plt.subplots_adjust(wspace=0.3)
plt.savefig("eps_fit_sample.png", dpi=150, bbox_inches="tight")
