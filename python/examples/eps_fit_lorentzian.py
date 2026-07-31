"""Fit a complex refractive index profile to a sum of Lorentzian (and,
optionally, Drude) susceptibility terms and build a `Medium` from the result.

The complex permittivity is modeled as

    ε(f) = ε_inf + Σ_n  A_n / (ω_n² − f² − i f γ_n)

and the parameters (A_n, ω_n, γ_n) are found by nonlinear least-squares fitting
to tabulated data using SciPy (scipy.optimize.least_squares). The fit is
repeated from several random starting points to avoid poor local minima.

Usage:
    python eps_fit_lorentzian.py mymaterial.csv --num-lorentzians 2 \\
        --eps-inf 1.0 --wl-min 400 --wl-max 700

The input CSV has three comma-separated columns and no header:
    wavelength (nm), real(n), imag(n)
"""
import argparse
from typing import Tuple

import matplotlib.pyplot as plt
import meep as mp
import numpy as np
from scipy.optimize import least_squares


def lorentzfunc(params: np.ndarray, freqs: np.ndarray) -> np.ndarray:
    """Returns the complex ε - ε_inf profile at ``freqs``.

    ``params`` is a flat array of Lorentzian parameters
    (A_0, ω_0, γ_0, A_1, ω_1, γ_1, ...).
    """
    num_terms = len(params) // 3
    eps = np.zeros(len(freqs), dtype=complex)
    for n in range(num_terms):
        a_n, w_n, g_n = params[3 * n : 3 * n + 3]
        eps += a_n / (np.square(w_n) - np.square(freqs) - 1j * freqs * g_n)
    return eps


def _residuals(params: np.ndarray, freqs: np.ndarray, eps: np.ndarray) -> np.ndarray:
    """Real-valued residual vector [real; imag] for ``least_squares``."""
    diff = lorentzfunc(params, freqs) - eps
    return np.concatenate((diff.real, diff.imag))


def _jacobian(params: np.ndarray, freqs: np.ndarray, eps: np.ndarray) -> np.ndarray:
    """Analytic Jacobian of ``_residuals`` with respect to ``params``."""
    num_terms = len(params) // 3
    jac = np.zeros((len(freqs), 3 * num_terms), dtype=complex)
    for n in range(num_terms):
        a_n, w_n, g_n = params[3 * n : 3 * n + 3]
        denom = np.square(w_n) - np.square(freqs) - 1j * freqs * g_n
        jac[:, 3 * n + 0] = 1 / denom
        jac[:, 3 * n + 1] = -2 * w_n * a_n / np.square(denom)
        jac[:, 3 * n + 2] = 1j * freqs * a_n / np.square(denom)
    return np.vstack((jac.real, jac.imag))


def lorentzfit(
    p0: np.ndarray,
    freqs: np.ndarray,
    eps: np.ndarray,
    tol: float = 1e-15,
    maxeval: int = 25000,
) -> Tuple[np.ndarray, float]:
    """Fit Lorentzian parameters to ``eps`` sampled at ``freqs``.

    Returns the optimal (non-negative) parameters and the residual sum of
    squares Σ|ε(p) − eps|².
    """
    result = least_squares(
        _residuals,
        p0,
        jac=_jacobian,
        args=(freqs, eps),
        bounds=(0, np.inf),
        method="trf",
        ftol=tol,
        xtol=tol,
        gtol=tol,
        max_nfev=maxeval,
    )
    # least_squares minimizes 0.5 Σ residual²; report the full Σ|Δ|².
    return result.x, 2 * result.cost


def fit_medium(
    freqs: np.ndarray,
    eps: np.ndarray,
    eps_inf: float,
    num_lorentzians: int,
    num_repeat: int,
    seed: int,
) -> Tuple[mp.Medium, np.ndarray, float]:
    """Fit ``num_lorentzians`` terms via ``num_repeat`` random restarts and
    return the resulting `Medium`, the optimal parameters, and the error."""
    rng = np.random.RandomState(seed)
    best_params = None
    best_err = np.inf
    for m in range(num_repeat):
        # Each of the three parameters per term is seeded in [1, 10).
        # Note: for a lossless material γ optimizes toward zero.
        # FIXME: the ranges for w_n and g_n should depend on the frequency range,
        #        and the range for a_n should depend on the scale of the ε data
        p_rand = 10 ** rng.random(3 * num_lorentzians)
        params, err = lorentzfit(p_rand, freqs, eps)
        params_str = "( " + ", ".join(f"{p:.4f}" for p in params) + " )"
        print(f"iteration: {m:3d}, {params_str}, {err:.6f}")
        if err < best_err:
            best_err = err
            best_params = params

    params_str = "( " + ", ".join(f"{p:.4f}" for p in best_params) + " )"
    print(f"optimal: {params_str}, {best_err:.6f}")

    susceptibilities = []
    for n in range(num_lorentzians):
        a_n, w_n, g_n = best_params[3 * n : 3 * n + 3]
        if w_n == 0:
            susceptibilities.append(
                mp.DrudeSusceptibility(frequency=1.0, gamma=g_n, sigma=a_n)
            )
        else:
            susceptibilities.append(
                mp.LorentzianSusceptibility(
                    frequency=w_n, gamma=g_n, sigma=a_n / np.square(w_n)
                )
            )
    medium = mp.Medium(epsilon=eps_inf, E_susceptibilities=susceptibilities)
    return medium, best_params, best_err


def plot_fit(
    medium: mp.Medium,
    wl: np.ndarray,
    freqs: np.ndarray,
    eps: np.ndarray,
    eps_inf: float,
    output: str,
) -> None:
    """Plot the fitted permittivity against the input data."""
    eps_fit = np.array([medium.epsilon(f)[0][0] for f in freqs])

    fig, ax = plt.subplots(ncols=2)
    ax[0].plot(wl, np.real(eps) + eps_inf, "bo-", label="actual")
    ax[0].plot(wl, np.real(eps_fit), "ro-", label="fit")
    ax[0].set_xlabel("wavelength (nm)")
    ax[0].set_ylabel(r"real($\epsilon$)")
    ax[0].legend()

    ax[1].plot(wl, np.imag(eps), "bo-", label="actual")
    ax[1].plot(wl, np.imag(eps_fit), "ro-", label="fit")
    ax[1].set_xlabel("wavelength (nm)")
    ax[1].set_ylabel(r"imag($\epsilon$)")
    ax[1].legend()

    fig.suptitle(
        "Comparison of Actual Material Data and Fit\n"
        "using Drude-Lorentzian Susceptibility"
    )
    fig.subplots_adjust(wspace=0.3)
    fig.savefig(output, dpi=150, bbox_inches="tight")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "csv",
        nargs="?",
        default="mymaterial.csv",
        help="input CSV: wavelength (nm), real(n), imag(n) (default: mymaterial.csv)",
    )
    parser.add_argument(
        "--eps-inf",
        type=float,
        default=1.0,
        help="instantaneous (infinite-frequency) permittivity. Must be ≥ 1.0 for "
        "stability with default Meep Courant factor (default: 1.0)",
    )
    parser.add_argument(
        "--num-lorentzians",
        type=int,
        default=2,
        help="number of Lorentzian terms in the fit (default: 2)",
    )
    parser.add_argument(
        "--wl-min",
        type=float,
        default=400.0,
        help="minimum wavelength to fit, in nm (default: 400)",
    )
    parser.add_argument(
        "--wl-max",
        type=float,
        default=700.0,
        help="maximum wavelength to fit, in nm (default: 700)",
    )
    parser.add_argument(
        "--num-repeat",
        type=int,
        default=30,
        help="number of random restarts of the local optimization (default: 30)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="seed for the random restarts, for reproducitiblity (default: 0)",
    )
    parser.add_argument(
        "--output",
        default="eps_fit_sample.png",
        help="filename for the comparison plot (default: eps_fit_sample.png)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Import the complex refractive index profile.
    data = np.genfromtxt(args.csv, delimiter=",")
    wl = data[:, 0]
    n = data[:, 1] + 1j * data[:, 2]
    eps = np.square(n) - args.eps_inf

    # Restrict the fit to [wl_min, wl_max]. The fitting variable is frequency
    # (units of 1/μm) rather than wavelength.
    mask = (wl >= args.wl_min) & (wl <= args.wl_max)
    wl_reduced = wl[mask]
    eps_reduced = eps[mask]
    freqs_reduced = 1000 / wl_reduced  # wl is in nm

    medium, _, _ = fit_medium(
        freqs_reduced,
        eps_reduced,
        args.eps_inf,
        args.num_lorentzians,
        args.num_repeat,
        args.seed,
    )

    plot_fit(medium, wl_reduced, freqs_reduced, eps_reduced, args.eps_inf, args.output)


if __name__ == "__main__":
    main()
