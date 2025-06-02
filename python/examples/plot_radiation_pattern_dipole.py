"""Auxiliary plotting routines for the radiation pattern of a dipole."""

import math
from typing import Tuple

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np


NUM_POLAR = 30
NUM_AZIMUTH = 50


def plot_radiation_pattern_3D(dipole_pol: str, radial_flux: np.ndarray):
    """Plots the 3D radiation pattern as a contour plot.

    Args:
        dipole_pol: the polarization the electric dipole. Either x or y.
        radial_flux: the radial flux in polar coordinates.
    """
    polar_rad = np.linspace(0, 0.5 * np.pi, NUM_POLAR)
    azimuth_rad = np.linspace(0, 2 * np.pi, NUM_AZIMUTH)
    x = np.sin(polar_rad[:, np.newaxis]) * np.cos(azimuth_rad)
    y = np.sin(polar_rad[:, np.newaxis]) * np.sin(azimuth_rad)
    normalized_radial_flux = radial_flux / np.max(radial_flux)

    if dipole_pol == "x":
        analytic_radial_flux = np.sin(np.arccos(x)) ** 2
    elif dipole_pol == "y":
        analytic_radial_flux = np.sin(np.arccos(y)) ** 2

    fig, ax = plt.subplots(ncols=2, figsize=(8.5, 4), constrained_layout=True)

    _ = ax[0].tricontourf(
        x.flatten(),
        y.flatten(),
        normalized_radial_flux.flatten(),
        levels=100,
        cmap="inferno",
    )
    ax[0].set_aspect("equal")
    ax[0].axis(False)
    ax[0].set_title("Meep")

    im = ax[1].tricontourf(
        x.flatten(),
        y.flatten(),
        analytic_radial_flux.flatten(),
        levels=100,
        cmap="inferno",
    )
    ax[1].set_aspect("equal")
    ax[1].axis(False)
    ax[1].set_title("analytic")

    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = fig.colorbar(im, cax=cax, orientation="vertical")
    cbar.set_ticks(list(np.arange(0, 1.1, 0.2)))
    cbar.set_ticklabels([f"{t:.1f}" for t in np.arange(0, 1.1, 0.2)])

    fig.suptitle(
        f"radiation pattern of an $E_{dipole_pol}$ dipole in vacuum", size="x-large"
    )

    fig.savefig(
        "dipole_radiation_pattern_3D.png",
        dpi=150,
        bbox_inches="tight",
    )


def plot_radiation_pattern(dipole_pol: str, radial_flux: np.ndarray):
    """Plots the radiation pattern for φ=0 (xz plane) in polar coordinates.

    The angles increase clockwise with zero in the +z direction (the "pole")
    and π/2 in the +r direction (the "equator").

    Args:
        dipole_pol: the polarization the electric dipole. Either x or y.
        radial_flux: the radial flux in polar coordinates.
    """
    zero_azimuth_idx = 0
    normalized_radial_flux = radial_flux[:, zero_azimuth_idx] / np.max(
        radial_flux[:, zero_azimuth_idx]
    )

    polar_rad = np.linspace(0, 0.5 * np.pi, NUM_POLAR)

    # Omit polar angles with zero radial flux.
    nonzero_mask = np.nonzero(normalized_radial_flux)
    polar_rad = polar_rad[nonzero_mask]
    normalized_radial_flux = normalized_radial_flux[nonzero_mask]

    if dipole_pol == "x":
        dipole_radial_flux = np.cos(polar_rad) ** 2
        dipole_radial_flux_label = r"$\cos^2θ$"
        dipole_name = "$E_x$"
    else:
        dipole_radial_flux = np.ones(nonzero_mask[0].shape[0])
        dipole_radial_flux_label = "constant (1.0)"
        dipole_name = "$E_y$"

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(6, 6))
    ax.plot(polar_rad, normalized_radial_flux, "b-", label="Meep")
    ax.plot(polar_rad, dipole_radial_flux, "r--", label=dipole_radial_flux_label)
    ax.legend()
    ax.set_theta_direction(-1)
    ax.set_theta_offset(0.5 * math.pi)
    ax.set_thetalim(0, 0.5 * math.pi)
    ax.set_rmax(1)
    ax.set_rticks([0, 0.5, 1])
    ax.grid(True)
    ax.set_rlabel_position(22)
    ax.set_ylabel("radial flux (a.u.)")
    ax.set_title(f"radiation pattern (φ = 0) of an {dipole_name} dipole in vacuum")

    fig.savefig(
        "dipole_radiation_pattern_phi0.png",
        dpi=150,
        bbox_inches="tight",
    )

    relative_error = np.linalg.norm(
        normalized_radial_flux - dipole_radial_flux
    ) / np.linalg.norm(dipole_radial_flux)
    print(f"relative error in radiation pattern (φ = 0):, {relative_error}")


if __name__ == "__main__":
    data = np.load("dipole_radiation_pattern.npz")
    NUM_POLAR = data["NUM_POLAR"]
    NUM_AZIMUTH = data["NUM_AZIMUTH"]
    dipole_pol = data["dipole_pol"]
    radial_flux = data["radial_flux"]

    plot_radiation_pattern(dipole_pol, radial_flux)
    plot_radiation_pattern_3D(dipole_pol, radial_flux)
