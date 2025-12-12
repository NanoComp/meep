"""Computes the extraction efficiency of a collection of dipoles in a disc.

tutorial reference:
https://meep.readthedocs.io/en/latest/Python_Tutorials/Near_to_Far_Field_Spectra/#extraction-efficiency-of-a-disc-in-cylindrical-coordinates
"""

import math
from typing import Tuple

import matplotlib.pyplot as plt
import meep as mp
import numpy as np


RESOLUTION_UM = 50
WAVELENGTH_UM = 1.0
N_DISC = 2.4
DISC_RADIUS_UM = 1.2
DISC_THICKNESS_UM = 0.7 * WAVELENGTH_UM / N_DISC
NUM_FARFIELD_PTS = 200
FARFIELD_RADIUS_UM = 1e6 * WAVELENGTH_UM
NUM_DIPOLES = 11
GREENCYL_TOL = 1e-6

farfield_angles = np.linspace(0, 0.5 * math.pi, NUM_FARFIELD_PTS)


def plot_radiation_pattern_polar(radial_flux: np.ndarray):
    """Plots the radiation pattern in polar coordinates.

    Args:
      radial_flux: radial flux of the far fields at each angle.
    """
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(6, 6))
    ax.plot(
        farfield_angles,
        radial_flux,
        "b-",
    )
    ax.set_theta_direction(-1)
    ax.set_theta_offset(0.5 * math.pi)
    ax.set_thetalim(0, 0.5 * math.pi)
    ax.grid(True)
    ax.set_rlabel_position(22)
    ax.set_ylabel("radial flux (a.u.)")
    ax.set_title("radiation pattern in polar coordinates")

    if mp.am_master():
        fig.savefig(
            "disc_radiation_pattern_polar.png",
            dpi=150,
            bbox_inches="tight",
        )


def plot_radiation_pattern_3d(radial_flux: np.ndarray):
    """Plots the radiation pattern in 3d Cartesian coordinates.

    Args:
      radial_flux: radial flux of the far fields at each angle.
    """
    phis = np.linspace(0, 2 * np.pi, NUM_FARFIELD_PTS)

    xs = np.zeros((NUM_FARFIELD_PTS, NUM_FARFIELD_PTS))
    ys = np.zeros((NUM_FARFIELD_PTS, NUM_FARFIELD_PTS))
    zs = np.zeros((NUM_FARFIELD_PTS, NUM_FARFIELD_PTS))

    for i, theta in enumerate(farfield_angles):
        for j, phi in enumerate(phis):
            xs[i, j] = radial_flux[i] * np.sin(theta) * np.cos(phi)
            ys[i, j] = radial_flux[i] * np.sin(theta) * np.sin(phi)
            zs[i, j] = radial_flux[i] * np.cos(theta)

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(6, 6))
    ax.plot_surface(xs, ys, zs, cmap="inferno")
    ax.set_title("radiation pattern in 3d")
    ax.set_box_aspect((np.amax(xs), np.amax(ys), np.amax(zs)))
    ax.set_zlabel("radial flux (a.u.)")
    ax.set(xticklabels=[], yticklabels=[])

    if mp.am_master():
        fig.savefig(
            "disc_radiation_pattern_3d.png",
            dpi=150,
            bbox_inches="tight",
        )


def radiation_pattern(sim: mp.Simulation, n2f_mon: mp.DftNear2Far) -> np.ndarray:
    """Computes the radiation pattern from the near fields.

    Args:
      sim: a `Simulation` object.
      n2f_mon: a `DftNear2Far` object returned by `Simulation.add_near2far`.

    Returns:
      The radiation pattern (radial flux at each angle) as a 1d array.
    """
    e_field = np.zeros((NUM_FARFIELD_PTS, 3), dtype=np.complex128)
    h_field = np.zeros((NUM_FARFIELD_PTS, 3), dtype=np.complex128)
    for n in range(NUM_FARFIELD_PTS):
        far_field = sim.get_farfield(
            n2f_mon,
            mp.Vector3(
                FARFIELD_RADIUS_UM * math.sin(farfield_angles[n]),
                0,
                FARFIELD_RADIUS_UM * math.cos(farfield_angles[n]),
            ),
            GREENCYL_TOL,
        )
        e_field[n, :] = [far_field[j] for j in range(3)]
        h_field[n, :] = [far_field[j + 3] for j in range(3)]

    flux_x = np.real(
        np.conj(e_field[:, 1]) * h_field[:, 2] - np.conj(e_field[:, 2]) * h_field[:, 1]
    )
    flux_z = np.real(
        np.conj(e_field[:, 0]) * h_field[:, 1] - np.conj(e_field[:, 1]) * h_field[:, 0]
    )
    flux_r = np.sqrt(np.square(flux_x) + np.square(flux_z))

    return flux_r


def radiation_pattern_flux(radial_flux: np.ndarray) -> float:
    """Computes the total flux from the radiation pattern.

    Based on integrating the radiation pattern over solid angles
    spanned by polar angles in the range of [0, π/2].

    Args:
      radial_flux: radial flux of the far fields at each angle.
    """
    total_flux = (
        2
        * math.pi
        * FARFIELD_RADIUS_UM**2
        * np.trapezoid(radial_flux * np.sin(farfield_angles), farfield_angles)
    )

    return total_flux


def dipole_in_disc(zpos: float, rpos_um: float, m: int) -> Tuple[float, np.ndarray]:
    """Computes the total flux and radiation pattern of a dipole in a disc.

    Args:
      zpos: height of dipole above ground plane as fraction of disc thickness.
      rpos_um: radial position of dipole.
      m: angular φ dependence of the fields exp(imφ).

    Returns:
      A 2-tuple of the total flux and the radiation pattern.
    """
    pml_um = 1.0  # thickness of PML
    padding_um = 1.0  # thickness of air padding above disc
    r_um = 4.0  # length of cell in r

    frequency = 1 / WAVELENGTH_UM  # center frequency of source/monitor

    # Runtime termination criteria.
    dft_decay_threshold = 1e-4

    size_r = r_um + pml_um
    size_z = DISC_THICKNESS_UM + padding_um + pml_um
    cell_size = mp.Vector3(size_r, 0, size_z)

    boundary_layers = [
        mp.PML(pml_um, direction=mp.R),
        mp.PML(pml_um, direction=mp.Z, side=mp.High),
    ]

    src_pt = mp.Vector3(rpos_um, 0, -0.5 * size_z + zpos * DISC_THICKNESS_UM)
    sources = [
        mp.Source(
            src=mp.GaussianSource(frequency, fwidth=0.05 * frequency),
            component=mp.Er,
            center=src_pt,
        )
    ]

    geometry = [
        mp.Block(
            material=mp.Medium(index=N_DISC),
            center=mp.Vector3(
                0.5 * DISC_RADIUS_UM, 0, -0.5 * size_z + 0.5 * DISC_THICKNESS_UM
            ),
            size=mp.Vector3(DISC_RADIUS_UM, mp.inf, DISC_THICKNESS_UM),
        )
    ]

    sim = mp.Simulation(
        resolution=RESOLUTION_UM,
        cell_size=cell_size,
        dimensions=mp.CYLINDRICAL,
        m=m,
        boundary_layers=boundary_layers,
        sources=sources,
        geometry=geometry,
        force_complex_fields=True,
    )

    n2f_mon = sim.add_near2far(
        frequency,
        0,
        1,
        mp.FluxRegion(
            center=mp.Vector3(0.5 * r_um, 0, 0.5 * size_z - pml_um),
            size=mp.Vector3(r_um, 0, 0),
        ),
        mp.FluxRegion(
            center=mp.Vector3(
                r_um, 0, 0.5 * size_z - pml_um - 0.5 * (padding_um + DISC_THICKNESS_UM)
            ),
            size=mp.Vector3(0, 0, padding_um + DISC_THICKNESS_UM),
        ),
    )

    sim.run(
        mp.dft_ldos(frequency, 0, 1),
        until_after_sources=mp.stop_when_dft_decayed(
            tol=dft_decay_threshold,
        ),
    )

    delta_vol = 2 * np.pi * rpos_um / (RESOLUTION_UM**2)
    dipole_flux = -np.real(sim.ldos_Fdata[0] * np.conj(sim.ldos_Jdata[0])) * delta_vol

    dipole_radiation_pattern = radiation_pattern(sim, n2f_mon)

    return dipole_flux, dipole_radiation_pattern


if __name__ == "__main__":
    dipole_height = 0.5
    dipole_rpos_um = np.linspace(0, DISC_RADIUS_UM, NUM_DIPOLES)
    delta_rpos_um = DISC_RADIUS_UM / (NUM_DIPOLES - 1)

    # 1. Er source at r = 0 requires a single simulation with m = ±1.

    # An Er source at r = 0 needs to be slighty offset due to a bug.
    # https://github.com/NanoComp/meep/issues/2704
    dipole_rpos_um[0] = 1.5 / RESOLUTION_UM

    m = -1
    dipole_flux, dipole_radiation_pattern = dipole_in_disc(
        dipole_height,
        dipole_rpos_um[0],
        m,
    )

    flux_total = dipole_flux * dipole_rpos_um[0] * delta_rpos_um
    radiation_pattern_total = (
        dipole_radiation_pattern * dipole_rpos_um[0] * delta_rpos_um
    )

    print(
        f"dipole:, {dipole_rpos_um[0]:.4f}, "
        f"{radiation_pattern_flux(dipole_radiation_pattern):.6f}"
    )

    # 2. Er source at r > 0 requires Fourier-series expansion of φ.

    # Threshold flux to determine when to truncate expansion.
    flux_decay_threshold = 1e-2

    for rpos_um in dipole_rpos_um[1:]:
        dipole_flux_total = 0
        dipole_radiation_pattern_total = np.zeros(NUM_FARFIELD_PTS)
        dipole_radiation_pattern_flux_max = 0
        m = 0
        while True:
            dipole_flux, dipole_radiation_pattern = dipole_in_disc(
                dipole_height, rpos_um, m
            )
            dipole_flux_total += dipole_flux * (1 if m == 0 else 2)
            dipole_radiation_pattern_total += dipole_radiation_pattern * (
                1 if m == 0 else 2
            )

            dipole_radiation_pattern_flux = radiation_pattern_flux(
                dipole_radiation_pattern
            )
            print(
                f"dipole:, {rpos_um:.4f}, {m}, " f"{dipole_radiation_pattern_flux:.6f}"
            )

            if dipole_radiation_pattern_flux > dipole_radiation_pattern_flux_max:
                dipole_radiation_pattern_flux_max = dipole_radiation_pattern_flux

            if (
                m > 0
                and (dipole_radiation_pattern_flux / dipole_radiation_pattern_flux_max)
                < flux_decay_threshold
            ):
                break
            else:
                m += 1

        dipole_position_scale_factor = 0.5 * (dipole_rpos_um[0] / rpos_um) ** 2
        flux_total += (
            dipole_flux_total * dipole_position_scale_factor * rpos_um * delta_rpos_um
        )
        radiation_pattern_total += (
            dipole_radiation_pattern_total
            * dipole_position_scale_factor
            * rpos_um
            * delta_rpos_um
        )

    radiation_pattern_total_flux = radiation_pattern_flux(radiation_pattern_total)
    extraction_efficiency = radiation_pattern_total_flux / flux_total
    print(f"exteff:, {extraction_efficiency:.6f}")

    radiation_pattern_scaled = radiation_pattern_total * FARFIELD_RADIUS_UM**2
    plot_radiation_pattern_polar(radiation_pattern_scaled)
    plot_radiation_pattern_3d(radiation_pattern_scaled)
