"""Computes the extraction efficiency of a collection of dipoles in a disc.

tutorial reference: https://meep.readthedocs.io/en/latest/Python_Tutorials/Near_to_Far_Field_Spectra/#extraction-efficiency-of-a-disc-in-cylindrical-coordinates
"""

import math
from typing import Tuple

import matplotlib
import meep as mp
import numpy as np

matplotlib.use("agg")
import matplotlib.pyplot as plt


N_DISC = 2.4  # refractive index of disc
DISC_RADIUS_UM = 1.2  # radius of disc
WAVELENGTH_UM = 1.0  # wavelength (in vacuum)

# radius of quarter circle for computing flux in far field
FARFIELD_RADIUS_UM = 1000 * WAVELENGTH_UM


def plot_radiation_pattern_polar(theta_rad: np.ndarray, radial_flux: np.ndarray):
    """Plots the radiation pattern in polar coordinates.

    Args:
      theta_rad: angles of the radiation pattern. The angles increase clockwise
        with zero at the pole (+z direction) and π/2 at the equator (+r
        direction).
      radial_flux: radial flux of the far fields in polar coordinates.
    """
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(6, 6))
    ax.plot(
        theta_rad,
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


def plot_radiation_pattern_3d(theta_rad: np.ndarray, radial_flux: np.ndarray):
    """Plots the radiation pattern in 3d Cartesian coordinates.

    Args:
      theta_rad: angles of the radiation pattern.
      radial_flux: radial flux of the far fields in polar coordinates.
    """
    phis = np.linspace(0, 2 * np.pi, num_angles)

    xs = np.zeros((len(theta_rad), len(phis)))
    ys = np.zeros((len(theta_rad), len(phis)))
    zs = np.zeros((len(theta_rad), len(phis)))

    for i, theta in enumerate(theta_rad):
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


def radiation_pattern(
    theta_rad: np.ndarray, sim: mp.Simulation, n2f_mon: mp.DftNear2Far
) -> np.ndarray:
    """Computes the radiation pattern from the near fields.

    Args:
      theta_rad: angles of the radiation pattern.
      sim: a `Simulation` object.
      n2f_mon: a `DftNear2Far` object returned by `Simulation.add_near2far`.
    """
    e_field = np.zeros((theta_rad.shape[0], 3), dtype=np.complex128)
    h_field = np.zeros((theta_rad.shape[0], 3), dtype=np.complex128)
    for n in range(num_angles):
        far_field = sim.get_farfield(
            n2f_mon,
            mp.Vector3(
                FARFIELD_RADIUS_UM * math.sin(theta_rad[n]),
                0,
                FARFIELD_RADIUS_UM * math.cos(theta_rad[n]),
            ),
        )
        e_field[n, :] = [np.conj(far_field[j]) for j in range(3)]
        h_field[n, :] = [far_field[j + 3] for j in range(3)]

    flux_r = np.real(e_field[:, 1] * h_field[:, 2] - e_field[:, 2] * h_field[:, 1])
    flux_z = np.real(e_field[:, 0] * h_field[:, 1] - e_field[:, 1] * h_field[:, 0])
    flux_rz = np.sqrt(np.square(flux_r) + np.square(flux_z))

    return flux_rz


def radiation_pattern_flux(theta_rad: np.ndarray, radial_flux: np.ndarray) -> float:
    """Computes the total flux from the radiation pattern.

    Based on integrating the radiation pattern over solid angles spanned by
    polar angles in the range of [0, π/2].

    Args:
      theta_rad: angles of the radiation pattern.
      radial_flux: radial flux of the far fields in polar coordinates.
    """
    dphi = 2 * math.pi
    dtheta_rad = theta_rad[1] - theta_rad[0]

    total_flux = (
        np.sum(radial_flux * np.sin(theta_rad))
        * FARFIELD_RADIUS_UM**2
        * dtheta_rad
        * dphi
    )

    return total_flux


def dipole_in_disc(
    t_disc_um: float, h_disc: float, rpos: float, m: int, theta_rad: np.ndarray
) -> Tuple[float, np.ndarray]:
    """Computes the total flux and radiation pattern of a dipole in a disc.

    Args:
      t_disc_um: thickness of disc.
      h_disc: height of dipole above ground plane as fraction of t_disc_um.
      rpos: radial position of dipole.
      m: angular φ dependence of the fields exp(imφ).
      theta_rad: angles of the radiation pattern.

    Returns:
      A 2-tuple of the total flux and the radiation pattern.
    """
    resolution = 50  # pixels/μm

    t_pml_um = 0.5  # thickness of PML
    t_air_um = 1.0  # thickness of air padding above disc
    length_r_um = 5.0  # length of cell in r

    frequency = 1 / WAVELENGTH_UM  # center frequency of source/monitor

    # field decay threshold for runtime termination criteria
    decay_tol = 1e-6

    size_r = length_r_um + t_pml_um
    size_z = t_disc_um + t_air_um + t_pml_um
    cell_size = mp.Vector3(size_r, 0, size_z)

    boundary_layers = [
        mp.PML(t_pml_um, direction=mp.R),
        mp.PML(t_pml_um, direction=mp.Z, side=mp.High),
    ]

    # An Er source at r=0 needs to be slighty offset.
    # https://github.com/NanoComp/meep/issues/2704
    if rpos == 0:
        rpos = 1.5 / resolution

    src_cmpt = mp.Er
    src_pt = mp.Vector3(rpos, 0, -0.5 * size_z + h_disc * t_disc_um)
    sources = [
        mp.Source(
            src=mp.GaussianSource(frequency, fwidth=0.1 * frequency),
            component=src_cmpt,
            center=src_pt,
        )
    ]

    geometry = [
        mp.Block(
            material=mp.Medium(index=N_DISC),
            center=mp.Vector3(0.5 * DISC_RADIUS_UM, 0, -0.5 * size_z + 0.5 * t_disc_um),
            size=mp.Vector3(DISC_RADIUS_UM, mp.inf, t_disc_um),
        )
    ]

    sim = mp.Simulation(
        resolution=resolution,
        cell_size=cell_size,
        dimensions=mp.CYLINDRICAL,
        m=m,
        boundary_layers=boundary_layers,
        sources=sources,
        geometry=geometry,
    )

    n2f_mon = sim.add_near2far(
        frequency,
        0,
        1,
        mp.FluxRegion(
            center=mp.Vector3(0.5 * length_r_um, 0, 0.5 * size_z - t_pml_um),
            size=mp.Vector3(length_r_um, 0, 0),
        ),
        mp.FluxRegion(
            center=mp.Vector3(
                length_r_um, 0, 0.5 * size_z - t_pml_um - 0.5 * (t_air_um + t_disc_um)
            ),
            size=mp.Vector3(0, 0, t_air_um + t_disc_um),
        ),
    )

    sim.run(
        mp.dft_ldos(frequency, 0, 1),
        until_after_sources=mp.stop_when_fields_decayed(
            50,
            src_cmpt,
            src_pt,
            decay_tol,
        ),
    )

    delta_vol = 2 * np.pi * rpos / (resolution**2)
    dipole_flux = -np.real(sim.ldos_Fdata[0] * np.conj(sim.ldos_Jdata[0])) * delta_vol

    dipole_radiation_pattern = radiation_pattern(theta_rad, sim, n2f_mon)

    return dipole_flux, dipole_radiation_pattern


if __name__ == "__main__":
    disc_thickness = 0.7 * WAVELENGTH_UM / N_DISC
    dipole_height = 0.5
    num_dipoles = 11
    dipole_rpos = np.linspace(0, DISC_RADIUS_UM, num_dipoles)

    delta_rpos = dipole_rpos[2] - dipole_rpos[1]

    # number of angular grid points in [0, π/2]
    num_angles = 100

    # grid of polar angles for computing radiated flux in far field
    theta_rad = np.linspace(0, 0.5 * math.pi, num_angles)

    # r = 0 source requires a single simulation with m = ±1.
    m = -1
    dipole_flux, dipole_radiation_pattern = dipole_in_disc(
        disc_thickness, dipole_height, dipole_rpos[0], m, theta_rad
    )

    flux_total = dipole_flux * dipole_rpos[0] * delta_rpos
    radiation_pattern_total = dipole_radiation_pattern * dipole_rpos[0] * delta_rpos

    print(
        f"dipole:, {dipole_rpos[0]:.4f}, "
        f"{radiation_pattern_flux(theta_rad, dipole_radiation_pattern):.6f}"
    )

    # r > 0 source requires Fourier-series expansion of φ.
    flux_tol = 1e-3  # threshold flux to determine when to truncate expansion
    for rpos in dipole_rpos[1:]:
        dipole_flux_total = 0
        dipole_radiation_pattern_total = np.zeros((num_angles,))
        dipole_radiation_pattern_flux_max = 0
        m = 0
        while True:
            dipole_flux, dipole_radiation_pattern = dipole_in_disc(
                disc_thickness, dipole_height, rpos, m, theta_rad
            )
            dipole_flux_total += dipole_flux if m == 0 else 2 * dipole_flux
            dipole_radiation_pattern_total += (
                dipole_radiation_pattern if m == 0 else 2 * dipole_radiation_pattern
            )

            dipole_radiation_pattern_flux = radiation_pattern_flux(
                theta_rad, dipole_radiation_pattern
            )
            if dipole_radiation_pattern_flux > dipole_radiation_pattern_flux_max:
                dipole_radiation_pattern_flux_max = dipole_radiation_pattern_flux

            if m > 0 and (
                (dipole_radiation_pattern_flux / dipole_radiation_pattern_flux_max)
                < flux_tol
            ):
                break

            print(
                f"dipole-m:, {rpos:.4f}, {m}, " f"{dipole_radiation_pattern_flux:.6f}"
            )
            m += 1

        scale_factor = 1 / (2 * (rpos / dipole_rpos[0]) ** 2)
        flux_total += dipole_flux_total * scale_factor * rpos * delta_rpos
        radiation_pattern_total += (
            dipole_radiation_pattern_total * scale_factor * rpos * delta_rpos
        )

        dipole_radiation_pattern_total_flux = radiation_pattern_flux(
            theta_rad, dipole_radiation_pattern_total
        )
        print(
            f"dipole:, {rpos:.4f}, {m}, " f"{dipole_radiation_pattern_total_flux:.6f}"
        )

    radiation_pattern_total_flux = radiation_pattern_flux(
        theta_rad, radiation_pattern_total
    )
    extraction_efficiency = radiation_pattern_total_flux / flux_total
    print(f"exteff:, {extraction_efficiency:.6f}")

    plot_radiation_pattern_polar(radiation_pattern_total * FARFIELD_RADIUS_UM**2)
    plot_radiation_pattern_3d(radiation_pattern_total * FARFIELD_RADIUS_UM**2)
