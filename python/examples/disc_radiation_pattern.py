"""Radiated flux of a lossless dielectric disc in cylindrical coordinates.

Tutorial Reference:

https://meep.readthedocs.io/en/latest/Python_Tutorials/Near_to_Far_Field_Spectra/#radiation-pattern-of-a-disc-in-cylindrical-coordinates
"""

import math
from typing import Tuple

import matplotlib.pyplot as plt
import meep as mp
import numpy as np


RESOLUTION_UM = 100
PML_UM = 0.5
AIR_UM = 1.0
DISC_RADIUS_UM = 1.2
N_DISC = 2.4
WAVELENGTH_UM = 1.0
NUM_AZIMUTH = 50
NUM_POLAR = 100
FARFIELD_RADIUS_UM = 1e6 * WAVELENGTH_UM
FIELD_DECAY_TOL = 1e-8
FIELD_DECAY_PERIOD = 50
GREENCYL_TOL = 1e-8

frequency = 1 / WAVELENGTH_UM
polar_rad = np.linspace(0, 0.5 * math.pi, NUM_POLAR)
azimuth_rad = np.linspace(0, 2 * math.pi, NUM_AZIMUTH)


def plot_radiation_pattern_polar(radial_flux: np.ndarray):
    """Plots the radiation pattern in polar coordinates.

    The angles increase clockwise with zero at the top (+z direction).

    Args:
        radial_flux: radial flux of the far fields in polar coordinates.
    """
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(6, 6))
    ax.plot(
        polar_rad,
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
        radial_flux: radial flux of the far fields in polar coordinates.
    """
    x_coord = np.zeros((NUM_POLAR, NUM_AZIMUTH))
    y_coord = np.zeros((NUM_POLAR, NUM_AZIMUTH))
    z_coord = np.zeros((NUM_POLAR, NUM_AZIMUTH))

    for i in range(NUM_POLAR):
        for j in range(NUM_AZIMUTH):
            x_coord[i, j] = (
                radial_flux[i] * np.sin(polar_rad[i]) * np.cos(azimuth_rad[j])
            )
            y_coord[i, j] = (
                radial_flux[i] * np.sin(polar_rad[i]) * np.sin(azimuth_rad[j])
            )
            z_coord[i, j] = radial_flux[i] * np.cos(polar_rad[i])

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(6, 6))
    ax.plot_surface(x_coord, y_coord, z_coord, cmap="inferno")
    ax.set_title("radiation pattern in 3d")
    ax.set_box_aspect((np.amax(x_coord), np.amax(y_coord), np.amax(z_coord)))
    ax.set_zlabel("radial flux (a.u.)")
    ax.set(xticklabels=[], yticklabels=[])

    if mp.am_master():
        fig.savefig(
            "disc_radiation_pattern_3d.png",
            dpi=150,
            bbox_inches="tight",
        )


def radiation_pattern(sim: mp.Simulation, n2f_mon: mp.DftNear2Far) -> np.ndarray:
    """Computes the radiation pattern from the far fields.

    Args:
        sim: a `Simulation` object.
        n2f_mon: a `DftNear2Far` object returned by `Simulation.add_near2far`.

    Returns:
        Array of radial Poynting flux, one for each point on the circumference of
        a quarter circle with angular range of [0, π/2] rad. 0 rad is the +z
        direction and π/2 is +r.
    """
    e_field = np.zeros((NUM_POLAR, 3), dtype=np.complex128)
    h_field = np.zeros((NUM_POLAR, 3), dtype=np.complex128)
    for i in range(NUM_POLAR):
        far_field = sim.get_farfield(
            n2f_mon,
            mp.Vector3(
                FARFIELD_RADIUS_UM * math.sin(polar_rad[i]),
                0,
                FARFIELD_RADIUS_UM * math.cos(polar_rad[i]),
            ),
            GREENCYL_TOL,
        )
        e_field[i, :] = [far_field[j] for j in range(3)]
        h_field[i, :] = [far_field[j + 3] for j in range(3)]

    flux_x = np.real(
        np.conj(e_field[:, 1]) * h_field[:, 2] - np.conj(e_field[:, 2]) * h_field[:, 1]
    )
    flux_z = np.real(
        np.conj(e_field[:, 0]) * h_field[:, 1] - np.conj(e_field[:, 1]) * h_field[:, 0]
    )
    flux_r = np.sqrt(np.square(flux_x) + np.square(flux_z))

    return flux_r


def disc_radiated_flux(disc_um: float, source_zpos: float) -> Tuple[float, float]:
    """Computes the  radiated flux from a "ring" current source within a disc.

    Args:
        disc_um: thickness of dielectric disc.
        source_zpos: height of the dipole source above the ground plane as
            a fraction of disc_um.

    Returns:
        A 2-tuple of the total flux computed using the near and far fields,
        respectively.
    """
    cell_r_um = 6.0
    sr = cell_r_um + PML_UM
    sz = disc_um + AIR_UM + PML_UM
    cell_size = mp.Vector3(sr, 0, sz)

    boundary_layers = [
        mp.PML(PML_UM, direction=mp.R),
        mp.PML(PML_UM, direction=mp.Z, side=mp.High),
    ]

    src_cmpt = mp.Er
    src_pt = mp.Vector3(0.5 * DISC_RADIUS_UM, 0, -0.5 * sz + source_zpos * disc_um)
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
            center=mp.Vector3(0.5 * DISC_RADIUS_UM, 0, -0.5 * sz + 0.5 * disc_um),
            size=mp.Vector3(DISC_RADIUS_UM, mp.inf, disc_um),
        )
    ]

    sim = mp.Simulation(
        resolution=RESOLUTION_UM,
        cell_size=cell_size,
        dimensions=mp.CYLINDRICAL,
        m=-1,
        boundary_layers=boundary_layers,
        sources=sources,
        geometry=geometry,
    )

    # flux monitor
    flux_mon = sim.add_flux(
        frequency,
        0,
        1,
        mp.FluxRegion(
            center=mp.Vector3(0.5 * cell_r_um, 0, 0.5 * sz - PML_UM),
            size=mp.Vector3(cell_r_um, 0, 0),
        ),
        mp.FluxRegion(
            center=mp.Vector3(
                cell_r_um, 0, 0.5 * sz - PML_UM - 0.5 * (AIR_UM + disc_um)
            ),
            size=mp.Vector3(0, 0, AIR_UM + disc_um),
        ),
    )

    # near-field monitor
    n2f_mon = sim.add_near2far(
        frequency,
        0,
        1,
        mp.FluxRegion(
            center=mp.Vector3(0.5 * cell_r_um, 0, 0.5 * sz - PML_UM),
            size=mp.Vector3(cell_r_um, 0, 0),
        ),
        mp.FluxRegion(
            center=mp.Vector3(
                cell_r_um, 0, 0.5 * sz - PML_UM - 0.5 * (AIR_UM + disc_um)
            ),
            size=mp.Vector3(0, 0, AIR_UM + disc_um),
        ),
    )

    fig, ax = plt.subplots()
    sim.plot2D(ax=ax)
    if mp.am_master():
        fig.savefig("disc_simulation_layout.png", dpi=150, bbox_inches="tight")

    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            FIELD_DECAY_PERIOD,
            src_cmpt,
            src_pt,
            FIELD_DECAY_TOL,
        ),
    )

    flux_near = mp.get_fluxes(flux_mon)[0]

    radial_flux = radiation_pattern(sim, n2f_mon)
    radial_flux_scaled = FARFIELD_RADIUS_UM * FARFIELD_RADIUS_UM * radial_flux
    plot_radiation_pattern_polar(radial_flux_scaled)
    plot_radiation_pattern_3d(radial_flux_scaled)

    flux_far = (
        2
        * math.pi
        * FARFIELD_RADIUS_UM**2
        * np.trapezoid(radial_flux * np.sin(polar_rad), polar_rad)
    )

    return flux_near, flux_far


if __name__ == "__main__":
    disc_thickness_um = 0.7 * WAVELENGTH_UM / N_DISC
    dipole_height = 0.5

    flux_near, flux_far = disc_radiated_flux(disc_thickness_um, dipole_height)

    err = abs(flux_near - flux_far) / flux_near
    print(
        f"total_flux:, {flux_near:.5f} (near), {flux_far:.5f} (far), "
        f"{100 * err:.2f}% (error)"
    )
