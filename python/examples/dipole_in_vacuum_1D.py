"""Radiation pattern of a dipole in vacuum using Brillouin-zone integration."""

from enum import Enum
from typing import Tuple
import math

import matplotlib.pyplot as plt
import meep as mp
import numpy as np


Current = Enum("Current", "Jx Jy")

RESOLUTION_UM = 50
WAVELENGTH_UM = 1.0
NUM_POLAR = 50
NUM_AZIMUTHAL = 50

frequency = 1 / WAVELENGTH_UM


def planewave_in_vacuum(
    current_component: Current, kx: float, ky: float, kz: float
) -> Tuple[float, float, float]:
    """
    Returns the Poynting flux of a linearly polarized planewave in vacuum.

    Args:
        current_component: the component of the dipole.
        kx, ky, kz: the wavevector components of the planewave.

    Returns:
        The Poynting flux in x, y, z as a 3-tuple.
    """
    pml_um = 1.0
    air_um = 10.0
    size_z_um = pml_um + air_um + pml_um
    cell_size = mp.Vector3(0, 0, size_z_um)
    pml_layers = [mp.PML(thickness=pml_um, direction=mp.Z)]
    k_point = mp.Vector3(kx, ky, kz)

    if current_component.name == "Jx":
        src_cmpt = mp.Ex
    elif current_component.name == "Jy":
        src_cmpt = mp.Ey

    src_pt = mp.Vector3(0, 0, -0.5 * air_um)
    sources = [
        mp.Source(
            src=mp.GaussianSource(frequency, fwidth=0.1 * frequency),
            component=src_cmpt,
            center=src_pt,
            size=mp.Vector3(),
            amplitude=kz / frequency,
        )
    ]

    sim = mp.Simulation(
        resolution=RESOLUTION_UM,
        cell_size=cell_size,
        sources=sources,
        boundary_layers=pml_layers,
        k_point=k_point,
    )

    dft_flux_x = sim.add_flux(
        frequency,
        0,
        1,
        mp.FluxRegion(
            center=mp.Vector3(0, 0, 0.5 * air_um), size=mp.Vector3(), direction=mp.X
        ),
    )

    dft_flux_y = sim.add_flux(
        frequency,
        0,
        1,
        mp.FluxRegion(
            center=mp.Vector3(0, 0, 0.5 * air_um), size=mp.Vector3(), direction=mp.Y
        ),
    )

    dft_flux_z = sim.add_flux(
        frequency,
        0,
        1,
        mp.FluxRegion(
            center=mp.Vector3(0, 0, 0.5 * air_um), size=mp.Vector3(), direction=mp.Z
        ),
    )

    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            25, src_cmpt, mp.Vector3(0, 0, 0.5 * air_um), 1e-6
        )
    )

    flux_x = mp.get_fluxes(dft_flux_x)[0]
    flux_y = mp.get_fluxes(dft_flux_y)[0]
    flux_z = mp.get_fluxes(dft_flux_z)[0]

    return flux_x, flux_y, flux_z


def spherical_to_cartesian(polar_rad, azimuthal_rad) -> Tuple[float, float, float]:
    """Converts a point on the unit sphere from spherical to Cartesian coords.

    Args:
        polar_rad: polar angle of the point. 0° is +z.
        azimuthal_rad: azimuthal angle of the point. 0° is +x.

    Returns:
        The x,y,z coordinates of the point as a 3-tuple.
    """
    x = np.sin(polar_rad) * np.cos(azimuthal_rad)
    y = np.sin(polar_rad) * np.sin(azimuthal_rad)
    z = np.cos(polar_rad)

    return (x, y, z)


if __name__ == "__main__":
    # Radial flux is defined on the surface of a hemisphere.
    polar_rad = np.linspace(0, 0.5 * np.pi, NUM_POLAR)
    azimuthal_rad = np.linspace(0, 2 * np.pi, NUM_AZIMUTHAL)
    radial_flux = np.zeros((NUM_POLAR, NUM_AZIMUTHAL))

    for i in range(NUM_POLAR):
        for j in range(NUM_AZIMUTHAL):
            rx, ry, rz = spherical_to_cartesian(polar_rad[i], azimuthal_rad[j])
            kx = frequency * rx
            ky = frequency * ry
            kz = frequency * rz

            # Skip wavevectors which are close to the light cone
            # due to poor absorption by PML.
            if np.sqrt(kx**2 + ky**2) > (0.95 * frequency):
                continue

            flux_x, flux_y, flux_z = planewave_in_vacuum(Current.Jx, kx, ky, kz)
            radial_flux[i, j] = rx * flux_x + ry * flux_y + rz * flux_z

    if mp.am_master():
        np.savez(
            "dipole_radiation_pattern.npz",
            RESOLUTION_UM=RESOLUTION_UM,
            WAVELENGTH_UM=WAVELENGTH_UM,
            NUM_POLAR=NUM_POLAR,
            NUM_AZIMUTHAL=NUM_AZIMUTHAL,
            polar_rad=polar_rad,
            azimuthal_rad=azimuthal_rad,
            radial_flux=radial_flux,
        )

    radiation_pattern = radial_flux[:, 0] / np.max(radial_flux[:, 0])

    # in plot of radiation pattern, omit polar angles with zero radial flux
    nonzero_mask = np.nonzero(radiation_pattern)
    polar_rad = polar_rad[nonzero_mask]
    radiation_pattern = radiation_pattern[nonzero_mask]

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(6, 6))
    ax.plot(polar_rad, radiation_pattern, "b-", label="Meep")
    ax.plot(polar_rad, np.cos(polar_rad) ** 2, "r--", label="cos$^2$(θ)")
    ax.set_theta_direction(-1)
    ax.set_theta_offset(0.5 * math.pi)
    ax.set_thetalim(0, 0.5 * math.pi)
    ax.set_rmax(1)
    ax.set_rticks([0, 0.5, 1])
    ax.grid(True)
    ax.set_rlabel_position(22)
    ax.legend()
    ax.set_title(r"radiation pattern for J$_x$ current source and $\phi = 0$")
    fig.savefig("dipole_radiation_pattern.png", dpi=150, bbox_inches="tight")
