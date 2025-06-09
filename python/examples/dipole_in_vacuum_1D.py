"""Radiation pattern of a dipole in vacuum using Brillouin-zone integration."""

import argparse
from enum import Enum
from typing import Tuple
import math

import matplotlib.pyplot as plt
import meep as mp
import numpy as np


RESOLUTION_UM = 50
WAVELENGTH_UM = 1.0
NUM_POLAR = 30
NUM_AZIMUTH = 50
FIELD_DECAY_THRESHOLD = 1e-6

frequency = 1 / WAVELENGTH_UM


def planewave_in_vacuum(dipole_pol: str, kx: float, ky: float, kz: float) -> float:
    """
    Returns the Poynting flux of a linearly polarized planewave in vacuum.

    Args:
        dipole_pol: the polarization of the electric dipole. Either x or y.
        kx, ky, kz: the wavevector components of the planewave.

    Returns:
        The Poynting flux in z.
    """
    pml_um = 1.0
    air_um = 10.0
    size_z_um = pml_um + air_um + pml_um
    cell_size = mp.Vector3(0, 0, size_z_um)
    pml_layers = [mp.PML(thickness=pml_um, direction=mp.Z)]
    k_point = mp.Vector3(kx, ky, kz)

    if dipole_pol == "x":
        src_cmpt = mp.Ex
    elif dipole_pol == "y":
        src_cmpt = mp.Ey

    src_pt = mp.Vector3(0, 0, -0.5 * air_um)
    sources = [
        mp.Source(
            src=mp.GaussianSource(frequency, fwidth=0.1 * frequency),
            component=src_cmpt,
            center=src_pt,
            size=mp.Vector3(),
        )
    ]

    sim = mp.Simulation(
        resolution=RESOLUTION_UM,
        cell_size=cell_size,
        sources=sources,
        boundary_layers=pml_layers,
        k_point=k_point,
    )

    mon_pt = mp.Vector3(0, 0, 0.5 * air_um)
    dft_flux_z = sim.add_flux(
        frequency,
        0,
        1,
        mp.FluxRegion(center=mon_pt, size=mp.Vector3(), direction=mp.Z),
    )

    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            25, src_cmpt, mon_pt, FIELD_DECAY_THRESHOLD
        )
    )

    flux_z = mp.get_fluxes(dft_flux_z)[0]

    return flux_z


def spherical_to_cartesian(polar_rad, azimuth_rad) -> Tuple[float, float, float]:
    """Converts a point on the unit sphere from spherical to Cartesian coords.

    Args:
        polar_rad: polar angle of the point. 0° is +z.
        azimuth_rad: azimuthal angle of the point. 0° is +x.

    Returns:
        The x,y,z coordinates of the point as a 3-tuple.
    """
    x = np.sin(polar_rad) * np.cos(azimuth_rad)
    y = np.sin(polar_rad) * np.sin(azimuth_rad)
    z = np.cos(polar_rad)

    return x, y, z


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "dipole_pol",
        type=str,
        choices=["x", "y"],
        help="polarization of the electric dipole (x or y)",
    )
    args = parser.parse_args()

    # Radial flux is defined on the surface of a hemisphere.
    polar_rad = np.linspace(0, 0.5 * np.pi, NUM_POLAR)
    azimuth_rad = np.linspace(0, 2 * np.pi, NUM_AZIMUTH)
    radial_flux = np.zeros((NUM_POLAR, NUM_AZIMUTH))

    for i in range(NUM_POLAR):
        for j in range(NUM_AZIMUTH):
            rx, ry, rz = spherical_to_cartesian(polar_rad[i], azimuth_rad[j])
            # Specify the components of the wavevector of the outgoing
            # planewave in vacuum.
            kx = frequency * rx
            ky = frequency * ry
            kz = frequency * rz

            # Skip wavevectors which are close to the light cone
            # due to poor absorption by PML (i.e. glancing-angle waves).
            if np.sqrt(kx**2 + ky**2) > (0.95 * frequency):
                continue

            flux_z = planewave_in_vacuum(args.dipole_pol, kx, ky, kz)
            radial_flux[i, j] = rz * flux_z

    if mp.am_master():
        np.savez(
            "dipole_radiation_pattern.npz",
            NUM_AZIMUTH=NUM_AZIMUTH,
            NUM_POLAR=NUM_POLAR,
            RESOLUTION_UM=RESOLUTION_UM,
            WAVELENGTH_UM=WAVELENGTH_UM,
            azimuth_rad=azimuth_rad,
            dipole_pol=args.dipole_pol,
            polar_rad=polar_rad,
            radial_flux=radial_flux,
        )
