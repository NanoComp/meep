"""Radiation pattern of a dipole antenna above a ground plane using Brillouin-zone integration.

Tutorial reference:

https://meep.readthedocs.io/en/latest/Python_Tutorials/Near_to_Far_Field_Spectra/#antenna-above-a-perfect-electric-conductor-ground-plane
"""

import math
from typing import Tuple

import matplotlib.pyplot as plt
import meep as mp
import numpy as np


RESOLUTION_UM = 50
PML_UM = 1.0
BULK_UM = 10.0
N_BACKGROUND = 1.2
ANTENNA_HEIGHT_UM = 1.25
WAVELENGTH_UM = 0.65
NUM_POLAR = 50
FIELD_DECAY_THRESHOLD = 1e-6
FIELD_DECAY_PERIOD = 25

frequency = 1 / WAVELENGTH_UM


def planewave_ground_plane(kx: float, ky: float, kz: float) -> float:
    """
    Returns the radiated flux of a linearly polarized planewave above ground plane.

    Args:
        kx, ky, kz: the wavevector components of the planewave.

    Returns:
        The Poynting flux in z.
    """
    size_z_um = BULK_UM + PML_UM
    pml_layers = [mp.PML(thickness=PML_UM, direction=mp.Z, side=mp.High)]
    src_pt = mp.Vector3(0, 0, -0.5 * size_z_um + ANTENNA_HEIGHT_UM)

    cell_size = mp.Vector3(0, 0, size_z_um)

    src_cmpt = mp.Ey
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
        default_material=mp.Medium(index=N_BACKGROUND),
        cell_size=cell_size,
        sources=sources,
        boundary_layers=pml_layers,
        k_point=mp.Vector3(kx, ky, kz),
    )

    sim.set_boundary(mp.Low, mp.Z, mp.Metallic)
    sim.set_boundary(mp.High, mp.Z, mp.Metallic)

    mon_pt = mp.Vector3(0, 0, 0.5 * size_z_um - PML_UM)
    dft_flux_z = sim.add_flux(
        frequency,
        0,
        1,
        mp.FluxRegion(center=mon_pt, size=mp.Vector3(), direction=mp.Z),
    )

    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            FIELD_DECAY_PERIOD, src_cmpt, mon_pt, FIELD_DECAY_THRESHOLD
        )
    )

    flux_z = mp.get_fluxes(dft_flux_z)[0]

    return flux_z


def planewave_wavevector(polar_rad: float) -> Tuple[float, float, float]:
    """Returns the wavevector of the outgoing planewave in free space.

    Args:
        polar_rad: polar angle of the wavevector for φ = 0 (xz plane). 0 is
            the +z direction.

    Returns:
        The wavevector of the outgoing planewave.
    """
    kx = N_BACKGROUND * frequency * np.sin(polar_rad)
    ky = 0
    kz = N_BACKGROUND * frequency * np.cos(polar_rad)

    return kx, ky, kz


if __name__ == "__main__":
    polar_rad = np.linspace(0, 0.5 * np.pi, NUM_POLAR)
    radial_flux_meep = np.zeros(NUM_POLAR)

    for i in range(NUM_POLAR):
        kx, ky, kz = planewave_wavevector(polar_rad[i])

        # Skip wavevectors which are close to the light cone
        # due to poor absorption by PML (i.e. glancing-angle waves).
        if kx > (0.95 * N_BACKGROUND * frequency):
            continue

        flux_z = planewave_ground_plane(kx, ky, kz)
        radial_flux_meep[i] = np.cos(polar_rad[i]) * flux_z

    # The radiation pattern of a two-element antenna
    # array is equivalent to the radiation pattern of
    # a single antenna multiplied by its array factor
    # as derived in Section 6.2 "Two-Element Array" of
    # Antenna Theory: Analysis and Design, Fourth Edition
    # (2016) by C.A. Balanis.
    k_free_space = 2 * math.pi / (WAVELENGTH_UM / N_BACKGROUND)
    radial_flux_analytic = (
        np.sin(k_free_space * ANTENNA_HEIGHT_UM * np.cos(polar_rad)) ** 2
    )

    # Find the smallest angle at which the radiation pattern is a maximum using
    # the formula for the radiation pattern. This angle is used to compute the
    # radial flux which is then used for normalization of the radiation pattern.
    polar_rad_max = math.acos(9 * math.pi / (2 * k_free_space * ANTENNA_HEIGHT_UM))
    kx, ky, kz = planewave_wavevector(polar_rad_max)
    flux_z = planewave_ground_plane(kx, ky, kz)
    radial_flux_meep_max = np.cos(polar_rad_max) * flux_z

    radial_flux_meep = radial_flux_meep / radial_flux_meep_max

    rel_err = np.linalg.norm(radial_flux_meep - radial_flux_analytic) / np.linalg.norm(
        radial_flux_analytic
    )
    print(f"error:, {rel_err:.6f}")

    fig, ax = plt.subplots()
    ax.plot(np.degrees(polar_rad), radial_flux_meep, "b-", label="Meep")
    ax.plot(np.degrees(polar_rad), radial_flux_analytic, "r-", label="theory")
    ax.set_xlabel("angle (degrees)")
    ax.set_ylabel("radiation pattern")
    ax.set_title("antenna above PEC ground plane\n using Brillouin-zone integration")
    ax.set_xlim(0, 90.0)
    ax.set_ylim(0, 1.0)
    ax.legend()
    if mp.am_master():
        fig.savefig("antenna_radiation_pattern_1D.png", dpi=150, bbox_inches="tight")
