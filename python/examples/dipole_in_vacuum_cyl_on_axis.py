"""Radiation pattern of an on-axis dipole in cylindrical coordinates.

Tutorial Reference:

https://meep.readthedocs.io/en/latest/Python_Tutorials/Near_to_Far_Field_Spectra/#radiation-pattern-of-an-antenna-in-cylindrical-coordinates
"""

import argparse
import cmath
import math
from typing import Tuple

import matplotlib.pyplot as plt
import meep as mp
import numpy as np


RESOLUTION_UM = 50
WAVELENGTH_UM = 1.0
PML_UM = 1.0
FARFIELD_RADIUS_UM = 1e6 * WAVELENGTH_UM
NUM_FARFIELD_PTS = 50
AZIMUTHAL_RAD = 0
POWER_DECAY_THRESHOLD = 1e-4

frequency = 1 / WAVELENGTH_UM
polar_rad = np.linspace(0, 0.5 * math.pi, NUM_FARFIELD_PTS)


def plot_radiation_pattern(dipole_pol: str, radial_flux: np.ndarray):
    """Plots the radiation pattern in polar coordinates.

    The angles increase clockwise with zero in the +z direction (the "pole")
    and π/2 in the +r direction (the "equator").

    Args:
        dipole_pol: the polarization the electric dipole. Either x or z.
        radial_flux: the radial flux in polar coordinates.
    """
    normalized_radial_flux = radial_flux / np.max(radial_flux)

    if dipole_pol == "x":
        dipole_radial_flux = np.square(np.cos(polar_rad))
        dipole_radial_flux_label = r"$\cos^2θ$"
        dipole_name = "$E_x$"
    else:
        dipole_radial_flux = np.square(np.sin(polar_rad))
        dipole_radial_flux_label = r"$\sin^2θ$"
        dipole_name = "$E_z$"

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
    ax.set_title(f"radiation pattern (φ = 0) of an on-axis {dipole_name} dipole")

    if mp.am_master():
        fig.savefig(
            "dipole_radiation_pattern_on_axis.png",
            dpi=150,
            bbox_inches="tight",
        )

    relative_error = np.linalg.norm(
        normalized_radial_flux - dipole_radial_flux
    ) / np.linalg.norm(dipole_radial_flux)
    print(f"relative error in radiation pattern:, {relative_error}")


def radiation_pattern(e_field: np.ndarray, h_field: np.ndarray) -> np.ndarray:
    """Computes the radiation pattern from the far fields.

    Args:
        e_field, h_field: the electric (Er, Ep, Ez) and magnetic (Hr, Hp, Hz)
          far fields, respectively.

    Returns:
        The radial Poynting flux as a 1D array. One element for each point on
        the circumference of a quarter circle with angular range of
        [0, π/2] rad. 0 radians is the +z direction (the "pole") and π/2 is
        the +r direction (the "equator").
    """
    flux_x = np.real(
        e_field[:, 1] * np.conj(h_field[:, 2]) - e_field[:, 2] * np.conj(h_field[:, 1])
    )
    flux_z = np.real(
        e_field[:, 0] * np.conj(h_field[:, 1]) - e_field[:, 1] * np.conj(h_field[:, 0])
    )
    flux_r = np.sqrt(np.square(flux_x) + np.square(flux_z))

    return flux_r


def get_farfields(
    sim: mp.Simulation, n2f_mon: mp.DftNear2Far
) -> Tuple[np.ndarray, np.ndarray]:
    """Computes the far fields from the near fields for φ = 0 (rz plane).

    Args:
        sim: a `Simulation` object.
        n2f_mon: a `DftNear2Far` object returned by `Simulation.add_near2far`.

    Returns:
        The electric (Er, Ep, Ez) and magnetic (Hr, Hp, Hz) far fields. One row
        for each point on the circumference of a quarter circle with angular
        range of [0, π/2] rad. Each row has six columns for the fields.
        0 radians is the +z direction (the "pole") and π/2 is the +r direction
        (the "equator").
    """
    e_field = np.zeros((NUM_FARFIELD_PTS, 3), dtype=np.complex128)
    h_field = np.zeros((NUM_FARFIELD_PTS, 3), dtype=np.complex128)
    for n in range(NUM_FARFIELD_PTS):
        far_field = sim.get_farfield(
            n2f_mon,
            mp.Vector3(
                FARFIELD_RADIUS_UM * math.sin(polar_rad[n]),
                0,
                FARFIELD_RADIUS_UM * math.cos(polar_rad[n]),
            ),
        )
        e_field[n, :] = [far_field[j] for j in range(3)]
        h_field[n, :] = [far_field[j + 3] for j in range(3)]

    return e_field, h_field


def dipole_in_vacuum(dipole_pol: str, m: int) -> Tuple[np.ndarray, np.ndarray]:
    """Computes the far fields of an on-axis point source.

    Args:
        dipole_pol: the polarization of the electric dipole. Either x or z.
        m: angular φ dependence of the fields exp(imφ).

    Returns:
        A 2-tuple containing the electric and magnetic far fields as 1D arrays.
    """
    sr = 2.0
    sz = 4.0
    cell_size = mp.Vector3(sr + PML_UM, 0, sz + 2 * PML_UM)

    boundary_layers = [mp.PML(thickness=PML_UM)]

    if dipole_pol == "x":
        # An Er source at r = 0 needs to be slighty offset due to a bug.
        # https://github.com/NanoComp/meep/issues/2704
        dipole_pos_r = 1.5 / RESOLUTION_UM
        sources = [
            mp.Source(
                src=mp.GaussianSource(frequency, fwidth=0.1 * frequency),
                component=mp.Er,
                center=mp.Vector3(dipole_pos_r, 0, 0),
            ),
            mp.Source(
                src=mp.GaussianSource(frequency, fwidth=0.1 * frequency),
                component=mp.Ep,
                center=mp.Vector3(dipole_pos_r, 0, 0),
                amplitude=1j if m == 1 else -1j,
            ),
        ]
    else:
        dipole_pos_r = 0
        sources = [
            mp.Source(
                src=mp.GaussianSource(frequency, fwidth=0.1 * frequency),
                component=mp.Ez,
                center=mp.Vector3(dipole_pos_r, 0, 0),
            )
        ]

    sim = mp.Simulation(
        resolution=RESOLUTION_UM,
        cell_size=cell_size,
        dimensions=mp.CYLINDRICAL,
        m=m,
        boundary_layers=boundary_layers,
        sources=sources,
        force_complex_fields=True,
    )

    nearfields_monitor = sim.add_near2far(
        frequency,
        0,
        1,
        mp.FluxRegion(
            center=mp.Vector3(0.5 * sr, 0, 0.5 * sz), size=mp.Vector3(sr, 0, 0)
        ),
        mp.FluxRegion(center=mp.Vector3(sr, 0, 0), size=mp.Vector3(0, 0, sz)),
        mp.FluxRegion(
            center=mp.Vector3(0.5 * sr, 0, -0.5 * sz),
            size=mp.Vector3(sr, 0, 0),
            weight=-1.0,
        ),
    )

    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            20.0,
            mp.Er if dipole_pol == "x" else mp.Ez,
            mp.Vector3(dipole_pos_r, 0, 0),
            1e-6,
        )
    )

    e_field, h_field = get_farfields(sim, nearfields_monitor)

    return e_field, h_field


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "dipole_pol",
        type=str,
        choices=["x", "z"],
        help="polarization of the electric dipole (x or z)",
    )
    args = parser.parse_args()

    e_field_total = np.zeros((NUM_FARFIELD_PTS, 3), dtype=np.complex128)
    h_field_total = np.zeros((NUM_FARFIELD_PTS, 3), dtype=np.complex128)

    if args.dipole_pol == "x":
        # An x-polarized dipole can be formed from the superposition
        # of left- and right-circularly polarized dipoles.

        e_field, h_field = dipole_in_vacuum("x", +1)
        e_field_total += 0.5 * e_field * cmath.exp(1j * AZIMUTHAL_RAD)
        h_field_total += 0.5 * h_field * cmath.exp(1j * AZIMUTHAL_RAD)

        e_field, h_field = dipole_in_vacuum("x", -1)
        e_field_total += 0.5 * e_field * cmath.exp(-1j * AZIMUTHAL_RAD)
        h_field_total += 0.5 * h_field * cmath.exp(-1j * AZIMUTHAL_RAD)
    else:
        e_field, h_field = dipole_in_vacuum("z", 0)
        e_field_total += e_field
        h_field_total += h_field

    dipole_radiation_pattern = radiation_pattern(e_field_total, h_field_total)
    dipole_radiation_pattern_scaled = dipole_radiation_pattern * FARFIELD_RADIUS_UM**2
    plot_radiation_pattern(args.dipole_pol, dipole_radiation_pattern_scaled)

    if mp.am_master():
        np.savez(
            "dipole_farfields_on_axis.npz",
            AZIMUTHAL_RAD=AZIMUTHAL_RAD,
            FARFIELD_RADIUS_UM=FARFIELD_RADIUS_UM,
            PML_UM=PML_UM,
            POWER_DECAY_THRESHOLD=POWER_DECAY_THRESHOLD,
            RESOLUTION_UM=RESOLUTION_UM,
            WAVELENGTH_UM=WAVELENGTH_UM,
            dipole_pol=args.dipole_pol,
            dipole_radiation_pattern=dipole_radiation_pattern,
            e_field_total=e_field_total,
            h_field_total=h_field_total,
            polar_rad=polar_rad,
        )
