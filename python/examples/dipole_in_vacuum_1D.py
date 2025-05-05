"""Radiation pattern of a dipole in vacuum obtained using a 1D calculation."""

from enum import Enum
from typing import Dict, Tuple
import math

import matplotlib.pyplot as plt
import meep as mp
import numpy as np


Current = Enum("Current", "Jx Jy Kx Ky")

RESOLUTION_UM = 25
WAVELENGTH_UM = 1.0
NUM_POLAR = 50
NUM_AZIMUTHAL = 50
DEBUG_OUTPUT = True

frequency = 1 / WAVELENGTH_UM


def planewave_in_vacuum(
    dipole_component: Current, kx: float, ky: float
) -> Tuple[complex, complex, complex, complex]:
    """
    Returns the near fields of a dipole in vacuum.

    Args:
        dipole_component: the component of the dipole.
        kx, ky: the wavevector components.

    Returns:
        The surface-tangential electric and magnetic near fields as a 4-tuple.
    """
    pml_um = 2.0
    air_um = 10.0
    size_z_um = pml_um + air_um + pml_um
    cell_size = mp.Vector3(0, 0, size_z_um)
    pml_layers = [mp.PML(thickness=pml_um, direction=mp.Z)]

    if dipole_component.name == "Jx":
        src_cmpt = mp.Ex
    elif dipole_component.name == "Jy":
        src_cmpt = mp.Ey
    elif dipole_component.name == "Kx":
        src_cmpt = mp.Hx
    elif dipole_component.name == "Ky":
        src_cmpt = mp.Hy

    sources = [
        mp.Source(
            src=mp.GaussianSource(frequency, fwidth=0.1 * frequency),
            component=src_cmpt,
            center=mp.Vector3(0, 0, -0.5 * air_um),
        )
    ]

    sim = mp.Simulation(
        resolution=RESOLUTION_UM,
        force_complex_fields=True,
        cell_size=cell_size,
        sources=sources,
        boundary_layers=pml_layers,
        k_point=mp.Vector3(kx, ky, 0),
    )

    dft_fields = sim.add_dft_fields(
        [mp.Ex, mp.Ey, mp.Hx, mp.Hy],
        frequency,
        0,
        1,
        center=mp.Vector3(0, 0, 0.5 * air_um),
        size=mp.Vector3(),
    )

    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            15, src_cmpt, mp.Vector3(0, 0, 0.5423), 1e-6
        )
    )

    ex_dft = sim.get_dft_array(dft_fields, mp.Ex, 0)
    ey_dft = sim.get_dft_array(dft_fields, mp.Ey, 0)
    hx_dft = sim.get_dft_array(dft_fields, mp.Hx, 0)
    hy_dft = sim.get_dft_array(dft_fields, mp.Hy, 0)

    return ex_dft, ey_dft, hx_dft, hy_dft


def equivalent_currents(
    ex_dft: complex, ey_dft: complex, hx_dft: complex, hy_dft: complex
) -> Tuple[Tuple[complex, complex], Tuple[complex, complex]]:
    """Computes the equivalent electric and magnetic currents on a surface.

    Args:
        ex_dft, ey_dft, hx_dft, hy_dft: the surface tangential DFT fields.

    Returns:
        A 2-tuple of the electric and magnetic sheet currents in x and y.
    """

    electric_current = (-hy_dft, hx_dft)
    magnetic_current = (ey_dft, -ex_dft)

    return (electric_current, magnetic_current)


def field_amplitudes_from_sheet_current(
    kx: float,
    kz: float,
    current_amplitude: complex,
    current_component: Current,
) -> Tuple[complex, complex, complex]:
    """Computes the S- or P-polarized field amplitudes from a sheet current.

    Assumes ky=0 such that wavevector is in xz plane.

    Args:
        kx, kz: wavevector of the outgoing planewave in the x,z direction. Units of 2π.
        current_amplitude: amplitude of the sheet current.
        current_component: component of the sheet current.

    Returns:
        A 3-tuple of the per-polarization electric and magnetic far fields.
    """
    if current_component.name == "Jx":
        # Jx --> (Ex, Hy, Ez) [P pol.]
        ex0 = kz * current_amplitude / (2 * frequency)
        hy0 = frequency * ex0 / kz
    elif current_component.name == "Jy":
        # Jy --> (Hx, Ey, Hz) [S pol.]
        ey0 = frequency * current_amplitude / (2 * kz)
        hx0 = -kz * ey0 / frequency
    elif current_component.name == "Kx":
        # Kx --> (Hx, Ey, Hz) [S pol.]
        hx0 = kz * current_amplitude / (2 * frequency)
        ey0 = -frequency * hx0 / kz
    elif current_component.name == "Ky":
        # Ky --> (Ex, Hy, Ez) [P pol.]
        hy0 = -frequency * current_amplitude / (2 * kz)
        ex0 = kz * hy0 / frequency

    if current_component.name == "Jx" or current_component.name == "Ky":
        ez0 = -kx * hy0 / frequency
        return (ex0, hy0, ez0)
    elif current_component.name == "Jy" or current_component.name == "Kx":
        hz0 = -kx * ey0 / frequency
        return (hx0, ey0, hz0)


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


def farfields_at_k_point(kx: float, ky: float) -> Dict[str, complex]:
    """Computes the farfields from a linearly polarized dipole in vacuum.

    Args:
        kx, ky: in-plane wavevector components.

    Returns:
        A dictionary with the six field components as keys and their farfields
        as values.
    """
    dipole_component = Current.Jx
    current_components = [Current.Jx, Current.Jy, Current.Kx, Current.Ky]
    farfield_components = ["Ex", "Ey", "Ez", "Hx", "Hy", "Hz"]

    kz = (frequency**2 - kx**2 - ky**2) ** 0.5

    ex_dft, ey_dft, hx_dft, hy_dft = planewave_in_vacuum(dipole_component, kx, ky)

    # Rotation angle around z axis to force ky=0. 0 radians is +x.
    rotation_rad = -math.atan2(ky, kx)

    if DEBUG_OUTPUT:
        print(
            f"rotation angle:, ({kx:.6f}, {ky:.6f}, {kz:.6f}), "
            f"{math.degrees(rotation_rad):.2f}°"
        )

    rotation_matrix = np.array(
        [
            [np.cos(rotation_rad), -np.sin(rotation_rad)],
            [np.sin(rotation_rad), np.cos(rotation_rad)],
        ]
    )

    k_rotated = rotation_matrix @ np.transpose(np.array([kx, ky]))
    if k_rotated[1] > 1e-10:
        raise ValueError(f"rotated ky is nonzero: {k_rotated[1]}")

    electric_current, magnetic_current = equivalent_currents(
        ex_dft, ey_dft, hx_dft, hy_dft
    )

    electric_current_rotated = rotation_matrix @ np.transpose(
        np.array([electric_current[0], electric_current[1]])
    )

    magnetic_current_rotated = rotation_matrix @ np.transpose(
        np.array([magnetic_current[0], magnetic_current[1]])
    )

    current_amplitudes = [
        electric_current_rotated[0],
        electric_current_rotated[1],
        magnetic_current_rotated[0],
        magnetic_current_rotated[1],
    ]

    # Obtain the far fields at a single point on the quarter circle (φ = 0)
    # given a sheet current with linear in-plane polarization.
    # Note: the phase is omitted because it is a constant (ωR).

    farfields = {}
    for component in farfield_components:
        farfields[component] = 0

    for current_component, current_amplitude in zip(
        current_components, current_amplitudes
    ):
        if abs(current_amplitude) == 0:
            continue

        farfield_amplitudes_pol = field_amplitudes_from_sheet_current(
            k_rotated[0],
            kz,
            current_amplitude,
            current_component,
        )

        farfield_pol = np.array(farfield_amplitudes_pol)

        if current_component.name == "Jx" or current_component.name == "Ky":
            # P polarization
            farfields["Ex"] += farfield_pol[0]
            farfields["Hy"] += farfield_pol[1]
            farfields["Ez"] += farfield_pol[2]
        elif current_component.name == "Jy" or current_component.name == "Kx":
            # S polarization
            farfields["Hx"] += farfield_pol[0]
            farfields["Ey"] += farfield_pol[1]
            farfields["Hz"] += farfield_pol[2]

    rotated_farfields = {}
    for component in farfield_components:
        rotated_farfields[component] = 0

    inverse_rotation_matrix = np.transpose(rotation_matrix)

    rotated_farfields["Ex"] = (
        inverse_rotation_matrix[0, 0] * farfields["Ex"]
        + inverse_rotation_matrix[0, 1] * farfields["Ey"]
    )
    rotated_farfields["Ey"] = (
        inverse_rotation_matrix[1, 0] * farfields["Ex"]
        + inverse_rotation_matrix[1, 1] * farfields["Ey"]
    )
    rotated_farfields["Ez"] = farfields["Ez"]

    rotated_farfields["Hx"] = (
        inverse_rotation_matrix[0, 0] * farfields["Hx"]
        + inverse_rotation_matrix[0, 1] * farfields["Hy"]
    )
    rotated_farfields["Hy"] = (
        inverse_rotation_matrix[1, 0] * farfields["Hx"]
        + inverse_rotation_matrix[1, 1] * farfields["Hy"]
    )
    rotated_farfields["Hz"] = farfields["Hz"]

    return rotated_farfields


if __name__ == "__main__":
    # Far fields are defined on the surface of a hemisphere.
    polar_rad = np.linspace(0, 0.5 * np.pi, NUM_POLAR)
    azimuthal_rad = np.linspace(0, 2 * np.pi, NUM_AZIMUTHAL)
    radial_flux = np.zeros((NUM_POLAR, NUM_AZIMUTHAL))

    for i in range(NUM_POLAR):
        for j in [0]:  # range(NUM_AZIMUTHAL):
            kx, ky, kz = spherical_to_cartesian(polar_rad[i], azimuthal_rad[j])
            kx *= frequency
            ky *= frequency
            kz *= frequency

            # Skip wavevectors which are outside the light cone.
            # maximum polar angle: math.degrees(math.asin(0.85)) = 58.2°
            if np.sqrt(kx**2 + ky**2) > (0.85 * frequency):
                if DEBUG_OUTPUT:
                    print(
                        f"skipping:, {math.degrees(polar_rad[i]):2.2f}°, "
                        f"({kx:.6f}, {ky:.6f}, {kz:.6f})"
                    )
                continue

            farfields = farfields_at_k_point(kx, ky)

            # (Ex, Hy, Ez) are nonzero
            flux_x = np.real(
                np.conj(farfields["Ey"]) * farfields["Hz"]
                - np.conj(farfields["Ez"]) * farfields["Hy"]
            )
            flux_y = np.real(
                np.conj(farfields["Ez"]) * farfields["Hx"]
                - np.conj(farfields["Ex"]) * farfields["Hz"]
            )
            flux_z = np.real(
                np.conj(farfields["Ex"]) * farfields["Hy"]
                - np.conj(farfields["Ey"]) * farfields["Hx"]
            )

            if DEBUG_OUTPUT:
                print(
                    f"directional_flux:, {math.degrees(polar_rad[i]):2.2f}°, "
                    f"{flux_x:.6g}, {flux_y:.6g}, {flux_z:.6g}"
                )

            rx, ry, rz = spherical_to_cartesian(polar_rad[i], azimuthal_rad[j])
            radial_flux[i, j] = rx * flux_x + ry * flux_y + rz * flux_z

            if DEBUG_OUTPUT:
                print(
                    f"radial_flux:, {math.degrees(polar_rad[i]):2.2f}°, "
                    f"({rx:.6f}, {ry:.6f}, {rz:.6f}), {radial_flux[i, j]:.6g}"
                )

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
    max_idx = np.argmin(radial_flux[:, 0])
    polar_rad = polar_rad[:max_idx]
    radiation_pattern = radiation_pattern[:max_idx]

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(6, 6))
    ax.plot(polar_rad, radiation_pattern, "b-", label="Meep")
    ax.plot(polar_rad, np.cos(polar_rad) ** 2, "r--", label="cos$^2$(θ)")
    ax.plot(polar_rad, np.cos(polar_rad) ** 4, "g--", label="cos$^4$(θ)")
    ax.set_theta_direction(-1)
    ax.set_theta_offset(0.5 * math.pi)
    ax.set_thetalim(0, 0.5 * math.pi)
    ax.set_rmax(1)
    ax.set_rticks([0, 0.5, 1])
    ax.grid(True)
    ax.set_rlabel_position(22)
    ax.legend()
    ax.set_title(r"radiation pattern for J$_x$ current source and $\phi = 0$")
    fig.savefig("radiation_pattern_1D.png", dpi=150, bbox_inches="tight")
