"""Radiation pattern of a dipole in vacuum obtained using a 1D calculation."""

from enum import Enum
from typing import Tuple
import math
import warnings

import meep as mp
import numpy as np


Current = Enum("Current", "Jx Jy Kx Ky")

RESOLUTION_UM = 25
WAVELENGTH_UM = 1.0
FARFIELD_RADIUS_UM = 1e6 * WAVELENGTH_UM
NUM_K = 25
NUM_POLAR = 80
NUM_AZIMUTHAL = 80
DEBUG_OUTPUT = True

frequency = 1 / WAVELENGTH_UM


def dipole_in_vacuum(
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
    pml_um = 1.0
    air_um = 10.0
    size_z_um = pml_um + air_um + pml_um
    cell_size = mp.Vector3(0, 0, size_z_um)
    pml_layers = [mp.PML(thickness=pml_um)]

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
            10, src_cmpt, mp.Vector3(0, 0, 0.5423), 1e-6
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


def farfield_amplitudes(
    kx: float,
    kz: float,
    rx: float,
    rz: float,
    current_amplitude: complex,
    current_component: Current,
) -> Tuple[complex, complex, complex]:
    """Computes the S- or P-polarized far-field amplitudes from a sheet current.                             
                                                                                                             
    Assumes ky=0 such that wavevector is in xz plane.                                                        
                                                                                                             
    Args:                                                                                                    
        kx, kz: wavevector of the outgoing planewave in the x,z direction. Units of 2π.                      
        rx, rz: x,z coordinate of a point on the surface of a hemicircle.                                    
        current_amplitude: amplitude of the sheet current.                                                   
        current_component: component of the sheet current.                                                   
                                                                                                             
    Returns:                                                                                                 
        A 3-tuple of the per-polarization electric and magnetic far fields.                                  
    """
    if current_component.name == "Jx":
        # Jx --> (Ex, Hy, Ez) [P pol.]                                                                       
        ex0 = frequency * current_amplitude / (2 * kz)
    elif current_component.name == "Jy":
        # Jy --> (Hx, Ey, Hz) [S pol.]                                                                       
        ey0 = -frequency * current_amplitude / (2 * kz)
    elif current_component.name == "Kx":
        # Kx --> (Hx, Ey, Hz) [S pol.]                                                                       
        ey0 = -current_amplitude / 2
    elif current_component.name == "Ky":
        # Ky --> (Ex, Hy, Ez) [P pol.]                                                                       
        ex0 = current_amplitude / 2

    if current_component.name == "Jx" or current_component.name == "Ky":
        hy0 = frequency * ex0 / kz
        ez0 = -kx * hy0 / frequency
        return (ex0, hy0, ez0)
    elif current_component.name == "Jy" or current_component.name == "Kx":
        hx0 = -kz * ey0 / frequency
        hz0 = kx * ey0 / frequency
        return (hx0, ey0, hz0)


def spherical_to_cartesian(polar_rad, azimuthal_rad) -> Tuple[float, float, float]:
    """Converts a far point in spherical to Cartesian coordinates.                                           
                                                                                                             
    Args:                                                                                                    
        polar_rad: polar angle of the point.                                                                 
        azimuthal_rad: azimuthal angle of the point.                                                         
                                                                                                             
    Returns:                                                                                                 
        The x,y,z coordinates of the far point as a 3-tuple.                                                 
    """
    x = FARFIELD_RADIUS_UM * np.sin(polar_rad) * np.cos(azimuthal_rad)
    y = FARFIELD_RADIUS_UM * np.sin(polar_rad) * np.sin(azimuthal_rad)
    z = FARFIELD_RADIUS_UM * np.cos(polar_rad)

    return (x, y, z)


if __name__ == "__main__":
    dipole_component = Current.Jx

    kxs = np.linspace(-frequency, frequency, NUM_K)
    kys = np.linspace(-frequency, frequency, NUM_K)

    # Far fields are defined on the surface of a hemisphere.                                                 
    polar_rad = np.linspace(0, 0.5 * np.pi, NUM_POLAR)
    azimuthal_rad = np.linspace(0, 2 * np.pi, NUM_AZIMUTHAL)
    delta_azimuthal_rad = 2 * np.pi / (NUM_AZIMUTHAL - 1)

    current_components = [Current.Jx, Current.Jy, Current.Kx, Current.Ky]
    farfield_components = ["Ex", "Ey", "Ez", "Hx", "Hy", "Hz"]

    total_farfields = {}
    for component in farfield_components:
        total_farfields[component] = np.zeros(
            (NUM_POLAR, NUM_AZIMUTHAL),
            dtype=np.complex128
        )

    # Brillouin zone integration over 2D grid of (kx, ky).                                                   
    for kx in kxs:
        for ky in kys:

            # Skip wavevectors which are outside the light cone.                                             
            if np.sqrt(kx**2 + ky**2) >= (0.95 * frequency):
                continue

            kz = (frequency**2 - kx**2 - ky**2) ** 0.5

            if DEBUG_OUTPUT:
                print(f"(kx, ky, kz) = ({kx:.6f}, {ky:.6f}, {kz:.6f})")

            ex_dft, ey_dft, hx_dft, hy_dft = dipole_in_vacuum(dipole_component, kx, ky)

            if DEBUG_OUTPUT:
                print(f"ex_dft:, {ex_dft}")
                print(f"ey_dft:, {ey_dft}")
                print(f"hx_dft:, {hx_dft}")
                print(f"hy_dft:, {hy_dft}")

            # Rotation angle around z axis to force ky=0. 0 radians is +x.                                   
            rotation_rad = -math.atan2(ky, kx)

            if DEBUG_OUTPUT:
                print(f"rotation angle:, {math.degrees(rotation_rad):.2f}")

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

            if DEBUG_OUTPUT:
                print(f"electric_current:, {electric_current}")
                print(f"magnetic_current:, {magnetic_current}")
                print(f"electric_current_rotated:, {electric_current_rotated}")
                print(f"magnetic_current_rotated:, {magnetic_current_rotated}")

            current_amplitudes = [
                electric_current_rotated[0],
                electric_current_rotated[1],
                magnetic_current_rotated[0],
                magnetic_current_rotated[1],
            ]

            # Obtain the far fields over a hemisphere given a sheet current                                  
            # with linear in-plane polarization.                                                             

            azimuthal_index_shift = int(abs(rotation_rad) // delta_azimuthal_rad)
            inverse_rotation_matrix = np.transpose(rotation_matrix)

            farfields = {}
            for component in farfield_components:
                farfields[component] = np.zeros(
                    (NUM_POLAR, NUM_AZIMUTHAL),
                    dtype=np.complex128
                )

            for i, polar in enumerate(polar_rad):
                for j, azimuthal in enumerate(azimuthal_rad):
                    rx, ry, rz = spherical_to_cartesian(polar, azimuthal)
                    r_rotated = rotation_matrix @ np.transpose(np.array([rx, ry]))
                    # Note: r_rotated[1] (ry) is not used because k_rotated[1]=0 (ky).                       
                    phase = np.exp(1j * 2 * np.pi * (k_rotated[0] * r_rotated[0] + kz * rz))

                    # Obtain the farfields for each type of sheet current                                    
                    # and keep a running sum.                                                                
                    for current_component, current_amplitude in zip(
                        current_components, current_amplitudes
                    ):
                        if abs(current_amplitude) == 0:
                            continue

                        farfield_amplitudes_pol = farfield_amplitudes(
                            k_rotated[0],
                            kz,
                            r_rotated[0],
                            rz,
                            current_amplitude,
                            current_component,
                        )

                        farfield_pol = np.array(farfield_amplitudes_pol) * phase

                        if (
                            current_component.name == "Jx" or
                            current_component.name == "Ky"
                        ):
                            # P polarization                                                                 
                            farfields["Ex"][i, j] += farfield_pol[0]
                            farfields["Hy"][i, j] += farfield_pol[1]
                            farfields["Ez"][i, j] += farfield_pol[2]
                        elif (
                            current_component.name == "Jy" or
                            current_component.name == "Kx"
                        ):
                            # S polarization                                                                 
                            farfields["Hx"][i, j] += farfield_pol[0]
                            farfields["Ey"][i, j] += farfield_pol[1]
                            farfields["Hz"][i, j] += farfield_pol[2]

            for i in range(NUM_POLAR):
                for j in range(NUM_AZIMUTHAL):
                    # Map the rotated point to the closest azimuthal angle.                                  
                    azimuthal_index = (j + azimuthal_index_shift) % NUM_AZIMUTHAL
                    total_farfields["Ex"][i, azimuthal_index] += (
                        inverse_rotation_matrix[0, 0] * farfields["Ex"][i, j] +
                        inverse_rotation_matrix[0, 1] * farfields["Ey"][i, j]
                    )
                    total_farfields["Ey"][i, azimuthal_index] += (
                        inverse_rotation_matrix[1, 0] * farfields["Ex"][i, j] +
                        inverse_rotation_matrix[1, 1] * farfields["Ey"][i, j]
                    )
                    total_farfields["Ez"][i, azimuthal_index] += (
                        farfields["Ez"][i, j]
                    )

                    total_farfields["Hx"][i, azimuthal_index] += (
                        inverse_rotation_matrix[0, 0] * farfields["Hx"][i, j] +
                        inverse_rotation_matrix[0, 1] * farfields["Hy"][i, j]
                    )
                    total_farfields["Hy"][i, azimuthal_index] += (
                        inverse_rotation_matrix[1, 0] * farfields["Hx"][i, j] +
                        inverse_rotation_matrix[1, 1] * farfields["Hy"][i, j]
                    )
                    total_farfields["Hz"][i, azimuthal_index] += (
                        farfields["Hz"][i, j]
                    )

    if mp.am_master():
        np.savez(
            "dipole_farfields.npz",
            **total_farfields,
            RESOLUTION_UM=RESOLUTION_UM,
            WAVELENGTH_UM=WAVELENGTH_UM,
            FARFIELD_RADIUS_UM=FARFIELD_RADIUS_UM,
            NUM_POLAR=NUM_POLAR,
            NUM_AZIMUTHAL=NUM_AZIMUTHAL,
            polar_rad=polar_rad,
            azimuthal_rad=azimuthal_rad,
            kxs=kxs,
            kys=kys,
        )

