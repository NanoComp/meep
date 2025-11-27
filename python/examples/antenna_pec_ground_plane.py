"""Radiation pattern of a dipole antenna above a ground plane.

Tutorial reference:

https://meep.readthedocs.io/en/latest/Python_Tutorials/Near_to_Far_Field_Spectra/#antenna-above-a-perfect-electric-conductor-ground-plane
"""

import math

import matplotlib.pyplot as plt
import meep as mp
import numpy as np


RESOLUTION_UM = 100
PML_UM = 0.5
BULK_UM = 6.0
N_BACKGROUND = 1.2
ANTENNA_HEIGHT_UM = 1.25
WAVELENGTH_UM = 0.65
NUM_POLAR = 50
FARFIELD_RADIUS_UM = 1000 * WAVELENGTH_UM
GREENCYL_TOL = 1e-8

frequency = 1 / WAVELENGTH_UM
polar_rad = np.linspace(0, 0.5 * math.pi, NUM_POLAR)


def radiation_pattern(sim: mp.Simulation, n2f_mon: mp.DftNear2Far) -> np.ndarray:
    """Computes the radiation pattern from the far fields.

    Args:
        sim: a `Simulation` object.
        n2f_mon: a `DftNear2Far` object returned by `Simulation.add_near2far`.

    Returns:
        Array of radial Poynting flux, one for each point on the circumference of
        a quarter circle with angular range of [0, π/2] rad. 0 rad is the +y
        direction and π/2 is +x.
    """
    e_field = np.zeros((NUM_POLAR, 3), dtype=np.complex128)
    h_field = np.zeros((NUM_POLAR, 3), dtype=np.complex128)
    for i in range(NUM_POLAR):
        far_field = sim.get_farfield(
            n2f_mon,
            mp.Vector3(
                FARFIELD_RADIUS_UM * math.sin(polar_rad[i]),
                FARFIELD_RADIUS_UM * math.cos(polar_rad[i]),
                0,
            ),
            GREENCYL_TOL,
        )
        e_field[i, :] = [far_field[j] for j in range(3)]
        h_field[i, :] = [far_field[j + 3] for j in range(3)]

    flux_x = np.real(
        np.conj(e_field[:, 1]) * h_field[:, 2] - np.conj(e_field[:, 2]) * h_field[:, 1]
    )
    flux_y = np.real(
        np.conj(e_field[:, 2]) * h_field[:, 0] - np.conj(e_field[:, 0]) * h_field[:, 2]
    )
    flux_r = np.sqrt(np.square(flux_x) + np.square(flux_y))

    return flux_r


def antenna_above_ground_plane(dipole_polarization: int) -> np.ndarray:
    """Returns the radiation pattern of an antenna above a ground plane.

    Args:
        dipole_polarization: polarization of the dipole antenna (Ez or Hz).

    Returns:
        The radiation pattern as a 1D array.
    """
    if dipole_polarization not in (mp.Ez, mp.Hz):
        raise ValueError("dipole_polarization must be either Ez or Hz")

    sxy = PML_UM + BULK_UM + PML_UM
    cell_size = mp.Vector3(sxy, sxy, 0)
    boundary_layers = [mp.PML(PML_UM)]

    # The near-to-far field transformation feature can only be applied
    # to structures which can be fully enclosed. It cannot be applied to
    # the ground plane which extends to infinity. As a workaround, we use
    # two antennas of opposite sign surrounded by a single near2far box
    # which encloses both antennas. We then use an odd mirror symmetry to
    # divide the computational cell in half which is effectively
    # equivalent to a PEC boundary condition on one side.
    # Note: This setup means that the radiation pattern can only
    # be measured in the top half above the dipole antenna.
    sources = [
        mp.Source(
            src=mp.GaussianSource(frequency, fwidth=0.2 * frequency),
            component=dipole_polarization,
            center=mp.Vector3(0, ANTENNA_HEIGHT_UM),
        ),
        mp.Source(
            src=mp.GaussianSource(frequency, fwidth=0.2 * frequency),
            component=dipole_polarization,
            center=mp.Vector3(0, -ANTENNA_HEIGHT_UM),
            amplitude=-1.0 if dipole_polarization == mp.Ez else +1.0,
        ),
    ]

    if dipole_polarization == mp.Hz:
        symmetries = [mp.Mirror(mp.X, phase=-1), mp.Mirror(mp.Y)]
    else:
        symmetries = [mp.Mirror(mp.X), mp.Mirror(mp.Y, phase=-1)]

    sim = mp.Simulation(
        resolution=RESOLUTION_UM,
        cell_size=cell_size,
        boundary_layers=boundary_layers,
        sources=sources,
        symmetries=symmetries,
        default_material=mp.Medium(index=N_BACKGROUND),
    )

    nearfield_box = sim.add_near2far(
        frequency,
        0,
        1,
        mp.Near2FarRegion(
            center=mp.Vector3(0, 0.5 * BULK_UM), size=mp.Vector3(BULK_UM, 0)
        ),
        mp.Near2FarRegion(
            center=mp.Vector3(0, -0.5 * BULK_UM), size=mp.Vector3(BULK_UM, 0), weight=-1
        ),
        mp.Near2FarRegion(
            center=mp.Vector3(0.5 * BULK_UM, 0), size=mp.Vector3(0, BULK_UM)
        ),
        mp.Near2FarRegion(
            center=mp.Vector3(-0.5 * BULK_UM, 0), size=mp.Vector3(0, BULK_UM), weight=-1
        ),
    )

    fig, ax = plt.subplots()
    sim.plot2D(ax=ax)
    if mp.am_master():
        fig.savefig(
            "antenna_above_ground_plane_layout.png", dpi=150, bbox_inches="tight"
        )

    sim.run(until_after_sources=mp.stop_when_dft_decayed())

    return radiation_pattern(sim, nearfield_box)


if __name__ == "__main__":
    dipole_polarization = mp.Ez
    radial_flux_meep = antenna_above_ground_plane(dipole_polarization)
    # Normalize radiation pattern to lie in range [0, 1].
    radial_flux_meep = radial_flux_meep / np.max(radial_flux_meep)

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

    rel_err = np.linalg.norm(radial_flux_meep - radial_flux_analytic) / np.linalg.norm(
        radial_flux_analytic
    )
    print(f"error:, {rel_err:.6f}")

    fig, ax = plt.subplots()
    ax.plot(np.degrees(polar_rad), radial_flux_meep, "b-", label="Meep")
    ax.plot(np.degrees(polar_rad), radial_flux_analytic, "r-", label="theory")
    ax.set_xlabel("angle (degrees)")
    ax.set_ylabel("radiation pattern")
    ax.set_title(
        f"antenna with {'E' if dipole_polarization==mp.Ez else 'H'}$_z$ "
        "polarization above PEC ground plane"
    )
    ax.set_xlim(0, 90.0)
    ax.set_ylim(0, 1.0)
    ax.legend()
    if mp.am_master():
        fig.savefig(
            "antenna_above_ground_plane_radiation_pattern.png",
            dpi=150,
            bbox_inches="tight",
        )
