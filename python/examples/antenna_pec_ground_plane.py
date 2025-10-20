"""Radiation pattern of a dipole antenna above a ground plane.

Tutorial reference:

https://meep.readthedocs.io/en/latest/Python_Tutorials/Near_to_Far_Field_Spectra/#antenna-above-a-perfect-electric-conductor-ground-plane
"""

import math

import matplotlib.pyplot as plt
import meep as mp
import numpy as np


RESOLUTION_UM = 200
PML_UM = 1.0
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


def antenna_radiation_pattern(dipole_polarization: int, free_space: bool) -> np.ndarray:
    """Returns the radiation pattern of an antenna.

    Args:
        dipole_polarization: polarization of the dipole antenna (Ez or Hz).
        free_space: whether to use free space or PEC ground plane.

    Returns:
        The radiation pattern as a 1D array.
    """
    if dipole_polarization not in (mp.Ez, mp.Hz):
        raise ValueError("dipole_polarization must be either Ez or Hz")

    cell_um = 8.0
    sxy = PML_UM + cell_um + PML_UM
    cell_size = mp.Vector3(sxy, sxy, 0)
    boundary_layers = [mp.PML(PML_UM)]

    if free_space:
        sources = [
            mp.Source(
                src=mp.GaussianSource(frequency, fwidth=0.2 * frequency),
                center=mp.Vector3(),
                component=dipole_polarization,
            )
        ]

        if dipole_polarization == mp.Hz:
            symmetries = [mp.Mirror(mp.X, phase=-1), mp.Mirror(mp.Y, phase=-1)]
        else:
            symmetries = [mp.Mirror(mp.X), mp.Mirror(mp.Y)]

        layout_filename = "free_space_layout"

    else:
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

        layout_filename = "ground_plane_layout"

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
            center=mp.Vector3(0, 0.5 * cell_um), size=mp.Vector3(cell_um, 0)
        ),
        mp.Near2FarRegion(
            center=mp.Vector3(0, -0.5 * cell_um), size=mp.Vector3(cell_um, 0), weight=-1
        ),
        mp.Near2FarRegion(
            center=mp.Vector3(0.5 * cell_um, 0), size=mp.Vector3(0, cell_um)
        ),
        mp.Near2FarRegion(
            center=mp.Vector3(-0.5 * cell_um, 0), size=mp.Vector3(0, cell_um), weight=-1
        ),
    )

    fig, ax = plt.subplots()
    sim.plot2D(ax=ax)
    if mp.am_master():
        fig.savefig(layout_filename + ".png", dpi=150, bbox_inches="tight")

    sim.run(until_after_sources=mp.stop_when_dft_decayed())

    return radiation_pattern(sim, nearfield_box)


if __name__ == "__main__":
    dipole_polarization = mp.Ez
    radial_flux_free_space = antenna_radiation_pattern(dipole_polarization, True)
    radial_flux_ground_plane = antenna_radiation_pattern(dipole_polarization, False)

    # The radiation pattern of a two-element antenna
    # array is equivalent to the radiation pattern of
    # a single antenna multiplied by its array factor
    # as derived in Section 6.2 "Two-Element Array" of
    # Antenna Theory: Analysis and Design, Fourth Edition
    # (2016) by C.A. Balanis.
    k_free_space = 2 * math.pi / (WAVELENGTH_UM / N_BACKGROUND)
    radial_flux_analytic = np.zeros(NUM_POLAR)
    for i in range(NUM_POLAR):
        radial_flux_analytic[i] = (
            radial_flux_free_space[i]
            * 2
            * math.sin(k_free_space * ANTENNA_HEIGHT_UM * math.cos(polar_rad[i]))
        )

    radial_flux_ground_plane_norm = radial_flux_ground_plane / np.max(
        radial_flux_ground_plane
    )
    radial_flux_analytic_norm = (radial_flux_analytic / max(radial_flux_analytic)) ** 2
    rel_err = np.linalg.norm(
        radial_flux_ground_plane_norm - radial_flux_analytic_norm
    ) / np.linalg.norm(radial_flux_analytic_norm)
    print(f"error:, {rel_err:.6f}")

    fig, ax = plt.subplots()
    ax.plot(np.degrees(polar_rad), radial_flux_ground_plane_norm, "b-", label="Meep")
    ax.plot(np.degrees(polar_rad), radial_flux_analytic_norm, "r-", label="theory")
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
        fig.savefig("antenna_radiation_pattern.png", dpi=150, bbox_inches="tight")
