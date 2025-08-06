# Computes the Purcell enhancement factor of an in-plane dipole in a planar
# dielectric cavity with lossless metallic walls. The result is computed in
# cylindrical and 3D coordinates and compared with the analytic result from:
# I. Abram et al., IEEE J. Quantum Electronics, Vol. 34, pp. 71-76 (1998).

# tutorial reference:
# https://meep.readthedocs.io/en/latest/Python_Tutorials/Local_Density_of_States/#planar-cavity-with-lossless-metallic-walls

from typing import Optional

import matplotlib.pyplot as plt
import meep as mp
import numpy as np


# Note: Meep may round the cell dimensions to an integer number of pixels which
# could modify the cavity structure.
RESOLUTION_UM = 71

PML_UM = 0.5
BULK_UM = 6.0
N_CAVITY = 2.4
WAVELENGTH_UM = 1.0
FIELD_DECAY_TOL = 1e-6
FIELD_DECAY_PERIOD = 20

frequency = 1 / WAVELENGTH_UM


def ldos_cyl(cavity_um: Optional[float] = None) -> float:
    """Computes the LDOS of a dipole in a cavity or bulk media in cyl. coords.

    Args:
        cavity_um: thickness of the cavity. If None, bulk media is used.

    Returns:
        The LDOS of the dipole.
    """
    if cavity_um is None:
        cell_z_um = BULK_UM + 2 * PML_UM
        pml_layers = [mp.PML(thickness=PML_UM)]
    else:
        cell_z_um = cavity_um
        pml_layers = [mp.PML(thickness=PML_UM, direction=mp.R)]

    cell_r_um = BULK_UM + PML_UM
    cell_size = mp.Vector3(cell_r_um, 0, cell_z_um)

    # An Er source at r = 0 and m=±1 needs to be slightly offset.
    # https://github.com/NanoComp/meep/issues/2704
    dipole_rpos_um = 1.5 / RESOLUTION_UM

    src_pt = mp.Vector3(dipole_rpos_um, 0, 0)
    sources = [
        mp.Source(
            src=mp.GaussianSource(frequency, fwidth=0.2 * frequency),
            component=mp.Er,
            center=src_pt,
        )
    ]

    sim = mp.Simulation(
        resolution=RESOLUTION_UM,
        cell_size=cell_size,
        boundary_layers=pml_layers,
        sources=sources,
        dimensions=mp.CYLINDRICAL,
        m=-1,
        default_material=mp.Medium(index=N_CAVITY),
    )

    sim.run(
        mp.dft_ldos(frequency, 0, 1),
        until_after_sources=mp.stop_when_fields_decayed(
            FIELD_DECAY_PERIOD, mp.Er, src_pt, FIELD_DECAY_TOL
        ),
    )

    return sim.ldos_data[0]


def ldos_3d(cavity_um: Optional[float] = None) -> float:
    """Computes the LDOS of a dipole in a cavity or bulk media in 3D coords.

    Args:
        cavity_um: thickness of the cavity. If None, bulk media is used.

    Returns:
        The LDOS of the dipole.
    """
    if cavity_um is None:
        size_z_um = BULK_UM + 2 * PML_UM
        pml_layers = [mp.PML(thickness=PML_UM)]
    else:
        size_z_um = cavity_um
        pml_layers = [
            mp.PML(thickness=PML_UM, direction=mp.X),
            mp.PML(thickness=PML_UM, direction=mp.Y),
        ]

    size_xy_um = BULK_UM + 2 * PML_UM
    cell_size = mp.Vector3(size_xy_um, size_xy_um, size_z_um)

    sources = [
        mp.Source(
            src=mp.GaussianSource(frequency, fwidth=0.2 * frequency),
            component=mp.Ex,
            center=mp.Vector3(),
        )
    ]

    symmetries = [
        mp.Mirror(direction=mp.X, phase=-1),
        mp.Mirror(direction=mp.Y),
        mp.Mirror(direction=mp.Z),
    ]

    sim = mp.Simulation(
        resolution=RESOLUTION_UM,
        cell_size=cell_size,
        boundary_layers=pml_layers,
        sources=sources,
        symmetries=symmetries,
        default_material=mp.Medium(index=N_CAVITY),
    )

    sim.run(
        mp.dft_ldos(frequency, 0, 1),
        until_after_sources=mp.stop_when_fields_decayed(
            FIELD_DECAY_PERIOD, mp.Ex, mp.Vector3(), FIELD_DECAY_TOL
        ),
    )

    return sim.ldos_data[0]


if __name__ == "__main__":
    ldos_bulk_cyl = ldos_cyl()
    ldos_bulk_3d = ldos_3d()

    cavity_um = np.arange(0.50, 2.55, 0.05)
    vacuum_cavity_um = cavity_um * WAVELENGTH_UM / N_CAVITY

    num_cavity_um = cavity_um.shape[0]
    ldos_cavity_cyl = np.zeros(num_cavity_um)
    ldos_cavity_3d = np.zeros(num_cavity_um)

    for j in range(num_cavity_um):
        ldos_cavity_cyl[j] = ldos_cyl(vacuum_cavity_um[j])
        ldos_cavity_3d[j] = ldos_3d(vacuum_cavity_um[j])
        purcell_cyl = ldos_cavity_cyl[j] / ldos_bulk_cyl
        purcell_3d = ldos_cavity_3d[j] / ldos_bulk_3d
        print(f"purcell:, {cavity_um[j]:.3f}, {purcell_cyl:.6f}, {purcell_3d:.6f}")

    # Purcell enhancement factor (relative to bulk medium)
    purcell_meep_cyl = ldos_cavity_cyl / ldos_bulk_cyl
    purcell_meep_3d = ldos_cavity_3d / ldos_bulk_3d

    # Equation 7 of 1998 reference.
    purcell_theory = 3 * np.fix(cavity_um + 0.5) / (4 * cavity_um) + (
        4 * np.power(np.fix(cavity_um + 0.5), 3) - np.fix(cavity_um + 0.5)
    ) / (16 * np.power(cavity_um, 3))

    if mp.am_master():
        fig, ax = plt.subplots()
        ax.plot(cavity_um, purcell_meep_3d, "b-", label="Meep (3d)")
        ax.plot(cavity_um, purcell_meep_cyl, "r-", label="Meep (cylin.)")
        ax.plot(cavity_um, purcell_theory, "g-", label="theory")
        ax.plot(cavity_um, np.ones(len(cavity_um)), "k--")
        ax.set_xlabel("cavity thickness (in media)")
        ax.set_ylabel("Purcell enhancement factor")
        ax.set_title(
            "in-plane dipole at λ=1.0 μm in a planar cavity\n"
            "with n=2.4 and lossless metallic walls"
        )
        ax.axis([0.5, 2.5, 0.4, 3.1])
        ax.legend()
        fig.savefig(
            "cavity_purcell_factor_vs_thickness.png", dpi=150, bbox_inches="tight"
        )
