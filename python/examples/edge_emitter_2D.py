"""Computes the extraction efficiency of a cleaved stack using a 2D sim."""

import argparse
from typing import Tuple

import meep as mp
import matplotlib.pyplot as plt
import numpy as np


RESOLUTION_UM = 25
WAVELENGTH_UM = 1.0
CLADDING_UM = 0.5
ACTIVE_UM = 0.3
SUBSTRATE_UM = 1.0
PML_UM = 0.5
AIR_UM = 1.0
SIDE_UM = 3.0
N_SUBSTRATE = 3.2
N_CLADDING = 1.5
N_ACTIVE = 2.3
FIELD_DECAY_PERIOD = 25.0
FIELD_DECAY_TOL = 1e-6
DELTA_KZ = 0.02
DIPOLE_FLUX_TOL = 0.02
DEBUG_OUTPUT = 0


def line_current_in_cleaved_stack(
    dipole_pol: str, dipole_pos_um: float, kz: float
) -> Tuple[float, float]:
    """Computes the flux of a line current in a cleaved multilayer stack.

    Args:
        dipole_pol: polarization of the dipole ("x", "y", or "z").
        dipole_pos_um: position of the dipole relative to the edge of the
          active region.
        kz: z-component of the Bloch-periodic wavevector.

    Returns:
        The flux emitted by a line current in the active region and into air as
          a 2-tuple.
    """
    if dipole_pol == "x":
        src_cmpt = mp.Ex
    elif dipole_pol == "y":
        src_cmpt = mp.Ey
    elif dipole_pol == "z":
        src_cmpt = mp.Ez
    else:
        raise ValueError("dipole_pol must be x, y, or z.")

    cell_x_um = PML_UM + SIDE_UM + AIR_UM + PML_UM
    cell_y_um = (
        PML_UM + SUBSTRATE_UM + CLADDING_UM + ACTIVE_UM + CLADDING_UM + AIR_UM + PML_UM
    )
    cell_size = mp.Vector3(cell_x_um, cell_y_um, 0)

    boundary_layers = [mp.PML(thickness=PML_UM)]

    frequency = 1 / WAVELENGTH_UM

    src_pos = mp.Vector3(
        -0.5 * cell_x_um + PML_UM + SIDE_UM - dipole_pos_um,
        -0.5 * cell_y_um + PML_UM + SUBSTRATE_UM + CLADDING_UM + 0.5 * ACTIVE_UM,
        0,
    )
    sources = [
        mp.Source(
            src=mp.GaussianSource(frequency, fwidth=0.1 * frequency),
            center=src_pos,
            component=src_cmpt,
        )
    ]

    geometry = [
        mp.Block(
            material=mp.Medium(index=N_SUBSTRATE),
            center=mp.Vector3(
                -0.5 * cell_x_um + 0.5 * (PML_UM + SIDE_UM),
                -0.5 * cell_y_um + 0.5 * (PML_UM + SUBSTRATE_UM),
                0,
            ),
            size=mp.Vector3(PML_UM + SIDE_UM, PML_UM + SUBSTRATE_UM, mp.inf),
        ),
        mp.Block(
            material=mp.Medium(index=N_CLADDING),
            center=mp.Vector3(
                -0.5 * cell_x_um + 0.5 * (PML_UM + SIDE_UM),
                -0.5 * cell_y_um + PML_UM + SUBSTRATE_UM + 0.5 * CLADDING_UM,
                0,
            ),
            size=mp.Vector3(PML_UM + SIDE_UM, CLADDING_UM, mp.inf),
        ),
        mp.Block(
            material=mp.Medium(index=N_ACTIVE),
            center=mp.Vector3(
                -0.5 * cell_x_um + 0.5 * (PML_UM + SIDE_UM),
                -0.5 * cell_y_um
                + PML_UM
                + SUBSTRATE_UM
                + CLADDING_UM
                + 0.5 * ACTIVE_UM,
                0,
            ),
            size=mp.Vector3(PML_UM + SIDE_UM, ACTIVE_UM, mp.inf),
        ),
        mp.Block(
            material=mp.Medium(index=N_CLADDING),
            center=mp.Vector3(
                -0.5 * cell_x_um + 0.5 * (PML_UM + SIDE_UM),
                -0.5 * cell_y_um
                + PML_UM
                + SUBSTRATE_UM
                + CLADDING_UM
                + ACTIVE_UM
                + 0.5 * CLADDING_UM,
                0,
            ),
            size=mp.Vector3(PML_UM + SIDE_UM, CLADDING_UM, mp.inf),
        ),
    ]

    sim = mp.Simulation(
        resolution=RESOLUTION_UM,
        cell_size=cell_size,
        k_point=mp.Vector3(0, 0, kz),
        kz_2d="real/imag",
        force_complex_fields=False,
        sources=sources,
        geometry=geometry,
        boundary_layers=boundary_layers,
    )

    flux_mon = sim.add_flux(
        frequency,
        0,
        1,
        mp.FluxRegion(
            center=mp.Vector3(
                -0.5 * cell_x_um + PML_UM + 0.5 * (cell_x_um - 2 * PML_UM),
                0.5 * cell_y_um - PML_UM,
                0,
            ),
            size=mp.Vector3(cell_x_um - 2 * PML_UM, 0, 0),
        ),
        mp.FluxRegion(
            center=mp.Vector3(
                0.5 * cell_x_um - PML_UM,
                0.5 * cell_y_um - PML_UM - 0.5 * (cell_y_um - 2 * PML_UM),
                0,
            ),
            size=mp.Vector3(0, cell_y_um - 2 * PML_UM, 0),
        ),
        mp.FluxRegion(
            center=mp.Vector3(
                0.5 * cell_x_um - PML_UM - +0.5 * AIR_UM, -0.5 * cell_y_um + PML_UM, 0
            ),
            size=mp.Vector3(AIR_UM, 0, 0),
            weight=-1.0,
        ),
    )

    if DEBUG_OUTPUT:
        if mp.am_master():
            fig, ax = plt.subplots()
            sim.plot2D(ax=ax)
            fig.savefig("edge_emitter_layout.png", dpi=150, bbox_inches="tight")

    sim.run(
        mp.dft_ldos(frequency, 0, 1),
        until_after_sources=mp.stop_when_fields_decayed(
            FIELD_DECAY_PERIOD,
            src_cmpt,
            src_pos,
            FIELD_DECAY_TOL,
        ),
    )

    pixel_area = (1 / RESOLUTION_UM) ** 2
    dipole_flux = -np.real(sim.ldos_Fdata[0] * np.conj(sim.ldos_Jdata[0])) * pixel_area

    air_flux = mp.get_fluxes(flux_mon)[0]

    return dipole_flux, air_flux


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "dipole_pol",
        type=str,
        choices=["x", "y", "z"],
        help="polarization of the electric dipole (x, y, z)",
    )
    parser.add_argument(
        "dipole_pos_um",
        type=float,
        help="position of the dipole relative to the edge of the cleaved facet",
    )
    args = parser.parse_args()

    num_kz = 0
    max_dipole_flux = 0
    dipole_flux = []
    air_flux = []

    # Brillouin-zone integration over kz.
    while True:
        kz = num_kz * DELTA_KZ
        dipole_flux_kz, air_flux_kz = line_current_in_cleaved_stack(
            args.dipole_pol, args.dipole_pos_um, kz
        )
        print(f"kz:, {num_kz}, {kz:.2f}, {dipole_flux_kz:.2f}, " f"{air_flux_kz:.2f}")

        dipole_flux.append(dipole_flux_kz)
        air_flux.append(air_flux_kz)
        num_kz += 1

        if dipole_flux_kz > max_dipole_flux:
            max_dipole_flux = dipole_flux_kz
        elif (dipole_flux_kz / max_dipole_flux) < DIPOLE_FLUX_TOL:
            break

    dipole_flux = np.array(dipole_flux)
    air_flux = np.array(air_flux)
    extraction_efficiency = (air_flux[0] + 2 * np.sum(air_flux[1:])) / (
        dipole_flux[0] + 2 * np.sum(dipole_flux[1:])
    )

    print(
        f"extraction_efficiency:, {args.dipole_pol}, {args.dipole_pos_um:.2f},"
        f" {100 * extraction_efficiency:.2f}%"
    )
