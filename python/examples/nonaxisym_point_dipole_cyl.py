"""
Verifies that the radiated flux from a point dipole in a dielectric layer
above a lossless ground plane computed in cylindrical coordinates
is the same for a dipole placed at either r=0 or r≠0.
"""

import argparse
import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
import meep as mp
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument(
    "-res", type=float, default=25.0, help="resolution (default: 25 pixels/um)"
)
parser.add_argument(
    "-L", type=float, default=5.0, help="length of non-PML region (default: 5.0 um)"
)
parser.add_argument(
    "-dpml", type=float, default=1.0, help="PML thickness in r and z (default: 1.0 um)"
)
args = parser.parse_args()

resolution = args.res
L = args.L
dpml = args.dpml

dair = 1.0  # thickness of air padding
n = 2.4  # refractive index of surrounding medium
wvl = 1.0  # wavelength (in vacuum)
fcen = 1 / wvl  # center frequency of source/monitor


def radiated_flux(dmat: float, h: float, rpos: float, m: int) -> float:
    """Computes the radiated flux of a point dipole embedded
       within a dielectric layer above a lossless ground plane in
       cylindrical coordinates.

    Args:
      dmat: thickness of dielectric layer.
      h: height of dipole above ground plane as fraction of dmat.
      rpos: position of source in r direction.
      m: angular φ dependence of the fields exp(imφ).
    """
    sr = L + dpml
    sz = dmat + dair + dpml
    cell_size = mp.Vector3(sr, 0, sz)

    boundary_layers = [
        mp.PML(dpml, direction=mp.R),
        mp.PML(dpml, direction=mp.Z, side=mp.High),
    ]

    src_cmpt = mp.Er
    src_pt = mp.Vector3(rpos, 0, -0.5 * sz + h * dmat)
    sources = [
        mp.Source(
            src=mp.GaussianSource(fcen, fwidth=0.1 * fcen),
            component=src_cmpt,
            center=src_pt,
        ),
    ]

    geometry = [
        mp.Block(
            material=mp.Medium(index=n),
            center=mp.Vector3(0, 0, -0.5 * sz + 0.5 * dmat),
            size=mp.Vector3(mp.inf, mp.inf, dmat),
        )
    ]

    sim = mp.Simulation(
        resolution=resolution,
        cell_size=cell_size,
        dimensions=mp.CYLINDRICAL,
        m=m,
        boundary_layers=boundary_layers,
        sources=sources,
        geometry=geometry,
    )

    flux_air_z = sim.add_flux(
        fcen,
        0,
        1,
        mp.FluxRegion(
            center=mp.Vector3(0.5 * L, 0, 0.5 * sz - dpml),
            size=mp.Vector3(L, 0, 0),
        ),
    )

    flux_air_r = sim.add_flux(
        fcen,
        0,
        1,
        mp.FluxRegion(
            center=mp.Vector3(L, 0, 0.5 * sz - dpml - 0.5 * dair),
            size=mp.Vector3(0, 0, dair),
        ),
    )

    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            50.0,
            src_cmpt,
            src_pt,
            1e-8,
        ),
    )

    flux_z = mp.get_fluxes(flux_air_z)[0]
    flux_r = mp.get_fluxes(flux_air_r)[0]
    flux_tot = flux_r + flux_z
    print(f"flux:, {flux_r:.10f}, {flux_z:.10f}, {flux_tot:.10f}")

    return flux_tot


if __name__ == "__main__":
    layer_thickness = 0.7 * wvl / n
    dipole_height = 0.5

    # r = 0
    rpos = 0
    m = +1
    ref_out_flux = radiated_flux(
        layer_thickness,
        dipole_height,
        rpos,
        m,
    )
    print(f"flux0-m:, {m}, {ref_out_flux}")

    # r ≠ 0
    rpos = 1.0
    # empirical
    cutoff_M = lambda rp: int(16 * rp + 8)
    # analytic
    # cutoff_M = int(2 * rp * 2 * np.pi * fcen * n) + 2
    ms = range(cutoff_M(rp) + 1)
    out_flux = []
    for m in ms:
        out_flux.append(
            radiated_flux(
                layer_thickness,
                dipole_height,
                rp,
                m,
            )
        )
        print(f"flux1-m:, {rp}, {m}, {out_flux[-1]}")

    flux_sum = out_flux[0] + 2 * sum(out_flux[1:])
    print(f"flux1:, {rp}, {flux_sum:.8f}")
