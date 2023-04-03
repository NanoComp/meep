"""Tutorial example for point-dipole sources in cylindrical coordinates.

This example demonstrates that the total and radiated flux from a point dipole
in a dielectric layer (a quantum well)) above a lossless ground plane (an LED)
computed in cylindrical coordinates as part of the calculation of the extraction
efficiency is independent of the dipole's position in the radial direction.

reference: https://meep.readthedocs.io/en/latest/Python_Tutorials/Cylindrical_Coordinates/#nonaxisymmetric-dipole-sources
"""

from typing import Tuple

import meep as mp
import numpy as np


resolution = 80  # pixels/μm
n = 2.4  # refractive index of dielectric layer
wvl = 1.0  # wavelength (in vacuum)
fcen = 1 / wvl  # center frequency of source/monitor


def led_flux(dmat: float, h: float, rpos: float, m: int) -> Tuple[float, float]:
    """Computes the radiated and total flux of a point source embedded
       within a dielectric layer above a lossless ground plane.

    Args:
       dmat: thickness of dielectric layer.
       h: height of dipole above ground plane as a fraction of dmat.
       rpos: position of source in radial direction.
       m: angular φ dependence of the fields exp(imφ).

    Returns:
       The radiated and total flux as a 2-Tuple.
    """
    L = 20  # length of non-PML region in radial direction
    dair = 1.0  # thickness of air padding
    dpml = 1.0  # PML thickness
    sr = L + dpml
    sz = dmat + dair + dpml
    cell_size = mp.Vector3(sr, 0, sz)

    boundary_layers = [
        mp.PML(dpml, direction=mp.R),
        mp.PML(dpml, direction=mp.Z, side=mp.High),
    ]

    src_pt = mp.Vector3(rpos, 0, -0.5 * sz + h * dmat)
    sources = [
        mp.Source(
            src=mp.GaussianSource(fcen, fwidth=0.1 * fcen),
            component=mp.Er,
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

    flux_air_mon = sim.add_flux(
        fcen,
        0,
        1,
        mp.FluxRegion(
            center=mp.Vector3(0.5 * L, 0, 0.5 * sz - dpml),
            size=mp.Vector3(L, 0, 0),
        ),
        mp.FluxRegion(
            center=mp.Vector3(L, 0, 0.5 * sz - dpml - 0.5 * dair),
            size=mp.Vector3(0, 0, dair),
        ),
    )

    sim.run(
        mp.dft_ldos(fcen, 0, 1),
        until_after_sources=mp.stop_when_fields_decayed(
            50.0,
            mp.Er,
            src_pt,
            1e-8,
        ),
    )

    flux_air = mp.get_fluxes(flux_air_mon)[0]

    if rpos == 0:
        dV = np.pi / (resolution**3)
    else:
        dV = 2 * np.pi * rpos / (resolution**2)

    # total flux from point source via LDOS
    flux_src = -np.real(sim.ldos_Fdata[0] * np.conj(sim.ldos_Jdata[0])) * dV

    print(f"flux-cyl:, {rpos:.2f}, {m:3d}, {flux_src:.6f}, {flux_air:.6f}")

    return flux_air, flux_src


if __name__ == "__main__":
    layer_thickness = 0.7 * wvl / n
    dipole_height = 0.5

    # r = 0 source requires a single simulation with m = ±1
    rpos = 0
    m = 1
    flux_air, flux_src = led_flux(
        layer_thickness,
        dipole_height,
        rpos,
        m,
    )
    ext_eff = flux_air / flux_src
    print(f"exteff:, {rpos}, {ext_eff:.6f}")

    # r > 0 source requires Fourier-series expansion of φ
    flux_tol = 1e-5  # threshold flux to determine when to truncate expansion
    rpos = [3.5, 6.7, 9.5]
    for rp in rpos:
        cutoff_M = int(rp * 2 * np.pi * fcen * n)  # analytic upper bound on m
        ms = range(cutoff_M + 1)
        flux_src_tot = 0
        flux_air_tot = 0
        flux_air_max = 0
        for m in ms:
            flux_air, flux_src = led_flux(
                layer_thickness,
                dipole_height,
                rp,
                m,
            )
            flux_air_tot += flux_air if m == 0 else 2 * flux_air
            flux_src_tot += flux_src if m == 0 else 2 * flux_src
            if flux_air > flux_air_max:
                flux_air_max = flux_air
            if m > 0 and (flux_air / flux_air_max) < flux_tol:
                break

        ext_eff = flux_air_tot / flux_src_tot
        print(f"exteff:, {rp}, {ext_eff:.6f}")
