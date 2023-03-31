"""Tutorial example for non-axisymmetric point sources in cylindrical coords.

This example demonstrates that the radiated flux from a point dipole in
a dielectric layer above a lossless ground plane computed in cylindrical
coordinates is the same for a dipole placed anywhere along the radial
direction.

Reference: https://meep.readthedocs.io/en/latest/Python_Tutorials/Cylindrical_Coordinates/#nonaxisymmetric-dipole-sources
"""

import meep as mp
import numpy as np


resolution = 30  # pixels/μm
n = 2.4  # refractive index of dielectric layer
wvl = 1.0  # wavelength (in vacuum)
fcen = 1 / wvl  # center frequency of source/monitor


def radiated_flux_cyl(dmat: float, h: float, rpos: float, m: int) -> float:
    """Computes the radiated flux of a point dipole embedded
       within a dielectric layer above a lossless ground plane in
       cylindrical coordinates.

    Args:
       dmat: thickness of dielectric layer.
       h: height of dipole above ground plane as a fraction of dmat.
       rpos: position of source in radial direction.
       m: angular φ dependence of the fields exp(imφ).

    Returns:
       The radiated flux in cylindrical coordinates.
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

    return flux_tot


def radiated_flux_3d(rpos: float) -> float:
    """Computes the radiated flux in 3d Cartesian coordinates using a
       point-dipole source in cylindrical coordinates.

    Args:
       rpos: position of source in radial direction.

    Returns:
       The radiated flux in 3d Cartesian coordinates.
    """
    layer_thickness = 0.7 * wvl / n
    dipole_height = 0.5

    if rpos == 0:
        # r = 0 source requires a single simulation with m = ±1
        m = 1
        flux_cyl = radiated_flux_cyl(
            layer_thickness,
            dipole_height,
            rpos,
            m,
        )

        flux_3d = flux_cyl * 2 * (resolution / np.pi) ** 2

    else:
        # r > 0 source requires Fourier-series expansion of φ
        flux_tol = 1e-6  # relative tolerance for flux for truncating expansion
        cutoff_M = int(2 * rpos * 2 * np.pi * fcen * n)
        ms = range(cutoff_M + 1)
        flux_cyl_tot = 0
        flux_cyl_max = 0
        for m in ms:
            flux_cyl = radiated_flux_cyl(
                layer_thickness,
                dipole_height,
                rpos,
                m,
            )
            print(f"flux-m:, {rpos}, {m}, {flux_cyl:.6f}")
            flux_cyl_tot += flux_cyl if m == 0 else 2 * flux_cyl
            if flux_cyl > flux_cyl_max:
                flux_cyl_max = flux_cyl
            if m > 0 and (flux_cyl / flux_cyl_max) < flux_tol:
                break

        flux_3d = (flux_cyl_tot / rpos**2) * (2 / np.pi**5)

    print(f"flux-3d:, {rpos:.2f}, {flux_3d:.6f}")

    return flux_3d


if __name__ == "__main__":
    # reference result for dipole at r = 0
    Pcyl_ref = radiated_flux_3d(0)

    rs = [3.3, 7.5, 12.1]
    for r in rs:
        Pcyl = radiated_flux_3d(r)
        err = abs(Pcyl - Pcyl_ref) / Pcyl_ref
        print(f"err:, {r}, {Pcyl:.6f}, {err:.6f}")
