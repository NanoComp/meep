# Verifies that the extraction efficiency of a point dipole in a
# dielectric layer above a lossless ground plane computed in
# cylindrical and 3D Cartesian coordinates agree.

import numpy as np
import meep as mp
import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt


resolution = 80  # pixels/μm
dpml = 0.5  # thickness of PML
dair = 1.0  # thickness of air padding
L = 6.0  # length of non-PML region
n = 2.4  # refractive index of surrounding medium
wvl = 1.0  # wavelength (in vacuum)

fcen = 1 / wvl  # center frequency of source/monitor

# runtime termination criteria
tol = 1e-8


def extraction_eff_cyl(dmat: float, h: float) -> float:
    """Computes the extraction efficiency of a point dipole embedded
    within a dielectric layer above a lossless ground plane in
    cylindrical coordinates.

    Args:
      dmat: thickness of dielectric layer.
      h: height of dipole above ground plane as fraction of dmat.
    """
    sr = L + dpml
    sz = dmat + dair + dpml
    cell_size = mp.Vector3(sr, 0, sz)

    boundary_layers = [
        mp.PML(dpml, direction=mp.R),
        mp.PML(dpml, direction=mp.Z, side=mp.High),
    ]

    src_cmpt = mp.Er
    src_pt = mp.Vector3(0, 0, -0.5 * sz + h * dmat)
    sources = [
        mp.Source(
            src=mp.GaussianSource(fcen, fwidth=0.1 * fcen),
            component=src_cmpt,
            center=src_pt,
        )
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
        m=-1,
        boundary_layers=boundary_layers,
        sources=sources,
        geometry=geometry,
    )

    flux_air = sim.add_flux(
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
        until_after_sources=mp.stop_when_fields_decayed(20, src_cmpt, src_pt, tol),
    )

    out_flux = mp.get_fluxes(flux_air)[0]
    dV = np.pi / (resolution**3)
    total_flux = -np.real(sim.ldos_Fdata[0] * np.conj(sim.ldos_Jdata[0])) * dV
    ext_eff = out_flux / total_flux
    print(f"extraction efficiency (cyl):, " f"{dmat:.4f}, {h:.4f}, {ext_eff:.6f}")

    return ext_eff


def extraction_eff_3D(dmat: float, h: float) -> float:
    """Computes the extraction efficiency of a point dipole embedded
    within a dielectric layer above a lossless ground plane in
    3D Cartesian coordinates.

    Args:
      dmat: thickness of dielectric layer.
      h: height of dipole above ground plane as fraction of dmat.
    """
    sxy = L + 2 * dpml
    sz = dmat + dair + dpml
    cell_size = mp.Vector3(sxy, sxy, sz)

    symmetries = [mp.Mirror(direction=mp.X, phase=-1), mp.Mirror(direction=mp.Y)]

    boundary_layers = [
        mp.PML(dpml, direction=mp.X),
        mp.PML(dpml, direction=mp.Y),
        mp.PML(dpml, direction=mp.Z, side=mp.High),
    ]

    src_cmpt = mp.Ex
    src_pt = mp.Vector3(0, 0, -0.5 * sz + h * dmat)
    sources = [
        mp.Source(
            src=mp.GaussianSource(fcen, fwidth=0.1 * fcen),
            component=src_cmpt,
            center=src_pt,
        )
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
        boundary_layers=boundary_layers,
        sources=sources,
        geometry=geometry,
        symmetries=symmetries,
    )

    flux_air = sim.add_flux(
        fcen,
        0,
        1,
        mp.FluxRegion(
            center=mp.Vector3(0, 0, 0.5 * sz - dpml),
            size=mp.Vector3(L, L, 0),
        ),
        mp.FluxRegion(
            center=mp.Vector3(0.5 * L, 0, 0.5 * sz - dpml - 0.5 * dair),
            size=mp.Vector3(0, L, dair),
        ),
        mp.FluxRegion(
            center=mp.Vector3(-0.5 * L, 0, 0.5 * sz - dpml - 0.5 * dair),
            size=mp.Vector3(0, L, dair),
            weight=-1.0,
        ),
        mp.FluxRegion(
            center=mp.Vector3(0, 0.5 * L, 0.5 * sz - dpml - 0.5 * dair),
            size=mp.Vector3(L, 0, dair),
        ),
        mp.FluxRegion(
            center=mp.Vector3(0, -0.5 * L, 0.5 * sz - dpml - 0.5 * dair),
            size=mp.Vector3(L, 0, dair),
            weight=-1.0,
        ),
    )

    sim.run(
        mp.dft_ldos(fcen, 0, 1),
        until_after_sources=mp.stop_when_fields_decayed(20, src_cmpt, src_pt, tol),
    )

    out_flux = mp.get_fluxes(flux_air)[0]
    dV = 1 / (resolution**3)
    total_flux = -np.real(sim.ldos_Fdata[0] * np.conj(sim.ldos_Jdata[0])) * dV
    ext_eff = out_flux / total_flux
    print(f"extraction efficiency (3D):, " f"{dmat:.4f}, {h:.4f}, {ext_eff:.6f}")

    return ext_eff


if __name__ == "__main__":
    layer_thickness = 0.7 * wvl / n
    dipole_height = np.linspace(0.1, 0.9, 21)

    exteff_cyl = np.zeros(len(dipole_height))
    exteff_3D = np.zeros(len(dipole_height))
    for j in range(len(dipole_height)):
        exteff_cyl[j] = extraction_eff_cyl(layer_thickness, dipole_height[j])
        exteff_3D[j] = extraction_eff_3D(layer_thickness, dipole_height[j])

    plt.plot(dipole_height, exteff_cyl, "bo-", label="cylindrical")
    plt.plot(dipole_height, exteff_3D, "ro-", label="3D Cartesian")
    plt.xlabel(f"height of dipole above ground plane " f"(fraction of layer thickness)")
    plt.ylabel("extraction efficiency")
    plt.legend()

    if mp.am_master():
        plt.savefig("extraction_eff_vs_dipole_height.png", dpi=150, bbox_inches="tight")
