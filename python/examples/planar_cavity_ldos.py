# Computes the Purcell enhancement factor of a parallel dipole in a planar
# dielectric cavity with lossless metallic walls. The result is computed in
# cylindrical and 3D coordinates and compared with the analytic theory from:
# I. Abram et al., IEEE J. Quantum Electronics, Vol. 34, pp. 71-76 (1998).
import matplotlib
import numpy as np

import meep as mp

matplotlib.use("agg")
import matplotlib.pyplot as plt

# important note:
# Meep may round the cell dimensions to an integer number
# of pixels which could modify the cavity structure.
resolution = 70  # pixels/μm


dpml = 0.5  # thickness of PML
L = 6.0  # length of non-PML region
n = 2.4  # refractive index of surrounding medium
wvl = 1.0  # wavelength (in vacuum)

fcen = 1 / wvl


def bulk_ldos_cyl():
    sr = L + dpml
    sz = L + 2 * dpml
    cell_size = mp.Vector3(sr, 0, sz)

    pml_layers = [mp.PML(dpml)]

    sources = [
        mp.Source(
            src=mp.GaussianSource(fcen, fwidth=0.2 * fcen),
            component=mp.Er,
            center=mp.Vector3(),
        )
    ]

    sim = mp.Simulation(
        resolution=resolution,
        cell_size=cell_size,
        boundary_layers=pml_layers,
        sources=sources,
        dimensions=mp.CYLINDRICAL,
        m=-1,
        default_material=mp.Medium(index=n),
    )

    sim.run(
        mp.dft_ldos(fcen, 0, 1),
        until_after_sources=mp.stop_when_fields_decayed(20, mp.Er, mp.Vector3(), 1e-6),
    )

    return sim.ldos_data[0]


def cavity_ldos_cyl(sz):
    sr = L + dpml
    cell_size = mp.Vector3(sr, 0, sz)

    pml_layers = [mp.PML(dpml, direction=mp.R)]

    sources = [
        mp.Source(
            src=mp.GaussianSource(fcen, fwidth=0.2 * fcen),
            component=mp.Er,
            center=mp.Vector3(),
        )
    ]

    sim = mp.Simulation(
        resolution=resolution,
        cell_size=cell_size,
        boundary_layers=pml_layers,
        sources=sources,
        dimensions=mp.CYLINDRICAL,
        m=-1,
        default_material=mp.Medium(index=n),
    )

    sim.run(
        mp.dft_ldos(fcen, 0, 1),
        until_after_sources=mp.stop_when_fields_decayed(20, mp.Er, mp.Vector3(), 1e-6),
    )

    return sim.ldos_data[0]


def bulk_ldos_3D():
    s = L + 2 * dpml
    cell_size = mp.Vector3(s, s, s)

    pml_layers = [mp.PML(dpml)]

    sources = [
        mp.Source(
            src=mp.GaussianSource(fcen, fwidth=0.2 * fcen),
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
        resolution=resolution,
        cell_size=cell_size,
        boundary_layers=pml_layers,
        sources=sources,
        symmetries=symmetries,
        default_material=mp.Medium(index=n),
    )

    sim.run(
        mp.dft_ldos(fcen, 0, 1),
        until_after_sources=mp.stop_when_fields_decayed(20, mp.Ex, mp.Vector3(), 1e-6),
    )

    return sim.ldos_data[0]


def cavity_ldos_3D(sz):
    sxy = L + 2 * dpml
    cell_size = mp.Vector3(sxy, sxy, sz)

    boundary_layers = [mp.PML(dpml, direction=mp.X), mp.PML(dpml, direction=mp.Y)]

    sources = [
        mp.Source(
            src=mp.GaussianSource(fcen, fwidth=0.2 * fcen),
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
        resolution=resolution,
        cell_size=cell_size,
        boundary_layers=boundary_layers,
        sources=sources,
        symmetries=symmetries,
        default_material=mp.Medium(index=n),
    )

    sim.run(
        mp.dft_ldos(fcen, 0, 1),
        until_after_sources=mp.stop_when_fields_decayed(20, mp.Ex, mp.Vector3(), 1e-6),
    )

    return sim.ldos_data[0]


if __name__ == "__main__":
    ldos_bulk_cyl = bulk_ldos_cyl()
    ldos_bulk_3D = bulk_ldos_3D()

    # units of wavelength in medium
    cavity_thickness = np.arange(0.50, 2.55, 0.05)

    gap = cavity_thickness * wvl / n

    ldos_cavity_cyl = np.zeros(len(cavity_thickness))
    ldos_cavity_3D = np.zeros(len(cavity_thickness))
    for idx, g in enumerate(gap):
        ldos_cavity_cyl[idx] = cavity_ldos_cyl(g)
        ldos_cavity_3D[idx] = cavity_ldos_3D(g)
        print(
            "purcell-enh:, {:.3f}, "
            "{:.6f} (cyl.), {:.6f} (3D)".format(
                cavity_thickness[idx],
                ldos_cavity_cyl[idx] / ldos_bulk_cyl,
                ldos_cavity_3D[idx] / ldos_bulk_3D,
            )
        )

    # Purcell enhancement factor (relative to bulk medium)
    pe_meep_cyl = ldos_cavity_cyl / ldos_bulk_cyl
    pe_meep_3D = ldos_cavity_3D / ldos_bulk_3D

    # equation 7 of reference
    pe_theory = 3 * np.fix(cavity_thickness + 0.5) / (4 * cavity_thickness) + (
        4 * np.power(np.fix(cavity_thickness + 0.5), 3) - np.fix(cavity_thickness + 0.5)
    ) / (16 * np.power(cavity_thickness, 3))

    if mp.am_master():
        plt.plot(cavity_thickness, pe_meep_3D, "b-", label="Meep (3D)")
        plt.plot(cavity_thickness, pe_meep_cyl, "r-", label="Meep (cylin.)")
        plt.plot(cavity_thickness, pe_theory, "g-", label="theory")
        plt.plot(cavity_thickness, np.ones(len(cavity_thickness)), "k--")
        plt.xlabel(r"cavity thickness, $nL/\lambda$")
        plt.ylabel("Purcell enhancement factor")
        plt.title(
            "planar point dipole at λ=1.0 μm in a planar cavity\n"
            "with n=2.4 and lossless metallic walls"
        )
        plt.axis([0.5, 2.5, 0.4, 3.1])
        plt.legend()
        plt.savefig("cavity_purcell_factor_vs_thickness.png", bbox_inches="tight")
