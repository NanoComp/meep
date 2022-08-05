# Computes the radiation pattern of a dipole antenna
# positioned a given height above a perfect-electric
# conductor (PEC) ground plane and compares the result
# to analytic theory.
import math

import matplotlib
import numpy as np

import meep as mp

matplotlib.use("agg")
import matplotlib.pyplot as plt

resolution = 200  # pixels/um
n = 1.2  # refractive index of surrounding medium
h = 1.25  # height of antenna (point dipole source) above ground plane
wvl = 0.65  # vacuum wavelength
r = 1000 * wvl  # radius of far-field circle
npts = 50  # number of points in [0,pi/2) range of angles

angles = 0.5 * math.pi / npts * np.arange(npts)


def radial_flux(sim, nearfield_box, r):
    E = np.zeros((npts, 3), dtype=np.complex128)
    H = np.zeros((npts, 3), dtype=np.complex128)

    for n in range(npts):
        ff = sim.get_farfield(
            nearfield_box, mp.Vector3(r * math.sin(angles[n]), r * math.cos(angles[n]))
        )
        E[n, :] = [np.conj(ff[j]) for j in range(3)]
        H[n, :] = [ff[j + 3] for j in range(3)]

    Px = np.real(E[:, 1] * H[:, 2] - E[:, 2] * H[:, 1])  # Ey*Hz-Ez*Hy
    Py = np.real(E[:, 2] * H[:, 0] - E[:, 0] * H[:, 2])  # Ez*Hx-Ex*Hz
    return np.sqrt(np.square(Px) + np.square(Py))


def free_space_radiation(src_cmpt):
    sxy = 4
    dpml = 1
    cell_size = mp.Vector3(sxy + 2 * dpml, sxy + 2 * dpml)
    pml_layers = [mp.PML(dpml)]

    fcen = 1 / wvl
    sources = [
        mp.Source(
            src=mp.GaussianSource(fcen, fwidth=0.2 * fcen),
            center=mp.Vector3(),
            component=src_cmpt,
        )
    ]

    if src_cmpt == mp.Hz:
        symmetries = [mp.Mirror(mp.X, phase=-1), mp.Mirror(mp.Y, phase=-1)]
    elif src_cmpt == mp.Ez:
        symmetries = [mp.Mirror(mp.X, phase=+1), mp.Mirror(mp.Y, phase=+1)]
    else:
        symmetries = []

    sim = mp.Simulation(
        cell_size=cell_size,
        resolution=resolution,
        sources=sources,
        symmetries=symmetries,
        boundary_layers=pml_layers,
        default_material=mp.Medium(index=n),
    )

    nearfield_box = sim.add_near2far(
        fcen,
        0,
        1,
        mp.Near2FarRegion(center=mp.Vector3(0, +0.5 * sxy), size=mp.Vector3(sxy, 0)),
        mp.Near2FarRegion(
            center=mp.Vector3(0, -0.5 * sxy), size=mp.Vector3(sxy, 0), weight=-1
        ),
        mp.Near2FarRegion(center=mp.Vector3(+0.5 * sxy, 0), size=mp.Vector3(0, sxy)),
        mp.Near2FarRegion(
            center=mp.Vector3(-0.5 * sxy, 0), size=mp.Vector3(0, sxy), weight=-1
        ),
    )

    sim.run(until_after_sources=mp.stop_when_dft_decayed())

    return radial_flux(sim, nearfield_box, r)


def pec_ground_plane_radiation(src_cmpt=mp.Hz):
    L = 8.0  # length of non-PML region
    dpml = 1.0  # thickness of PML
    sxy = dpml + L + dpml
    cell_size = mp.Vector3(sxy, sxy, 0)
    boundary_layers = [mp.PML(dpml)]

    fcen = 1 / wvl

    # The near-to-far field transformation feature only supports
    # homogeneous media which means it cannot explicitly take into
    # account the ground plane. As a workaround, we use two antennas
    # of opposite sign surrounded by a single near2far box which
    # encloses both antennas. We then use an odd mirror symmetry to
    # divide the computational cell in half which is effectively
    # equivalent to a PEC boundary condition on one side.
    # Note: This setup means that the radiation pattern can only
    # be measured in the top half above the dipole.
    sources = [
        mp.Source(
            src=mp.GaussianSource(fcen, fwidth=0.2 * fcen),
            component=src_cmpt,
            center=mp.Vector3(0, +h),
        ),
        mp.Source(
            src=mp.GaussianSource(fcen, fwidth=0.2 * fcen),
            component=src_cmpt,
            center=mp.Vector3(0, -h),
            amplitude=-1 if src_cmpt == mp.Ez else +1,
        ),
    ]

    symmetries = [
        mp.Mirror(direction=mp.X, phase=+1 if src_cmpt == mp.Ez else -1),
        mp.Mirror(direction=mp.Y, phase=-1 if src_cmpt == mp.Ez else +1),
    ]

    sim = mp.Simulation(
        resolution=resolution,
        cell_size=cell_size,
        boundary_layers=boundary_layers,
        sources=sources,
        symmetries=symmetries,
        default_material=mp.Medium(index=n),
    )

    nearfield_box = sim.add_near2far(
        fcen,
        0,
        1,
        mp.Near2FarRegion(
            center=mp.Vector3(0, 2 * h), size=mp.Vector3(2 * h, 0), weight=+1
        ),
        mp.Near2FarRegion(
            center=mp.Vector3(0, -2 * h), size=mp.Vector3(2 * h, 0), weight=-1
        ),
        mp.Near2FarRegion(
            center=mp.Vector3(h, 0), size=mp.Vector3(0, 4 * h), weight=+1
        ),
        mp.Near2FarRegion(
            center=mp.Vector3(-h, 0), size=mp.Vector3(0, 4 * h), weight=-1
        ),
    )

    sim.plot2D()
    plt.savefig("antenna_pec_ground_plane.png", bbox_inches="tight")

    sim.run(until_after_sources=mp.stop_when_dft_decayed())

    return radial_flux(sim, nearfield_box, r)


if __name__ == "__main__":
    src_cmpt = mp.Ez  # TM/P: Hz or TE/S: Ez
    Pr_fsp = free_space_radiation(src_cmpt)
    Pr_pec = pec_ground_plane_radiation(src_cmpt)

    # The radiation pattern of a two-element antenna
    # array is equivalent to the radiation pattern of
    # a single antenna multiplied by its array factor
    # as derived in Section 6.2 "Two-Element Array" of
    # Antenna Theory: Analysis and Design, Fourth Edition
    # (2016) by C.A. Balanis.
    k = 2 * np.pi / (wvl / n)  # wavevector in free space
    Pr_theory = np.zeros(
        npts,
    )
    for i, ang in enumerate(angles):
        Pr_theory[i] = Pr_fsp[i] * 2 * np.sin(k * h * np.cos(ang))

    Pr_pec_norm = Pr_pec / np.max(Pr_pec)
    Pr_theory_norm = (Pr_theory / max(Pr_theory)) ** 2

    plt.figure()
    plt.plot(np.degrees(angles), Pr_pec_norm, "b-", label="Meep")
    plt.plot(np.degrees(angles), Pr_theory_norm, "r-", label="theory")
    plt.xlabel("angle (degrees)")
    plt.ylabel("radial flux (normalized by maximum flux)")
    plt.title(
        f"antenna with {'E' if src_cmpt==mp.Ez else r'H'}$_z$ polarization above PEC ground plane"
    )

    plt.axis([0, 90, 0, 1.0])
    plt.legend()
    plt.savefig("radiation_pattern.png", bbox_inches="tight")

    print(f"norm:, {np.linalg.norm(Pr_pec_norm - Pr_theory_norm):.6f}")
