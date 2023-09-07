# Verifies that the total flux from a light-emitting diode (LED) computed
# in cylindrical coordinates using its near fields is equivalent to using
# its far fields via the radiation pattern obtained using a near-to-far field
# transformation.

# tutorial reference: https://meep.readthedocs.io/en/latest/Python_Tutorials/Near_to_Far_Field_Spectra/#radiation-pattern-of-a-light-emitting-diode-led

import math
from typing import Tuple

import numpy as np
import meep as mp
import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt


resolution = 100  # pixels/μm
dpml = 0.5  # thickness of PML
dair = 1.0  # thickness of air padding
L = 10.0  # length of non-PML region
n = 2.4  # refractive index of surrounding medium
wvl = 1.0  # wavelength (in vacuum)

fcen = 1 / wvl  # center frequency of source/monitor

# field decay threshold for runtime termination criteria
tol = 1e-8

# number of angular grid points
npts = 100

# grid of polar angles for evaluating radiated flux in far field
thetas = np.linspace(0, 0.5 * math.pi, npts)

# radius of quarter circle for computing flux in far field
r = 1000 * wvl


def plot_radiation_pattern_polar(Ptheta: np.ndarray):
    """Plots the radiation pattern in polar coordinates.

    Args:
        Ptheta: radial flux of the far fields in polar coordinates.
    """
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    ax.plot(
        thetas,
        Ptheta,
        "b-",
    )
    ax.set_thetalim(0, 0.5 * math.pi)
    ax.grid(True)
    ax.set_rlabel_position(22)
    ax.set_title("radiation pattern in polar coordinates")

    if mp.am_master():
        fig.savefig("led_radpattern_polar.png", dpi=150, bbox_inches="tight")


def plot_radiation_pattern_3d(Ptheta: np.ndarray):
    """Plots the radiation pattern in spherical coordinates.

    Args:
        Ptheta: radial flux of the far fields in polar coordinates.
    """
    phis = np.linspace(0, 2 * np.pi, npts)

    xs = np.zeros((len(thetas), len(phis)))
    ys = np.zeros((len(thetas), len(phis)))
    zs = np.zeros((len(thetas), len(phis)))

    for i, theta in enumerate(thetas):
        for j, phi in enumerate(phis):
            xs[i, j] = Ptheta[i] * np.cos(theta) * np.cos(phi)
            ys[i, j] = Ptheta[i] * np.cos(theta) * np.sin(phi)
            zs[i, j] = Ptheta[i] * np.sin(theta)

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.plot_surface(xs, ys, zs, cmap="inferno")
    ax.set_title("radiation pattern in 3d")

    if mp.am_master():
        fig.savefig(
            "led_radpattern_3d.png",
            dpi=150,
            bbox_inches="tight",
        )


def radiation_pattern(sim: mp.Simulation, n2f_mon: mp.DftNear2Far) -> np.ndarray:
    """Computes the radiation pattern from the far fields.

    Args:
        sim: a `Simulation` object.
        n2f_mon: a `DftNear2Far` object returned by `Simulation.add_near2far`.

    Returns:
        Array of radial Poynting flux, one for each point on the circumference of
        a quarter circle with angular range of [0, π/2] rad. 0 rad is the +r
        direction and π/2 is +z.
    """
    E = np.zeros((npts, 3), dtype=np.complex128)
    H = np.zeros((npts, 3), dtype=np.complex128)
    for n in range(npts):
        ff = sim.get_farfield(
            n2f_mon, mp.Vector3(r * math.cos(thetas[n]), 0, r * math.sin(thetas[n]))
        )
        E[n, :] = [np.conj(ff[j]) for j in range(3)]
        H[n, :] = [ff[j + 3] for j in range(3)]

    Pr = np.real(E[:, 1] * H[:, 2] - E[:, 2] * H[:, 1])
    Pz = np.real(E[:, 0] * H[:, 1] - E[:, 1] * H[:, 0])
    Prz = np.sqrt(np.square(Pr) + np.square(Pz))

    return Prz


def led_total_flux(dmat: float, h: float) -> Tuple[float, float]:
    """Computes the total radiated flux from a point dipole embedded
    within a dielectric layer above a lossless ground plane using
    its near and far fields as separate calculations.

    Args:
        dmat: thickness of dielectric layer.
        h: height of dipole above ground plane as fraction of dmat.

    Returns:
        A 2-tuple of the total flux computed using the near and far fields,
        respectively.
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

    # flux monitor
    flux_mon = sim.add_flux(
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

    # near-field monitor
    n2f_mon = sim.add_near2far(
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

    fig, ax = plt.subplots()
    sim.plot2D(ax=ax)
    if mp.am_master():
        fig.savefig("led_plot2D.png", dpi=150, bbox_inches="tight")

    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            50,
            src_cmpt,
            src_pt,
            tol,
        ),
    )

    flux_near = mp.get_fluxes(flux_mon)[0]

    Ptheta = radiation_pattern(sim, n2f_mon)
    plot_radiation_pattern_polar(r * r * Ptheta)
    plot_radiation_pattern_3d(r * r * Ptheta)

    dtheta = 0.5 * math.pi / npts
    dphi = 2 * math.pi
    flux_far = np.sum(Ptheta * np.cos(thetas)) * r * r * dtheta * dphi

    err = abs(flux_near - flux_far) / flux_near
    print(
        f"total_flux:, {flux_near:.5f} (near), {flux_far:.5f} (far), "
        f"{err:.5f} (error)"
    )

    return flux_near, flux_far


if __name__ == "__main__":
    layer_thickness = 0.7 * wvl / n
    dipole_height = 0.5

    near_flux, far_flux = led_total_flux(layer_thickness, dipole_height)
