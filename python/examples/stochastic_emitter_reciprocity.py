"""A demonstration of computing the flux emitted in a single direction from a
   collection of dipole emitters in an LED-like structure using reciprocity.
"""

from typing import List
import numpy as np
import matplotlib.pyplot as plt
import meep as mp
from meep.materials import Ag

resolution = 200  # pixels/μm

nfreq = 100  # number of frequencies
ndipole = 10  # number of point dipoles in forward simulation

fcen = 1.0  # center frequency of Gaussian source/monitors
df = 0.2  # frequency bandwidth of source/monitors

dpml = 1.0  # PML thickness
dair = 2.0  # air padding thickness
hrod = 0.7  # grating height
wrod = 0.5  # graing width
dsub = 5.0  # substrate thickness
dAg = 0.5  # Ag back reflecter thickness

sx = 1.1
sy = dpml + dair + hrod + dsub + dAg

cell_size = mp.Vector3(sx, sy)

pml_layers = [mp.PML(direction=mp.Y, thickness=dpml, side=mp.High)]


def substrate_geometry(is_textured: bool):
    """Returns the geometry of the LED-like structure.

    Args:
      is_textured: whether the substrate is textured or not.
    """
    geometry = [
        mp.Block(
            material=mp.Medium(index=3.45),
            center=mp.Vector3(0, 0.5 * sy - dpml - dair - hrod - 0.5 * dsub),
            size=mp.Vector3(mp.inf, dsub, mp.inf),
        ),
        mp.Block(
            material=Ag,
            center=mp.Vector3(0, -0.5 * sy + 0.5 * dAg),
            size=mp.Vector3(mp.inf, dAg, mp.inf),
        ),
    ]

    if is_textured:
        geometry.append(
            mp.Block(
                material=mp.Medium(index=3.45),
                center=mp.Vector3(0, 0.5 * sy - dpml - dair - 0.5 * hrod),
                size=mp.Vector3(wrod, hrod, mp.inf),
            )
        )

    return geometry


def forward(n: int, rt: int, is_textured: bool) -> [List, np.ndarray]:
    """Computes the Poynting flux in the +y direction in air
    given a point dipole source positioned somewhere along a
    line in the middle of the high-index substrate.

    Args:
      n: n'th position along a line of equally spaced dipoles.
      rt: runtime of simulation after the source has turned off
          in units of nfreq/df.
      is_textured: whether the substrate is textured or not.

    Returns:
      The frequency and Poynting flux spectra.
    """
    sources = [
        mp.Source(
            mp.GaussianSource(fcen, fwidth=df),
            component=mp.Ez,
            center=mp.Vector3(
                sx * (-0.5 + n / ndipole),
                -0.5 * sy + dAg + 0.5 * dsub,
            ),
        )
    ]

    geometry = substrate_geometry(is_textured)

    sim = mp.Simulation(
        cell_size=cell_size,
        resolution=resolution,
        k_point=mp.Vector3(),
        boundary_layers=pml_layers,
        geometry=geometry,
        sources=sources,
    )

    flux_mon = sim.add_flux(
        fcen,
        df,
        nfreq,
        mp.FluxRegion(center=mp.Vector3(0, 0.5 * sy - dpml), size=mp.Vector3(sx)),
    )

    run_time = rt * nfreq / df
    sim.run(until_after_sources=run_time)

    res = sim.get_eigenmode_coefficients(flux_mon, [1], eig_parity=mp.ODD_Z)

    flux = np.abs(res.alpha[0, :, 0]) ** 2
    freqs = mp.get_flux_freqs(flux_mon)

    return freqs, flux


def backward(rt: int, is_textured: bool) -> [List, np.ndarray]:
    """Computes the Poynting flux spectrum of the dipole emission using
       an overlap integral of the DFT fields from a line monitor in the
       high-index substrate given a planewave source in air propagating
       in the -y direction.

    Args:
      rt: runtime of simulation after the source has turned off
          in units of nfreq/df.
      is_textured: whether the substrate is textured or not.

    Returns:
      The frequency and Poynting flux spectra.
    """
    sources = [
        mp.Source(
            mp.GaussianSource(fcen, fwidth=df),
            component=mp.Ez,
            center=mp.Vector3(0, 0.5 * sy - dpml),
            size=mp.Vector3(sx, 0),
        )
    ]

    geometry = substrate_geometry(is_textured)

    sim = mp.Simulation(
        cell_size=cell_size,
        resolution=resolution,
        k_point=mp.Vector3(),
        boundary_layers=pml_layers,
        geometry=geometry,
        sources=sources,
    )

    dft_mon = sim.add_dft_fields(
        [mp.Ez],
        fcen,
        df,
        nfreq,
        center=mp.Vector3(0, -0.5 * sy + dAg + 0.5 * dsub),
        size=mp.Vector3(sx),
    )

    run_time = rt * nfreq / df
    sim.run(until_after_sources=run_time)

    freqs = mp.get_flux_freqs(dft_mon)

    abs_flux = np.zeros(nfreq)
    for nf in range(nfreq):
        dft_ez = sim.get_dft_array(dft_mon, mp.Ez, nf)
        abs_flux[nf] = np.sum(np.abs(dft_ez) ** 2)

    return freqs, abs_flux


if __name__ == "__main__":
    fwd_flat_flux = np.zeros((nfreq, ndipole))
    fwd_text_flux = np.zeros((nfreq, ndipole))
    for d in range(ndipole):
        fwd_freqs, fwd_flat_flux[:, d] = forward(d, 2, False)
        _, fwd_text_flux[:, d] = forward(d, 4, True)

    fwd_norm_flux = np.mean(fwd_text_flux, axis=1) / np.mean(fwd_flat_flux, axis=1)

    bwd_freqs, bwd_flat_flux = backward(2, False)
    _, bwd_text_flux = backward(4, True)
    bwd_norm_flux = bwd_text_flux / bwd_flat_flux

    plt.figure()
    plt.semilogy(fwd_freqs, fwd_norm_flux, "b-", label="forward")
    plt.semilogy(bwd_freqs, bwd_norm_flux, "r-", label="backward")
    plt.xlabel("frequency")
    plt.ylabel("normalized flux")
    plt.legend()

    if mp.am_master():
        plt.savefig(
            "forward_vs_backward_flux_spectrum.png",
            bbox_inches="tight",
            dpi=150,
        )
