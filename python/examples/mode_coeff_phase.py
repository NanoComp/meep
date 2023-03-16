"""
Verifies that the magnitude and phase of the reflection coefficient of a
total internal reflected (TIR) mode of a flat interface of two lossless
materials given an incident planewave at oblique incidence computed using
the mode-decomposition feature matches the Fresnel equations.

ref: https://meep.readthedocs.io/en/latest/Python_Tutorials/Mode_Decomposition/#phase-of-a-total-internal-reflected-mode
"""

import cmath
from enum import Enum
import math

import meep as mp
import numpy as np

Polarization = Enum("Polarization", "S P")

# refractive indices of materials 1 and 2
n1 = 1.5
n2 = 1.0


def refl_coeff_meep(pol: Polarization, theta: float, L: float) -> complex:
    """Computes the reflection coefficient of a TIR mode using mode
       decomposition.

    Args:
        pol: polarization of the incident planewave (S or P).
        theta: angle of the incident planewave (radians).
        L: position of the mode monitor relative to the flat interface. 0 is
            interface position.
    """
    resolution = 50.0

    # cell size is arbitrary
    sx = 7.0
    sy = 3.0
    dpml = 2.0
    cell_size = mp.Vector3(sx + 2 * dpml, sy, 0)
    pml_layers = [mp.PML(dpml, direction=mp.X)]

    fcen = 1.0  # center frequency
    df = 0.1 * fcen

    # k (in source medium) with correct length
    # plane of incidence is xy
    k = mp.Vector3(n1 * fcen, 0, 0).rotate(mp.Vector3(0, 0, 1), theta)

    # planewave amplitude function (for source)
    def pw_amp(k, x0):
        def _pw_amp(x):
            return cmath.exp(1j * 2 * math.pi * k.dot(x + x0))

        return _pw_amp

    src_pt = mp.Vector3(-0.5 * sx, 0, 0)

    if pol.name == "S":
        src_cmpt = mp.Ez
        eig_parity = mp.ODD_Z
    elif pol.name == "P":
        src_cmpt = mp.Hz
        eig_parity = mp.EVEN_Z
    else:
        raise ValueError("pol must be S or P, only.")

    sources = [
        mp.Source(
            mp.GaussianSource(fcen, fwidth=df),
            component=src_cmpt,
            center=src_pt,
            size=mp.Vector3(0, cell_size.y, 0),
            amp_func=pw_amp(k, src_pt),
        ),
    ]

    sim = mp.Simulation(
        resolution=resolution,
        cell_size=cell_size,
        default_material=mp.Medium(index=n1),
        boundary_layers=pml_layers,
        k_point=k,
        sources=sources,
    )

    # DFT monitor for incident fields
    mode_mon = sim.add_mode_monitor(
        fcen,
        0,
        1,
        mp.FluxRegion(
            center=mp.Vector3(-L, 0, 0),
            size=mp.Vector3(0, cell_size.y, 0),
        ),
    )

    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            50,
            src_cmpt,
            mp.Vector3(-L, 0, 0),
            1e-6,
        ),
    )

    res = sim.get_eigenmode_coefficients(
        mode_mon,
        bands=[1],
        eig_parity=eig_parity,
        kpoint_func=lambda *not_used: k,
        direction=mp.NO_DIRECTION,
    )

    input_mode_coeff = res.alpha[0, 0, 0]
    input_flux_data = sim.get_flux_data(mode_mon)

    sim.reset_meep()

    geometry = [
        mp.Block(
            material=mp.Medium(index=n1),
            center=mp.Vector3(-0.25 * (sx + 2 * dpml), 0, 0),
            size=mp.Vector3(0.5 * (sx + 2 * dpml), mp.inf, mp.inf),
        ),
        mp.Block(
            material=mp.Medium(index=n2),
            center=mp.Vector3(0.25 * (sx + 2 * dpml), 0, 0),
            size=mp.Vector3(0.5 * (sx + 2 * dpml), mp.inf, mp.inf),
        ),
    ]

    sim = mp.Simulation(
        resolution=resolution,
        cell_size=cell_size,
        boundary_layers=pml_layers,
        k_point=k,
        sources=sources,
        geometry=geometry,
    )

    # DFT monitor for reflected fields
    mode_mon = sim.add_mode_monitor(
        fcen,
        0,
        1,
        mp.FluxRegion(
            center=mp.Vector3(-L, 0, 0),
            size=mp.Vector3(0, cell_size.y, 0),
        ),
    )

    sim.load_minus_flux_data(mode_mon, input_flux_data)

    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            50,
            mp.Ez,
            mp.Vector3(-L, 0, 0),
            1e-6,
        ),
    )

    res = sim.get_eigenmode_coefficients(
        mode_mon,
        bands=[1],
        eig_parity=eig_parity,
        kpoint_func=lambda *not_used: k,
        direction=mp.NO_DIRECTION,
    )

    # mode coefficient of reflected planewave
    refl_mode_coeff = res.alpha[0, 0, 1]

    # reflection coefficient
    refl_coeff = refl_mode_coeff / input_mode_coeff

    # apply phase correction factor
    refl_coeff /= cmath.exp(1j * k.x * 2 * math.pi * 2 * L)

    return refl_coeff


def refl_coeff_Fresnel(pol: Polarization, theta: float) -> complex:
    """Computes the reflection coefficient of a TIR mode using the Fresnel
       equations.

    Args:
        pol: polarization of the incident planewave (S or P).
        theta: angle of the incident planewave (degrees).
    """
    if pol.name == "S":
        refl_coeff = (
            math.cos(theta) - ((n2 / n1) ** 2 - math.sin(theta) ** 2) ** 0.5
        ) / (math.cos(theta) + ((n2 / n1) ** 2 - math.sin(theta) ** 2) ** 0.5)
    else:
        refl_coeff = (
            -((n2 / n1) ** 2) * math.cos(theta)
            + ((n2 / n1) ** 2 - math.sin(theta) ** 2) ** 0.5
        ) / (
            (n2 / n1) ** 2 * math.cos(theta)
            + ((n2 / n1) ** 2 - math.sin(theta) ** 2) ** 0.5
        )

    return refl_coeff


if __name__ == "__main__":
    thetas = [54.3, 48.5]  # angle of incident planewave (degrees)
    Ls = [0.4, 1.2]  # position of mode monitor relative to flat interface
    pols = [Polarization.S, Polarization.P]  # polarization of incident planewave

    for pol, theta, L in zip(pols, thetas, Ls):
        theta_rad = np.radians(theta)
        rc_m = refl_coeff_meep(pol, theta_rad, L)
        rc_f = refl_coeff_Fresnel(pol, theta_rad)

        rc_m_str = f"{rc_m.real:.5f}{rc_m.imag:+.5f}j"
        rc_f_str = f"{rc_f.real:.5f}{rc_f.imag:+.5f}j"
        print(
            f"refl-coeff:, {pol.name}, {theta}, {rc_m_str} (Meep), "
            f"{rc_f_str} (Fresnel)"
        )

        mag_m = abs(rc_m)
        mag_f = abs(rc_f)
        err_mag = abs(mag_m - mag_f) / mag_f
        print(
            f"magnitude:, {mag_m:.5f} (Meep), {mag_f:.5f} (Fresnel), "
            f"{err_mag:.5f} (error)"
        )

        phase_m = cmath.phase(rc_m)
        phase_f = cmath.phase(rc_f)
        err_phase = abs(phase_m - phase_f) / abs(phase_f)
        print(
            f"phase:, {phase_m:.5f} (Meep), {phase_f:.5f} (Fresnel), "
            f"{err_phase:.5f} (error)"
        )
