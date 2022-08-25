### compute the transmitted diffraction orders of a binary grating using mode decomposition
### based on two different methods: (1) MPB eigensolver and (2) DiffractedPlanewave object.
### Also, verify that the total power in all the orders is equivalent to the Poynting flux.
### for normal incidence, compute only positive diff. orders (total transmittance <= 0.50)
### for oblique incidence, compute ALL diff. orders (total transmittance <= 1.00)
import cmath
import math

import numpy as np

import meep as mp


def binary_grating_diffraction(gp, gh, gdc, theta):

    resolution = 50  # pixels/μm

    dpml = 1.0  # PML thickness
    dsub = 3.0  # substrate thickness
    dpad = 3.0  # length of padding between grating and PML

    sx = dpml + dsub + gh + dpad + dpml
    sy = gp

    cell_size = mp.Vector3(sx, sy, 0)
    pml_layers = [mp.PML(thickness=dpml, direction=mp.X)]

    wvl = 0.5  # center wavelength
    fcen = 1 / wvl  # center frequency
    df = 0.05 * fcen  # frequency width

    ng = 1.5
    glass = mp.Medium(index=ng)

    # rotation angle of incident planewave; counter clockwise (CCW) about Z axis, 0 degrees along +X axis
    theta_in = math.radians(theta)

    eig_parity = mp.EVEN_Z

    # k (in source medium) with correct length (plane of incidence: XY)
    k = mp.Vector3(fcen * ng).rotate(mp.Vector3(z=1), theta_in)

    symmetries = []
    if theta_in == 0:
        k = mp.Vector3()
        eig_parity += mp.ODD_Y
        symmetries = [mp.Mirror(direction=mp.Y, phase=-1)]

    def pw_amp(k, x0):
        def _pw_amp(x):
            return cmath.exp(1j * 2 * math.pi * k.dot(x + x0))

        return _pw_amp

    src_pt = mp.Vector3(-0.5 * sx + dpml, 0, 0)
    sources = [
        mp.Source(
            mp.GaussianSource(fcen, fwidth=df),
            component=mp.Hz,
            center=src_pt,
            size=mp.Vector3(0, sy, 0),
            amp_func=pw_amp(k, src_pt),
        )
    ]

    sim = mp.Simulation(
        resolution=resolution,
        cell_size=cell_size,
        boundary_layers=pml_layers,
        k_point=k,
        default_material=glass,
        sources=sources,
        symmetries=symmetries,
    )

    tran_pt = mp.Vector3(0.5 * sx - dpml, 0, 0)
    tran_mon = sim.add_flux(
        fcen, 0, 1, mp.FluxRegion(center=tran_pt, size=mp.Vector3(0, sy, 0))
    )

    sim.run(until_after_sources=50)

    input_flux = mp.get_fluxes(tran_mon)

    sim.reset_meep()

    geometry = [
        mp.Block(
            material=glass,
            size=mp.Vector3(dpml + dsub, mp.inf, mp.inf),
            center=mp.Vector3(-0.5 * sx + 0.5 * (dpml + dsub), 0, 0),
        ),
        mp.Block(
            material=glass,
            size=mp.Vector3(gh, gdc * gp, mp.inf),
            center=mp.Vector3(-0.5 * sx + dpml + dsub + 0.5 * gh, 0, 0),
        ),
    ]

    sim = mp.Simulation(
        resolution=resolution,
        cell_size=cell_size,
        boundary_layers=pml_layers,
        geometry=geometry,
        k_point=k,
        sources=sources,
        symmetries=symmetries,
    )

    tran_mon = sim.add_mode_monitor(
        fcen, 0, 1, mp.FluxRegion(center=tran_pt, size=mp.Vector3(0, sy, 0))
    )

    sim.run(until_after_sources=100)

    # number of (non-evanescent) transmitted orders
    nm_t = np.floor((fcen - k.y) * gp) - np.ceil((-fcen - k.y) * gp)
    if theta_in == 0:
        nm_t = nm_t / 2
    nm_t = int(nm_t) + 1

    bands = range(1, nm_t + 1)

    if theta_in == 0:
        orders = range(0, nm_t)
    else:
        orders = range(
            int(np.ceil((-fcen - k.y) * gp)), int(np.floor((fcen - k.y) * gp)) + 1
        )

    eig_sum = 0
    dp_sum = 0

    for band, order in zip(bands, orders):
        res = sim.get_eigenmode_coefficients(tran_mon, [band], eig_parity=eig_parity)
        if res is not None:
            tran_eig = abs(res.alpha[0, 0, 0]) ** 2 / input_flux[0]
            if theta_in == 0:
                tran_eig = 0.5 * tran_eig
        else:
            tran_eig = 0
        eig_sum += tran_eig

        res = sim.get_eigenmode_coefficients(
            tran_mon, mp.DiffractedPlanewave((0, order, 0), mp.Vector3(0, 1, 0), 0, 1)
        )
        if res is not None:
            tran_dp = abs(res.alpha[0, 0, 0]) ** 2 / input_flux[0]
            if (theta_in == 0) and (order == 0):
                tran_dp = 0.5 * tran_dp
        else:
            tran_dp = 0
        dp_sum += tran_dp

        if theta_in == 0:
            err = abs(tran_eig - tran_dp) / tran_eig
            print(
                "tran:, {:2d}, {:.8f}, {:2d}, {:.8f}, {:.8f}".format(
                    band, tran_eig, order, tran_dp, err
                )
            )
        else:
            print(
                "tran:, {:2d}, {:.8f}, {:2d}, {:.8f}".format(
                    band, tran_eig, order, tran_dp
                )
            )

    flux = mp.get_fluxes(tran_mon)
    t_flux = flux[0] / input_flux[0]
    if theta_in == 0:
        t_flux = 0.5 * t_flux

    err = abs(dp_sum - t_flux) / t_flux
    print(f"flux:, {eig_sum:.8f}, {dp_sum:.8f}, {t_flux:.8f}, {err:.8f}")


if __name__ == "__main__":
    binary_grating_diffraction(2.6, 0.4, 0.3, 0)
    binary_grating_diffraction(3.7, 0.6, 0.4, 13.5)
