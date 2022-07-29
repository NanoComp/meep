import argparse
import math

import meep as mp


def main(args):
    resolution = args.res

    dpml = 1.0  # PML thickness
    sz = 10 + 2 * dpml  # size of computational cell (without PMLs)
    cell_size = mp.Vector3(0, 0, sz)
    pml_layers = [mp.PML(dpml)]

    wvl_min = 0.4  # min wavelength
    wvl_max = 0.8  # max wavelength
    fmin = 1 / wvl_max  # min frequency
    fmax = 1 / wvl_min  # max frequency
    fcen = 0.5 * (fmin + fmax)  # center frequency
    df = fmax - fmin  # frequency width
    nfreq = 50  # number of frequency bins

    # rotation angle (in degrees) of source: CCW around Y axis, 0 degrees along +Z axis
    theta_r = math.radians(args.theta)

    # plane of incidence is XZ
    k = mp.Vector3(math.sin(theta_r), 0, math.cos(theta_r)).scale(fmin)

    # if normal incidence, force number of dimensions to be 1
    dimensions = 1 if theta_r == 0 else 3
    sources = [
        mp.Source(
            mp.GaussianSource(fcen, fwidth=df),
            component=mp.Ex,
            center=mp.Vector3(0, 0, -0.5 * sz + dpml),
        )
    ]

    sim = mp.Simulation(
        cell_size=cell_size,
        boundary_layers=pml_layers,
        sources=sources,
        k_point=k,
        dimensions=dimensions,
        resolution=resolution,
    )

    refl_fr = mp.FluxRegion(center=mp.Vector3(0, 0, -0.25 * sz))
    refl = sim.add_flux(fcen, df, nfreq, refl_fr)

    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            50, mp.Ex, mp.Vector3(0, 0, -0.5 * sz + dpml), 1e-9
        )
    )

    empty_flux = mp.get_fluxes(refl)
    empty_data = sim.get_flux_data(refl)
    sim.reset_meep()

    # add a block with n=3.5 for the air-dielectric interface
    geometry = [
        mp.Block(
            size=mp.Vector3(mp.inf, mp.inf, 0.5 * sz),
            center=mp.Vector3(0, 0, 0.25 * sz),
            material=mp.Medium(index=3.5),
        )
    ]

    sim = mp.Simulation(
        cell_size=cell_size,
        geometry=geometry,
        boundary_layers=pml_layers,
        sources=sources,
        k_point=k,
        dimensions=dimensions,
        resolution=resolution,
    )

    refl = sim.add_flux(fcen, df, nfreq, refl_fr)
    sim.load_minus_flux_data(refl, empty_data)

    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            50, mp.Ex, mp.Vector3(0, 0, -0.5 * sz + dpml), 1e-9
        )
    )

    refl_flux = mp.get_fluxes(refl)
    freqs = mp.get_flux_freqs(refl)

    for i in range(nfreq):
        print(
            f"refl:, {k.x}, {1 / freqs[i]}, {math.degrees(math.asin(k.x/freqs[i]))}, {-refl_flux[i] / empty_flux[i]}"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-res", type=int, default=200, help="resolution (default: 200 pixels/um)"
    )
    parser.add_argument(
        "-theta",
        type=float,
        default=0,
        help="angle of incident planewave (default: 0 degrees)",
    )
    args = parser.parse_args()
    main(args)
