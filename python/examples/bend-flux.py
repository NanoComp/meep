# From the Meep tutorial: transmission around a 90-degree waveguide bend in 2d.
from __future__ import division

import argparse
import meep as mp


def main(args):
    sx = 16  # size of cell in X direction
    sy = 32  # size of cell in Y direction
    cell = mp.Vector3(sx, sy, 0)

    pad = 4  # padding distance between waveguide and cell edge
    w = 1  # width of waveguide

    wvg_ycen = -0.5 * (sy - w - (2 * pad))  # y center of horiz. wvg
    wvg_xcen = 0.5 * (sx - w - (2 * pad))  # x center of vert. wvg

    if args.no_bend:
        geometry = [mp.Block(mp.Vector3(1e20, w, 1e20), center=mp.Vector3(0, wvg_ycen),
                             material=mp.Medium(epsilon=12))]
    else:
        geometry = [mp.Block(mp.Vector3(sx - pad, w, 1e20), center=mp.Vector3(-0.5 * pad, wvg_ycen),
                             material=mp.Medium(epsilon=12)),
                    mp.Block(mp.Vector3(w, sy - pad, 1e20), center=mp.Vector3(wvg_xcen, 0.5 * pad),
                             material=mp.Medium(epsilon=12))]

    fcen = 0.15  # pulse center frequency
    df = 0.1  # pulse width (in frequency)

    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ez,
                         center=mp.Vector3(1 + (-0.5 * sx), wvg_ycen), size=mp.Vector3(0, w))]

    pml_layers = [mp.PML(1.0)]
    resolution = 10

    nfreq = 100  # number of frequencies at which to compute flux

    sim = mp.Simulation(cell_size=cell,
                        boundary_layers=pml_layers,
                        geometry=geometry,
                        sources=sources,
                        resolution=resolution)

    if args.no_bend:
        fr = mp.FluxRegion(center=mp.Vector3((sx / 2) - 1.5, wvg_ycen), size=mp.Vector3(0, w * 2))
    else:
        fr = mp.FluxRegion(center=mp.Vector3(wvg_xcen, (sy / 2) - 1.5), size=mp.Vector3(w * 2, 0))

    trans = sim.add_flux(fcen, df, nfreq, fr)

    refl_fr = mp.FluxRegion(center=mp.Vector3((-0.5 * sx) + 1.5, wvg_ycen),
                            size=mp.Vector3(0, w * 2))

    # reflected flux
    refl = sim.add_flux(fcen, df, nfreq, refl_fr)

    # for normal run, load negated fields to subtract incident from refl. fields
    if not args.no_bend:
        sim.load_minus_flux('refl-flux', refl)

    if args.no_bend:
        pt = mp.Vector3((sx / 2) - 1.5, wvg_ycen)
    else:
        pt = mp.Vector3(wvg_xcen, (sy / 2) - 1.5)

    sim.run(mp.at_beginning(mp.output_epsilon),
            until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, pt, 1e-3))

    # for normalization run, save flux fields for refl. plane
    if args.no_bend:
        sim.save_flux('refl-flux', refl)

    sim.display_fluxes(trans, refl)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--no_bend', action='store_true', default=False,
                        help="Straight waveguide without bend.")
    args = parser.parse_args()
    main(args)
