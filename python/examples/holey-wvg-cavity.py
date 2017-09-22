# Meep Tutorial: TE transmission and reflection through a cavity
# formed by a periodic sequence of holes in a dielectric waveguide,
# with a defect formed by a larger spacing between one pair of holes.

# This structure is based on one analyzed in:
#    S. Fan, J. N. Winn, A. Devenyi, J. C. Chen, R. D. Meade, and
#    J. D. Joannopoulos, "Guided and defect modes in periodic dielectric
#    waveguides," J. Opt. Soc. Am. B, 12 (7), 1267-1272 (1995).
from __future__ import division

import argparse
import meep as mp


def main(args):
    # Some parameters to describe the geometry:
    eps = 13  # dielectric constant of waveguide
    w = 1.2  # width of waveguide
    r = 0.36  # radius of holes
    d = 1.4  # defect spacing (ordinary spacing = 1)
    N = 3  # number of holes on either side of defect

    # The cell dimensions
    sy = 6  # size of cell in y direction (perpendicular to wvg.)
    pad = 2  # padding between last hole and PML edge
    dpml = 1  # PML thickness

    sx = (2 * (pad + dpml + N)) + d - 1  # size of cell in x direction

    cell = mp.Vector3(sx, sy, 0)

    blk = mp.Block(size=mp.Vector3(1e20, w, 1e20),
                   material=mp.Medium(epsilon=eps))

    geometry = [blk]

    for i in range(3):
        geometry.append(mp.Cylinder(r, center=mp.Vector3(d / 2 + i)))

    for i in range(3):
        geometry.append(mp.Cylinder(r, center=mp.Vector3(d / -2 - i)))

    fcen = 0.25  # pulse center frequency
    df = 0.2  # pulse width (in frequency)

    nfreq = 500  # number of frequencies at which to compute flux

    sim = mp.Simulation(cell_size=cell,
                        geometry=geometry,
                        sources=[],
                        boundary_layers=[mp.Pml(dpml)],
                        resolution=20)

    if args.resonant_modes:
        sim.sources.append(mp.Source(mp.GaussianSource(fcen, fwidth=df), mp.Hz, mp.Vector3()))

        sim.symmetries.append(mp.Mirror(mp.Y, phase=-1))
        sim.symmetries.append(mp.Mirror(mp.X, phase=-1))

        sim.run(mp.at_beginning(mp.output_epsilon),
                mp.after_sources(mp.Harminv(mp.Hz, mp.Vector3(), fcen, df)),
                until_after_sources=400)

        sim.run(mp.at_every(1 / fcen / 20, mp.output_hfield_z), until=1 / fcen)

    else:
        sim.sources.append(mp.Source(mp.GaussianSource(fcen, fwidth=df), mp.Ey,
                           mp.Vector3(dpml + (-0.5 * sx)), size=mp.Vector3(0, w)))

        sim.symmetries.append(mp.Mirror(mp.Y, phase=-1))

        freg = mp.FluxRegion(center=mp.Vector3((0.5 * sx) - dpml - 0.5),
                             size=mp.Vector3(0, 2 * w))

        # transmitted flux
        trans = sim.add_flux(fcen, df, nfreq, freg)

        vol = mp.Volume(mp.Vector3(), size=mp.Vector3(sx))

        sim.run(
            mp.at_beginning(mp.output_epsilon),
            mp.during_sources(
                mp.in_volume(vol, mp.to_appended("hz-slice", mp.at_every(0.4, mp.output_hfield_z)))),
            until_after_sources=mp.stop_when_fields_decayed(50, mp.Ey, mp.Vector3((0.5 * sx) - dpml - 0.5, 0), 1e-3)
        )

        sim.display_fluxes(trans)  # print out the flux spectrum


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--resonant_modes', action='store_true', default=False,
                        help="Compute resonant modes. Default is transmission spectrum.")
    args = parser.parse_args()
    main(args)
