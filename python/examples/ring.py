# Calculating 2d ring-resonator modes, from the Meep tutorial.
from __future__ import division

import sys
import meep as mp


def main(args):
    n = 3.4  # index of waveguide
    w = 1  # width of waveguide
    r = 1  # inner radius of ring
    pad = 4  # padding between waveguide and edge of PML
    dpml = 2  # thickness of PML
    sxy = 2 * (r + w + pad + dpml)  # cell size

    # Create a ring waveguide by two overlapping cylinders - later objects
    # take precedence over earlier objects, so we put the outer cylinder first.
    # and the inner (air) cylinder second.
    dielectric = mp.Medium(epsilon_diag=mp.epsilon(n * n))
    air = mp.Medium(epsilon_diag=mp.Vector3(1, 1, 1))

    c1 = mp.Cylinder(r + w, material=dielectric)
    c2 = mp.Cylinder(r, material=air)

    # If we don't want to excite a specific mode symmetry, we can just
    # put a single point source at some arbitrary place, pointing in some
    # arbitrary direction.  We will only look for TM modes (E out of the plane).

    fcen = 0.15  # pulse center frequency
    df = 0.1  # pulse width (in frequency)

    src = mp.Source(mp.GaussianSource(fcen, df), mp.Ez, mp.Vector3(r + 0.1))

    sim = mp.Simulation(cell_size=mp.Vector3(sxy, sxy),
                        geometry=[c1, c2],
                        sources=[src],
                        resolution=10,
                        symmetries=[mp.Mirror(mp.Y)],
                        boundary_layers=[mp.Pml(dpml)])

    sim.run(
        sim.at_beginning(sim.output_epsilon),
        sim.after_sources(sim.harminv(mp.Ez, mp.Vector3(r + 0.1), fcen, df)),
        until=300, sources=True
    )

    # Output fields for one period at the end.  (If we output
    # at a single time, we might accidentally catch the Ez field when it is
    # almost zero and get a distorted view.)
    sim.run(sim.at_every((1 / fcen / 20), sim.output_efield_z), until=(1 / fcen))

if __name__ == '__main__':
    sys.exit(main(sys.argv))
