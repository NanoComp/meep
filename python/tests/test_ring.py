# Python port of meep/examples/ring.ctl
# Calculating 2d ring-resonator modes, from the Meep tutorial.
from __future__ import division

import sys

import meep as mp
import meep.geom as gm
from meep.source import GaussianSource


# dummy material function needed to pass to structure( )
# constructor as a placeholder before we can call
# set_materials_from_geometry

def dummy_eps(vec):
    return 1.0


def main(args):

    n = 3.4  # index of waveguide
    w = 1.0  # width of waveguide
    r = 1.0  # inner radius of ring

    pad = 4  # padding between waveguide and edge of PML
    dpml = 2  # thickness of PML

    sxy = 2.0 * (r + w + pad + dpml)  # cell size
    resolution = 10.0

    gv = mp.voltwo(sxy, sxy, resolution)
    gv.center_origin()

    sym = mp.mirror(mp.Y, gv)

    # exploit the mirror symmetry in structure+source:
    the_structure = mp.structure(gv, dummy_eps, mp.pml(dpml), sym)

    # Create a ring waveguide by two overlapping cylinders - later objects
    # take precedence over earlier objects, so we put the outer cylinder first.
    # and the inner (air) cylinder second.

    objects = []
    n2 = n * n
    dielectric = gm.Medium(epsilon_diag=gm.Vector3(n2, n2, n2))
    objects.append(gm.Cylinder(r + w, material=dielectric))
    objects.append(gm.Cylinder(r))

    mp.set_materials_from_geometry(the_structure, objects)
    f = mp.fields(the_structure)

    # If we don't want to excite a specific mode symmetry, we can just
    # put a single point source at some arbitrary place, pointing in some
    # arbitrary direction.  We will only look for TM modes (E out of the plane).
    fcen = 0.15  # pulse center frequency
    df = 0.1
    src = GaussianSource(fcen, df)
    v = mp.volume(mp.vec(r + 0.1, 0.0), mp.vec(0.0, 0.0))
    f.add_volume_source(mp.Ez, src.swigobj, v)

    T = 300.0
    stop_time = f.last_source_time() + T
    while f.round_time() < stop_time:
        f.step()

    # TODO: translate call to harminv
    # int bands = do_harminv (... Ez, vec3(r+0.1), fcen, df)

    # Output fields for one period at the end.  (If we output
    # at a single time, we might accidentally catch the Ez field
    # when it is almost zero and get a distorted view.)
    DeltaT = 1.0 / (20 * fcen)
    NextOutputTime = f.round_time() + DeltaT
    while f.round_time() < 1.0 / fcen:
        f.step()
        if f.round_time() >= NextOutputTime:
            f.output_hdf5(mp.Ez, f.total_volume())
            NextOutputTime += DeltaT


if __name__ == '__main__':
    sys.exit(main(sys.argv))
