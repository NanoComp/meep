from __future__ import division

import argparse
import sys

import meep as mp
import meep.geom as gm
from meep.source import GaussianSource


# Simple test for libmeepgeom, modeled after meep_test.ctl

# Dummy material function needed to pass to structure( ) constructor as a
# placeholder before we can call set_materials_from_geometry
def dummy_eps(vec):
    return 1.0


def main(args):

    src_cmpt = mp.Ez

    if args.polarization in 'Pp':
        src_cmpt = mp.Hz
        print("Using P-polarization")
    else:
        print("Using S-polarization")

    resolution = 100.0

    gv = mp.voltwo(10.0, 10.0, resolution)
    gv.center_origin()

    if src_cmpt == mp.Ez:
        sym = mp.mirror(mp.X, gv) + mp.mirror(mp.Y, gv)
    else:
        sym = -mp.mirror(mp.X, gv) - mp.mirror(mp.Y, gv)

    the_structure = mp.structure(gv, dummy_eps, mp.pml(1.0), sym)

    n = 3.5  # index of refraction
    nsqr = n * n
    dielectric = gm.Medium(epsilon_diag=gm.Vector3(nsqr, nsqr, nsqr))
    objects = []
    radius = 3.0
    height = float('inf')
    size = gm.Vector3(1.0, 2.0, float('inf'))
    objects.append(gm.Cylinder(material=dielectric, radius=radius, height=height))
    objects.append(gm.Ellipsoid(size=size))
    mp.set_materials_from_geometry(the_structure, objects)

    f = mp.fields(the_structure)
    fcen = 1.0
    df = 0.1
    src = GaussianSource(fcen, df)
    src_point = mp.vec(0.0, 0.0)
    src_size = mp.vec(10.0, 10.0)
    f.add_volume_source(src_cmpt, src.swigobj, mp.volume(src_point, src_size))

    f.output_hdf5(mp.Dielectric, f.total_volume())
    stop_time = 23.0

    while f.round_time() < stop_time:
        f.step()

    eval_pt = mp.vec(4.13, 3.75)
    out_field = f.get_field(src_cmpt, eval_pt)

    print("field: {} + i{}".format(out_field.real, out_field.imag))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--polarization', default='S', help="'S' for TE polarization or 'P' for TM polarization")
    args = parser.parse_args()

    sys.exit(main(args))
