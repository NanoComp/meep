from __future__ import division

# Simple test for libmeepgeom, modeled after meep_test.ctl
import unittest

import meep as mp
from meep.geom import Cylinder, Ellipsoid, Medium, Vector3
from meep.source import GaussianSource


# Dummy material function needed to pass to structure( ) constructor as a
# placeholder before we can call set_materials_from_geometry
def dummy_eps(vec):
    return 1.0


def set_materials(structure):
    n = 3.5  # index of refraction
    nsqr = n * n
    dielectric = Medium(epsilon_diag=Vector3(nsqr, nsqr, nsqr))
    objects = []
    radius = 3.0
    height = float(1.0e20)
    size = Vector3(1.0, 2.0, 1.0e20)
    objects.append(Cylinder(material=dielectric, radius=radius, height=height))
    objects.append(Ellipsoid(size=size))
    mp.set_materials_from_geometry(structure, objects)


def add_source(fields, src_cmpt):
    fcen = 1.0
    df = 0.1
    src = GaussianSource(fcen, df)
    src_point = mp.vec(0.0, 0.0)
    fields.add_point_source(src_cmpt, src, src_point)


def compare_hdf5_datasets(test_file, ref_file):
    # TODO
    return True


class TestCylEllipsoid(unittest.TestCase):

    def setUp(self):

        # TODO(chogan): Need abs path to this file
        self.eps_ref_file = 'cyl-ellipsoid-eps-ref.h5'
        self.src_cmpt = mp.Hz
        resolution = 100.0

        gv = mp.voltwo(10.0, 10.0, resolution)
        gv.center_origin()

        sym = mp.mirror(mp.X, gv) + mp.mirror(mp.Y, gv)
        the_structure = mp.structure(gv, dummy_eps, mp.pml(1.0), sym)

        set_materials(the_structure)

        self.f = mp.fields(the_structure)
        add_source(self.f, self.src_cmpt)

    def test_hdf5_output(self):
        self.f.output_hdf5(mp.Dielectric, self.f.total_volume())
        # TODO(chogan): Get hdf5 filename programatically? Use abs path.
        status = compare_hdf5_datasets('eps-000000000.h5', self.eps_ref_file)
        self.assertTrue(status, "HDF5 output differs")

    def test_fields(self):
        duration = 23.0
        start_time = self.f.round_time()
        stop_time = start_time + duration

        while self.f.round_time() < stop_time:
            self.f.step()

        ref_ez = -8.29555720049629e-5
        ref_hz = -4.5623185899766e-5
        ref_out_field = ref_ez if self.src_cmpt == mp.Ez else ref_hz

        out_field = self.f.get_field(self.src_cmpt, mp.vec(4.13, 3.75)).real
        diff = abs(out_field - ref_out_field)

        self.assertTrue(abs(diff) <= 0.05 * abs(ref_out_field), "Field output differs")

        print("field: {} + i{}".format(out_field.real, out_field.imag))


if __name__ == '__main__':
    unittest.main()


# def main(args):

#     if args.polarization in 'Pp':
#         src_cmpt = mp.Hz
#         print("Using P-polarization")
#     else:
#         src_cmpt = mp.Ez
#         print("Using S-polarization")

#     resolution = 100.0

#     # mp.geometry_lattice.size.x = 10.0
#     # mp.geometry_lattice.size.y = 10.0
#     # mp.geometry_lattice.size.z = 0.0

#     gv = mp.voltwo(10.0, 10.0, resolution)
#     gv.center_origin()

#     if src_cmpt == mp.Ez:
#         sym = mp.mirror(mp.X, gv) + mp.mirror(mp.Y, gv)
#     else:
#         sym = -mp.mirror(mp.X, gv) - mp.mirror(mp.Y, gv)

#     the_structure = mp.structure(gv, dummy_eps, mp.pml(1.0), sym)

#     set_materials(the_structure)

#     f = mp.fields(the_structure)
#     add_source(f, src_cmpt)

#     f.output_hdf5(mp.Dielectric, f.total_volume())

#     # TODO(chogan): Get hdf5 filename programatically?
#     status = compare_hdf5_datasets('eps-000000000.h5', args.eps_ref_file)

#     assert status, "HDF5 output differs"

#     duration = 23.0
#     start_time = f.round_time()
#     stop_time = start_time + duration

#     while f.round_time() < stop_time:
#         f.step()

#     ref_ez = -8.29555720049629e-5
#     ref_hz = -4.5623185899766e-5
#     ref_out_field = ref_ez if src_cmpt == mp.Ez else ref_hz

#     out_field = f.get_field(src_cmpt, mp.vec(4.13, 3.75)).real
#     diff = abs(out_field - ref_out_field)

#     assert abs(diff) <= 0.05 * abs(ref_out_field), "Field output differs"

#     print("field: {} + i{}".format(out_field.real, out_field.imag))


# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#     parser.add_argument(
#         '-p',
#         '--polarization',
#         default='S',
#         help="'S' for TE polarization or 'P' for TM polarization"
#     )
#     parser.add_argument(
#         '-f',
#         '--eps-ref-file',
#         default='cyl-ellipsoid-eps-ref.h5',
#         help='The reference h5 file for results comparison'
#     )
#     args = parser.parse_args()

#     sys.exit(main(args))
