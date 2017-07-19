from __future__ import division

# Simple test for libmeepgeom, modeled after meep_test.ctl
import os
import subprocess
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


def compare_hdf5_datasets(file1, file2):

    ret = subprocess.check_call(['h5diff', file1, file2])

    return not ret

    # f1 = mp.h5file(file1, mp.h5file.READONLY, False)
    # data1, rank1, dims1 = f1.read(name1)
    # if not data1:
    #     return False

    # f2 = mp.h5file(file2, mp.h5file.READONLY, False)
    # data2, rank2, dims2 = f2.read(name2)
    # if not data2:
    #     return False

    # if len(dims1) != expected_rank or len(dims2) != expected_rank:
    #     return False

    # size = 1
    # for r in range(expected_rank):
    #     if dims1[r] != dims2[r]:
    #         return False
    #     size *= dims1[r]

    # for n in range(size):
    #     d1 = data1[n]
    #     d2 = data2[n]
    #     diff = abs(d1 - d2)
    #     maximum = max(abs(d1), abs(d2))
    #     if diff > abs_tol or diff > maximum * rel_tol:
    #         return False

    # return True


class TestCylEllipsoid(unittest.TestCase):

    def setUp(self):

        test_dir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))
        libmeepgeom_dir = os.path.join(test_dir, '..', '..', 'libmeepgeom')
        self.eps_ref_file = os.path.join(libmeepgeom_dir, 'cyl-ellipsoid-eps-ref.h5')
        self.src_cmpt = mp.Ez
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
        # TODO(chogan): Get hdf5 filename programatically?
        status = compare_hdf5_datasets('eps-000000000.h5', self.eps_ref_file,)
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
