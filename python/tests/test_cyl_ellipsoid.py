from __future__ import division

import unittest

import meep as mp
from meep.geom import Cylinder, Ellipsoid, Medium, Vector3
from meep.source import GaussianSource


# Simple test for libmeepgeom, modeled after meep_test.ctl

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
    src = GaussianSource(fcen, fwidth=df)
    src_point = mp.vec(0.0, 0.0)
    fields.add_point_source(src_cmpt, src, src_point)


class TestCylEllipsoid(unittest.TestCase):

    ref_Ez = -8.29555720049629e-5
    ref_Hz = -4.5623185899766e-5

    def init(self):
        resolution = 100.0
        gv = mp.voltwo(10.0, 10.0, resolution)
        gv.center_origin()

        if self.src_cmpt == mp.Ez:
            sym = mp.mirror(mp.X, gv) + mp.mirror(mp.Y, gv)
        elif self.src_cmpt == mp.Hz:
            sym = -mp.mirror(mp.X, gv) - mp.mirror(mp.Y, gv)

        the_structure = mp.structure(gv, dummy_eps, mp.pml(1.0), sym)

        set_materials(the_structure)

        self.f = mp.fields(the_structure)
        self.duration = 23.0
        self.start_time = self.f.round_time()
        self.stop_time = self.start_time + self.duration

    def run_simulation(self):
        add_source(self.f, self.src_cmpt)

        while self.f.round_time() < self.stop_time:
            self.f.step()

        ref_out_field = self.ref_Ez if self.src_cmpt == mp.Ez else self.ref_Hz
        out_field = self.f.get_field(self.src_cmpt, mp.vec(4.13, 3.75)).real
        diff = abs(out_field - ref_out_field)

        self.assertTrue(abs(diff) <= 0.05 * abs(ref_out_field), "Field output differs")

    def test_ez_field(self):

        self.src_cmpt = mp.Ez
        self.init()
        self.run_simulation()

    def test_hz_field(self):

        self.src_cmpt = mp.Hz
        self.init()
        self.run_simulation()


if __name__ == '__main__':
    unittest.main()
    mp.all_wait()
