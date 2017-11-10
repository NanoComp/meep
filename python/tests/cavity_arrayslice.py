from __future__ import division

import unittest

import meep as mp
import numpy as np


class TestCavityArraySlice(unittest.TestCase):

    def setUp(self):

        r = 0.36
        d = 1.4
        sy = 6
        pad = 2
        dpml = 1
        sx = (2 * (pad + dpml + 3)) + d - 1

        cell = mp.Vector3(sx, sy, 0)

        blk = mp.Block(size=mp.Vector3(1e20, 1.2, 1e20),
                       material=mp.Medium(epsilon=13))

        geometry = [blk]

        for i in range(3):
            geometry.append(mp.Cylinder(r, center=mp.Vector3(d / 2 + i)))

        for i in range(3):
            geometry.append(mp.Cylinder(r, center=mp.Vector3(d / -2 - i)))

        sources = [mp.Source(mp.GaussianSource(0.25, fwidth=0.2), mp.Hz, mp.Vector3())]

        self.sim = mp.Simulation(
            cell_size=cell,
            geometry=geometry,
            sources=sources,
            boundary_layers=[mp.PML(dpml)],
            resolution=20
        )

        self.x_min = -0.25 * sx
        self.x_max = +0.25 * sx
        self.y_min = -0.15 * sy
        self.y_max = +0.15 * sy

    def test_1d_slice(self):
        self.sim.run(until_after_sources=0)

        # Low level 1D slice of Hz data
        dims1d = np.zeros(3, dtype=np.int32)
        v1d = mp.volume(mp.vec(self.x_min, 0.0), mp.vec(self.x_max, 0.0))
        self.sim.fields.get_array_slice_dimensions(v1d, dims1d)
        NX1 = dims1d[0]
        slice1d = np.zeros(NX1, dtype=np.double)
        self.sim.fields.get_array_slice(v1d, mp.Hz, slice1d)

        # High level equivalent
        size = mp.Vector3(self.x_max - self.x_min)
        center = mp.Vector3((self.x_min + self.x_max) / 2)
        hl_slice1d = self.sim.get_array(center=center, size=size, component=mp.Hz)

        np.testing.assert_allclose(slice1d, hl_slice1d)

    def test_2d_slice(self):
        self.sim.run(until_after_sources=0)

        # Low level 2D slice of Hz data
        dims2d = np.zeros(3, dtype=np.int32)
        v2d = mp.volume(mp.vec(self.x_min, self.y_min), mp.vec(self.x_max, self.y_max))
        self.sim.fields.get_array_slice_dimensions(v2d, dims2d)
        NX2 = dims2d[0]
        NY2 = dims2d[1]
        N2 = NX2 * NY2
        slice2d = np.zeros(N2, dtype=np.double)
        self.sim.fields.get_array_slice(v2d, mp.Hz, slice2d)

        # High level equivalent
        size = mp.Vector3(self.x_max - self.x_min, self.y_max - self.y_min)
        center = mp.Vector3((self.x_min + self.x_max) / 2, (self.y_min + self.y_max) / 2)
        hl_slice2d = self.sim.get_array(center=center, size=size, component=mp.Hz)

        np.testing.assert_allclose(slice2d, hl_slice2d)


if __name__ == '__main__':
    unittest.main()
