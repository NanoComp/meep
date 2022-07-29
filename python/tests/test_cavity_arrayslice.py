import os
import unittest

import numpy as np
from utils import ApproxComparisonTestCase

import meep as mp


class TestCavityArraySlice(ApproxComparisonTestCase):

    data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "data"))
    expected_1d = np.load(os.path.join(data_dir, "cavity_arrayslice_1d.npy"))
    expected_2d = np.load(os.path.join(data_dir, "cavity_arrayslice_2d.npy"))

    def setUp(self):

        r = 0.36
        d = 1.4
        sy = 6
        pad = 2
        dpml = 1
        sx = (2 * (pad + dpml + 3)) + d - 1

        cell = mp.Vector3(sx, sy, 0)

        blk = mp.Block(
            size=mp.Vector3(mp.inf, 1.2, mp.inf), material=mp.Medium(epsilon=13)
        )

        geometry = [blk]

        geometry.extend(mp.Cylinder(r, center=mp.Vector3(d / 2 + i)) for i in range(3))
        geometry.extend(mp.Cylinder(r, center=mp.Vector3(d / -2 - i)) for i in range(3))

        sources = [mp.Source(mp.GaussianSource(0.25, fwidth=0.2), mp.Hz, mp.Vector3())]

        self.sim = mp.Simulation(
            cell_size=cell,
            geometry=geometry,
            sources=sources,
            boundary_layers=[mp.PML(dpml)],
            resolution=20,
        )

        self.x_min = -0.25 * sx
        self.x_max = +0.25 * sx
        self.y_min = -0.15 * sy
        self.y_max = +0.15 * sy

        self.size_1d = mp.Vector3(self.x_max - self.x_min)
        self.center_1d = mp.Vector3((self.x_min + self.x_max) / 2)

        self.size_2d = mp.Vector3(self.x_max - self.x_min, self.y_max - self.y_min)
        self.center_2d = mp.Vector3(
            (self.x_min + self.x_max) / 2, (self.y_min + self.y_max) / 2
        )

    def test_1d_slice(self):
        self.sim.run(until_after_sources=0)
        vol = mp.Volume(center=self.center_1d, size=self.size_1d)
        hl_slice1d = self.sim.get_array(mp.Hz, vol)
        tol = 1e-5 if mp.is_single_precision() else 1e-8
        self.assertClose(self.expected_1d, hl_slice1d, epsilon=tol)

    def test_2d_slice(self):
        self.sim.run(until_after_sources=0)
        vol = mp.Volume(center=self.center_2d, size=self.size_2d)
        hl_slice2d = self.sim.get_array(mp.Hz, vol)
        tol = 1e-5 if mp.is_single_precision() else 1e-8
        self.assertClose(self.expected_2d, hl_slice2d, epsilon=tol)

    def test_1d_slice_user_array(self):
        self.sim.run(until_after_sources=0)
        arr = np.zeros(
            126, dtype=np.float32 if mp.is_single_precision() else np.float64
        )
        vol = mp.Volume(center=self.center_1d, size=self.size_1d)
        self.sim.get_array(mp.Hz, vol, arr=arr)
        tol = 1e-5 if mp.is_single_precision() else 1e-8
        self.assertClose(self.expected_1d, arr, epsilon=tol)

    def test_2d_slice_user_array(self):
        self.sim.run(until_after_sources=0)
        arr = np.zeros(
            (126, 38), dtype=np.float32 if mp.is_single_precision() else np.float64
        )
        vol = mp.Volume(center=self.center_2d, size=self.size_2d)
        self.sim.get_array(mp.Hz, vol, arr=arr)
        tol = 1e-5 if mp.is_single_precision() else 1e-8
        self.assertClose(self.expected_2d, arr, epsilon=tol)

    def test_illegal_user_array(self):
        self.sim.run(until_after_sources=0)

        with self.assertRaises(ValueError):
            arr = np.zeros(128)
            vol = mp.Volume(center=self.center_1d, size=self.size_1d)
            self.sim.get_array(mp.Hz, vol, arr=arr)

        with self.assertRaises(ValueError):
            arr = np.zeros((126, 39))
            vol = mp.Volume(center=self.center_2d, size=self.size_2d)
            self.sim.get_array(mp.Hz, vol, arr=arr)

        with self.assertRaises(ValueError):
            arr = np.zeros((126, 38))
            vol = mp.Volume(center=self.center_2d, size=self.size_2d)
            self.sim.get_array(mp.Hz, vol, cmplx=True, arr=arr)

    def test_1d_complex_slice(self):
        self.sim.run(until_after_sources=0)
        vol = mp.Volume(center=self.center_1d, size=self.size_1d)
        hl_slice1d = self.sim.get_array(mp.Hz, vol, cmplx=True)
        self.assertTrue(
            hl_slice1d.dtype == np.complex64
            if mp.is_single_precision()
            else np.complex128
        )
        self.assertTrue(hl_slice1d.shape[0] == 126)

    def test_2d_complex_slice(self):
        self.sim.run(until_after_sources=0)
        vol = mp.Volume(center=self.center_2d, size=self.size_2d)
        hl_slice2d = self.sim.get_array(mp.Hz, vol, cmplx=True)
        self.assertTrue(
            hl_slice2d.dtype == np.complex64
            if mp.is_single_precision()
            else np.complex128
        )
        self.assertTrue(hl_slice2d.shape[0] == 126 and hl_slice2d.shape[1] == 38)


if __name__ == "__main__":
    unittest.main()
