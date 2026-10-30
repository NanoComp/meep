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

    def test_user_array_wrong_rank(self):
        # A flat array must not be accepted for a 2d slice: zip()-based shape
        # checking used to truncate to the shorter sequence and let this
        # through, after which C++ wrote past the end of the buffer.
        self.sim.run(until_after_sources=0)
        dtype = np.float32 if mp.is_single_precision() else np.float64
        vol = mp.Volume(center=self.center_2d, size=self.size_2d)

        with self.assertRaises(ValueError):
            self.sim.get_array(mp.Hz, vol, arr=np.zeros(126, dtype=dtype))
        with self.assertRaises(ValueError):
            self.sim.get_array(mp.Hz, vol, arr=np.zeros((126, 38, 1), dtype=dtype))

    def test_user_array_wrong_dtype(self):
        # The SWIG typemap is a bare pointer cast, so a mismatched float width
        # would make C++ write 2x the bytes the buffer can hold.
        self.sim.run(until_after_sources=0)
        wrong = np.float64 if mp.is_single_precision() else np.float32
        vol = mp.Volume(center=self.center_1d, size=self.size_1d)

        with self.assertRaises(ValueError):
            self.sim.get_array(mp.Hz, vol, arr=np.zeros(126, dtype=wrong))
        with self.assertRaises(ValueError):
            # complex buffer for a real slice is equally unusable
            self.sim.get_array(mp.Hz, vol, arr=np.zeros(126, dtype=np.complex128))

    def test_user_array_not_contiguous(self):
        # np.require() used to silently substitute a copy here, so the caller's
        # array was never actually written to.
        self.sim.run(until_after_sources=0)
        dtype = np.float32 if mp.is_single_precision() else np.float64
        vol = mp.Volume(center=self.center_2d, size=self.size_2d)

        with self.assertRaises(ValueError):
            self.sim.get_array(mp.Hz, vol, arr=np.zeros((126, 76), dtype=dtype)[:, ::2])

        readonly = np.zeros((126, 38), dtype=dtype)
        readonly.flags.writeable = False
        with self.assertRaises(ValueError):
            self.sim.get_array(mp.Hz, vol, arr=readonly)

    def test_slice_outside_cell(self):
        # A subvolume that intersects no chunk used to leave slice_size
        # uninitialized in do_get_array_slice, so the allocation was sized from
        # stack garbage.
        self.sim.run(until_after_sources=0)
        arr = self.sim.get_array(
            mp.Hz, center=mp.Vector3(1000, 1000), size=mp.Vector3(1, 1)
        )
        self.assertEqual(np.count_nonzero(arr), 0)

        # The same path also left an am_now_working_on(FieldOutput) push on the
        # timing stack, which pauses the parent sink for the rest of the run.
        before = sum(self.sim.time_spent_on(mp.Stepping))
        self.sim.run(until=50)
        self.assertGreater(sum(self.sim.time_spent_on(mp.Stepping)), before)

        # ...and a normal slice must still work afterwards.
        vol = mp.Volume(center=self.center_2d, size=self.size_2d)
        self.assertEqual(self.sim.get_array(mp.Hz, vol).shape, (126, 38))

    def test_invalid_component(self):
        self.sim.run(until_after_sources=0)
        # out of range -> a bad value
        with self.assertRaises(ValueError):
            self.sim.get_array(component=-3)
        # wrong type -> a bad type
        with self.assertRaises(TypeError):
            self.sim.get_array(component="Hz")
        with self.assertRaises(TypeError):
            self.sim.get_array(mp.Volume(center=self.center_1d, size=self.size_1d))
        # component has no default any more
        with self.assertRaises(TypeError):
            self.sim.get_array()

    def test_volume_and_center_are_exclusive(self):
        self.sim.run(until_after_sources=0)
        vol = mp.Volume(center=self.center_2d, size=self.size_2d)
        with self.assertRaises(ValueError):
            self.sim.get_array(mp.Hz, vol=vol, center=self.center_2d)
        with self.assertRaises(ValueError):
            self.sim.get_array(mp.Hz, vol=vol, size=self.size_2d)

    def test_slice_dimensions_named_tuple(self):
        self.sim.run(until_after_sources=0)
        vol = mp.Volume(center=self.center_2d, size=self.size_2d)
        dims = self.sim.get_array_slice_dimensions(mp.Hz, vol=vol)

        # still unpacks positionally
        dim_sizes, min_corner, max_corner = dims
        self.assertIs(dim_sizes, dims.dim_sizes)
        self.assertIs(min_corner, dims.min_corner)
        self.assertIs(max_corner, dims.max_corner)
        self.assertIsInstance(dims.min_corner, mp.Vector3)

    def test_slice_dimensions_match_get_array(self):
        """dim_sizes must equal the shape get_array actually returns."""
        self.sim.run(until_after_sources=0)
        regions = {
            "2d": mp.Volume(center=self.center_2d, size=self.size_2d),
            "line": mp.Volume(center=self.center_1d, size=self.size_1d),
            "point": mp.Volume(center=self.center_1d, size=mp.Vector3()),
        }
        for name, vol in regions.items():
            # get_array always interpolates onto the centered grid, so the
            # default component (mp.Centered) must predict every component.
            dims, _, _ = self.sim.get_array_slice_dimensions(vol=vol)
            for component in (mp.Hz, mp.Ex, mp.Ey, mp.Dielectric):
                for snap in (False, True):
                    with self.subTest(region=name, component=component, snap=snap):
                        actual = self.sim.get_array(component, vol=vol, snap=snap).shape
                        self.assertEqual(dims, actual)

        # whole cell, and the component argument is genuinely optional
        self.assertEqual(
            self.sim.get_array_slice_dimensions().dim_sizes,
            self.sim.get_array(mp.Hz).shape,
        )

    def test_slice_dimensions_component_validation(self):
        self.sim.run(until_after_sources=0)
        vol = mp.Volume(center=self.center_2d, size=self.size_2d)
        default = self.sim.get_array_slice_dimensions(vol=vol).dim_sizes

        # aliases for the centered grid
        for centered in (mp.Centered, mp.Dielectric, mp.Permeability):
            self.assertEqual(
                self.sim.get_array_slice_dimensions(centered, vol=vol).dim_sizes,
                default,
            )
        # Er does not exist in 2d Cartesian: this used to reach meep::abort
        # inside loop_in_chunks, which takes down every rank under MPI.
        # NO_COMPONENT used to return a meaningless shape, and a derived
        # component an opaque SWIG "wrong number or type of arguments" error.
        for bad in (mp.Er, mp.NO_COMPONENT, mp.Sz):
            with self.assertRaises(ValueError):
                self.sim.get_array_slice_dimensions(bad, vol=vol)
        with self.assertRaises(TypeError):
            self.sim.get_array_slice_dimensions("Ez", vol=vol)

    def test_slice_dimensions_component_absent_from_grid(self):
        # Ez has no coordinate mismatch in 1d but the grid carries no Ez, so
        # there cannot be Yee-grid data for it; Ex in cylindrical is a mismatch.
        src = mp.GaussianSource(0.15, fwidth=0.1)
        for cell, dims, bad, ok in (
            (mp.Vector3(0, 0, 6), 1, mp.Ez, mp.Ex),
            (mp.Vector3(3, 0, 4), mp.CYLINDRICAL, mp.Ex, mp.Er),
        ):
            sim = mp.Simulation(
                cell_size=cell,
                resolution=10,
                dimensions=dims,
                sources=[mp.Source(src, ok, mp.Vector3())],
            )
            sim.init_sim()
            with self.subTest(dimensions=dims):
                with self.assertRaises(ValueError):
                    sim.get_array_slice_dimensions(bad)
                sim.get_array_slice_dimensions(ok)  # no raise

    def test_keyword_only_options(self):
        """cmplx/arr/frequency/snap may no longer be passed positionally."""
        self.sim.run(until_after_sources=0)
        vol = mp.Volume(center=self.center_1d, size=self.size_1d)
        with self.assertRaises(TypeError):
            self.sim.get_array(mp.Hz, vol, None, None, True)
        # the region arguments stay positional
        self.assertEqual(self.sim.get_array(mp.Hz, vol).shape, (126,))

    def test_frequency_only_for_materials(self):
        # `frequency` is only consumed when evaluating chi1inv, so silently
        # accepting it for a field component hid a mistake in the caller.
        self.sim.run(until_after_sources=0)
        with self.assertRaises(ValueError):
            self.sim.get_array(mp.Hz, frequency=0.25)
        self.sim.get_array(mp.Dielectric, frequency=0.25)  # no raise

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
