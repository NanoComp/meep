import unittest
import warnings

import meep as mp
import numpy as np
from utils import ApproxComparisonTestCase


class TestArrayMetadata(ApproxComparisonTestCase):
    def test_array_metadata(self):
        """
        Verifies that the CW fields via the modal volume of a ring resonator
        are the same for the CW solver and time stepping using a pulsed source.
        """
        resolution = 25

        n = 3.4
        w = 1
        r = 1
        pad = 4
        dpml = 2

        sxy = 2 * (r + w + pad + dpml)
        cell_size = mp.Vector3(sxy, sxy)

        nonpml_vol = mp.Volume(
            mp.Vector3(), size=mp.Vector3(sxy - 2 * dpml, sxy - 2 * dpml)
        )

        geometry = [
            mp.Cylinder(radius=r + w, material=mp.Medium(index=n)),
            mp.Cylinder(radius=r),
        ]

        fcen = 0.118
        df = 0.08

        symmetries = [mp.Mirror(mp.X, phase=-1), mp.Mirror(mp.Y, phase=+1)]

        pml_layers = [mp.PML(dpml)]

        # CW source
        src = [
            mp.Source(
                mp.ContinuousSource(fcen, fwidth=df),
                mp.Ez,
                mp.Vector3(r + 0.1),
            ),
            mp.Source(
                mp.ContinuousSource(fcen, fwidth=df),
                mp.Ez,
                mp.Vector3(-(r + 0.1)),
                amplitude=-1,
            ),
        ]

        sim = mp.Simulation(
            cell_size=cell_size,
            geometry=geometry,
            sources=src,
            resolution=resolution,
            force_complex_fields=True,
            symmetries=symmetries,
            boundary_layers=pml_layers,
        )

        sim.init_sim()
        # The convergence properties of the CW solver's biCGSTAB-L algorithm
        # is sensitive to the floating-point precision of the fields. This
        # requires using a smaller L for single compared to double precision.
        sim.solve_cw(
            1e-4 if mp.is_single_precision() else 1e-6,  # tol
            1000,  # maxiters
            7 if mp.is_single_precision() else 10,  # L
        )

        def electric_energy(r, ez, eps):
            return np.real(eps * np.conj(ez) * ez)

        def vec_func(r):
            return r.x**2 + 2 * r.y**2

        electric_energy_total = sim.integrate_field_function(
            [mp.Ez, mp.Dielectric], electric_energy, nonpml_vol
        )
        electric_energy_max = sim.max_abs_field_function(
            [mp.Ez, mp.Dielectric], electric_energy, nonpml_vol
        )
        vec_func_total = sim.integrate_field_function([], vec_func, nonpml_vol)
        cw_modal_volume = (electric_energy_total / electric_energy_max) * vec_func_total

        sim.reset_meep()

        # pulsed source
        src = [
            mp.Source(
                mp.GaussianSource(fcen, fwidth=df),
                mp.Ez,
                mp.Vector3(r + 0.1),
            ),
            mp.Source(
                mp.GaussianSource(fcen, fwidth=df),
                mp.Ez,
                mp.Vector3(-(r + 0.1)),
                amplitude=-1,
            ),
        ]

        sim = mp.Simulation(
            cell_size=cell_size,
            geometry=geometry,
            k_point=mp.Vector3(),
            sources=src,
            resolution=resolution,
            symmetries=symmetries,
            boundary_layers=pml_layers,
        )

        dft_obj = sim.add_dft_fields([mp.Ez], fcen, 0, 1, where=nonpml_vol)
        sim.run(until_after_sources=100)

        Ez = sim.get_dft_array(dft_obj, mp.Ez, 0)
        (X, Y, Z, W) = sim.get_array_metadata(dft_cell=dft_obj)
        Eps = sim.get_array(vol=nonpml_vol, component=mp.Dielectric)
        EpsE2 = np.real(Eps * np.conj(Ez) * Ez)
        # W is indexed (x, y) so the meshgrid must "ij" indexing; the
        # default "xy" would transpose it relative to W.
        xm, ym = np.meshgrid(X, Y, indexing="ij")
        vec_func_sum = np.sum(W * (xm**2 + 2 * ym**2))
        pulse_modal_volume = np.sum(W * EpsE2) / np.max(EpsE2) * vec_func_sum

        tol = 0.05 if mp.is_single_precision() else 0.01
        self.assertClose(
            cw_modal_volume / pulse_modal_volume,
            1.0,
            epsilon=tol,
        )

    def test_metadata_types_and_shapes(self):
        """x/y/z must be NumPy arrays, and w must match get_array's shape."""
        sim = mp.Simulation(
            cell_size=mp.Vector3(6, 4),
            resolution=10,
            sources=[
                mp.Source(mp.GaussianSource(0.15, fwidth=0.1), mp.Ez, mp.Vector3())
            ],
        )
        sim.init_sim()

        # a non-square region, so a transposed `w` cannot go unnoticed
        vols = [
            mp.Volume(center=mp.Vector3(), size=mp.Vector3(3, 2)),  # 2d
            mp.Volume(center=mp.Vector3(), size=mp.Vector3(3, 0)),  # zero-thickness
            mp.Volume(center=mp.Vector3(), size=mp.Vector3(0, 0)),  # single point
        ]
        for vol in vols:
            x, y, z, w = sim.get_array_metadata(vol=vol)
            for tics in (x, y, z):
                self.assertIsInstance(tics, np.ndarray)
                self.assertEqual(tics.dtype, np.float64)
            self.assertEqual(w.shape, sim.get_array(mp.Ez, vol=vol).shape)

    def test_metadata_named_tuple(self):
        sim = mp.Simulation(
            cell_size=mp.Vector3(6, 4),
            resolution=10,
            sources=[
                mp.Source(mp.GaussianSource(0.15, fwidth=0.1), mp.Ez, mp.Vector3())
            ],
        )
        sim.init_sim()
        box = mp.Volume(center=mp.Vector3(), size=mp.Vector3(3, 2))

        meta = sim.get_array_metadata(vol=box)

        # still unpacks positionally, and the names refer to the same objects
        x, y, z, w = meta
        self.assertIs(x, meta.x)
        self.assertIs(y, meta.y)
        self.assertIs(z, meta.z)
        self.assertIs(w, meta.w)

        # `points` is computed on access and shaped like `w`
        self.assertEqual(meta.points.shape, meta.w.shape)
        self.assertIsInstance(meta.points[0, 0], mp.Vector3)
        self.assertEqual(meta.points[0, 0], mp.Vector3(meta.x[0], meta.y[0], 0))

        # a weighted integral over the region reads naturally
        total = np.sum(meta.w)
        self.assertAlmostEqual(total, 3 * 2, places=6)

    def test_metadata_return_pw_deprecated(self):
        sim = mp.Simulation(
            cell_size=mp.Vector3(6, 4),
            resolution=10,
            sources=[
                mp.Source(mp.GaussianSource(0.15, fwidth=0.1), mp.Ez, mp.Vector3())
            ],
        )
        sim.init_sim()
        box = mp.Volume(center=mp.Vector3(), size=mp.Vector3(3, 2))

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            p, w = sim.get_array_metadata(vol=box, return_pw=True)
        self.assertTrue(
            any(issubclass(c.category, DeprecationWarning) for c in caught),
            "return_pw should emit a DeprecationWarning",
        )

        meta = sim.get_array_metadata(vol=box)
        self.assertEqual(p.shape, w.shape)
        self.assertEqual(p.shape, meta.points.shape)
        np.testing.assert_array_equal(w, meta.w)

    def test_metadata_cylindrical_raises(self):
        sim = mp.Simulation(
            cell_size=mp.Vector3(4, 0, 4),
            dimensions=mp.CYLINDRICAL,
            resolution=10,
            sources=[
                mp.Source(mp.GaussianSource(0.15, fwidth=0.1), mp.Er, mp.Vector3(1))
            ],
        )
        sim.init_sim()
        # C++ calls meep::abort() here; Python should raise instead.
        with self.assertRaises(ValueError):
            sim.get_array_metadata()


if __name__ == "__main__":
    unittest.main()
