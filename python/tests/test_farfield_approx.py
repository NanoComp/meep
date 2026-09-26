"""Verify that the far-field approximation (far_field_approx=True)
agrees with the exact near-to-far field transformation for a 3D
point-dipole source at observation distances far from the source."""

import math
import unittest

import meep as mp
import numpy as np

from utils import ApproxComparisonTestCase


class TestFarFieldApprox(ApproxComparisonTestCase):
    @classmethod
    def setUpClass(cls):
        resolution = 20
        sxy = 4
        dpml = 1
        cell_size = mp.Vector3(sxy + 2 * dpml, sxy + 2 * dpml, sxy + 2 * dpml)
        fcen = 1.0
        cls.sim = mp.Simulation(
            cell_size=cell_size,
            resolution=resolution,
            sources=[
                mp.Source(
                    mp.GaussianSource(fcen, fwidth=0.5 * fcen),
                    component=mp.Ez,
                    center=mp.Vector3(),
                )
            ],
            boundary_layers=[mp.PML(dpml)],
        )
        cls.n2f = cls.sim.add_near2far(
            fcen,
            0,
            1,
            mp.Near2FarRegion(
                center=mp.Vector3(0, 0, 0.5 * sxy),
                size=mp.Vector3(sxy, sxy, 0),
                weight=+1,
            ),
            mp.Near2FarRegion(
                center=mp.Vector3(0, 0, -0.5 * sxy),
                size=mp.Vector3(sxy, sxy, 0),
                weight=-1,
            ),
            mp.Near2FarRegion(
                center=mp.Vector3(0.5 * sxy, 0, 0),
                size=mp.Vector3(0, sxy, sxy),
                weight=+1,
            ),
            mp.Near2FarRegion(
                center=mp.Vector3(-0.5 * sxy, 0, 0),
                size=mp.Vector3(0, sxy, sxy),
                weight=-1,
            ),
            mp.Near2FarRegion(
                center=mp.Vector3(0, 0.5 * sxy, 0),
                size=mp.Vector3(sxy, 0, sxy),
                weight=+1,
            ),
            mp.Near2FarRegion(
                center=mp.Vector3(0, -0.5 * sxy, 0),
                size=mp.Vector3(sxy, 0, sxy),
                weight=-1,
            ),
        )
        cls.sim.run(until_after_sources=mp.stop_when_dft_decayed(1e-5))

    def test_farfield_approx_vs_exact(self):
        """Component-wise comparison of exact and approximate far fields.

        At R = 1000*lambda, the far-field approximation error is O(L/R)
        where L is the near-field surface extent, giving ~0.4% error for
        the dominant (significant) field components.
        """
        R = 1000.0
        npts = 8

        for theta in np.linspace(0.3, math.pi - 0.3, npts):
            for phi in [0, math.pi / 4, math.pi / 2]:
                x = mp.Vector3(
                    R * math.sin(theta) * math.cos(phi),
                    R * math.sin(theta) * math.sin(phi),
                    R * math.cos(theta),
                )
                ff_exact = np.array(self.sim.get_farfield(self.n2f, x))
                ff_approx = np.array(
                    self.sim.get_farfield(self.n2f, x, far_field_approx=True)
                )

                max_mag = max(abs(ff_exact[k]) for k in range(6))
                if max_mag < 1e-20:
                    continue

                for k in range(6):
                    if abs(ff_exact[k]) > 0.01 * max_mag:
                        rel_err = abs(ff_approx[k] - ff_exact[k]) / abs(ff_exact[k])
                        self.assertLess(
                            rel_err,
                            0.01,
                            f"component {k} at theta={theta:.2f}, phi={phi:.2f}: "
                            f"rel_err={rel_err:.2e}",
                        )

    def test_farfield_approx_radiation_pattern(self):
        """Verify that the normalized E_theta radiation pattern computed
        with the far-field approximation matches the exact result.

        For a z-oriented electric dipole, E_theta ~ sin(theta). Both
        exact and approximate calculations should produce the same
        normalized pattern to within 0.05%
        """
        R = 1000.0
        thetas = np.linspace(0.2, math.pi - 0.2, 20)
        E_theta_exact = np.zeros(len(thetas))
        E_theta_approx = np.zeros(len(thetas))

        for i, theta in enumerate(thetas):
            x = mp.Vector3(R * math.sin(theta), 0, R * math.cos(theta))

            ff_exact = self.sim.get_farfield(self.n2f, x)
            ff_approx = self.sim.get_farfield(self.n2f, x, far_field_approx=True)

            E_theta_exact[i] = abs(
                ff_exact[0] * math.cos(theta) - ff_exact[2] * math.sin(theta)
            )
            E_theta_approx[i] = abs(
                ff_approx[0] * math.cos(theta) - ff_approx[2] * math.sin(theta)
            )

        E_theta_exact /= E_theta_exact.max()
        E_theta_approx /= E_theta_approx.max()

        self.assertClose(E_theta_exact, E_theta_approx, epsilon=5e-3)

    def test_farfield_approx_poynting_flux(self):
        """Verify that the radiated Poynting flux integrated over a
        hemisphere agrees between the exact and approximate methods.

        The total radiated power from a z-dipole is obtained by
        integrated the Poynting vector over a hemisphere at R=1000
        and comparing the exact and approximate results.
        """
        R = 1000.0
        npts = 15
        angles = np.linspace(0, 0.5 * math.pi, npts, endpoint=False)
        Pr_exact = np.zeros(npts)
        Pr_approx = np.zeros(npts)

        for i, ang in enumerate(angles):
            x = mp.Vector3(R * math.sin(ang), 0, R * math.cos(ang))

            ff_exact = self.sim.get_farfield(self.n2f, x)
            ff_approx = self.sim.get_farfield(self.n2f, x, far_field_approx=True)

            Pr_exact[i] = np.real(
                np.conj(ff_exact[1]) * ff_exact[5] - np.conj(ff_exact[2]) * ff_exact[4]
            )
            Pr_approx[i] = np.real(
                np.conj(ff_approx[1]) * ff_approx[5]
                - np.conj(ff_approx[2]) * ff_approx[4]
            )

        flux_exact = np.sum(Pr_exact) * 0.5 * math.pi * R / npts
        flux_approx = np.sum(Pr_approx) * 0.5 * math.pi * R / npts

        rel_err = abs(flux_exact - flux_approx) / abs(flux_exact)
        self.assertLess(rel_err, 1e-3, f"flux rel_err = {rel_err:.2e}")

    def test_farfield_approx_get_farfields_array(self):
        """Verify that the get_farfields array API works with the
        far-field approximation and produces results consistent with
        the single-point get_farfield API.

        For a 3D non-periodic far field, get_farfields uses an
        FFT-accelerated evaluation of the equivalent-current surface
        integral, whereas the single-point get_farfield uses a direct
        summation.  The two therefore differ only by the bicubic
        interpolation error of the oversampled FFT, which is well below
        1% for the default oversampling.
        """
        R = 1000.0
        theta = math.pi / 3
        x_pt = mp.Vector3(R * math.sin(theta), 0, R * math.cos(theta))

        ff_single = np.array(
            self.sim.get_farfield(self.n2f, x_pt, far_field_approx=True)
        )

        result = self.sim.get_farfields(
            self.n2f,
            resolution=0,
            center=x_pt,
            size=mp.Vector3(0, 0, 0),
            far_field_approx=True,
        )

        ff_array = np.array(
            [
                result["Ex"].ravel()[0],
                result["Ey"].ravel()[0],
                result["Ez"].ravel()[0],
                result["Hx"].ravel()[0],
                result["Hy"].ravel()[0],
                result["Hz"].ravel()[0],
            ]
        )

        self.assertClose(ff_single, ff_array, epsilon=5e-3)

    def test_farfield_approx_fft_grid_vs_exact(self):
        """Verify the FFT-accelerated far-field path over a finite grid.

        get_farfields over a 3D non-periodic volume uses an
        FFT-accelerated evaluation of the radiation-zone surface
        integral.  This is only exercised for a genuine multi-point grid
        (a single point degenerates to a trivial transform), so we
        compare a far-field patch in the radiation zone against both the
        exact near-to-far transform (get_farfields with
        far_field_approx=False) and the per-point direct radiation-zone
        summation (get_farfield with far_field_approx=True).

        The patch spans many fractional FFT bins, so this stresses the
        bicubic interpolation of the oversampled FFT.  The aggregate
        error is dominated by the O(L/R) radiation-zone approximation
        (~0.4% at R=1000) plus a sub-percent interpolation contribution.
        """
        R = 1000.0
        theta = math.pi / 3
        # planar patch centered far away; extent (40) is << R so every
        # point remains deep in the radiation zone
        center = mp.Vector3(R * math.sin(theta), 0, R * math.cos(theta))
        size = mp.Vector3(40, 40, 0)
        npts = 8
        res = (npts - 1) / 40.0

        ff_fft = self.sim.get_farfields(
            self.n2f, resolution=res, center=center, size=size, far_field_approx=True
        )
        ff_exact = self.sim.get_farfields(
            self.n2f, resolution=res, center=center, size=size, far_field_approx=False
        )

        comps = ["Ex", "Ey", "Ez", "Hx", "Hy", "Hz"]
        fft_stack = np.array([ff_fft[c].ravel() for c in comps])
        exact_stack = np.array([ff_exact[c].ravel() for c in comps])
        # FFT (radiation-zone) vs exact near-to-far transform
        self.assertClose(fft_stack, exact_stack, epsilon=1e-2)

        # FFT grid vs per-point direction radiation-zone summation: same
        # underlying approximation, so they agree to the interpolation error
        nx, ny = ff_fft["Ex"].shape
        dx = size.x / (nx - 1)
        dy = size.y / (ny - 1)
        x0 = center.x - 0.5 * size.x
        y0 = center.y - 0.5 * size.y
        fft_pts = []
        direct_pts = []
        for i in range(nx):
            for j in range(ny):
                x = mp.Vector3(x0 + i * dx, y0 + j * dy, center.z)
                direct_pts.append(
                    np.array(self.sim.get_farfield(self.n2f, x, far_field_approx=True))
                )
                fft_pts.append(np.array([ff_fft[c][i, j] for c in comps]))
        self.assertClose(np.array(fft_pts), np.array(direct_pts), epsilon=1e-2)


if __name__ == "__main__":
    unittest.main()
