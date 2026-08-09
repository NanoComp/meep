"""Regression coverage for MaterialGrid subpixel smoothing in 3d.

Every adjoint test upstream runs in 2d or cylindrical coordinates, which left two
3d-only defects unnoticed:

  1. get_uproj_w()'s D3 branch normalized by `4 / 3 * pi * rad^3`. `4 / 3` is
     integer division, so the denominator was `pi * rad^3` and the smoothing
     kernel integrated to 4/3 instead of 1, scaling epsilon everywhere the
     material grid had a nonzero gradient.

  2. set_materials_from_geom_epsilon() stored the raw `maxeval` on the
     geom_epsilon regardless of `use_anisotropic_averaging`, while
     structure_chunk::set_chi1inv() zeroes it when averaging is off. With
     eps_averaging=False the forward solve therefore saw an unsmoothed operator
     while the adjoint differentiated a smoothed one.
"""

import unittest

import numpy as np

import meep as mp


# Allow finite-difference and MPI chunking noise while remaining far below the
# 70-100% disagreement caused by the regression covered by these tests.
FD_REL_TOL = 1e-2


def _epsilon(dim, do_averaging, n1, n2, resolution=20, n=20):
    """Dielectric array over a MaterialGrid carrying a linear ramp in u."""
    design, pad = 1.0, 0.6
    sz = 0 if dim == 2 else design + 2 * pad

    ramp = np.linspace(0.25, 0.75, n)
    weights = np.repeat(ramp[:, None], n, axis=1)[:, :, None]

    grid = mp.MaterialGrid(
        mp.Vector3(n, n, 1),
        mp.Medium(index=n1),
        mp.Medium(index=n2),
        weights=weights,
        do_averaging=do_averaging,
        beta=0,
    )
    sim = mp.Simulation(
        cell_size=mp.Vector3(design + 2 * pad, design + 2 * pad, sz),
        resolution=resolution,
        default_material=mp.Medium(index=n1),
        geometry=[
            mp.Block(
                center=mp.Vector3(),
                size=mp.Vector3(design, design, mp.inf if dim == 2 else design),
                material=grid,
            )
        ],
        dimensions=dim,
        eps_averaging=True,
    )
    sim.init_sim()
    extent = design * 0.8
    return np.asarray(
        sim.get_array(
            component=mp.Dielectric,
            center=mp.Vector3(),
            size=mp.Vector3(extent, extent, 0 if dim == 2 else extent),
        )
    )


class TestSubpixelKernelNormalization(unittest.TestCase):
    """The smoothing kernel must integrate to 1 in every dimensionality.

    Measuring that through eps needs care: smoothing produces an anisotropic
    tensor whose normal component is a harmonic mean and whose transverse
    components are arithmetic means. With real contrast those differ, and the
    resulting (legitimate) shift swamps the normalization signal. Removing the
    contrast collapses both means onto the same value, leaving

        eps(do_averaging=True) / eps(do_averaging=False) == integral(w)

    so the ratio reads the normalization off directly.
    """

    def test_kernel_normalized_2d(self):
        ratio = _epsilon(2, True, 1.44, 1.4414) / _epsilon(2, False, 1.44, 1.4414)
        self.assertAlmostEqual(float(np.median(ratio)), 1.0, places=6)

    def test_kernel_normalized_3d(self):
        ratio = _epsilon(3, True, 1.44, 1.4414) / _epsilon(3, False, 1.44, 1.4414)
        # Before the fix this was 18/17 = 1.058823..., the closed form of
        # 3 / (K + 2/K) at K = 4/3.
        self.assertAlmostEqual(float(np.median(ratio)), 1.0, places=6)


def _directional_fd(do_averaging, eps_averaging, resolution=12, n=6, dp=1e-3, seed=0):
    """Return (fd, adjoint) directional derivatives along the adjoint gradient."""
    import autograd.numpy as npa
    import meep.adjoint as mpa

    si, clad = mp.Medium(index=3.48), mp.Medium(index=1.44)
    pml, port_pad, side_pad = 0.5, 0.8, 0.5
    design, thickness = 1.0, 0.22

    sx = 2 * pml + 2 * port_pad + design
    sy = 2 * pml + 2 * side_pad + design
    sz = 2 * pml + 2 * side_pad + thickness

    rng = np.random.default_rng(seed)
    weights = rng.uniform(0.2, 0.8, size=n * n)

    grid = mp.MaterialGrid(
        mp.Vector3(n, n, 1),
        clad,
        si,
        weights=weights.reshape(n, n, 1),
        do_averaging=do_averaging,
    )
    region = mpa.DesignRegion(
        grid,
        volume=mp.Volume(
            center=mp.Vector3(), size=mp.Vector3(design, design, thickness)
        ),
    )

    fcen = 1 / 1.55
    port = mp.Vector3(0, sy - 2 * pml, sz - 2 * pml)
    sim = mp.Simulation(
        cell_size=mp.Vector3(sx, sy, sz),
        resolution=resolution,
        boundary_layers=[mp.PML(pml)],
        default_material=clad,
        geometry=[
            mp.Block(
                center=mp.Vector3(),
                size=mp.Vector3(mp.inf, 0.5, thickness),
                material=si,
            ),
            mp.Block(center=region.center, size=region.size, material=grid),
        ],
        sources=[
            mp.EigenModeSource(
                mp.GaussianSource(fcen, fwidth=0.1 * fcen),
                center=mp.Vector3(-(design / 2 + port_pad / 2)),
                size=port,
                eig_band=1,
            )
        ],
        eps_averaging=eps_averaging,
    )
    monitor = mpa.EigenmodeCoefficient(
        sim,
        mp.Volume(center=mp.Vector3(design / 2 + port_pad / 2), size=port),
        mode=1,
    )
    opt = mpa.OptimizationProblem(
        simulation=sim,
        objective_functions=[lambda c: npa.abs(c) ** 2],
        objective_arguments=[monitor],
        design_regions=[region],
        frequencies=[fcen],
        decay_by=1e-6,
    )

    _, gradient = opt([weights], need_gradient=True)
    gradient = np.asarray(np.real(np.squeeze(gradient)), dtype=np.float64).reshape(-1)

    # Along the gradient the directional derivative is as large as it gets, which
    # keeps the difference quotient clear of the forward solve's convergence noise.
    direction = gradient / np.linalg.norm(gradient)
    f_plus, _ = opt([weights + dp * direction], need_gradient=False)
    f_minus, _ = opt([weights - dp * direction], need_gradient=False)
    fd = (float(np.squeeze(f_plus)) - float(np.squeeze(f_minus))) / (2 * dp)
    return fd, float(gradient @ direction)


class TestAdjointGradient3D(unittest.TestCase):
    def test_gradient_matches_fd_without_smoothing(self):
        fd, adj = _directional_fd(do_averaging=False, eps_averaging=False)
        self.assertAlmostEqual(fd / adj, 1.0, delta=FD_REL_TOL)

    def test_gradient_matches_fd_with_smoothing(self):
        fd, adj = _directional_fd(do_averaging=True, eps_averaging=True)
        self.assertAlmostEqual(fd / adj, 1.0, delta=FD_REL_TOL)

    def test_do_averaging_ignored_when_eps_averaging_off(self):
        """do_averaging=True with eps_averaging=False must not change the gradient.

        The forward solve discards MaterialGrid smoothing when eps_averaging is
        off, so the gradient has to discard it too. Before the fix the adjoint
        kept smoothing and this pair disagreed by ~70-100%.
        """
        fd, adj = _directional_fd(do_averaging=True, eps_averaging=False)
        self.assertAlmostEqual(fd / adj, 1.0, delta=FD_REL_TOL)


if __name__ == "__main__":
    unittest.main()
