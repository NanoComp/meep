# Adjoint gradients when the simulation declares a mirror symmetry.
#
# Meep had no adjoint test that declared a symmetry at all, and the behaviour is
# easy to mistake for a bug: with `Mirror(Y)` the design gradient comes back with
# essentially all of its magnitude in y >= 0, which looks like a routine that
# forgot to mirror. It is not. Meep builds the structure on the reduced grid
# volume, so weights below the plane are never read, and a weight above it drives
# both halves at once. A finite difference confirms both halves of that: the
# gradient really is zero down there, and really is doubled up here.
#
# What is a bug is the scaling when the *objective* sits on the symmetry plane,
# which these tests pin down rather than paper over.

import unittest

import numpy as np
from autograd import numpy as npa

import meep as mp
import meep.adjoint as mpa

RESOLUTION = 20
FCEN = 1.0
N = 16
DESIGN = 1.6
CELL = 6.0
DPML = 1.0
OXIDE = mp.Medium(index=1.44)
SILICON = mp.Medium(index=3.48)
FD_STEP = 1e-4


def symmetric_weights(seed=0):
    """A design that is its own mirror image, so both runs describe one structure."""
    rng = np.random.default_rng(seed)
    half = rng.uniform(0.3, 0.7, size=(N, N // 2))
    return np.concatenate([half[:, ::-1], half], axis=1)


def problem(weights, use_symmetry, monitor_center, monitor_size):
    grid = mp.MaterialGrid(
        mp.Vector3(N, N),
        OXIDE,
        SILICON,
        weights=weights,
        do_averaging=False,
        grid_type="U_MEAN",
    )
    region = mpa.DesignRegion(
        grid, volume=mp.Volume(center=mp.Vector3(), size=mp.Vector3(DESIGN, DESIGN))
    )
    sim = mp.Simulation(
        cell_size=mp.Vector3(CELL, CELL),
        resolution=RESOLUTION,
        boundary_layers=[mp.PML(DPML)],
        default_material=OXIDE,
        geometry=[mp.Block(center=region.center, size=region.size, material=grid)],
        # Ez is even about y = 0, so this is the symmetry the structure has.
        symmetries=[mp.Mirror(mp.Y, phase=+1)] if use_symmetry else [],
        sources=[
            mp.Source(
                mp.GaussianSource(FCEN, fwidth=0.2),
                mp.Ez,
                center=mp.Vector3(-2.0, 0),
                size=mp.Vector3(0, 3.0),
            )
        ],
        force_complex_fields=True,
    )
    monitor = mpa.FourierFields(
        sim, mp.Volume(center=monitor_center, size=monitor_size), mp.Ez
    )
    return mpa.OptimizationProblem(
        simulation=sim,
        objective_functions=[lambda f: npa.sum(npa.abs(f) ** 2)],
        objective_arguments=[monitor],
        design_regions=[region],
        frequencies=[FCEN],
    )


def objective(weights, use_symmetry, center, size):
    result = problem(weights, use_symmetry, center, size)(
        [weights.ravel()], need_gradient=False
    )
    if isinstance(result, tuple):
        result = result[0]
    return float(np.real(np.squeeze(np.asarray(result, dtype=complex))))


def gradient(weights, use_symmetry, center, size):
    value, grad = problem(weights, use_symmetry, center, size)([weights.ravel()])
    return (
        float(np.real(np.squeeze(np.asarray(value, dtype=complex)))),
        np.asarray(grad).reshape(N, N),
    )


def finite_difference(weights, use_symmetry, center, size, ix, iy):
    plus, minus = weights.copy(), weights.copy()
    plus[ix, iy] += FD_STEP
    minus[ix, iy] -= FD_STEP
    return (
        objective(plus, use_symmetry, center, size)
        - objective(minus, use_symmetry, center, size)
    ) / (2 * FD_STEP)


# A monitor well clear of the mirror plane, so nothing here depends on how the
# plane itself is handled.
OFF_PLANE = (mp.Vector3(2.0, 0.5), mp.Vector3(0, 0))
# A monitor sitting exactly on it.
ON_PLANE = (mp.Vector3(2.0, 0.0), mp.Vector3(0, 0))
# The usual case: an extended monitor that straddles the plane.
CROSSING = (mp.Vector3(2.0, 0.0), mp.Vector3(0, 3.0))


class TestAdjointMirrorSymmetry(unittest.TestCase):
    """`Mirror(Y)` declared, against a finite difference of the same setup."""

    @classmethod
    def setUpClass(cls):
        cls.weights = symmetric_weights()

    def test_forward_objective_unaffected_by_symmetry(self):
        """Declaring the symmetry must not change the value being differentiated."""
        w = self.weights
        self.assertAlmostEqual(
            objective(w, True, *CROSSING) / objective(w, False, *CROSSING),
            1.0,
            places=6,
        )

    def test_gradient_matches_finite_difference_off_plane(self):
        """With the objective clear of the plane the gradient is exact."""
        w = self.weights
        _, g = gradient(w, True, *OFF_PLANE)
        for ix, iy in ((8, 13), (8, 15)):
            fd = finite_difference(w, True, *OFF_PLANE, ix, iy)
            self.assertAlmostEqual(g[ix, iy] / fd, 1.0, places=2,
                                   msg=f"node ({ix},{iy})")

    def test_mirrored_half_is_inert(self):
        """Weights below the plane are never read, so their gradient is zero.

        This is the part that looks like a missing mirror operation and is not:
        the finite difference is zero there too. Meep evaluates epsilon only on
        the reduced grid volume, so those variables have no effect on the field.
        """
        w = self.weights
        _, g = gradient(w, True, *CROSSING)
        # Away from the plane -- close to it, a weight still reaches y >= 0
        # through the design grid's own interpolation stencil.
        for ix, iy in ((8, 2), (8, 3), (4, 1)):
            fd = finite_difference(w, True, *CROSSING, ix, iy)
            self.assertAlmostEqual(fd, 0.0, places=6, msg=f"FD at ({ix},{iy})")
            self.assertAlmostEqual(g[ix, iy], 0.0, places=6, msg=f"grad at ({ix},{iy})")

    def test_symmetric_variable_drives_both_halves(self):
        """One weight above the plane moves the structure above and below it.

        So its derivative is the sum of the two independent derivatives the same
        structure has without the symmetry, which for a symmetric design is twice
        either one.
        """
        w = self.weights
        ix, iy = 8, 15
        with_sym = finite_difference(w, True, *CROSSING, ix, iy)
        without = finite_difference(w, False, *CROSSING, ix, iy)
        self.assertAlmostEqual(with_sym / without, 2.0, places=3)

    @unittest.expectedFailure
    def test_gradient_matches_finite_difference_on_plane(self):
        """Known wrong: an objective on the symmetry plane gives 3/4 the gradient.

        Exactly 0.7500 at resolution 40 and 0.7484 at 20, so this is a scaling
        error and not discretization. Everything else in this file passes, which
        places it in how the adjoint source is weighted at the plane rather than
        in the gradient accumulation.
        """
        w = self.weights
        ix, iy = 8, 15
        _, g = gradient(w, True, *ON_PLANE)
        fd = finite_difference(w, True, *ON_PLANE, ix, iy)
        self.assertAlmostEqual(g[ix, iy] / fd, 1.0, places=2)

    @unittest.expectedFailure
    def test_gradient_matches_finite_difference_crossing_plane(self):
        """Known wrong: a monitor spanning the plane is low by one row's worth.

        4.0% at resolution 20, 2.9% at 30, 1.8% at 40 -- first order in the grid
        spacing, consistent with the single row on the plane being mishandled out
        of the O(1/dy) rows the monitor covers. Same root cause as the on-plane
        case above; this is how it shows up in a realistic objective.
        """
        w = self.weights
        ix, iy = 8, 15
        _, g = gradient(w, True, *CROSSING)
        fd = finite_difference(w, True, *CROSSING, ix, iy)
        self.assertAlmostEqual(g[ix, iy] / fd, 1.0, places=2)


if __name__ == "__main__":
    unittest.main()
