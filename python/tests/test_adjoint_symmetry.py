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
# What was a bug is the scaling when the *objective* sits on or spans the
# symmetry plane. `fourier_sourcedata` divided each quarter of an interpolated
# adjoint source by the multiplicity of the voxel *centre* rather than of the
# Yee site it was written to. Those differ exactly on a symmetry plane, where a
# point is its own image and so has multiplicity 1 while the centre half a pixel
# away has 2, leaving the source short there. The on- and crossing-plane tests
# cover it.
#
# The converse case was wrong too. A component half a pixel off the plane -- Hz
# or Ey under Mirror(Y) -- is stored on *both* sides of it, each side visited by
# one pass of the symmetry loop, yet each was still divided by the full
# multiplicity, so the row next to the plane got half its source: half the
# gradient for a point objective there. The TE tests cover that.

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


# The field each polarization is driven and measured in, and its eigenvalue
# under Mirror(Y): Ez sits on the plane's Yee rows, Hz half a pixel off them.
POLARIZATIONS = {"TM": (mp.Ez, +1), "TE": (mp.Hz, -1)}


def problem(weights, use_symmetry, monitor_center, monitor_size, polarization="TM"):
    component, phase = POLARIZATIONS[polarization]
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
        # Ez and Hz are both even about y = 0; the phase is that of Hz, a
        # pseudovector, under the mirror.
        symmetries=[mp.Mirror(mp.Y, phase=phase)] if use_symmetry else [],
        sources=[
            mp.Source(
                mp.GaussianSource(FCEN, fwidth=0.2),
                component,
                center=mp.Vector3(-2.0, 0),
                size=mp.Vector3(0, 3.0),
            )
        ],
        force_complex_fields=True,
    )
    monitor = mpa.FourierFields(
        sim, mp.Volume(center=monitor_center, size=monitor_size), component
    )
    return mpa.OptimizationProblem(
        simulation=sim,
        objective_functions=[lambda f: npa.sum(npa.abs(f) ** 2)],
        objective_arguments=[monitor],
        design_regions=[region],
        frequencies=[FCEN],
    )


def objective(weights, use_symmetry, center, size, polarization="TM"):
    result = problem(weights, use_symmetry, center, size, polarization)(
        [weights.ravel()], need_gradient=False
    )
    if isinstance(result, tuple):
        result = result[0]
    return float(np.real(np.squeeze(np.asarray(result, dtype=complex))))


def gradient(weights, use_symmetry, center, size, polarization="TM"):
    value, grad = problem(weights, use_symmetry, center, size, polarization)(
        [weights.ravel()]
    )
    return (
        float(np.real(np.squeeze(np.asarray(value, dtype=complex)))),
        np.asarray(grad).reshape(N, N),
    )


def finite_difference(weights, use_symmetry, center, size, ix, iy, polarization="TM"):
    plus, minus = weights.copy(), weights.copy()
    plus[ix, iy] += FD_STEP
    minus[ix, iy] -= FD_STEP
    return (
        objective(plus, use_symmetry, center, size, polarization)
        - objective(minus, use_symmetry, center, size, polarization)
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
            self.assertAlmostEqual(
                g[ix, iy] / fd, 1.0, places=2, msg=f"node ({ix},{iy})"
            )

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

    def test_gradient_matches_finite_difference_on_plane(self):
        """An objective sitting on the symmetry plane.

        This gave exactly 3/4 of the correct gradient while the adjoint source
        was divided by the voxel centre's multiplicity: of the source's three Yee
        rows, the one on the plane carries half the weight and was the one being
        halved, so 1/4 of the total went missing regardless of resolution.
        """
        w = self.weights
        ix, iy = 8, 15
        _, g = gradient(w, True, *ON_PLANE)
        fd = finite_difference(w, True, *ON_PLANE, ix, iy)
        self.assertAlmostEqual(g[ix, iy] / fd, 1.0, places=2)

    def test_gradient_matches_finite_difference_crossing_plane(self):
        """An extended monitor spanning the plane -- the usual arrangement.

        Same cause as the on-plane case, diluted over the rows the monitor
        covers, so it used to show up as a first-order error: 4.0% at resolution
        20, 2.9% at 30, 1.8% at 40.
        """
        w = self.weights
        ix, iy = 8, 15
        _, g = gradient(w, True, *CROSSING)
        fd = finite_difference(w, True, *CROSSING, ix, iy)
        self.assertAlmostEqual(g[ix, iy] / fd, 1.0, places=2)

    def test_te_gradient_matches_finite_difference(self):
        """Hz, half a pixel off the plane, for an objective on it and across it.

        On the plane this was exactly half; across it, short by the one row of
        the monitor next to the plane.
        """
        w = self.weights
        ix, iy = 8, 15
        for name, place in (("on plane", ON_PLANE), ("crossing", CROSSING)):
            with self.subTest(monitor=name):
                _, g = gradient(w, True, *place, polarization="TE")
                fd = finite_difference(w, True, *place, ix, iy, polarization="TE")
                self.assertAlmostEqual(g[ix, iy] / fd, 1.0, places=2)


if __name__ == "__main__":
    unittest.main()
