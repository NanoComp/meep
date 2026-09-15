import unittest
import numpy as np
from autograd import numpy as npa
import meep as mp
import meep.adjoint as mpa

RESOLUTION = 20
FCEN = 1 / 1.55
N = 12
DESIGN = 1.2
CELL_X, CELL_Y, DPML = 8.0, 6.0, 1.0
WG_W = 0.5
OXIDE = mp.Medium(index=1.44)
SILICON = mp.Medium(index=3.48)
FD_STEP = 1e-4


def symmetric_weights(seed=0):
    rng = np.random.default_rng(seed)
    half = rng.uniform(0.3, 0.7, size=(N, N // 2))
    return np.concatenate([half[:, ::-1], half], axis=1)


def problem(weights, use_symmetry):
    grid = mp.MaterialGrid(mp.Vector3(N, N), OXIDE, SILICON, weights=weights,
                           do_averaging=False, grid_type="U_MEAN")
    region = mpa.DesignRegion(
        grid, volume=mp.Volume(center=mp.Vector3(), size=mp.Vector3(DESIGN, DESIGN)))
    sim = mp.Simulation(
        cell_size=mp.Vector3(CELL_X, CELL_Y), resolution=RESOLUTION,
        boundary_layers=[mp.PML(DPML)], default_material=OXIDE,
        geometry=[mp.Block(center=mp.Vector3(), size=mp.Vector3(mp.inf, WG_W),
                           material=SILICON),
                  mp.Block(center=region.center, size=region.size, material=grid)],
        # Ez is even about y = 0, so +1 is this mode's mirror eigenvalue.
        symmetries=[mp.Mirror(mp.Y, phase=+1)] if use_symmetry else [],
        sources=[mp.EigenModeSource(mp.GaussianSource(FCEN, fwidth=0.2 * FCEN),
                                    center=mp.Vector3(-2.5, 0),
                                    size=mp.Vector3(0, CELL_Y - 2 * DPML),
                                    eig_band=1, eig_parity=mp.ODD_Z + mp.EVEN_Y)])
    mon = mpa.EigenmodeCoefficient(
        sim, mp.Volume(center=mp.Vector3(2.5, 0), size=mp.Vector3(0, CELL_Y - 2 * DPML)),
        mode=1)
    return mpa.OptimizationProblem(
        simulation=sim, objective_functions=[lambda a: npa.abs(a) ** 2],
        objective_arguments=[mon], design_regions=[region], frequencies=[FCEN])


def value_and_grad(w, sym):
    v, g = problem(w, sym)([w.ravel()])
    return float(np.real(np.squeeze(v))), np.asarray(g).reshape(N, N)


def objective(w, sym):
    r = problem(w, sym)([w.ravel()], need_gradient=False)
    if isinstance(r, tuple):
        r = r[0]
    return float(np.real(np.squeeze(np.asarray(r, dtype=complex))))


class TestEigenmodeCoefficientSymmetry(unittest.TestCase):
    """A mode monitor under Mirror, which FourierFields coverage does not reach.

    The distinction matters: a mirror in y flips y-components, so a mode's
    eigenvalue is not always +1, and declaring the wrong one silently projects
    onto the orthogonal sector rather than failing.  Here Ez is even, so +1 is
    correct, and everything must be invariant to declaring it.
    """

    @classmethod
    def setUpClass(cls):
        cls.w = symmetric_weights()

    def test_objective_invariant_to_declaring_the_symmetry(self):
        w = self.w
        self.assertAlmostEqual(objective(w, True) / objective(w, False), 1.0, places=4)

    def test_folded_gradient_invariant_to_declaring_the_symmetry(self):
        """The *folded* gradient is what must agree, not the raw one.

        Under the symmetry Meep reads epsilon only on y >= 0, so weights below
        the plane are inert and a weight above it drives both halves: the raw
        gradient comes back doubled on the upper half and zero on the lower.
        Folding the two halves together removes that difference and leaves the
        derivative with respect to a variable that owns both mirrored cells.
        """
        def fold(g):
            return g[:, N // 2:] + g[:, : N // 2][:, ::-1]

        w = self.w
        _, gs = value_and_grad(w, True)
        _, gn = value_and_grad(w, False)
        fs, fn = fold(gs), fold(gn)
        denom = np.abs(fn).max()
        self.assertGreater(denom, 0.0)
        self.assertLess(np.abs(fs - fn).max() / denom, 5e-2)

    def test_gradient_matches_finite_difference_under_symmetry(self):
        w = self.w
        _, g = value_and_grad(w, True)
        for ix, iy in ((6, 8), (4, 9)):
            p, m = w.copy(), w.copy()
            p[ix, iy] += FD_STEP
            m[ix, iy] -= FD_STEP
            fd = (objective(p, True) - objective(m, True)) / (2 * FD_STEP)
            if abs(fd) < 1e-9:
                continue
            self.assertAlmostEqual(g[ix, iy] / fd, 1.0, places=1, msg=f"({ix},{iy})")


if __name__ == "__main__":
    unittest.main()
