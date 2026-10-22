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

# The fundamental mode's parity and its eigenvalue under Mirror(Y).  A mirror
# in y also flips the sign of y-components, so the TE mode, whose Ey is even in
# y, has eigenvalue -1, while the TM mode (Ez even) has +1.
POLARIZATIONS = {
    "TM": dict(parity=mp.ODD_Z + mp.EVEN_Y, phase=+1),
    "TE": dict(parity=mp.EVEN_Z + mp.ODD_Y, phase=-1),
}


def symmetric_weights(seed=0):
    rng = np.random.default_rng(seed)
    half = rng.uniform(0.3, 0.7, size=(N, N // 2))
    return np.concatenate([half[:, ::-1], half], axis=1)


def problem(weights, polarization, phase=None):
    """`phase` is the declared mirror eigenvalue, or None for no symmetry."""
    parity = POLARIZATIONS[polarization]["parity"]
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
        cell_size=mp.Vector3(CELL_X, CELL_Y),
        resolution=RESOLUTION,
        boundary_layers=[mp.PML(DPML)],
        default_material=OXIDE,
        geometry=[
            mp.Block(
                center=mp.Vector3(), size=mp.Vector3(mp.inf, WG_W), material=SILICON
            ),
            mp.Block(center=region.center, size=region.size, material=grid),
        ],
        symmetries=[] if phase is None else [mp.Mirror(mp.Y, phase=phase)],
        sources=[
            mp.EigenModeSource(
                mp.GaussianSource(FCEN, fwidth=0.2 * FCEN),
                center=mp.Vector3(-2.5, 0),
                size=mp.Vector3(0, CELL_Y - 2 * DPML),
                eig_band=1,
                eig_parity=parity,
            )
        ],
    )
    mon = mpa.EigenmodeCoefficient(
        sim,
        mp.Volume(center=mp.Vector3(2.5, 0), size=mp.Vector3(0, CELL_Y - 2 * DPML)),
        mode=1,
        eig_parity=parity,
    )
    return mpa.OptimizationProblem(
        simulation=sim,
        objective_functions=[lambda a: npa.abs(a) ** 2],
        objective_arguments=[mon],
        design_regions=[region],
        frequencies=[FCEN],
    )


def value_and_grad(w, polarization, phase=None):
    v, g = problem(w, polarization, phase)([w.ravel()])
    return float(np.real(np.squeeze(v))), np.asarray(g).reshape(N, N)


def objective(w, polarization, phase=None):
    v, _ = problem(w, polarization, phase)([w.ravel()], need_gradient=False)
    return float(np.real(np.squeeze(v)))


def fold(g):
    """Combine each pair of mirrored cells (y is the second axis)."""
    return g[:, N // 2 :] + g[:, : N // 2][:, ::-1]


class TestEigenmodeCoefficientSymmetry(unittest.TestCase):
    """A mode monitor under Mirror, which FourierFields coverage does not reach.

    A mode's mirror eigenvalue is not always +1, and declaring the wrong one
    silently projects onto the orthogonal sector rather than failing.  With the
    right one, everything must be invariant to declaring it.
    """

    @classmethod
    def setUpClass(cls):
        cls.w = symmetric_weights()

    def test_objective_invariant_to_declaring_the_symmetry(self):
        for polarization, p in POLARIZATIONS.items():
            with self.subTest(polarization=polarization):
                self.assertAlmostEqual(
                    objective(self.w, polarization, p["phase"])
                    / objective(self.w, polarization),
                    1.0,
                    places=8,
                )

    def test_wrong_eigenvalue_changes_the_objective(self):
        # The trap: no error, just a different answer.
        right = objective(self.w, "TE")
        wrong = objective(self.w, "TE", phase=+1)
        self.assertGreater(abs(wrong / right - 1), 0.1)

    def test_folded_gradient_invariant_to_declaring_the_symmetry(self):
        """The *folded* gradient is what must agree, not the raw one.

        Under the symmetry Meep reads epsilon only on y >= 0, so weights below
        the plane are inert and a weight above it drives both halves: the raw
        gradient comes back doubled on the upper half and zero on the lower.
        Folding the two halves together removes that difference and leaves the
        derivative with respect to a variable that owns both mirrored cells.
        """
        for polarization, p in POLARIZATIONS.items():
            with self.subTest(polarization=polarization):
                _, gs = value_and_grad(self.w, polarization, p["phase"])
                _, gn = value_and_grad(self.w, polarization)
                fs, fn = fold(gs), fold(gn)
                denom = np.abs(fn).max()
                self.assertGreater(denom, 0.0)
                self.assertLess(np.abs(fs - fn).max() / denom, 1e-6)

    def test_gradient_matches_finite_difference_under_symmetry(self):
        # cells above the plane, which are the ones Meep reads
        for polarization, p in POLARIZATIONS.items():
            _, g = value_and_grad(self.w, polarization, p["phase"])
            for ix, iy in ((6, 8), (4, 9)):
                with self.subTest(polarization=polarization, cell=(ix, iy)):
                    plus, minus = self.w.copy(), self.w.copy()
                    plus[ix, iy] += FD_STEP
                    minus[ix, iy] -= FD_STEP
                    fd = (
                        objective(plus, polarization, p["phase"])
                        - objective(minus, polarization, p["phase"])
                    ) / (2 * FD_STEP)
                    self.assertGreater(abs(fd), 1e-6, "cell must affect the objective")
                    self.assertAlmostEqual(g[ix, iy] / fd, 1.0, places=3)


if __name__ == "__main__":
    unittest.main()
