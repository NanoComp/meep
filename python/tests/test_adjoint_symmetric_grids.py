"""A MaterialGrid shared by several transformed copies of itself.

Placing one MaterialGrid in several overlapping Blocks, combined with U_MEAN, is
how a design is constrained to a symmetry: a design variable then backs every
copy at once. The adjoint gradient must account for all of them, so it has to
agree with a finite difference taken with respect to that shared variable.
"""

import unittest

import numpy as np
from autograd import numpy as npa

import meep as mp
import meep.adjoint as mpa

RESOLUTION = 20
DPML = 1.0
CELL = 6.0
LX, LY = 1.0, 1.0
N = 8
FCEN = 1 / 1.55
SILICON = mp.Medium(index=3.4)
OXIDE = mp.Medium(index=1.44)


def _problem(transform):
    """One design grid, optionally overlapped with a transformed copy of itself.

    `transform` is None, "rot90" or "mirror". The copy reuses the same
    MaterialGrid object, which is what makes it the same design variable.
    """
    grid = mp.MaterialGrid(
        mp.Vector3(N, N),
        OXIDE,
        SILICON,
        weights=np.zeros((N, N)),
        do_averaging=False,
        beta=0,
        eta=0.5,
        grid_type="U_MEAN",
    )
    volume = mp.Volume(center=mp.Vector3(), size=mp.Vector3(LX, LY))
    blocks = [mp.Block(center=volume.center, size=volume.size, material=grid)]
    if transform == "rot90":
        blocks.append(
            mp.Block(
                center=volume.center,
                size=volume.size,
                material=grid,
                e1=mp.Vector3(1).rotate(mp.Vector3(z=1), np.pi / 2),
                e2=mp.Vector3(y=1).rotate(mp.Vector3(z=1), np.pi / 2),
            )
        )
    elif transform == "mirror":
        blocks.append(
            mp.Block(
                center=volume.center,
                size=volume.size,
                material=grid,
                e1=mp.Vector3(-1, 0, 0),
                e2=mp.Vector3(0, 1, 0),
            )
        )

    sim = mp.Simulation(
        cell_size=mp.Vector3(CELL, CELL),
        resolution=RESOLUTION,
        boundary_layers=[mp.PML(DPML)],
        geometry=blocks,
        sources=[
            mp.Source(
                mp.GaussianSource(FCEN, fwidth=0.1),
                component=mp.Ez,
                center=mp.Vector3(-CELL / 2 + DPML + 0.3, 0),
                size=mp.Vector3(0, CELL - 2 * DPML),
            )
        ],
        default_material=OXIDE,
        dimensions=2,
    )
    monitor = mp.Volume(
        center=mp.Vector3(CELL / 2 - DPML - 0.3, 0), size=mp.Vector3(0, 2.0)
    )
    return mpa.OptimizationProblem(
        simulation=sim,
        objective_functions=[lambda ez: npa.sum(npa.abs(ez) ** 2, axis=1)],
        objective_arguments=[mpa.FourierFields(sim, monitor, mp.Ez, yee_grid=False)],
        design_regions=[mpa.DesignRegion(grid, volume=volume)],
        frequencies=[FCEN],
        decay_by=1e-9,
    )


def _directional_ratio(opt, seed=0, step=3e-3):
    """adjoint / finite difference along one random direction.

    Directional rather than single-pixel: perturbing one design variable moves
    the objective by about as much as the DFT convergence floor.
    """
    rng = np.random.default_rng(seed)
    weights = rng.uniform(0.3, 0.7, N * N)
    _, gradient = opt(rho_vector=[weights.copy()])
    gradient = np.asarray(gradient).ravel()

    def objective(w):
        value, _ = opt(rho_vector=[w.copy()], need_gradient=False)
        return float(np.ravel(value)[0])

    direction = rng.standard_normal(weights.size)
    direction /= np.linalg.norm(direction)
    adjoint = float(np.dot(gradient, direction))
    finite_difference = (
        objective(weights + step * direction) - objective(weights - step * direction)
    ) / (2 * step)
    return adjoint / finite_difference


class TestSymmetricDesignGrids(unittest.TestCase):
    # What the field solve reproduces across the two runs of the central
    # difference, not what the gradient assembly contributes.
    tol = 2e-3

    def test_no_transform(self):
        """A single block, with no overlapping copy."""
        self.assertAlmostEqual(_directional_ratio(_problem(None)), 1.0, delta=self.tol)

    def test_rotated_copy(self):
        """Overlapped with a 90-degree rotation of itself."""
        self.assertAlmostEqual(
            _directional_ratio(_problem("rot90")), 1.0, delta=self.tol
        )

    def test_mirrored_copy(self):
        """Overlapped with a reflection of itself."""
        self.assertAlmostEqual(
            _directional_ratio(_problem("mirror")), 1.0, delta=self.tol
        )


if __name__ == "__main__":
    unittest.main()
