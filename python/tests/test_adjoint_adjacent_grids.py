"""Each design region's gradient must come only from its own MaterialGrid.

Design regions may be adjacent, separated, or overlapping; in every arrangement
the adjoint gradient of one region must agree with a finite difference of the
objective with respect to that region's own design variables.
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
FCEN = 1 / 1.55
SILICON = mp.Medium(index=3.4)
OXIDE = mp.Medium(index=1.44)


def _problem(grid_shapes, gap=0.0, overlap=False, grid_type="U_DEFAULT"):
    """Design region(s) stacked in y, optionally separated by `gap`.

    With `overlap`, the grids share one volume and only the first is a design
    region.
    """
    if overlap:
        boxes = [(-LY / 2, LY / 2)] * len(grid_shapes)
    else:
        if len(grid_shapes) == 1:
            boxes = [(-LY / 2, LY / 2)]
        else:
            boxes = [(-LY / 2, -gap / 2), (gap / 2, LY / 2)]

    grids = [
        mp.MaterialGrid(
            mp.Vector3(*s),
            OXIDE,
            SILICON,
            # Non-design grids are filled at mid density: an almost-empty cell
            # scatters too weakly for the gradient to converge at this decay_by.
            weights=np.full(s, 0.5) if overlap else np.zeros(s),
            do_averaging=False,
            beta=0,
            eta=0.5,
            grid_type=grid_type,
        )
        for s in grid_shapes
    ]
    vols = [
        mp.Volume(center=mp.Vector3(0, 0.5 * (a + b)), size=mp.Vector3(LX, b - a))
        for a, b in boxes
    ]
    sim = mp.Simulation(
        cell_size=mp.Vector3(CELL, CELL),
        resolution=RESOLUTION,
        boundary_layers=[mp.PML(DPML)],
        geometry=[
            mp.Block(center=v.center, size=v.size, material=g)
            for v, g in zip(vols, grids)
        ],
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
    opt = mpa.OptimizationProblem(
        simulation=sim,
        objective_functions=[lambda ez: npa.sum(npa.abs(ez) ** 2, axis=1)],
        objective_arguments=[mpa.FourierFields(sim, monitor, mp.Ez, yee_grid=False)],
        design_regions=(
            [mpa.DesignRegion(grids[0], volume=vols[0])]
            if overlap
            else [mpa.DesignRegion(g, volume=v) for g, v in zip(grids, vols)]
        ),
        frequencies=[FCEN],
        decay_by=1e-9,
    )
    return opt


def _directional_ratio(opt, grid_shapes, seed=0, step=3e-3, n_regions=None):
    """adjoint / finite difference along one random direction.

    Directional rather than single-pixel: perturbing one design variable moves
    the objective by about as much as the DFT convergence floor.
    """
    rng = np.random.default_rng(seed)
    if n_regions is None:
        n_regions = len(grid_shapes)
    weights = [rng.uniform(0.3, 0.7, int(np.prod(s))) for s in grid_shapes[:n_regions]]
    _, gradient = opt(rho_vector=[w.copy() for w in weights])
    if n_regions == 1:
        gradient = [gradient]
    gradient = [np.asarray(g).ravel() for g in gradient]

    def objective(ws):
        value, _ = opt(rho_vector=[w.copy() for w in ws], need_gradient=False)
        return float(np.ravel(value)[0])

    direction = [rng.standard_normal(w.size) for w in weights]
    norm = np.sqrt(sum(float(np.sum(d**2)) for d in direction))
    direction = [d / norm for d in direction]

    adjoint = sum(float(np.dot(gradient[k], direction[k])) for k in range(len(weights)))
    finite_difference = (
        objective([weights[k] + step * direction[k] for k in range(len(weights))])
        - objective([weights[k] - step * direction[k] for k in range(len(weights))])
    ) / (2 * step)
    return adjoint / finite_difference


class TestAdjacentDesignGrids(unittest.TestCase):
    # What the field solve reproduces across the two runs of the central
    # difference, not what the gradient assembly contributes.
    tol = 2e-3

    def test_single_design_region(self):
        """A single design region covering the whole area."""
        shapes = [(8, 16)]
        ratio = _directional_ratio(_problem(shapes), shapes)
        self.assertAlmostEqual(ratio, 1.0, delta=self.tol)

    def test_separated_design_regions(self):
        """Two design regions separated by a pixel of background."""
        shapes = [(8, 8), (8, 13)]
        ratio = _directional_ratio(_problem(shapes, gap=1.0 / RESOLUTION), shapes)
        self.assertAlmostEqual(ratio, 1.0, delta=self.tol)

    def test_adjacent_design_regions(self):
        """Two design regions sharing a boundary, with equal grid sizes."""
        shapes = [(8, 8), (8, 8)]
        ratio = _directional_ratio(_problem(shapes, gap=0.0), shapes)
        self.assertAlmostEqual(ratio, 1.0, delta=self.tol)

    def test_adjacent_design_regions_unequal_sizes(self):
        """Two design regions sharing a boundary, with unequal grid sizes.

        The larger neighbour's grid_size must not be used to index this
        region's smaller gradient array.
        """
        shapes = [(8, 8), (8, 13)]
        ratio = _directional_ratio(_problem(shapes, gap=0.0), shapes)
        self.assertAlmostEqual(ratio, 1.0, delta=self.tol)

    def test_overlapping_grids_are_unaffected(self):
        """Two grids overlapping in one volume, combined with U_MEAN."""
        shapes = [(8, 8), (8, 8)]
        ratio = _directional_ratio(
            _problem(shapes, overlap=True, grid_type="U_MEAN"), shapes, n_regions=1
        )
        self.assertAlmostEqual(ratio, 1.0, delta=self.tol)

    def test_overlapping_grids_unequal_sizes(self):
        """An overlapping grid must not supply the owner's interpolation size."""
        # Keep the owner larger than the top grid so the regression is a wrong
        # derivative rather than a test that relies on an out-of-bounds access.
        shapes = [(8, 13), (8, 8)]
        ratio = _directional_ratio(
            _problem(shapes, overlap=True, grid_type="U_MEAN"), shapes, n_regions=1
        )
        self.assertAlmostEqual(ratio, 1.0, delta=self.tol)


if __name__ == "__main__":
    unittest.main()
