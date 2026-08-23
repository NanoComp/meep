"""Adjacent MaterialGrid design regions must each get their own gradient.

`material_grids_addgradient_point` resolves the material grid by a spatial
lookup at each point (`geom_tree_search`) and then interpolates into the
supplied `v[]` using *that* grid's `grid_size` and box coordinates. A design
region's DFT monitor is snapped outward to enclosing grid nodes, so when two
material grids abut, one region's monitor covers a row of nodes that lie in the
neighbouring grid's object. Those nodes resolve to the neighbour, and their
contribution was interpolated in the neighbour's coordinates and summed into
this region's array.

Instrumented on the pre-fix build, the lower of two abutting 21x11-node regions
touched 210 nodes of its own grid and 21 of the neighbour's -- exactly the
shared boundary row -- and its directional derivative came out 35% high. With
the two grids sized differently the stray writes also ran off the end of the
smaller region's array, corrupting the heap.

The adjacent case is covered twice: with equal grid sizes, where a regression is
a wrong value and fails as an assertion, and with unequal sizes, where it is also
an out-of-bounds write.
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

    With `overlap`, the grids are instead placed on top of each other in the same
    volume -- the configuration used to impose symmetries -- and only the first
    is a design region.
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
            weights=np.zeros(s),
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

    A directional derivative rather than a single-pixel difference: perturbing
    one design variable moves the objective by about as much as the DFT
    convergence floor, whereas moving all of them is reproducible to five
    digits.
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
    # What the field solve reproduces across the two extra runs of the central
    # difference, not what the gradient assembly contributes. The pre-fix error
    # was 35%, three orders of magnitude outside this.
    tol = 2e-3

    def test_single_design_region(self):
        """Baseline: one grid over the whole area, no neighbour to confuse it."""
        shapes = [(8, 16)]
        ratio = _directional_ratio(_problem(shapes), shapes)
        self.assertAlmostEqual(ratio, 1.0, delta=self.tol)

    def test_separated_design_regions(self):
        """Two grids with a pixel of background between them: never affected."""
        shapes = [(8, 8), (8, 13)]
        ratio = _directional_ratio(_problem(shapes, gap=1.0 / RESOLUTION), shapes)
        self.assertAlmostEqual(ratio, 1.0, delta=self.tol)

    def test_adjacent_design_regions(self):
        """Two grids sharing a boundary. This is the regression.

        Equal grid sizes, so a regression shows up as a wrong value (1.35 on the
        pre-fix build) rather than as an out-of-bounds write -- keep the failure
        mode a clean assertion so the rest of the file still runs.
        """
        shapes = [(8, 8), (8, 8)]
        ratio = _directional_ratio(_problem(shapes, gap=0.0), shapes)
        self.assertAlmostEqual(ratio, 1.0, delta=self.tol)

    def test_adjacent_design_regions_unequal_sizes(self):
        """Same, with grids of different sizes.

        The stray contributions are interpolated with the neighbour's
        `grid_size`, so when the neighbour is larger they also run off the end of
        this region's array. Whether that corruption is fatal depends on the
        allocator's state: on the pre-fix build this was seen both as a wrong
        value (1.06) and as an abort in malloc.
        """
        shapes = [(8, 8), (8, 13)]
        ratio = _directional_ratio(_problem(shapes, gap=0.0), shapes)
        self.assertAlmostEqual(ratio, 1.0, delta=self.tol)

    def test_overlapping_grids_are_unaffected(self):
        """Deliberately overlapping grids must behave exactly as before.

        Overlapping a grid with transformed copies of itself, combined with
        U_MEAN, is the documented way to impose a symmetry. The owner filter is
        therefore scoped to U_DEFAULT, where "topmost wins" makes a foreign
        node's derivative genuinely zero. This pins that scoping: without it, an
        earlier revision of the fix zeroed this gradient entirely, and a second
        revision halved it.
        """
        shapes = [(8, 8), (8, 8)]
        ratio = _directional_ratio(
            _problem(shapes, overlap=True, grid_type="U_MEAN"), shapes, n_regions=1
        )
        # Looser than self.tol: overlapping U_MEAN gradients carry a ~0.3%
        # inaccuracy that predates this change and is step-independent, so it is
        # not finite-difference error. It is untouched here -- the guard
        # short-circuits for every kind other than U_DEFAULT, leaving that path
        # byte-identical to before. What this test is for is the 50-100% errors
        # that a wrongly-scoped owner filter produces.
        self.assertAlmostEqual(ratio, 1.0, delta=1e-2)


if __name__ == "__main__":
    unittest.main()
