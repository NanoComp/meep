"""The 3d adjoint gradient must not depend on how the cell is split into chunks.

`material_grids_addgradient` assembles the gradient by walking the forward and
adjoint DFT chunks. Before the chunk-independent rewrite it (a) paired forward
and adjoint chunks by list position, which is not a spatial correspondence once
a component is missing from some chunk, (b) skipped cross terms when the two
lists had different lengths, and (c) looked up the neighboring forward node with
chunk-local index arithmetic guarded only against the flat index leaving
[0, N) -- a 3d point outside the chunk can still have an in-range flat index, so
the lookup either aliased to an unrelated voxel or substituted zero.

All three only bite where a chunk boundary crosses the design region, so the
error appears with MPI but is *not* an MPI bug: forcing the same split in a
serial run reproduces it exactly. That is what this test does, so it runs in
every CI job rather than only the MPI ones.

A directional finite-difference check compresses the gradient onto a single
number, and this error survives that: measured on the pre-fix build, probing
along g/|g| gave fd/adjoint = 1.0008 while the gradient vector itself was 44%
off in L2. Comparing the gradient vectors directly is what catches it.

The tolerance is therefore what the *field solve* can reproduce across a
different split, not what the gradient assembly contributes; see the comment on
`tol` below.
"""

import unittest

import numpy as np
from autograd import numpy as npa

import meep as mp
import meep.adjoint as mpa

RESOLUTION = 12
DESIGN_N = 6
SEED = 0


def _gradient(num_chunks):
    """Adjoint gradient of |alpha|^2 for a fixed problem at a forced chunk count."""
    si, clad = mp.Medium(index=3.48), mp.Medium(index=1.44)
    pml, port_pad, side_pad = 0.5, 0.8, 0.5
    design, thickness = 1.0, 0.22

    sx = 2 * pml + 2 * port_pad + design
    sy = 2 * pml + 2 * side_pad + design
    sz = 2 * pml + 2 * side_pad + thickness

    weights = np.random.default_rng(SEED).uniform(0.2, 0.8, size=DESIGN_N * DESIGN_N)
    grid = mp.MaterialGrid(
        mp.Vector3(DESIGN_N, DESIGN_N, 1),
        clad,
        si,
        weights=weights.reshape(DESIGN_N, DESIGN_N, 1),
        do_averaging=False,
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
        resolution=RESOLUTION,
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
        eps_averaging=False,
        num_chunks=num_chunks,
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
    f0, gradient = opt([weights], need_gradient=True)
    return (
        float(np.squeeze(f0)),
        np.asarray(np.real(np.squeeze(gradient)), dtype=np.float64).reshape(-1),
    )


class TestAdjointChunks(unittest.TestCase):
    def test_gradient_independent_of_chunk_division(self):
        # Both counts are >= the process count, so this is a genuine change of
        # chunk division in serial and under MPI alike.
        np_ = mp.count_processors()
        f0_a, g_a = _gradient(num_chunks=np_)
        f0_b, g_b = _gradient(num_chunks=3 * np_)

        # A single-precision build cannot resolve either check to 1e-9. The
        # objective alone -- which no gradient code touches -- already moves by
        # one float ulp (1.2e-7) when the cell is split differently under MPI,
        # and the gradient inherits that. The looser bound is still four orders
        # of magnitude below the defect this test guards against.
        tol = 1e-5 if mp.is_single_precision() else 1e-9

        # the forward solve was never chunk-dependent; if this trips, the test
        # problem itself changed rather than the gradient assembly
        self.assertAlmostEqual(f0_a / f0_b, 1.0, delta=tol)

        rel = np.linalg.norm(g_a - g_b) / np.linalg.norm(g_a)
        # Pre-fix this comparison read 4.7e-1. Post-fix the only difference is
        # summation order in the collective reductions: measured 3e-16 (serial,
        # where there are none) to 2e-12 (np 8) in double precision, 1.5e-7 in
        # single.
        self.assertLess(
            rel, tol, f"gradient changed by {rel:.3e} under a different chunk split"
        )


class TestAdjointSourcePlacement(unittest.TestCase):
    """`create_adjoint_sources` must tolerate a rank that places nothing.

    `fourier_sourcedata` returns one entry per *local* chunk intersecting the
    monitor and `IndexedSource` is inherently per-chunk, so under MPI every rank
    but one places nothing for a point monitor. A bare `assert adjoint_sources`
    therefore made a point objective fail outright on more than one process.

    This is a unit test rather than a simulation because the condition cannot be
    reproduced serially: forcing many chunks does not help, since one rank still
    owns all of them. The contract is what matters -- an empty *local* list is
    legitimate, while a globally zero cotangent is not.
    """

    class _Monitor:
        def __init__(self, sources):
            self._sources = sources

        def place_adjoint_source(self, dJ):
            return list(self._sources)

    def test_empty_local_placement_is_allowed(self):
        monitor = self._Monitor([])
        out = mpa.utils.create_adjoint_sources([monitor], [np.ones(3)])
        self.assertEqual(out, [])

    def test_placements_from_several_monitors_are_concatenated(self):
        a, b = self._Monitor(["a"]), self._Monitor([])
        out = mpa.utils.create_adjoint_sources([a, b], [np.ones(2), np.ones(2)])
        self.assertEqual(out, ["a"])

    def test_globally_zero_cotangent_still_raises(self):
        # The condition that is a real user error is unchanged, and it is global:
        # the cotangents are identical on every rank.
        with self.assertRaisesRegex(RuntimeError, "gradient of all monitor values"):
            mpa.utils.create_adjoint_sources([self._Monitor(["a"])], [np.zeros(3)])


if __name__ == "__main__":
    unittest.main()
