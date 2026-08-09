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

Note the usual directional finite-difference check is nearly blind to this: it
probes along g/|g|, i.e. along the erroneous gradient itself, where the error
largely cancels. Measured on the pre-fix build, the gradient vector was 44% off
in L2 while fd/adjoint still read 1.0008. Comparing the gradient vectors
directly is what catches it.
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

        # the forward solve was never chunk-dependent; if this trips, the test
        # problem itself changed rather than the gradient assembly
        self.assertAlmostEqual(f0_a / f0_b, 1.0, places=9)

        rel = np.linalg.norm(g_a - g_b) / np.linalg.norm(g_a)
        # Pre-fix this was 1.2e-1 to 4.4e-1 depending on where the split landed.
        # Post-fix the only difference is MPI summation order.
        self.assertLess(
            rel, 1e-9, f"gradient changed by {rel:.3e} under a different chunk split"
        )


if __name__ == "__main__":
    unittest.main()
