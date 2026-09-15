import unittest

import autograd.numpy as npa
import numpy as np

import meep as mp
import meep.adjoint as mpa


class TestPadeAdjoint(unittest.TestCase):
    def _problem(self, weights, pade_samples, run_time, objective_functions=None):
        fcen = 0.5
        grid = mp.MaterialGrid(
            mp.Vector3(1, 1, len(weights)),
            mp.air,
            mp.Medium(epsilon=12),
            weights=np.asarray(weights).reshape(1, 1, -1),
            do_averaging=False,
        )
        region = mpa.DesignRegion(
            grid,
            volume=mp.Volume(center=mp.Vector3(), size=mp.Vector3(z=1)),
        )
        sim = mp.Simulation(
            cell_size=mp.Vector3(z=6),
            dimensions=1,
            resolution=20,
            boundary_layers=[mp.PML(1)],
            geometry=[
                mp.Block(
                    center=region.center,
                    size=region.size,
                    material=grid,
                )
            ],
            sources=[
                mp.Source(
                    mp.GaussianSource(fcen, fwidth=0.3),
                    component=mp.Ex,
                    center=mp.Vector3(z=-1.5),
                )
            ],
        )
        field = mpa.FourierFields(
            sim,
            mp.Volume(center=mp.Vector3(z=1.2), size=mp.Vector3(z=0.1)),
            mp.Ex,
        )
        if objective_functions is None:
            objective_functions = lambda value: npa.abs(value[0, 0]) ** 2
        return mpa.OptimizationProblem(
            simulation=sim,
            objective_functions=objective_functions,
            objective_arguments=[field],
            design_regions=[region],
            frequencies=[fcen],
            decay_by=1e-30,
            minimum_run_time=run_time,
            maximum_run_time=run_time,
            decimation_factor=1,
            pade_samples=pade_samples,
        )

    def _run(
        self,
        weights,
        pade_samples,
        run_time,
        objective_functions=None,
        need_gradient=True,
    ):
        problem = self._problem(
            weights, pade_samples, run_time, objective_functions=objective_functions
        )
        value, gradient = problem([np.asarray(weights)], need_gradient=need_gradient)
        return value, gradient, problem

    def test_gradient_matches_long_run_and_finite_difference(self):
        """Exercise corrected forward/adjoint design fields in the native gradient."""
        weights = np.array([0.35, 0.55, 0.7])
        short_value, short_gradient, short_problem = self._run(
            weights, pade_samples=12, run_time=1
        )
        # The Padé path normalizes its synthesized adjoint source over the
        # source's complete support. Run the raw reference past that support so
        # the comparison measures extrapolation rather than normalization time.
        long_value, long_gradient, long_problem = self._run(
            weights, pade_samples=0, run_time=220
        )

        tolerance = 0.02 if mp.is_single_precision() else 0.005
        np.testing.assert_allclose(short_value, long_value, rtol=tolerance)
        np.testing.assert_allclose(
            short_gradient, long_gradient, rtol=tolerance, atol=1e-8
        )
        self.assertLess(short_problem.sim.meep_time(), long_problem.sim.meep_time())

        direction = np.array([1.0, -0.5, 0.25])
        direction /= np.linalg.norm(direction)
        step = 1e-3
        plus, _, _ = self._run(
            weights + step * direction,
            pade_samples=0,
            run_time=220,
            need_gradient=False,
        )
        minus, _, _ = self._run(
            weights - step * direction,
            pade_samples=0,
            run_time=220,
            need_gradient=False,
        )
        finite_difference = (float(np.squeeze(plus)) - float(np.squeeze(minus))) / (
            2 * step
        )
        adjoint_derivative = float(np.asarray(short_gradient) @ direction)
        self.assertAlmostEqual(
            adjoint_derivative / finite_difference,
            1.0,
            delta=0.03 if mp.is_single_precision() else 0.01,
        )

    def test_streamed_multi_objective_matches_long_run(self):
        """Streamed Padé gradients match the legacy retained-monitor result."""
        weights = np.array([0.35, 0.55, 0.7])

        def left(value):
            return npa.abs(value[0, 0]) ** 2

        def right(value):
            return npa.abs(value[0, -1]) ** 2

        short_value, short_gradient, short_problem = self._run(
            weights,
            pade_samples=12,
            run_time=1,
            objective_functions=[left, right],
        )
        long_value, long_gradient, _ = self._run(
            weights,
            pade_samples=0,
            run_time=220,
            objective_functions=[left, right],
        )

        tolerance = 0.02 if mp.is_single_precision() else 0.005
        np.testing.assert_allclose(short_value, long_value, rtol=tolerance)
        np.testing.assert_allclose(
            short_gradient, long_gradient, rtol=tolerance, atol=1e-8
        )
        self.assertEqual(len(short_problem._streamed_gradients), 2)
        self.assertEqual(short_problem.adjoint_design_region_monitors, [None, None])


if __name__ == "__main__":
    unittest.main()
