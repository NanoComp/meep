import unittest
import warnings
import tempfile
from unittest import mock

import numpy as np

import meep as mp
import meep.adjoint as mpa
from meep import simulation


class PadeMonitorApiTest(unittest.TestCase):
    def setUp(self):
        self.sim = mp.Simulation(
            cell_size=mp.Vector3(1, 1),
            resolution=10,
            sources=[
                mp.Source(
                    mp.GaussianSource(1.0, fwidth=0.5),
                    component=mp.Ez,
                    center=mp.Vector3(),
                )
            ],
        )
        self.frequency = [1.0]
        self.volume = mp.Volume(center=mp.Vector3(), size=mp.Vector3(0, 1))

    def test_pade_samples_validation(self):
        for invalid in (True, 1.5, -1, 1, 2, 3):
            with self.subTest(invalid=invalid), self.assertRaises(
                (TypeError, ValueError)
            ):
                self.sim.add_dft_fields([mp.Ez], self.frequency, pade_samples=invalid)

    def test_pade_samples_are_deferred_with_each_monitor(self):
        monitors = [
            self.sim.add_dft_fields(
                [mp.Ez], self.frequency, where=self.volume, pade_samples=8
            ),
            self.sim.add_flux(
                self.frequency, mp.FluxRegion(volume=self.volume), pade_samples=8
            ),
            self.sim.add_mode_monitor(
                self.frequency, mp.ModeRegion(volume=self.volume), pade_samples=8
            ),
            self.sim.add_force(
                self.frequency, mp.ForceRegion(volume=self.volume), pade_samples=8
            ),
            self.sim.add_energy(
                self.frequency, mp.EnergyRegion(volume=self.volume), pade_samples=8
            ),
            self.sim.add_near2far(
                self.frequency, mp.Near2FarRegion(volume=self.volume), pade_samples=8
            ),
        ]
        self.assertTrue(all(monitor.args[-1] == 8 for monitor in monitors))
        self.assertTrue(all(monitor.pade_samples == 8 for monitor in monitors))
        self.sim._evaluate_dft_objects()
        self.assertTrue(all(monitor.swigobj is not None for monitor in monitors))

    def test_error_conversion_and_warning_suppression(self):
        monitor = simulation.DftObj(lambda: object(), [])
        status = ([0.1], [0.2], [0.3], 8, True, 1, 7, [1.0])
        with mock.patch.object(
            mp, "_get_dft_pade_error", return_value=status, create=True
        ):
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                first = monitor.get_pade_error()
                second = monitor.get_pade_error()
        self.assertIsInstance(first, simulation.PadeError)
        np.testing.assert_allclose(first.relative_estimate, [0.2])
        self.assertEqual(first.generation, second.generation)
        self.assertEqual(len(caught), 1)

    def test_small_dft_monitor_reports_native_status(self):
        sim = mp.Simulation(
            cell_size=mp.Vector3(z=4),
            dimensions=1,
            resolution=10,
            boundary_layers=[mp.PML(1)],
            sources=[
                mp.Source(
                    mp.GaussianSource(0.5, fwidth=0.4),
                    component=mp.Ex,
                    center=mp.Vector3(),
                )
            ],
        )
        frequencies = [0.4, 0.5, 0.6]
        monitor = sim.add_dft_fields(
            [mp.Ex], frequencies, where=mp.Volume(center=mp.Vector3()), pade_samples=4
        )
        sim.run(until_after_sources=2)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            status = monitor.get_pade_error()
            values = [sim.get_dft_array(monitor, mp.Ex, i) for i in range(3)]
        self.assertEqual(status.samples, 4)
        self.assertTrue(status.ready)
        self.assertEqual(status.relative_estimate.shape, (3,))
        self.assertTrue(np.all(np.isfinite(values)))

    def test_loading_dft_snapshot_clears_pade_history(self):
        monitor = self.sim.add_dft_fields(
            [mp.Ez], self.frequency, where=self.volume, pade_samples=4
        )
        self.sim.run(until_after_sources=2)
        self.assertEqual(monitor.get_pade_error().samples, 4)

        snapshot = self.sim.get_dft_data(monitor.chunks)
        mp._load_dft_data(monitor.chunks, snapshot)
        status = monitor.get_pade_error()
        self.assertEqual(status.samples, 0)
        self.assertFalse(status.ready)

    def test_full_checkpoint_rejects_pade_monitor(self):
        self.sim.add_dft_fields(
            [mp.Ez], self.frequency, where=self.volume, pade_samples=4
        )
        self.sim.run(until=1)
        with tempfile.TemporaryDirectory() as dirname, self.assertRaisesRegex(
            RuntimeError, "does not yet support resumable Padé"
        ):
            self.sim.dump(dirname, dump_structure=False, dump_fields=True)

    def test_short_corrected_dft_matches_long_run(self):
        sim = mp.Simulation(
            cell_size=mp.Vector3(z=6),
            dimensions=1,
            resolution=20,
            boundary_layers=[mp.PML(1)],
            geometry=[
                mp.Block(
                    center=mp.Vector3(),
                    size=mp.Vector3(z=1),
                    material=mp.Medium(epsilon=12),
                )
            ],
            sources=[
                mp.Source(
                    mp.GaussianSource(0.5, fwidth=0.3),
                    component=mp.Ex,
                    center=mp.Vector3(),
                )
            ],
        )
        where = mp.Volume(center=mp.Vector3(z=0.3))
        raw = sim.add_dft_fields([mp.Ex], [0.5], where=where)
        pade = sim.add_dft_fields([mp.Ex], [0.5], where=where, pade_samples=12)

        sim.run(until_after_sources=1)
        short_raw = complex(sim.get_dft_array(raw, mp.Ex, 0))
        short_corrected = complex(sim.get_dft_array(pade, mp.Ex, 0))
        sim.run(until_after_sources=40)
        long_reference = complex(sim.get_dft_array(raw, mp.Ex, 0))

        raw_error = abs(short_raw - long_reference)
        corrected_error = abs(short_corrected - long_reference)
        self.assertLess(corrected_error, 0.05 * raw_error)


class PadeStoppingTest(unittest.TestCase):
    class _Fields:
        t = 100
        dt = 1.0

        def __init__(self, decimation=1):
            self.norm = 0
            self.decimation = decimation

        def max_decimation(self):
            return self.decimation

        def last_source_time(self):
            return 0

        def dft_norm(self):
            self.norm += 1
            return self.norm

    class _Simulation:
        def __init__(self, decimation=1):
            self.fields = PadeStoppingTest._Fields(decimation)

        def round_time(self):
            return 10

        def meep_time(self):
            return 10

        def timestep(self):
            return self.fields.t

        def advance(self, timesteps=1):
            self.fields.t += timesteps

    class _Monitor:
        def __init__(self):
            self.generation = 2
            self.status_calls = 0
            self.result_calls = 0

        def get_pade_error(self):
            self.status_calls += 1
            self.generation += 2
            return simulation.PadeError(
                np.array([0.0]),
                np.array([1e-8]),
                np.array([1e-8]),
                4,
                True,
                0,
                self.generation,
                np.array([1.0]),
            )

        def get_result(self):
            self.result_calls += 1

    def test_requires_two_separated_passes(self):
        sim = self._Simulation()
        condition = simulation.stop_when_dft_pade_converged(
            [self._Monitor()], 1e-4, 4, raw_decay_tol=1e-12
        )
        self.assertFalse(condition(sim))
        sim.advance(2)
        self.assertTrue(condition(sim))
        self.assertEqual(condition.finalize(sim)["exit_reason"], "pade_converged")

    def test_two_level_confirmation(self):
        sim = self._Simulation()
        calls = []

        def confirmation():
            calls.append(1)
            return np.array([2.0, -1.0])

        condition = simulation.stop_when_dft_pade_converged(
            [self._Monitor()],
            1e-4,
            4,
            raw_decay_tol=1e-12,
            confirmation=confirmation,
        )
        self.assertFalse(condition(sim))
        sim.advance(2)
        self.assertTrue(condition(sim))
        self.assertEqual(len(calls), 2)

    def test_status_is_not_polled_between_accepted_checkpoints(self):
        sim = self._Simulation(decimation=2)
        monitor = self._Monitor()
        condition = simulation.stop_when_dft_pade_converged(
            [monitor], 1e-4, 4, raw_decay_tol=1e-12
        )

        self.assertFalse(condition(sim))
        self.assertEqual(monitor.status_calls, 1)
        for _ in range(3):
            sim.advance()
            monitor.get_result()
            self.assertFalse(condition(sim))
            self.assertEqual(monitor.status_calls, 1)
        sim.advance()
        monitor.get_result()
        self.assertTrue(condition(sim))
        self.assertEqual(monitor.status_calls, 2)
        self.assertEqual(monitor.result_calls, 4)

    def test_raw_convergence_wins_maximum_time_tie(self):
        sim = self._Simulation()
        with mock.patch.object(
            simulation, "stop_when_dft_decayed", return_value=lambda _sim: True
        ):
            condition = simulation.stop_when_dft_pade_converged(
                [self._Monitor()], 1e-4, 4, maximum_run_time=10
            )
        self.assertTrue(condition(sim))
        self.assertEqual(condition.finalize(sim)["exit_reason"], "raw_decay_fallback")

    def test_maximum_time_remains_strict_without_convergence(self):
        sim = self._Simulation()
        with mock.patch.object(
            simulation, "stop_when_dft_decayed", return_value=lambda _sim: False
        ):
            condition = simulation.stop_when_dft_pade_converged(
                [self._Monitor()], 1e-4, 4, maximum_run_time=10
            )
        self.assertTrue(condition(sim))
        self.assertEqual(condition.finalize(sim)["exit_reason"], "maximum_time")


class AdjointPadeOptionsTest(unittest.TestCase):
    class _Objective(mpa.ObjectiveQuantity):
        def __call__(self):
            return None

        def register_monitors(self, frequencies, pade_samples=0):
            self._frequencies = np.asarray(frequencies)
            self._pade_samples = pade_samples

        def place_adjoint_source(self, dJ):
            return []

    class _Simulation:
        resolution = 10

        class _Fields:
            dt = 0.1

        fields = _Fields()

        def __init__(self):
            self.time = 1.0
            self.sources = []

        def meep_time(self):
            return self.time

        def _infer_dimensions(self, unused):
            return 2

        def using_real_fields(self):
            return False

        k_point = None

    def test_option_validation(self):
        self.assertEqual(mpa.utils.validate_pade_options(0, None), 0)
        self.assertEqual(mpa.utils.validate_pade_options(np.int64(8), 1e-4), 8)
        with self.assertRaises(ValueError):
            mpa.utils.validate_pade_options(0, 1e-4)
        with self.assertRaises(TypeError):
            mpa.utils.validate_pade_options(True, None)

    def test_continuous_source_rejected_for_automatic_stopping(self):
        source = mp.Source(
            mp.ContinuousSource(1.0), component=mp.Ez, center=mp.Vector3()
        )
        with self.assertRaisesRegex(ValueError, "finite-duration"):
            mpa.utils.validate_finite_sources([source])

    def test_nonfinite_gaussian_end_time_is_rejected(self):
        source = mp.Source(
            mp.GaussianSource(1.0, width=np.inf),
            component=mp.Ez,
            center=mp.Vector3(),
        )
        with self.assertRaisesRegex(ValueError, "finite-duration"):
            mpa.utils.validate_finite_sources([source])

    def test_ldos_rejects_pade_samples_without_automatic_stopping(self):
        sim = self._Simulation()
        objective = mpa.LDOS(sim)
        with self.assertRaisesRegex(ValueError, "does not support LDOS"):
            objective.register_monitors([1.0], pade_samples=8)
        with self.assertRaisesRegex(ValueError, "does not support LDOS"):
            mpa.OptimizationProblem(
                simulation=sim,
                objective_functions=lambda value: value,
                objective_arguments=[objective],
                design_regions=[],
                frequencies=np.array([1.0]),
                pade_samples=8,
            )

    def test_adjoint_source_scale_is_independent_of_shortened_run(self):
        sim = self._Simulation()
        objective = self._Objective(sim)
        objective.register_monitors([1.0], pade_samples=8)
        short_scale = objective._adj_src_scale()
        sim.time = 100.0
        long_scale = objective._adj_src_scale()
        np.testing.assert_allclose(short_scale, long_scale, rtol=0, atol=0)

    def test_failed_adjoint_run_restores_transforms_and_removes_monitors(self):
        class FakeSimulation:
            is_cylindrical = True
            dimensions = mp.CYLINDRICAL

            def __init__(self):
                self.m = 3
                self.k_point = mp.Vector3(0.25, 0, 0)

            def change_m(self, value):
                self.m = value

            def change_k_point(self, value):
                self.k_point = value

            def restart_fields(self):
                pass

            def clear_dft_monitors(self):
                pass

            def change_sources(self, sources):
                pass

            def _evaluate_dft_objects(self):
                pass

        monitor = mock.Mock()
        problem = mpa.OptimizationProblem.__new__(mpa.OptimizationProblem)
        problem.sim = FakeSimulation()
        problem.objective_functions = [lambda: None]
        problem.design_regions = []
        problem.frequencies = [1.0]
        problem.decimation_factor = 1
        problem.pade_samples = 8
        problem.forward_design_region_monitors = []
        problem.finite_difference_step = 1e-3
        # adjoint_run also installs monitors for differentiable geometry and
        # sources; this test builds the object attribute by attribute, so it
        # has to name them even when there are none.
        problem.differentiable_objects = []
        problem.differentiable_sources = []
        problem.prepare_adjoint_run = mock.Mock(
            side_effect=lambda: setattr(problem, "adjoint_sources", [[]])
        )
        problem._run_with_stopping = mock.Mock(side_effect=RuntimeError("timeout"))
        original_k = problem.sim.k_point

        with mock.patch.object(
            mpa.utils,
            "install_design_region_monitors",
            return_value=[[monitor]],
        ), self.assertRaisesRegex(RuntimeError, "timeout"):
            problem.adjoint_run()

        self.assertEqual(problem.sim.m, 3)
        self.assertEqual(problem.sim.k_point, original_k)
        monitor.remove.assert_called_once_with()


if __name__ == "__main__":
    unittest.main()
