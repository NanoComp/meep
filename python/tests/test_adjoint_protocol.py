"""Tests for how `meep.adjoint` differentiates objective functions.

These cover the contract between `OptimizationProblem` and `ObjectiveQuantity`
rather than the physics: that a single vector-Jacobian product reproduces the
frequency contraction of the Jacobian that Meep used to build explicitly, that
`place_adjoint_source` validates the cotangent it is handed, and that an
objective function written in JAX differentiates identically to the same
function written in autograd.

Everything except `TestJaxObjectiveFunctionGradient` runs without an FDTD simulation.
"""
import os
import re
import unittest
import warnings

import numpy as np
from autograd import jacobian
from autograd import numpy as npa
from autograd.extend import defvjp, primitive

import meep as mp

try:
    import meep.adjoint as mpa
except ImportError:
    import adjoint as mpa

from meep.adjoint import utils
from utils import ApproxComparisonTestCase

try:
    import jax

    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp
except ImportError:
    jax = None


NUM_FREQ = 3


class _StubQuantity(mpa.ObjectiveQuantity):
    """A minimal `ObjectiveQuantity` that records the cotangent it receives."""

    def __init__(self, value):
        self._eval = np.asarray(value)
        self._frequencies = np.arange(NUM_FREQ, dtype=float)
        self.received = None

    def __call__(self):
        return self._eval

    def register_monitors(self, frequencies):
        self._frequencies = np.asarray(frequencies)
        return None

    def place_adjoint_source(self, dJ):
        self.received = self._cotangent(dJ)
        return []


def _random(shape, seed) -> np.ndarray:
    rng = np.random.RandomState(seed)
    return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)


_reverse_passes = []


@primitive
def _count_reverse_passes(x):
    """Identity whose VJP records that a reverse pass traversed it."""
    return x


defvjp(
    _count_reverse_passes,
    lambda ans, x: lambda g: (_reverse_passes.append(None), g)[1],
)


class TestObjectiveVJP(ApproxComparisonTestCase):
    """`utils.objective_vjp` versus the Jacobian it replaced."""

    def test_matches_legacy_jacobian_contraction(self):
        """One VJP reproduces the contracted Jacobian, for every value shape.

        Meep used to call `autograd.jacobian` per objective argument and let each
        `place_adjoint_source` sum away one of the two frequency axes. Those two
        routes agree exactly whenever the objective is block diagonal in
        frequency, which is the assumption the whole single-adjoint-simulation
        formulation rests on. The shapes below are those of
        `EigenmodeCoefficient`/`LDOS`, `FourierFields`, and `Near2FarFields`.
        """
        cases = {
            "per_frequency": (NUM_FREQ,),
            "plane": (NUM_FREQ, 4, 5),
            "near2far": (7, NUM_FREQ, 6),
        }
        for name, shape in cases.items():
            with self.subTest(shape=name):
                value = _random(shape, seed=hash(name) % 2**31)
                freq_axis = 1 if name == "near2far" else 0
                weights = _random(shape, seed=17)

                def objective(x):
                    # Block diagonal in frequency: entry f touches only the
                    # frequency-f slice of x.
                    return npa.sum(
                        npa.abs(x * weights) ** 2,
                        axis=tuple(i for i in range(len(shape)) if i != freq_axis),
                    )

                # Legacy: (nfreq,) + shape Jacobian, contracted along the
                # objective's output axis.
                legacy = np.sum(jacobian(objective)(value), axis=0)

                (cotangent,) = utils.objective_vjp(
                    objective, (value,), np.ones(NUM_FREQ)
                )

                self.assertEqual(cotangent.shape, shape)
                self.assertClose(
                    legacy.ravel(), cotangent.ravel(), epsilon=1e-13, msg=name
                )

    def test_multi_argument_single_pass(self):
        """Cotangents for every argument from one call, matching separate VJPs."""
        a = _random((NUM_FREQ,), seed=1)
        b = _random((NUM_FREQ, 6), seed=2)
        c = _random((NUM_FREQ,), seed=3)

        def objective(a, b, c):
            return npa.abs(a + npa.sum(b, axis=1) + c) ** 2

        together = utils.objective_vjp(objective, (a, b, c), np.ones(NUM_FREQ))
        self.assertEqual(len(together), 3)

        for i, arg in enumerate((a, b, c)):
            separate = np.sum(jacobian(objective, i)(a, b, c), axis=0)
            self.assertEqual(together[i].shape, arg.shape)
            self.assertClose(
                separate.ravel(), together[i].ravel(), epsilon=1e-13, msg=f"arg{i}"
            )

    def test_single_reverse_pass(self):
        """One forward trace and one reverse pass, not one per argument per frequency.

        Guards the performance claim. `autograd.jacobian` traces the forward pass
        once and then replays the tape once per output component, so the route
        this replaced -- a `jacobian` call per objective argument -- cost
        (num arguments) forward traces and (num arguments) x (num frequencies)
        reverse passes.
        """
        args = (_random((NUM_FREQ,), 4), _random((NUM_FREQ, 3), 5))
        forward_traces = []

        def objective(a, b):
            forward_traces.append(None)
            # `_count_reverse_passes` is the identity, but its VJP is recorded on
            # every tape that reaches the output, so it fires exactly once per
            # reverse pass.
            return _count_reverse_passes(npa.abs(a + npa.sum(b, axis=1)) ** 2)

        _reverse_passes.clear()
        utils.objective_vjp(objective, args, np.ones(NUM_FREQ))
        self.assertEqual(len(forward_traces), 1)
        self.assertEqual(len(_reverse_passes), 1)

        # For contrast, the route this replaced.
        forward_traces.clear()
        _reverse_passes.clear()
        for i in range(len(args)):
            jacobian(objective, i)(*args)
        self.assertEqual(len(forward_traces), len(args))
        self.assertEqual(len(_reverse_passes), len(args) * NUM_FREQ)

    def test_scalar_objective(self):
        """A scalar-valued objective needs a scalar seed, not a vector of ones."""
        value = _random((NUM_FREQ, 4), seed=6)

        def objective(x):
            return npa.sum(npa.abs(x) ** 2)

        (cotangent,) = utils.objective_vjp(objective, (value,), 1.0)
        self.assertEqual(cotangent.shape, value.shape)
        self.assertClose(
            jacobian(objective)(value).ravel(), cotangent.ravel(), epsilon=1e-13
        )


class TestCotangentContract(unittest.TestCase):
    """`place_adjoint_source` validates the shape of what it is handed."""

    def test_accepts_value_shape(self):
        quantity = _StubQuantity(_random((NUM_FREQ, 5), seed=7))
        cotangent = _random((NUM_FREQ, 5), seed=8)
        with warnings.catch_warnings():
            warnings.simplefilter("error")  # no deprecation for the new contract
            quantity.place_adjoint_source(cotangent)
        np.testing.assert_array_equal(quantity.received, cotangent)

    def test_accepts_legacy_jacobian_with_warning(self):
        """The pre-VJP Jacobian shape still works, contracted, with a warning."""
        quantity = _StubQuantity(_random((NUM_FREQ, 5), seed=9))
        legacy = _random((NUM_FREQ, NUM_FREQ, 5), seed=10)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            quantity.place_adjoint_source(legacy)
        deprecations = [w for w in caught if w.category is DeprecationWarning]
        self.assertEqual(len(deprecations), 1)
        np.testing.assert_allclose(quantity.received, np.sum(legacy, axis=0))

    def test_scalar_cotangent_for_single_frequency(self):
        """A bare scalar is accepted where the value is a single frequency."""
        quantity = _StubQuantity(np.array([1.0 + 2.0j]))
        quantity.place_adjoint_source(3.0 + 4.0j)
        np.testing.assert_array_equal(quantity.received, np.array([3.0 + 4.0j]))

    def test_rejects_mismatched_shape(self):
        quantity = _StubQuantity(_random((NUM_FREQ, 5), seed=11))
        with self.assertRaisesRegex(ValueError, "expected a cotangent of shape"):
            quantity.place_adjoint_source(_random((NUM_FREQ, 4), seed=12))


class TestGatherMonitorValues(unittest.TestCase):
    """`utils.gather_monitor_values` stacks what it can and tuples the rest."""

    def test_stacks_homogeneous_per_frequency_values(self):
        """Monitors that all yield one value per frequency stack as before.

        The rank-2 (monitor, frequency) array is what `MeepJaxWrapper` has always
        returned, and user code indexes it as `monitor_values[0, :]`.
        """
        monitors = [_StubQuantity(_random((NUM_FREQ,), seed=s)) for s in (13, 14)]
        values = utils.gather_monitor_values(monitors)
        self.assertIsInstance(values, np.ndarray)
        self.assertEqual(values.shape, (2, NUM_FREQ))

    def test_stacks_scalar_values(self):
        monitors = [_StubQuantity(np.array(1.0 + 1j)) for _ in range(3)]
        values = utils.gather_monitor_values(monitors)
        self.assertIsInstance(values, np.ndarray)
        self.assertEqual(values.shape, (3, 1))

    def test_tuples_heterogeneous_values(self):
        monitors = [
            _StubQuantity(_random((NUM_FREQ,), seed=15)),
            _StubQuantity(_random((NUM_FREQ, 4, 5), seed=16)),
        ]
        values = utils.gather_monitor_values(monitors)
        self.assertIsInstance(values, tuple)
        self.assertEqual([v.shape for v in values], [(NUM_FREQ,), (NUM_FREQ, 4, 5)])


class TestJaxIsOptional(unittest.TestCase):
    """JAX must stay an optional dependency of `meep.adjoint`."""

    def test_autograd_is_the_default_backend(self):
        """An autograd objective is unaffected by any registered backend."""
        value = _random((NUM_FREQ, 4), seed=24)

        def objective(x):
            return npa.sum(npa.abs(x) ** 2, axis=1)

        (cotangent,) = utils.objective_vjp(
            objective, (value,), np.ones(NUM_FREQ), value=objective(value)
        )
        np.testing.assert_allclose(cotangent, 2 * value.conj())

    def test_no_jax_import_outside_wrapper(self):
        """Only `wrapper.py` may import JAX, and it is guarded in `__init__.py`.

        Without this, adding `import jax` to any other module of the package
        would break every user who does not have JAX installed, and the failure
        would only show up in an environment nobody runs tests in.
        """
        adjoint_dir = os.path.dirname(os.path.abspath(mpa.__file__))
        pattern = re.compile(r"^\s*(?:import\s+jax|from\s+jax)", re.MULTILINE)
        offenders = []
        for name in sorted(os.listdir(adjoint_dir)):
            if not name.endswith(".py") or name == "wrapper.py":
                continue
            with open(os.path.join(adjoint_dir, name)) as f:
                if pattern.search(f.read()):
                    offenders.append(name)
        self.assertEqual(
            offenders,
            [],
            "JAX may only be imported from meep/adjoint/wrapper.py, which is "
            "imported behind a ModuleNotFoundError guard in __init__.py.",
        )


@unittest.skipIf(jax is None, "jax is not installed")
class TestJaxObjectiveFunction(ApproxComparisonTestCase):
    """An objective written in JAX differentiates like the autograd version.

    A JAX objective is a plain function passed straight to `OptimizationProblem`:
    there is nothing to wrap or declare. It is recognized because it returns a
    JAX array, and `wrapper` registers JAX as a way to differentiate such
    functions when it is imported.
    """

    def test_dispatches_on_return_type(self):
        """A bare jax.numpy function is differentiated by JAX, not autograd."""
        value = _random((NUM_FREQ, 4), seed=30)

        def objective(x):
            # autograd cannot trace this: jnp rejects autograd's boxes.
            return jnp.sum(jnp.abs(x) ** 2, axis=1)

        cotangent = np.ones(NUM_FREQ)
        (got,) = utils.objective_vjp(
            objective, (value,), cotangent, value=objective(value)
        )
        self.assertEqual(got.shape, value.shape)
        self.assertClose(2 * value.conj().ravel(), got.ravel(), epsilon=1e-13)

        # Without the return value there is nothing to dispatch on, so autograd
        # is used -- and fails, rather than silently doing something else.
        with self.assertRaises(Exception):
            utils.objective_vjp(objective, (value,), cotangent)

    def test_agrees_with_autograd(self):
        """The two frameworks' complex cotangent conventions must coincide.

        A real-valued objective of complex monitor values is the case that
        matters, and it is the one where a conjugation-convention mismatch would
        silently produce a gradient of the right magnitude pointing the wrong
        way.
        """
        a = _random((NUM_FREQ,), seed=18)
        b = _random((NUM_FREQ, 6), seed=19)
        weights = np.random.RandomState(20).standard_normal(6)

        def build(xp):
            def objective(a, b):
                overlap = xp.sum(b * weights, axis=1)
                return xp.abs(overlap) ** 2 / xp.abs(a) ** 2 + xp.real(
                    a * xp.conj(overlap)
                )

            return objective

        autograd_objective = build(npa)
        jax_objective = build(jnp)

        self.assertClose(
            autograd_objective(a, b),
            np.asarray(jax_objective(a, b)),
            epsilon=1e-13,
        )

        cotangent = np.ones(NUM_FREQ)
        expected = utils.objective_vjp(
            autograd_objective, (a, b), cotangent, value=autograd_objective(a, b)
        )
        actual = utils.objective_vjp(
            jax_objective, (a, b), cotangent, value=jax_objective(a, b)
        )

        self.assertEqual(len(actual), len(expected))
        for i, (want, got) in enumerate(zip(expected, actual)):
            self.assertEqual(got.shape, want.shape)
            self.assertClose(want.ravel(), got.ravel(), epsilon=1e-12, msg=f"arg{i}")
            # A conjugated convention would be a much larger error than the
            # tolerance above; assert it explicitly so the test names the failure.
            self.assertGreater(
                np.abs(want - np.conj(got)).max(),
                1e-6 * np.abs(want).max(),
                msg=f"arg{i}: cotangent is real, so the conjugation check is vacuous",
            )

    def test_scalar_objective(self):
        value = _random((NUM_FREQ, 4), seed=21)

        def objective(x):
            return jnp.sum(jnp.abs(x) ** 2)

        (cotangent,) = utils.objective_vjp(
            objective, (value,), 1.0, value=objective(value)
        )
        self.assertEqual(cotangent.shape, value.shape)
        self.assertClose(2 * value.conj().ravel(), cotangent.ravel(), epsilon=1e-13)

    def test_rejects_non_array_output(self):
        def objective(x):
            return (jnp.sum(x), jnp.sum(x))

        # A tuple return is not dispatched to JAX (it is not a jax.Array), so
        # force the backend directly to check its own diagnostic.
        from meep.adjoint import wrapper

        with self.assertRaisesRegex(TypeError, "single array or scalar"):
            wrapper._jax_objective_vjp(objective, (_random((NUM_FREQ,), 22),), 1.0)

    def test_explicit_vjp_takes_precedence(self):
        """A hand-written pullback wins over any registered backend.

        This is the escape hatch for an analytic derivative, or for a framework
        with no registered backend.
        """
        value = _random((NUM_FREQ, 4), seed=23)
        sentinel = np.full(value.shape, 7.0 + 0j)

        def objective(x):
            return jnp.sum(jnp.abs(x) ** 2, axis=1)

        objective.vjp = lambda cotangent, *args: (sentinel,)

        (cotangent,) = utils.objective_vjp(
            objective, (value,), np.ones(NUM_FREQ), value=objective(value)
        )
        np.testing.assert_array_equal(cotangent, sentinel)


@unittest.skipIf(jax is None, "jax is not installed")
class TestJaxObjectiveFunctionGradient(ApproxComparisonTestCase):
    """End-to-end: a JAX objective driving `OptimizationProblem`.

    The two objective arguments are `FourierFields` monitors of different rank,
    so this also exercises passing cotangents of heterogeneous shape through
    `utils.create_adjoint_sources`.

    Deliberately no `EigenmodeCoefficient`: MPB warm-starts from the previous
    solve, so repeating an otherwise identical run reproduces a mode coefficient
    only to about 1e-8. That is far below every physics tolerance in the suite,
    but it would put a floor under the framework-equivalence comparison below,
    which should be exact to roundoff.
    """

    @classmethod
    def setUpClass(cls):
        cls.resolution = 20
        cls.silicon = mp.Medium(epsilon=12)
        cls.sxy = 5.0
        cls.dpml = 1.0
        cls.fcen = 1 / 1.55
        cls.frequencies = [cls.fcen * 0.98, cls.fcen, cls.fcen * 1.02]
        cls.design_region_size = mp.Vector3(1.0, 1.0)
        cls.Nx = int(cls.design_region_size.x * cls.resolution) + 1
        cls.Ny = int(cls.design_region_size.y * cls.resolution) + 1

        rng = np.random.RandomState(4321)
        cls.p = 0.5 * rng.rand(cls.Nx * cls.Ny)
        cls.dp = 1e-4 * rng.rand(cls.Nx * cls.Ny)

    def _solve(self, design_params, use_jax: bool, need_gradient: bool = True):
        """Runs one forward (and optionally adjoint) solve of a small problem."""
        matgrid = mp.MaterialGrid(
            mp.Vector3(self.Nx, self.Ny),
            mp.air,
            self.silicon,
            weights=np.ones((self.Nx, self.Ny)),
        )
        design_region = mpa.DesignRegion(
            matgrid,
            volume=mp.Volume(
                center=mp.Vector3(),
                size=mp.Vector3(self.design_region_size.x, self.design_region_size.y),
            ),
        )
        sim = mp.Simulation(
            resolution=self.resolution,
            cell_size=mp.Vector3(self.sxy, self.sxy),
            boundary_layers=[mp.PML(self.dpml)],
            sources=[
                mp.Source(
                    mp.GaussianSource(self.fcen, fwidth=0.2 * self.fcen),
                    component=mp.Ez,
                    center=mp.Vector3(-1.25, 0),
                    size=mp.Vector3(0, 2.0),
                )
            ],
            geometry=[
                mp.Block(
                    center=design_region.center,
                    size=design_region.size,
                    material=matgrid,
                )
            ],
        )

        # Two objective arguments of different rank: a line of Ez, shaped
        # (nfreq, ny), and a patch of Hy, shaped (nfreq, nx, ny).
        obj_args = [
            mpa.FourierFields(
                sim,
                mp.Volume(center=mp.Vector3(1.0, 0), size=mp.Vector3(0, 0.6)),
                mp.Ez,
            ),
            mpa.FourierFields(
                sim,
                mp.Volume(center=mp.Vector3(-1.0, 0.2), size=mp.Vector3(0.3, 0.4)),
                mp.Hy,
            ),
        ]

        def build(xp):
            def objective(ez_line, hy_patch):
                return xp.sum(xp.abs(ez_line) ** 2, axis=1) + xp.sum(
                    xp.abs(hy_patch) ** 2, axis=(1, 2)
                )

            return objective

        # Nothing distinguishes the two at the call site: a JAX objective is a
        # plain function, recognized by the array type it returns.
        objective = build(jnp) if use_jax else build(npa)

        opt = mpa.OptimizationProblem(
            simulation=sim,
            objective_functions=objective,
            objective_arguments=obj_args,
            design_regions=[design_region],
            frequencies=self.frequencies,
        )
        return opt([design_params], need_gradient=need_gradient)

    def test_matches_autograd_gradient(self):
        """Same objective, two frameworks, same value and same gradient."""
        jax_value, jax_grad = self._solve(self.p, use_jax=True)
        autograd_value, autograd_grad = self._solve(self.p, use_jax=False)

        self.assertClose(jax_value, autograd_value, epsilon=1e-12)
        self.assertClose(
            np.asarray(jax_grad).ravel(),
            np.asarray(autograd_grad).ravel(),
            epsilon=1e-10,
        )

    def test_matches_finite_difference(self):
        """The JAX-derived adjoint gradient is the real gradient."""
        value, grad = self._solve(self.p, use_jax=True)
        perturbed, _ = self._solve(self.p + self.dp, use_jax=True, need_gradient=False)

        grad = np.atleast_2d(np.asarray(grad))
        adjoint_dd = (self.dp[None, :] @ grad).ravel()
        finite_difference_dd = np.asarray(perturbed) - np.asarray(value)

        tol = 0.075 if mp.is_single_precision() else 0.005
        self.assertClose(adjoint_dd, finite_difference_dd, epsilon=tol)


if __name__ == "__main__":
    unittest.main()
