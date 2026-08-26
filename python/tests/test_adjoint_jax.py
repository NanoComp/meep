import unittest
from typing import NamedTuple

import jax
import jax.numpy as jnp
import meep.adjoint as mpa
import numpy as onp
import parameterized
from utils import ApproxComparisonTestCase

import meep as mp

# The calculation of finite-difference gradients
# requires that JAX be operated with double precision
jax.config.update("jax_enable_x64", True)

# The step size for the finite-difference
# gradient calculation
_FD_STEP = 1e-4

# The tolerance for the adjoint and finite-difference
# gradient comparison
_TOL = 0.1 if mp.is_single_precision() else 0.025

# We expect 3 design region monitor pointers
# (one for each field component)
_NUM_DES_REG_MON = 3

mp.verbosity(0)


class _PostProcessing(NamedTuple):
    """Stands in for parameters that never reach Meep, e.g. a layer stack."""

    weight: float
    bias: onp.ndarray


def build_straight_wg_simulation(
    wg_width=0.5,
    wg_padding=1.0,
    wg_length=1.0,
    pml_width=1.0,
    source_to_pml=0.5,
    source_to_monitor=0.1,
    frequencies=[1 / 1.55],
    gaussian_rel_width=0.2,
    sim_resolution=20,
    design_region_resolution=20,
):
    """Builds a simulation of a straight waveguide with a design region segment."""
    design_region_shape = (1.0, wg_width)

    # Simulation domain size
    sx = 2 * pml_width + 2 * wg_length + design_region_shape[0]
    sy = (
        2 * pml_width
        + 2 * wg_padding
        + max(
            wg_width,
            design_region_shape[1],
        )
    )

    # Mean / center frequency
    fmean = onp.mean(frequencies)

    si = mp.Medium(index=3.4)
    sio2 = mp.Medium(index=1.44)

    sources = [
        mp.EigenModeSource(
            mp.GaussianSource(frequency=fmean, fwidth=fmean * gaussian_rel_width),
            eig_band=1,
            direction=mp.NO_DIRECTION,
            eig_kpoint=mp.Vector3(1, 0, 0),
            size=mp.Vector3(0, wg_width + 2 * wg_padding, 0),
            center=[-sx / 2 + pml_width + source_to_pml, 0, 0],
        ),
        mp.EigenModeSource(
            mp.GaussianSource(frequency=fmean, fwidth=fmean * gaussian_rel_width),
            eig_band=1,
            direction=mp.NO_DIRECTION,
            eig_kpoint=mp.Vector3(-1, 0, 0),
            size=mp.Vector3(0, wg_width + 2 * wg_padding, 0),
            center=[sx / 2 - pml_width - source_to_pml, 0, 0],
        ),
    ]
    nx, ny = int(design_region_shape[0] * design_region_resolution), int(
        design_region_shape[1] * design_region_resolution
    )
    mat_grid = mp.MaterialGrid(
        mp.Vector3(nx, ny),
        sio2,
        si,
        grid_type="U_DEFAULT",
    )

    design_regions = [
        mpa.DesignRegion(
            mat_grid,
            volume=mp.Volume(
                center=mp.Vector3(),
                size=mp.Vector3(
                    design_region_shape[0],
                    design_region_shape[1],
                    0,
                ),
            ),
        )
    ]

    geometry = [
        mp.Block(
            center=mp.Vector3(
                x=-design_region_shape[0] / 2 - wg_length / 2 - pml_width / 2
            ),
            material=si,
            size=mp.Vector3(wg_length + pml_width, wg_width, 0),
        ),  # left wg
        mp.Block(
            center=mp.Vector3(
                x=+design_region_shape[0] / 2 + wg_length / 2 + pml_width / 2
            ),
            material=si,
            size=mp.Vector3(wg_length + pml_width, wg_width, 0),
        ),  # right wg
        mp.Block(
            center=design_regions[0].center,
            size=design_regions[0].size,
            material=mat_grid,
        ),  # design region
    ]

    simulation = mp.Simulation(
        cell_size=mp.Vector3(sx, sy),
        boundary_layers=[mp.PML(pml_width)],
        geometry=geometry,
        sources=sources,
        resolution=sim_resolution,
    )

    monitor_centers = [
        mp.Vector3(-sx / 2 + pml_width + source_to_pml + source_to_monitor),
        mp.Vector3(sx / 2 - pml_width - source_to_pml - source_to_monitor),
    ]
    monitor_size = mp.Vector3(y=wg_width + 2 * wg_padding)

    monitors = [
        mpa.EigenmodeCoefficient(
            simulation,
            mp.Volume(center=center, size=monitor_size),
            mode=1,
            forward=forward,
        )
        for center in monitor_centers
        for forward in [True, False]
    ]
    return simulation, sources, monitors, design_regions, frequencies


class UtilsTest(unittest.TestCase):
    def setUp(self):
        super().setUp()
        (
            self.simulation,
            self.sources,
            self.monitors,
            self.design_regions,
            self.frequencies,
        ) = build_straight_wg_simulation()

    def test_mode_monitor_helpers(self):
        mpa.utils.register_monitors(self.monitors, self.frequencies)
        self.simulation.run(until=100)
        monitor_values = mpa.utils.gather_monitor_values(self.monitors)
        self.assertEqual(monitor_values.dtype, onp.complex128)
        self.assertEqual(
            monitor_values.shape, (len(self.monitors), len(self.frequencies))
        )

    def test_dist_dft_pointers(self):
        fwd_design_region_monitors = mpa.utils.install_design_region_monitors(
            self.simulation,
            self.design_regions,
            self.frequencies,
        )
        self.assertEqual(len(fwd_design_region_monitors[0]), _NUM_DES_REG_MON)


class WrapperTest(ApproxComparisonTestCase):
    @parameterized.parameterized.expand(
        [
            (
                "1550_singlefreq_gaussian_port1",
                [1 / 1.55],
                0.1,
                0.5,
                0,
            ),
            (
                "1550_singlefreq_gaussian_port2",
                [1 / 1.55],
                0.1,
                0.5,
                1,
            ),
            (
                "1500_1550bw_01relative_gaussian_port1",
                onp.linspace(1 / 1.50, 1 / 1.55, 3).tolist(),
                0.1,
                0.5,
                0,
            ),
            (
                "1550_1600bw_02relative_gaussian_port1",
                onp.linspace(1 / 1.55, 1 / 1.60, 3).tolist(),
                0.2,
                0.5,
                0,
            ),
            (
                "1500_1600bw_03relative_gaussian_port1",
                onp.linspace(1 / 1.50, 1 / 1.60, 4).tolist(),
                0.3,
                0.5,
                0,
            ),
            (
                "1500_1550bw_01relative_gaussian_port2",
                onp.linspace(1 / 1.50, 1 / 1.55, 3).tolist(),
                0.1,
                0.5,
                1,
            ),
            (
                "1550_1600bw_02relative_gaussian_port2",
                onp.linspace(1 / 1.55, 1 / 1.60, 3).tolist(),
                0.2,
                0.5,
                1,
            ),
            (
                "1500_1600bw_03relative_gaussian_port2",
                onp.linspace(1 / 1.50, 1 / 1.60, 4).tolist(),
                0.3,
                0.5,
                1,
            ),
        ]
    )
    def test_wrapper_gradients(
        self,
        _,
        frequencies,
        gaussian_rel_width,
        design_variable_fill_value,
        excite_port_idx,
    ):
        """Tests gradient from the JAX-Meep wrapper against finite differences."""
        (
            simulation,
            sources,
            monitors,
            design_regions,
            frequencies,
        ) = build_straight_wg_simulation(
            frequencies=frequencies, gaussian_rel_width=gaussian_rel_width
        )

        design_shape = tuple(
            int(i) for i in design_regions[0].design_parameters.grid_size
        )[:2]
        x = onp.ones(design_shape) * design_variable_fill_value

        # Define a loss function
        def loss_fn(x, excite_port_idx=0):
            wrapped_meep = mpa.MeepJaxWrapper(
                simulation,
                [sources[excite_port_idx]],
                monitors,
                design_regions,
                frequencies,
            )
            monitor_values = wrapped_meep([x])
            s1p, s1m, s2p, s2m = monitor_values
            t = s2p / s1p if excite_port_idx == 0 else s1m / s2m
            return jnp.mean(jnp.square(jnp.abs(t)))

        value, adjoint_grad = jax.value_and_grad(loss_fn)(
            x, excite_port_idx=excite_port_idx
        )

        projection = []
        fd_projection = []

        # Project along 5 random directions in the design parameter space.
        for seed in range(5):
            # Create dp
            random_perturbation_vector = _FD_STEP * jax.random.normal(
                jax.random.PRNGKey(seed),
                x.shape,
            )

            # Calculate p + dp
            x_perturbed = x + random_perturbation_vector

            # Calculate T(p + dp)
            value_perturbed = loss_fn(x_perturbed, excite_port_idx=excite_port_idx)

            projection.append(
                onp.dot(
                    random_perturbation_vector.ravel(),
                    adjoint_grad.ravel(),
                )
            )
            fd_projection.append(value_perturbed - value)

        projection = onp.stack(projection)
        fd_projection = onp.stack(fd_projection)

        # Check that dp . ∇T ~ T(p + dp) - T(p)
        self.assertClose(
            projection,
            fd_projection,
            epsilon=_TOL,
        )

    def _minimax_setup(self, num_frequencies=3):
        """A loss returning one value per frequency, with non-Meep parameters."""
        frequencies = onp.linspace(1 / 1.50, 1 / 1.60, num_frequencies).tolist()
        (
            simulation,
            sources,
            monitors,
            design_regions,
            frequencies,
        ) = build_straight_wg_simulation(frequencies=frequencies)
        wrapped_meep = mpa.MeepJaxWrapper(
            simulation, [sources[0]], monitors, design_regions, frequencies
        )
        design_shape = tuple(
            int(i) for i in design_regions[0].design_parameters.grid_size
        )[:2]

        def loss(rho, weight, bias):
            s1p, _, s2p, _ = wrapped_meep([rho])
            # `weight` and `bias` never reach Meep, so they exercise the path
            # that is differentiated directly rather than adjointly.
            return weight * jnp.abs(s2p / s1p) ** 2 + bias

        return wrapped_meep, loss, onp.ones(design_shape) * 0.5, frequencies

    def _count_simulations(self, wrapped_meep):
        """Patches the wrapper to count forward and adjoint timestepping runs."""
        counts = {"forward": 0, "adjoint": 0}
        run_forward = wrapped_meep._run_fwd_simulation
        run_adjoint = wrapped_meep._run_adjoint_simulation

        def counting_forward(*args, **kwargs):
            counts["forward"] += 1
            return run_forward(*args, **kwargs)

        def counting_adjoint(*args, **kwargs):
            counts["adjoint"] += 1
            return run_adjoint(*args, **kwargs)

        wrapped_meep._run_fwd_simulation = counting_forward
        wrapped_meep._run_adjoint_simulation = counting_adjoint
        return counts

    def test_value_and_jacobian(self):
        """A gradient per frequency, for Meep and non-Meep parameters alike.

        The reference for row `f` is `jax.grad` of component `f`, which costs one
        adjoint run per row. `value_and_jacobian` produces every row from one.
        """
        wrapped_meep, loss, x, frequencies = self._minimax_setup()
        num_frequencies = len(frequencies)
        weight, bias = 2.5, 0.75

        counts = self._count_simulations(wrapped_meep)
        values, (d_rho, d_weight, d_bias) = mpa.value_and_jacobian(
            loss, argnums=(0, 1, 2)
        )(x, weight, bias)

        self.assertEqual(counts, {"forward": 1, "adjoint": 1})
        self.assertEqual(values.shape, (num_frequencies,))
        self.assertEqual(d_rho.shape, (num_frequencies,) + x.shape)
        self.assertEqual(d_weight.shape, (num_frequencies,))
        self.assertEqual(d_bias.shape, (num_frequencies,))

        for f in range(num_frequencies):
            reference = jax.grad(lambda *a: loss(*a)[f], argnums=(0, 1, 2))(
                x, weight, bias
            )
            # Parameters that bypass Meep are differentiated directly, so they
            # agree exactly.
            self.assertClose(
                onp.asarray([d_weight[f], d_bias[f]]),
                onp.asarray([reference[1], reference[2]]),
                epsilon=1e-13,
                msg=f"row {f} non-Meep parameters",
            )
            # The design row agrees to the accuracy of the reference route, not
            # of this one: seeding a single frequency asks `FilteredSource` to
            # synthesize exact zeros at the others, which it can only
            # approximate. Cross-checked against `OptimizationProblem`, which
            # drives the adjoint the same way this does, the rows agree to
            # roundoff.
            tol = 1e-2 if mp.is_single_precision() else 5e-3
            self.assertClose(
                onp.asarray(d_rho[f]).ravel(),
                onp.asarray(reference[0]).ravel(),
                epsilon=tol,
                msg=f"row {f} design",
            )

        rows = onp.asarray(d_rho).reshape(num_frequencies, -1)
        spread = onp.abs(rows[0] - rows[-1]).max() / onp.abs(rows).max()
        self.assertGreater(spread, 1e-3, "rows are identical; test is vacuous")

    def test_value_and_jacobian_pytree(self):
        """Parameters may be an arbitrary pytree; the Jacobian mirrors it."""
        wrapped_meep, scalar_loss, x, frequencies = self._minimax_setup()
        num_frequencies = len(frequencies)

        # A NamedTuple nested inside a dict, since both are common and users
        # reach for NamedTuples and dataclasses for anything realistic.
        params = {
            "design": {"rho": x},
            "post": _PostProcessing(weight=2.5, bias=onp.array([0.75, 0.25])),
        }

        def loss(p):
            post = p["post"]
            s1p, _, s2p, _ = wrapped_meep([p["design"]["rho"]])
            return post.weight * jnp.abs(s2p / s1p) ** 2 + jnp.sum(post.bias)

        counts = self._count_simulations(wrapped_meep)
        values, jacobian = mpa.value_and_jacobian(loss)(params)

        self.assertEqual(counts, {"forward": 1, "adjoint": 1})
        self.assertEqual(
            jax.tree_util.tree_structure(jacobian),
            jax.tree_util.tree_structure(params),
        )
        shapes = jax.tree_util.tree_map(lambda leaf: leaf.shape, jacobian)
        self.assertIsInstance(shapes["post"], _PostProcessing)
        self.assertEqual(shapes["design"]["rho"], (num_frequencies,) + x.shape)
        self.assertEqual(shapes["post"].weight, (num_frequencies,))
        self.assertEqual(shapes["post"].bias, (num_frequencies, 2))

    def test_value_and_jacobian_parameter_on_both_paths(self):
        """A parameter that both shapes the design and is used downstream.

        `rho` reaches the loss twice here: through the simulation, which needs
        the adjoint run, and directly through a regularizer, which does not. The
        two contributions have to be summed, and nothing else in this file
        exercises that.
        """
        wrapped_meep, _, x, frequencies = self._minimax_setup()
        num_frequencies = len(frequencies)
        regularization = 5.0

        def loss(rho):
            s1p, _, s2p, _ = wrapped_meep([rho])
            return jnp.abs(s2p / s1p) ** 2 + regularization * jnp.sum(rho**2)

        counts = self._count_simulations(wrapped_meep)
        values, jacobian = mpa.value_and_jacobian(loss)(x)
        self.assertEqual(counts, {"forward": 1, "adjoint": 1})

        # Guard against the test being unable to see the direct term at all.
        direct = 2 * regularization * x
        self.assertGreater(
            onp.abs(direct).max() / onp.abs(onp.asarray(jacobian)).max(),
            0.1,
            "the direct contribution is negligible here, so dropping it "
            "entirely would still pass the comparison below",
        )

        tol = 1e-2 if mp.is_single_precision() else 5e-3
        for f in range(num_frequencies):
            reference = jax.grad(lambda r: loss(r)[f])(x)
            self.assertClose(
                onp.asarray(jacobian[f]).ravel(),
                onp.asarray(reference).ravel(),
                epsilon=tol,
                msg=f"row {f}",
            )

    def test_value_and_jacobian_errors(self):
        wrapped_meep, loss, x, frequencies = self._minimax_setup()

        with self.assertRaisesRegex(ValueError, "one per frequency"):
            mpa.value_and_jacobian(lambda r: loss(r, 1.0, 0.0)[:1])(x)

        with self.assertRaisesRegex(TypeError, "1-D JAX array"):
            mpa.value_and_jacobian(lambda r: jnp.sum(loss(r, 1.0, 0.0)))(x)

        with self.assertRaisesRegex(ValueError, "never called a MeepJaxWrapper"):
            mpa.value_and_jacobian(lambda r: jnp.zeros(len(frequencies)) + r.sum())(x)

        with self.assertRaisesRegex(ValueError, "exactly one call"):
            mpa.value_and_jacobian(lambda r: loss(r, 1.0, 0.0) + loss(r, 1.0, 0.0))(x)

    def test_monitor_values_are_stacked_when_homogeneous(self):
        """Mode monitors alone still come back as one (monitor, frequency) array.

        User code indexes the return value as `monitor_values[i, :]`, so monitors
        that all yield a single value per frequency must keep stacking even now
        that heterogeneous monitors are supported.
        """
        frequencies = onp.linspace(1 / 1.50, 1 / 1.60, 3).tolist()
        (
            simulation,
            sources,
            monitors,
            design_regions,
            frequencies,
        ) = build_straight_wg_simulation(frequencies=frequencies)

        wrapped_meep = mpa.MeepJaxWrapper(
            simulation, [sources[0]], monitors, design_regions, frequencies
        )
        design_shape = tuple(
            int(i) for i in design_regions[0].design_parameters.grid_size
        )[:2]
        monitor_values = wrapped_meep([onp.ones(design_shape) * 0.5])

        self.assertIsInstance(monitor_values, jnp.ndarray)
        self.assertEqual(monitor_values.shape, (len(monitors), len(frequencies)))
        self.assertEqual(monitor_values[0, :].shape, (len(frequencies),))

    def test_heterogeneous_monitor_gradients(self):
        """A FourierFields monitor alongside mode monitors.

        `FourierFields` contributes a whole plane of values rather than one per
        frequency, so the wrapper returns a tuple and `custom_vjp` carries a
        matching tuple of cotangents back to `place_adjoint_source`.
        """
        frequencies = onp.linspace(1 / 1.50, 1 / 1.60, 3).tolist()
        (
            simulation,
            sources,
            monitors,
            design_regions,
            frequencies,
        ) = build_straight_wg_simulation(frequencies=frequencies)

        # Inside the non-PML region, past the edge of the design region.
        monitors = monitors + [
            mpa.FourierFields(
                simulation,
                mp.Volume(center=mp.Vector3(0.75, 0), size=mp.Vector3(0, 0.5)),
                mp.Ez,
            )
        ]

        design_shape = tuple(
            int(i) for i in design_regions[0].design_parameters.grid_size
        )[:2]
        x = onp.ones(design_shape) * 0.5

        def loss_fn(x):
            wrapped_meep = mpa.MeepJaxWrapper(
                simulation, [sources[0]], monitors, design_regions, frequencies
            )
            monitor_values = wrapped_meep([x])
            # A tuple, not an array: the last entry has a spatial dimension.
            *mode_coeffs, dft_fields = monitor_values
            transmission = jnp.abs(mode_coeffs[2] / mode_coeffs[0]) ** 2
            intensity = jnp.sum(jnp.abs(dft_fields) ** 2, axis=1)
            return jnp.mean(transmission + intensity)

        value, adjoint_grad = jax.value_and_grad(loss_fn)(x)

        projection = []
        fd_projection = []
        for seed in range(3):
            perturbation = _FD_STEP * jax.random.normal(
                jax.random.PRNGKey(seed), x.shape
            )
            projection.append(onp.dot(perturbation.ravel(), adjoint_grad.ravel()))
            fd_projection.append(loss_fn(x + perturbation) - value)

        self.assertClose(onp.stack(projection), onp.stack(fd_projection), epsilon=_TOL)


if __name__ == "__main__":
    unittest.main()
