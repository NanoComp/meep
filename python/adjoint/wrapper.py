"""JAX interoperability for `meep.adjoint`.

This module holds everything in `meep.adjoint` that needs JAX, so that JAX stays
an optional dependency: it is imported behind a guard in `meep/adjoint/__init__.py`
and nothing else in the package imports it.

It does two things.

Importing it registers JAX with `utils.register_vjp_backend`, so an objective
function written with `jax.numpy` and handed to `OptimizationProblem` is
differentiated by JAX rather than autograd. Nothing needs to be declared or
wrapped -- writing a JAX objective differs from writing an autograd objective
only in which numpy it imports:

```
def objective(mode_coeff, dft_fields):
    return jnp.abs(mode_coeff)**2 + jnp.sum(jnp.abs(dft_fields)**2, axis=1)

opt = mpa.OptimizationProblem(
    simulation=simulation,
    objective_functions=[objective],
    objective_arguments=[mpa.EigenmodeCoefficient(...), mpa.FourierFields(...)],
    design_regions=[mpa.DesignRegion(...)],
    frequencies=frequencies,
)
```

`MeepJaxWrapper` is the other direction: it turns a whole Meep simulation into a
JAX-differentiable callable, so the design variables, the objective, and anything
else in the loss are all traced by JAX. Usage example:
```
import jax.numpy as jnp
import meep as mp
import meep.adjoint as mpa

sources = [
  mp.EigenModeSource(...)
]

monitors = [
  mpa.EigenmodeCoefficient(...),
  mpa.EigenmodeCoefficient(...),
]

design_regions = [
  mpa.DesignRegion(...)
]

frequencies = [1/1.55, 1/1.60, 1/1.65, ...]

simulation = mp.Simulation(...)

wrapped_meep = MeepJaxWrapper(
    simulation,
    sources,
    monitors,
    design_regions,
    frequencies,
    measurement_interval = 50.0,
    dft_field_components = (mp.Ez,),
    dft_threshold = 1e-6,
    minimum_run_time = 0,
    maximum_run_time = onp.inf,
    until_after_sources = True
)

def loss(x):
    monitor_values = wrapped_meep([x])
    t = monitor_values[0,:] / monitor_values[1,:]
    # Mean transmission vs wavelength
    return jnp.mean(jnp.abs(t))

value, grad = jax.value_and_grad(loss)(x)
```

For worst-case (minimax) optimization the loss returns one value per frequency
rather than a scalar, and `value_and_jacobian` returns a gradient for each of
them -- still from a single adjoint simulation. Parameters that never reach Meep
are differentiated alongside the design weights:

```
def loss(rho, thickness, tilt):
    (dft,) = wrapped_meep([rho])
    return postprocess(dft, thickness, tilt)   # (num frequencies,)

values, grads = mpa.value_and_jacobian(loss, argnums=(0, 1, 2))(rho, 1.8, 8.0)
```
"""
import contextlib
import warnings
from typing import Callable, Iterable, List, Optional, Sequence, Tuple, Union

import jax
import jax.numpy as jnp
import numpy as onp

import meep as mp
from meep.simulation import stop_when_dft_pade_converged

from . import DesignRegion, ObjectiveQuantity, utils
from . import source_gradient

_warned_about_precision = False


def _warn_if_precision_mismatch() -> None:
    """Warns when JAX would silently narrow Meep's double-precision fields.

    JAX defaults to 32-bit dtypes, so a double-precision Meep build feeding
    complex128 monitor values into JAX loses half its digits without any error.
    The usual symptom is a finite-difference gradient check that agrees to only
    three or four digits, which reads as a physics problem rather than a dtype
    problem.
    """
    global _warned_about_precision
    if _warned_about_precision:
        return
    if not jax.config.jax_enable_x64 and not mp.is_single_precision():
        _warned_about_precision = True
        warnings.warn(
            "Meep was built in double precision but JAX's x64 mode is off, so "
            "JAX will narrow complex128 monitor values to complex64 and the "
            "gradients will carry roughly single-precision accuracy. Add\n\n"
            "    jax.config.update('jax_enable_x64', True)\n\n"
            "before constructing the objective to keep double precision.",
            RuntimeWarning,
            stacklevel=3,
        )


def _returns_jax_array(value) -> bool:
    """Whether an objective function's return value came from JAX."""
    return isinstance(value, jax.Array)


def _jax_objective_vjp(
    objective_function: Callable, args: Tuple, cotangent
) -> Tuple[onp.ndarray, ...]:
    """Differentiates an objective function written with `jax.numpy`.

    Registered with `utils.register_vjp_backend`, so an objective function that
    returns a JAX array is differentiated by JAX automatically -- writing one is
    no different from writing an autograd objective, apart from which numpy it
    imports.
    """
    _warn_if_precision_mismatch()
    value, pullback = jax.vjp(objective_function, *[jnp.asarray(arg) for arg in args])
    if not isinstance(value, jax.Array):
        raise TypeError(
            "An objective function must return a single array or scalar, but "
            f"{objective_function} returned {type(value)}."
        )
    # `_objective_cotangent` builds the seed in NumPy, which is float64 even when
    # JAX is running in its default 32-bit mode.
    cotangent = jnp.asarray(cotangent).astype(value.dtype).reshape(value.shape)
    return tuple(onp.asarray(grad) for grad in pullback(cotangent))


utils.register_vjp_backend(_returns_jax_array, _jax_objective_vjp)


class MeepJaxWrapper:
    """Wraps a Meep simulation object into a JAX-differentiable callable.

    Attributes:
        simulation: the pre-configured Meep simulation object.
        sources: a list of Meep sources for the forward simulation.
        monitors: a list of objective quantities from the `meep.adjoint` module,
          such as `EigenmodeCoefficient` or `FourierFields`.
        design_regions: a list of design regions from the `meep.adjoint` module.
        frequencies: the list of frequencies, in normalized Meep units.
        measurement_interval: the time interval between DFT field convergence
          measurements, in Meep time units. The default value is 50.
        dft_field_components: a list of Meep field components, such as `mp.Ex`,
          `mp.Hy`, etc, whose DFT will be monitored for convergence to stop the
          simulation. The default is `mp.Ez`.
        dft_threshold: the threshold for DFT field convergence. Once the norm of the
          change in the fields (the maximum over all design regions and field
          components) is less than this value, the simulation will be stopped. The
          default value is 1e-6.
        minimum_run_time: the minimum run time of the simulation, in Meep time
          units. The default value is 0.
        maximum_run_time: the maximum run time of the simulation, in Meep time
          units. The default value is infinity.
        until_after_sources: whether `maximum_run_time` should be ignored until the
          sources have turned off. This parameter specifies whether `until` or
          `until_after_sources` is used. See
          https://meep.readthedocs.io/en/latest/Python_User_Interface/#Simulation
          for more information. The default is true.
        pade_samples: number of recent DFT samples retained for bounded-memory Padé
          extrapolation. Zero disables extrapolation.
        pade_tolerance: optional experimental monitor-level convergence tolerance.
          Requires `pade_samples >= 4` and finite-duration sources.
    """

    def __init__(
        self,
        simulation: mp.Simulation,
        sources: List[mp.Source],
        monitors: List[ObjectiveQuantity],
        design_regions: List[DesignRegion],
        frequencies: List[float],
        dft_threshold: float = 1e-11,
        minimum_run_time: float = 0,
        maximum_run_time: float = onp.inf,
        until_after_sources: bool = True,
        finite_difference_step: float = utils.FD_DEFAULT,
        pade_samples: int = 0,
        pade_tolerance: float = None,
    ):
        _warn_if_precision_mismatch()
        self.simulation = simulation
        self.sources = sources
        self.monitors = monitors
        self.design_regions = design_regions
        self.frequencies = frequencies
        self.dft_threshold = dft_threshold
        self.minimum_run_time = minimum_run_time
        self.maximum_run_time = maximum_run_time
        self.until_after_sources = until_after_sources
        self.finite_difference_step = finite_difference_step
        self.pade_samples = utils.validate_pade_options(pade_samples, pade_tolerance)
        self.pade_tolerance = pade_tolerance
        self.pade_diagnostics = []

        # A source whose amplitudes come from JAX is differentiated
        # automatically: the parameters that produced them live upstream in
        # JAX, so there is nothing for Meep to name. `amp_data` is the array
        # Meep interpolates onto the grid, so that is the cut point.
        self._source_shapes = []
        self.differentiable_sources = [
            s for s in sources if getattr(s, "amp_data", None) is not None
        ]
        for s in self.differentiable_sources:
            if "amp_data" not in getattr(s, "differentiable", ()):
                s.differentiable = tuple(getattr(s, "differentiable", ())) + (
                    "amp_data",
                )

        self._simulate_fn = self._initialize_callable()

    def _run_simulation(self, monitors, leg):
        run_keyword = "until_after_sources" if self.until_after_sources else "until"
        if self.pade_tolerance is None:
            condition = mp.stop_when_dft_decayed(
                self.dft_threshold, self.minimum_run_time, self.maximum_run_time
            )
        else:
            utils.validate_finite_sources(self.simulation.sources)
            condition = stop_when_dft_pade_converged(
                monitors,
                self.pade_tolerance,
                self.pade_samples,
                raw_decay_tol=self.dft_threshold,
                minimum_run_time=self.minimum_run_time,
                maximum_run_time=self.maximum_run_time,
                monitor_only=True,
            )
        self.simulation.run(**{run_keyword: condition})
        if self.pade_tolerance is not None:
            diagnostics = condition.finalize(self.simulation)
            diagnostics["leg"] = leg
            self.pade_diagnostics.append(diagnostics)
            if diagnostics["exit_reason"] == "maximum_time":
                raise RuntimeError(
                    "automatic Padé convergence was not reached before maximum_run_time"
                )

    def __call__(
        self, designs: List[jnp.ndarray], sources: Optional[List[jnp.ndarray]] = None
    ) -> Union[jnp.ndarray, Tuple[jnp.ndarray, ...]]:
        """Performs a Meep simulation, taking designs and returning monitor values.

        Args:
          sources: an `amp_data` array for each source given to the constructor
            that has one, in that order. These are differentiated automatically
            -- there is no need to flag them, because the parameters that
            produced them live upstream in JAX rather than in Meep. Omit when
            there are none.
          designs: a list of design variables as 1D, 2D, or 3D JAX arrays. Valid shapes for
          design variables are (Nx, Ny, Nz) where Nx{y,z} match the elements of the
          `grid_size` constructor argument of Meep's `MaterialGrid` used for the
          corresponding design region. Singleton dimensions of the `grid_size` may be
          omitted from the corresponding design variable. For example, a design variable
          with a shape of either (10, 20) or (10, 20, 1) would be compatible with a
          `grid_size` of (10, 20, 1). Similarly, a design variable with shapes of (25,),
          (25, 1), or (25, 1, 1) would be compatible with a `grid_size` of (25, 1, 1).

        Returns:
          The differentiable values of the monitors. If every monitor yields one
          value per frequency -- as `EigenmodeCoefficient` and `LDOS` do -- this
          is a single complex-valued JAX ndarray with a shape of (num monitors,
          num frequencies). If the monitors have heterogeneous shapes, for
          example when a `FourierFields` monitor contributes a whole plane of
          values, it is a tuple with one array per monitor, in the order the
          monitors were given.
        """
        if _active_recorder is not None:
            return _active_recorder(self, designs, sources)
        return self._simulate_fn(designs, sources)

    def _update_sources(self, source_variables) -> None:
        """Push JAX-supplied amplitudes into the differentiable sources."""
        if source_variables is None:
            return
        if len(source_variables) != len(self.differentiable_sources):
            raise ValueError(
                f"Got {len(source_variables)} source arrays but "
                f"{len(self.differentiable_sources)} differentiable sources."
            )
        # Meep wants amp_data as a 3D array, but the caller may hand in any
        # shape with the same number of entries; the cotangent has to go back
        # in the shape they used, so remember it.
        self._source_shapes = [onp.shape(a) for a in source_variables]
        for src, amps in zip(self.differentiable_sources, source_variables):
            src.amp_data = onp.asarray(amps, dtype=onp.complex128).reshape(
                onp.shape(src.amp_data)
            )

    def _run_fwd_simulation(
        self,
        design_variables: Iterable[onp.ndarray],
        source_variables: Optional[Iterable[onp.ndarray]] = None,
    ) -> Tuple[Union[jnp.ndarray, Tuple[jnp.ndarray, ...]], List[List[mp.DftFields]]]:
        """Runs forward simulation, returning monitor values and design region fields."""
        utils.validate_and_update_design(self.design_regions, design_variables)
        self._update_sources(source_variables)
        self.simulation.reset_meep()
        self.simulation.change_sources(self.sources)
        utils.register_monitors(
            self.monitors, self.frequencies, pade_samples=self.pade_samples
        )
        fwd_design_region_monitors = utils.install_design_region_monitors(
            self.simulation,
            self.design_regions,
            self.frequencies,
            pade_samples=self.pade_samples,
        )
        self.simulation.init_sim()
        self._run_simulation(
            [
                [monitor._monitor for monitor in self.monitors],
                fwd_design_region_monitors,
            ],
            "forward_monitor_only",
        )

        monitor_values = utils.gather_monitor_values(self.monitors)
        # A tuple when the monitors have heterogeneous shapes; `custom_vjp`
        # handles either, and the cotangent arrives with a matching structure.
        monitor_values = jax.tree_util.tree_map(jnp.asarray, monitor_values)
        return (monitor_values, fwd_design_region_monitors)

    def _run_adjoint_simulation(
        self, monitor_values_grad: Union[onp.ndarray, Sequence[onp.ndarray]]
    ) -> List[List[mp.DftFields]]:
        """Runs adjoint simulation, returning design region fields."""
        if not self.design_regions:
            raise RuntimeError(
                "An adjoint simulation was attempted when no design "
                "regions are present."
            )
        adjoint_sources = utils.create_adjoint_sources(
            self.monitors, monitor_values_grad
        )
        # TODO refactor with optimization_problem.py #
        self.simulation.restart_fields()
        self.simulation.clear_dft_monitors()
        self.simulation.change_sources(adjoint_sources)
        #                                            #
        adj_design_region_monitors = utils.install_design_region_monitors(
            self.simulation,
            self.design_regions,
            self.frequencies,
            pade_samples=self.pade_samples,
        )
        self.adj_source_monitors = source_gradient.install_source_gradient_monitors(
            self.simulation,
            self.differentiable_sources,
            self.frequencies,
        )
        self.simulation.init_sim()
        self._run_simulation(adj_design_region_monitors, "adjoint_monitor_only")

        return adj_design_region_monitors

    def _source_vjps(self, sum_freq_partials: bool = True) -> List[onp.ndarray]:
        """Cotangents with respect to each differentiable source's currents.

        With `sum_freq_partials` false the frequency axis is kept, which is what
        `value_and_jacobian` needs to keep its rows separable.
        """
        if not self.differentiable_sources:
            return []
        phase = self.monitors[0]._adj_src_phase()
        out = []
        for src, group in zip(self.differentiable_sources, self.adj_source_monitors):
            # one monitor per component the source drives; an amp_data source
            # drives exactly one
            monitor = group[0]
            scale = source_gradient.source_grad_scale(
                self.simulation, src, self.frequencies, phase
            )
            currents = monitor.gather(self.simulation, scale)
            positions = monitor.positions(self.simulation)
            shape = self._source_shapes[len(out)]
            # one row per frequency, each carried back through the trilinear
            # interpolation Meep applies to amp_data
            rows = onp.stack(
                [
                    source_gradient.amp_data_gradient(
                        self.simulation, src, positions, currents[f_index]
                    )
                    for f_index in range(currents.shape[0])
                ]
            )
            if sum_freq_partials:
                # the amplitudes are shared across the band, as design weights are
                out.append(onp.sum(rows, axis=0).reshape(shape))
            else:
                out.append(rows.reshape((rows.shape[0], *shape)))
        return out

    def _calculate_vjps(
        self,
        fwd_fields: List[List[mp.DftFields]],
        adj_fields: List[List[mp.DftFields]],
        design_variable_shapes: List[Tuple[int, ...]],
        sum_freq_partials: bool = True,
    ) -> List[onp.ndarray]:
        """Calculates the VJP for a given set of forward and adjoint fields."""
        return utils.calculate_vjps(
            self.simulation,
            self.design_regions,
            self.frequencies,
            fwd_fields,
            adj_fields,
            design_variable_shapes,
            sum_freq_partials=sum_freq_partials,
            finite_difference_step=self.finite_difference_step,
        )

    def _initialize_callable(self) -> Callable:
        """Initializes the callable JAX function and registers its VJP."""

        @jax.custom_vjp
        def simulate(design_variables: List[jnp.ndarray], source_variables):
            monitor_values, _ = self._run_fwd_simulation(
                design_variables, source_variables
            )
            return monitor_values

        def _simulate_fwd(design_variables, source_variables):
            """Runs forward simulation, returning monitor values and fields."""
            monitor_values, self.fwd_design_region_monitors = self._run_fwd_simulation(
                design_variables, source_variables
            )
            design_variable_shapes = [x.shape for x in design_variables]
            return monitor_values, (design_variable_shapes)

        def _simulate_rev(res, monitor_values_grad):
            """Runs adjoint simulation, returning VJP of design wrt monitor values."""
            design_variable_shapes = res
            self.adj_design_region_monitors = self._run_adjoint_simulation(
                monitor_values_grad
            )
            vjps = self._calculate_vjps(
                self.fwd_design_region_monitors,
                self.adj_design_region_monitors,
                design_variable_shapes,
            )
            source_vjps = [jnp.asarray(v) for v in self._source_vjps()]
            return (
                [jnp.asarray(vjp) for vjp in vjps],
                source_vjps if self.differentiable_sources else None,
            )

        simulate.defvjp(_simulate_fwd, _simulate_rev)

        return simulate


class _WrapperCallRecorder:
    """Splits a loss function apart at its `MeepJaxWrapper` call.

    `value_and_jacobian` needs to see the user's loss function as three pieces:
    whatever maps the parameters to design weights before the simulation, the
    simulation itself, and the post-processing after it. The user writes a single
    function, so the split is recovered by having the wrapper cooperate. The first
    pass runs the forward simulation and records what crossed the boundary; later
    passes replay that recording instead of simulating, which leaves the pieces
    before and after the call as ordinary differentiable JAX.
    """

    def __init__(self):
        self.wrapper = None
        self.designs = None
        self.design_shapes = None
        self.monitor_values = None
        self.fwd_monitors = None
        self.sources = None
        self.num_calls = 0
        self.substitute = None

    def __call__(self, wrapper: "MeepJaxWrapper", designs, sources=None):
        self.num_calls += 1
        if self.wrapper is None:
            self.wrapper = wrapper
        elif self.wrapper is not wrapper:
            raise ValueError(
                "value_and_jacobian supports a loss function that calls a single "
                "MeepJaxWrapper; this one called two or more. Each simulation "
                "needs its own adjoint run, so compute their Jacobians "
                "separately and combine the rows yourself."
            )
        # Overwritten on every pass. On the first it holds concrete arrays; on
        # the traced passes it holds tracers, which is what the design mapping's
        # pullback is recovered from.
        self.designs = list(designs)
        self.sources = None if sources is None else list(sources)
        if self.monitor_values is None:
            self.design_shapes = [onp.shape(design) for design in designs]
            self.monitor_values, self.fwd_monitors = wrapper._run_fwd_simulation(
                designs, sources
            )
        return self.monitor_values if self.substitute is None else self.substitute


_active_recorder = None


@contextlib.contextmanager
def _recording(recorder: _WrapperCallRecorder):
    global _active_recorder
    previous, _active_recorder = _active_recorder, recorder
    try:
        yield
    finally:
        _active_recorder = previous


def _split_args(args: Tuple, argnums):
    """Partitions positional arguments into differentiated and held-fixed."""
    indices = (argnums,) if isinstance(argnums, int) else tuple(argnums)
    differentiated = tuple(args[i] for i in indices)

    def rebuild(new_values):
        merged = list(args)
        for position, index in enumerate(indices):
            merged[index] = new_values[position]
        return tuple(merged)

    return indices, differentiated, rebuild


def value_and_jacobian(loss: Callable, argnums=0) -> Callable:
    """Transforms a per-frequency loss function into one that also returns its Jacobian.

    Worst-case (minimax) optimization over a bandwidth needs a separate gradient
    for each frequency rather than the gradient of a scalar reduction over them.
    This is the counterpart of `jax.value_and_grad` for that case:

        def loss(rho, thickness, tilt):
            (dft,) = wrapped_meep([rho])              # a MeepJaxWrapper call
            return postprocess(dft, thickness, tilt)  # one value per frequency

        values, grads = mpa.value_and_jacobian(loss, argnums=(0, 1, 2))(
            rho, thickness, tilt)

    Every leaf of `grads` is the corresponding parameter's shape with a leading
    frequency axis, so `grads[1][f]` is the gradient of `values[f]` with respect
    to `thickness`. Parameters that Meep knows nothing about -- a layer stack fed
    to a propagator, a fiber tilt -- are differentiated alongside the design
    weights, and arbitrary pytrees are supported throughout.

    The whole Jacobian costs one forward simulation and one adjoint simulation,
    independent of the number of frequencies. `jax.jacrev` cannot achieve this:
    it evaluates the reverse pass once per output component, and here each
    evaluation would be a full timestepping run.

    As in `OptimizationProblem`, the dependence of the objective on the *monitor
    values* must be block diagonal in frequency -- entry `f` may only read the
    monitor values at frequency `f` -- because a single adjoint field cannot
    carry independent sensitivities for different frequencies. Dependence on
    parameters that bypass the simulation is unrestricted, since those are
    differentiated directly.

    Args:
      loss: a function whose positional arguments are the parameters and whose
        return value is a 1-D array with one entry per frequency. It must call a
        `MeepJaxWrapper` exactly once.
      argnums: which positional arguments to differentiate, following the
        convention of `jax.value_and_grad`.

    Returns:
      A function with the same signature as `loss`, returning
      `(values, jacobian)`.
    """

    def evaluate(*args):
        recorder = _WrapperCallRecorder()
        indices, differentiated, rebuild = _split_args(args, argnums)

        with _recording(recorder):
            # Pass 1: the only forward simulation.
            values = loss(*args)
            if recorder.num_calls == 0:
                raise ValueError(
                    "The loss function never called a MeepJaxWrapper, so there "
                    "is no simulation to differentiate. Use jax.jacrev for a "
                    "function that does not involve Meep."
                )
            if recorder.num_calls != 1:
                raise ValueError(
                    f"The loss function called a MeepJaxWrapper "
                    f"{recorder.num_calls} times; value_and_jacobian supports "
                    "exactly one call, since it drives a single adjoint "
                    "simulation."
                )
            wrapper = recorder.wrapper
            num_frequencies = onp.asarray(wrapper.frequencies).size
            if not isinstance(values, jax.Array) or jnp.ndim(values) != 1:
                raise TypeError(
                    "The loss function must return a 1-D JAX array with one "
                    f"entry per frequency, but it returned {type(values)} with "
                    f"shape {jnp.shape(values)}."
                )
            if values.shape[0] != num_frequencies:
                raise ValueError(
                    f"The loss function returned {values.shape[0]} values; "
                    f"value_and_jacobian expects one per frequency, i.e. "
                    f"({num_frequencies},). Use jax.value_and_grad for a scalar "
                    "loss."
                )

            # From here on the recorder replays the recorded monitor values, so
            # nothing below runs another simulation.

            # (1) How the parameters reach the loss without passing through the
            # simulation -- a layer stack, a fiber tilt, a regularizer.
            explicit = jax.jacrev(loss, argnums=argnums)(*args)

            # (2) How the loss depends on the monitor values. Seeding with ones
            # extracts the frequency diagonal of that Jacobian, which is what the
            # adjoint sources need; see `utils.objective_vjp`.
            def evaluate_at_monitor_values(monitor_values):
                recorder.substitute = monitor_values
                try:
                    return loss(*args)
                finally:
                    recorder.substitute = None

            _, pull_monitor_values = jax.vjp(
                evaluate_at_monitor_values, recorder.monitor_values
            )
            cotangent = pull_monitor_values(jnp.ones_like(values))[0]

            # (3) The only adjoint simulation. Keeping the frequency axis of the
            # design-region contraction is what makes the rows separable.
            adjoint_monitors = wrapper._run_adjoint_simulation(cotangent)
            rows = wrapper._calculate_vjps(
                recorder.fwd_monitors,
                adjoint_monitors,
                recorder.design_shapes,
                sum_freq_partials=False,
            )
            rows = [jnp.moveaxis(jnp.asarray(row), -1, 0) for row in rows]

            # (4) Carry each row back through whatever produced the design
            # weights, and the source amplitudes if any came from JAX. Pure
            # JAX, so `vmap` over the frequency axis is free.
            source_rows = [
                jnp.asarray(v) for v in wrapper._source_vjps(sum_freq_partials=False)
            ]

            def inputs_from(*values_to_differentiate):
                loss(*rebuild(values_to_differentiate))
                return recorder.designs, recorder.sources

            _, pull_inputs = jax.vjp(inputs_from, *differentiated)
            implicit = jax.vmap(pull_inputs)(
                (rows, source_rows if recorder.sources is not None else None)
            )

        if isinstance(argnums, int):
            implicit = implicit[0]
        jacobian = jax.tree_util.tree_map(lambda a, b: a + b, explicit, implicit)
        return values, jacobian

    return evaluate
