from collections import namedtuple
from typing import Callable, List, Union, Optional, Tuple

import numpy as np

import meep as mp
from meep.simulation import stop_when_dft_pade_converged

from . import LDOS, DesignRegion, utils, ObjectiveQuantity
from . import source_gradient
from . import geometry_gradient


class OptimizationProblem:
    """Top-level class in the Meep adjoint module.

    Intended to be instantiated from user scripts with mandatory constructor
    input arguments specifying the data required to define an adjoint-based
    optimization.

    The class knows how to do one basic thing: Given an input vector
    of design variables, compute the objective function value (forward
    calculation) and optionally its gradient (adjoint calculation).
    This is done using the __call__ method.
    """

    def __init__(
        self,
        simulation: mp.Simulation,
        objective_functions: List[Callable],
        objective_arguments: List[ObjectiveQuantity],
        design_regions: List[DesignRegion],
        frequencies: Optional[Union[np.ndarray, List[float]]] = None,
        fcen: Optional[float] = None,
        df: Optional[float] = None,
        nf: Optional[int] = None,
        decay_by: Optional[float] = 1e-11,
        decimation_factor: Optional[int] = 0,
        minimum_run_time: Optional[float] = 0,
        maximum_run_time: Optional[float] = None,
        finite_difference_step: Optional[float] = utils.FD_DEFAULT,
        step_funcs: Optional[List[Callable]] = None,
        pade_samples: int = 0,
        pade_tolerance: Optional[float] = None,
    ):
        """Initialize an instance of OptimizationProblem.

        Args:
          simulation: the `meep.Simulation` object that describes the problem
            (i.e., defining sources, geometry, boundary layers, etc).
          objective_functions: list of differentiable functions (callable
            objects) whose arguments are given by objective_arguments. The
            functions should take all objective_arguments as arguments even if
            not all of them are used by each function. For example, if we are
            interested in functions f(A,B) and g(B,C) of quantities A, B, C,
            then objective_functions must be [f1, g1] where
            f1 = lambda A, B, C: f(A,B) and g1 = lambda A, B, C: g(B,C), and
            objective_arguments must be [A, B, C].
          objective_arguments: list of ObjectiveQuantity objects passed as
            arguments to the objective_functions.
          design_regions: list of DesignRegion objects to be optimized.
          frequencies: array/list of frequencies of interest in the problem.
            If not specified then the list of frequencies will be created from
            fcen, df, and nf as a list of size nf that goes from fcen-df/2 to
            fcen+df/2.
          fcen: the center frequency.
          df: the frequency width (i.e., maximum frequency - minimum frequency).
          nf: number of frequencies.
          decay_by: the threshold value by which all field components at each
            frequency of every DFT object have to decay relative to their
            maximum before the simulation is terminated. Default is 1e-11.
          decimation_factor: an integer used to specify the number of timesteps
            between updates of the DFT fields. The default is 0, at which the
            value is automatically determined from the Nyquist rate of the
            bandwidth-limited sources and the DFT monitors. Can be disabled by
            setting it to 1.
          minimum_run_time: the minimum runtime for each simulation. Default
            is 0.
          maximum_run_time: the maximum runtime for each simulation.
          finite_difference_step: the step size for calculation of the
            finite-difference gradients.
          step_funcs: list of step functions to be called at each timestep.
          pade_samples: number of recent DFT samples retained for bounded-memory
            Padé extrapolation. Zero (the default) disables extrapolation.
          pade_tolerance: optional convergence tolerance for automatic Padé stopping.
            Requires `pade_samples >= 4` and finite-duration sources.
        """

        self.step_funcs = step_funcs if step_funcs is not None else []
        self.sim = simulation

        if isinstance(objective_functions, list):
            self.objective_functions = objective_functions
        else:
            self.objective_functions = [objective_functions]
        self.objective_arguments = objective_arguments
        self.f_bank = []  # objective function evaluation history

        if isinstance(design_regions, list):
            self.design_regions = design_regions
        else:
            self.design_regions = [design_regions]

        self.num_design_params = [ni.num_design_params for ni in self.design_regions]
        self.num_design_regions = len(self.design_regions)

        if frequencies is not None:
            if isinstance(frequencies, (np.ndarray, List)):
                self.frequencies = frequencies
                if isinstance(frequencies, np.ndarray):
                    self.nf = frequencies.size
                else:
                    self.nf = len(frequencies)
            else:
                raise TypeError(
                    "frequencies argument of OptimizationProblem "
                    "constructor must be a Numpy array or List."
                )
        elif nf == 1:
            self.nf = nf
            self.frequencies = [fcen]
        else:
            fmax = fcen + 0.5 * df
            fmin = fcen - 0.5 * df
            dfreq = (fmax - fmin) / (nf - 1)
            self.frequencies = np.linspace(
                fmin,
                fmin + dfreq * nf,
                num=nf,
                endpoint=False,
            )
            self.nf = nf

        if self.nf == 1:
            self.fcen_idx = 0
        else:
            self.fcen_idx = int(
                np.argmin(
                    np.abs(
                        np.asarray(self.frequencies)
                        - np.mean(np.asarray(self.frequencies))
                    )
                    ** 2
                )
            )  # index of center frequency

        self.decay_by = decay_by
        self.decimation_factor = decimation_factor
        self.minimum_run_time = minimum_run_time
        self.maximum_run_time = maximum_run_time
        self.pade_samples = utils.validate_pade_options(pade_samples, pade_tolerance)
        self.pade_tolerance = pade_tolerance
        if self.pade_samples and any(
            isinstance(m, LDOS) for m in self.objective_arguments
        ):
            raise ValueError("Padé extrapolation does not support LDOS objectives")
        self.pade_diagnostics = []
        self.finite_difference_step = (
            finite_difference_step  # step size used in Aᵤ computation
        )

        # store sources for finite difference estimations
        self.forward_sources = self.sim.sources

        # Sources flagged with `differentiable=[...]` are a second category of
        # differentiable input alongside the design regions. Their gradient is
        # the adjoint field sampled over the source's own support, so it needs
        # no simulation beyond the adjoint run already being performed.
        self.differentiable_sources = source_gradient.differentiable_sources(self.sim)
        self.source_gradient = {}

        # Geometric objects flagged with `differentiable=[...]`: the shape
        # derivative of the structure itself, as opposed to the density inside
        # a design region.
        self.differentiable_objects = geometry_gradient.differentiable_objects(self.sim)
        self.geometry_gradient = {}

        # The optimizer has three allowable states : "INIT", "FWD", and "ADJ".
        #    INIT - The optimizer is initialized and ready to run a forward simulation
        #    FWD  - The optimizer has already run a forward simulation
        #    ADJ  - The optimizer has already run an adjoint simulation (but not yet calculated the gradient)
        self.current_state = "INIT"

        self.gradient = []
        self._objective_values = []

    def __call__(
        self,
        rho_vector: List[List[float]] = None,
        need_value: bool = True,
        need_gradient: bool = True,
        beta: float = None,
    ) -> Tuple[List[np.ndarray], List[List[np.ndarray]]]:
        """Evaluate value and/or gradient of objective function.

        Args:
          rho_vector: list of design weights (which is itself a list). Each
            list in the list represents the design weights for one design
            region. The design weights are updated to the specified values.
            The objective functions and their gradients are then evaluated
            using these design weights.
          need_value: whether forward simulations for evaluating the objective
            functions are necessary. Default is True.
          need_gradient: whether adjoint simulations for evaluating the
            gradients are necessary. Default is True.
          beta: the strength (or "bias") of projecting the design weights in
            rho_vector using a hyperbolic tangent function. Default is None.

        Returns:
          A 2-tuple (f0, gradient) for which:
            f0 is the list of values of the objective functions.
            gradient is a list (over objective functions) of lists (over design
            regions) of 2d arrays (design weights by frequencies) of
            derivatives. If there is only a single objective function, the
            outer 1-element list is replaced by just that element, and
            similarly if there is only one design region then those 1-element
            list are replaced by just those elements. In addition, if there is
            only one frequency then the innermost array is squeezed to a 1d
            array. For example, if there is only a single objective function,
            a single design region, and a single frequency, then gradient is
            simply a 1d array of the derivatives.
        """
        if rho_vector:
            self.update_design(rho_vector=rho_vector, beta=beta)

        # Run forward run if requested
        if need_value and self.current_state == "INIT":
            print("Starting forward run...")
            self.forward_run()

        # Run adjoint simulation and calculate gradient if requested
        if need_gradient:
            if self.current_state == "INIT":
                # we need to run a forward run before an adjoint run
                print("Starting forward run...")
                self.forward_run()
                print("Starting adjoint run...")
                self.adjoint_run()
                print("Calculating gradient...")
                self.calculate_gradient()
            elif self.current_state == "FWD":
                print("Starting adjoint run...")
                self.adjoint_run()
                print("Calculating gradient...")
                self.calculate_gradient()
            else:
                raise ValueError(
                    f"Incorrect solver state detected: {self.current_state}"
                )

        if self.differentiable_sources or self.differentiable_objects:
            # Only change the return shape when something extra was flagged;
            # otherwise this is exactly as before.
            return self.f0, {
                "design": self.gradient,
                **self.source_gradient,
                **self.geometry_gradient,
            }

        return self.f0, self.gradient

    def get_fdf_funcs(self) -> Tuple[Callable, Callable]:
        """Construct callable functions for objective functions and gradients.

        Returns
        -------
        2-tuple (f_func, df_func) of standalone (non-class-method) callables, where
            f_func(beta) = objective function value for design variables beta
           df_func(beta) = objective function gradient for design variables beta
        """

        def _f(x=None):
            (fq, _) = self.__call__(rho_vector=x, need_gradient=False)
            return fq

        def _df(x=None):
            (_, df) = self.__call__(need_value=False)
            return df

        return _f, _df

    def prepare_forward_run(self):
        # prepare forward run
        self.sim.reset_meep()

        # add forward sources
        self.sim.change_sources(self.forward_sources)

        # register user specified monitors
        self.forward_monitors = [
            (
                m.register_monitors(self.frequencies, pade_samples=self.pade_samples)
                if self.pade_samples
                else m.register_monitors(self.frequencies)
            )
            for m in self.objective_arguments
        ]

        # register design region
        self.forward_design_region_monitors = utils.install_design_region_monitors(
            self.sim,
            self.design_regions,
            self.frequencies,
            self.decimation_factor,
            self.pade_samples,
        )
        # an ordinary object has no monitor of its own, so its fields have to be
        # recorded in both runs
        self.forward_geometry_monitors = geometry_gradient.install_geometry_monitors(
            self.sim,
            self.differentiable_objects,
            self.frequencies,
            self.decimation_factor,
        )

    @staticmethod
    def _flatten_confirmation_values(values):
        flattened = []
        for value in values:
            array = np.asarray(value).ravel()
            flattened.extend(np.real(array))
            flattened.extend(np.imag(array))
        return np.asarray(flattened, dtype=float)

    def _forward_confirmation(self):
        """Values the Padé criterion requires to be stable, not just the fields.

        Uses the same backend dispatch as the adjoint itself rather than
        autograd's `jacobian` directly.  Hardcoding autograd here silently
        undid the framework-agnostic objective support: an objective written
        with `jax.numpy` got traced with autograd's ArrayBox, which JAX then
        refuses as an argument.  Going through `objective_vjp` also confirms
        stability of exactly the cotangents that drive the adjoint sources,
        which is the quantity the gradient actually depends on.
        """
        results = [monitor() for monitor in self.objective_arguments]
        values = list(results)
        for ar, objective in enumerate(self.objective_functions):
            value = objective(*results)
            values.append(value)
            values.extend(
                utils.objective_vjp(
                    objective,
                    results,
                    self._cotangent_for(value, ar),
                    value=value,
                )
            )
        return self._flatten_confirmation_values(values)

    def _run_with_stopping(
        self, monitors, confirmation=None, monitor_only=False, leg=None
    ):
        if self.pade_tolerance is None:
            self.sim.run(
                *self.step_funcs,
                until_after_sources=mp.stop_when_dft_decayed(
                    self.decay_by, self.minimum_run_time, self.maximum_run_time
                ),
            )
            return None

        utils.validate_finite_sources(self.sim.sources)
        condition = stop_when_dft_pade_converged(
            monitors,
            self.pade_tolerance,
            self.pade_samples,
            raw_decay_tol=self.decay_by,
            minimum_run_time=self.minimum_run_time,
            maximum_run_time=self.maximum_run_time,
            confirmation=confirmation,
            monitor_only=monitor_only,
        )
        self.sim.run(*self.step_funcs, until_after_sources=condition)
        diagnostics = condition.finalize(self.sim)
        diagnostics["leg"] = leg
        self.pade_diagnostics.append(diagnostics)
        if diagnostics["exit_reason"] == "maximum_time":
            raise RuntimeError(
                "automatic Padé convergence was not reached before maximum_run_time"
            )
        return diagnostics

    def forward_run(self):
        # set up monitors
        self.prepare_forward_run()

        # Forward run
        if any(isinstance(m, LDOS) for m in self.objective_arguments):
            self.sim.run(
                mp.dft_ldos(self.frequencies),
                *self.step_funcs,
                until_after_sources=mp.stop_when_dft_decayed(
                    self.decay_by, self.minimum_run_time, self.maximum_run_time
                ),
            )
        else:
            self._run_with_stopping(
                [self.forward_monitors, self.forward_design_region_monitors],
                confirmation=(
                    self._forward_confirmation
                    if self.pade_tolerance is not None
                    else None
                ),
                leg="forward",
            )

        # record objective quantities from user specified monitors
        self.results_list = [m() for m in self.objective_arguments]
        # evaluate objectives
        self._objective_values = [
            fi(*self.results_list) for fi in self.objective_functions
        ]
        self.f0 = list(self._objective_values)
        if len(self.f0) == 1:
            self.f0 = self.f0[0]

        # store objective function evaluation in memory
        self.f_bank.append(self.f0)

        # update solver's current state
        self.current_state = "FWD"

    def _objective_cotangent(self, ar: int) -> np.ndarray:
        return self._cotangent_for(self._objective_values[ar], ar)

    def _cotangent_for(self, value, ar: int) -> np.ndarray:
        """Builds the reverse-pass seed for one objective function.

        Meep runs a single adjoint simulation per objective function and recovers
        a gradient at every frequency from it, which is only valid when each
        component of a frequency-vector-valued objective depends on the monitor
        values at its own frequency alone. Under that assumption the quantity the
        adjoint sources need -- the frequency diagonal of the Jacobian -- is the
        pullback of a cotangent of ones.
        """
        if np.ndim(value) == 0:
            return np.ones((), dtype=np.asarray(value).dtype)
        value = np.asarray(value)
        if value.size not in (1, self.nf):
            raise ValueError(
                f"objective_functions[{ar}] returned {value.size} values. Meep's "
                "single-adjoint-simulation formulation requires each objective "
                f"function to be scalar or to return one value per frequency "
                f"({self.nf}). Pass independent scalar objectives as separate "
                "entries of `objective_functions` instead."
            )
        return np.ones(value.shape, dtype=value.dtype)

    def prepare_adjoint_run(self):
        # Compute adjoint sources. One reverse pass through each objective
        # function yields the cotangents for all of its objective arguments.
        self.adjoint_sources = [[] for _ in range(len(self.objective_functions))]
        for ar in range(len(self.objective_functions)):
            dJ = utils.objective_vjp(
                self.objective_functions[ar],
                self.results_list,
                self._objective_cotangent(ar),
                value=self._objective_values[ar],
            )
            self.adjoint_sources[ar] = utils.create_adjoint_sources(
                self.objective_arguments, dJ
            )

    def adjoint_run(self):
        # set up adjoint sources and monitors
        self.prepare_adjoint_run()

        is_cylindrical = utils._check_if_cylindrical(self.sim)
        original_m = self.sim.m if is_cylindrical else None
        original_k_point = self.sim.k_point
        current_adjoint_monitors = None
        completed = False
        self.adjoint_design_region_monitors = []
        self._streamed_gradients = []
        self.adjoint_source_monitors = []
        self.adjoint_geometry_monitors = []
        try:
            if is_cylindrical:
                self.sim.change_m(-original_m)
            if original_k_point:
                self.sim.change_k_point(-1 * original_k_point)

            for ar in range(len(self.objective_functions)):
                # Reset the fields
                self.sim.restart_fields()
                self.sim.clear_dft_monitors()
                current_adjoint_monitors = None

                # Update the sources
                self.sim.change_sources(self.adjoint_sources[ar])

                # register design dft fields
                current_adjoint_monitors = utils.install_design_region_monitors(
                    self.sim,
                    self.design_regions,
                    self.frequencies,
                    self.decimation_factor,
                    self.pade_samples,
                )
                self.adjoint_design_region_monitors.append(current_adjoint_monitors)

                # A monitor over each differentiable object's support, and over
                # each differentiable source's. These carry the shape and source
                # derivatives, so they get the same Padé treatment as the design
                # region -- and they are handed to the stopping condition below,
                # or the run could stop before *they* have converged.
                geometry_monitors = geometry_gradient.install_geometry_monitors(
                    self.sim,
                    self.differentiable_objects,
                    self.frequencies,
                    self.decimation_factor,
                    self.pade_samples,
                )
                self.adjoint_geometry_monitors.append(geometry_monitors)

                source_monitors = source_gradient.install_source_gradient_monitors(
                    self.sim,
                    self.differentiable_sources,
                    self.frequencies,
                    self.decimation_factor,
                    self.pade_samples,
                )
                self.adjoint_source_monitors.append(source_monitors)
                self.sim._evaluate_dft_objects()

                # Adjoint run
                self._run_with_stopping(
                    current_adjoint_monitors
                    + list(geometry_monitors)
                    + list(source_monitors),
                    leg=f"adjoint_{ar}",
                    monitor_only=True,
                )

                if self.pade_samples:
                    self._streamed_gradients.append(
                        [
                            dr.get_gradient(
                                self.sim,
                                current_adjoint_monitors[dri],
                                self.forward_design_region_monitors[dri],
                                self.frequencies,
                                self.finite_difference_step,
                            )
                            for dri, dr in enumerate(self.design_regions)
                        ]
                    )
                    # Only the design-region monitors can be released here; the
                    # geometry and source monitors are read after the loop.
                    for monitor_set in current_adjoint_monitors:
                        for monitor in monitor_set:
                            monitor.remove()
                    self.adjoint_design_region_monitors[-1] = None
                    current_adjoint_monitors = None
            completed = True
        finally:
            try:
                if not completed and current_adjoint_monitors is not None:
                    for monitor_set in current_adjoint_monitors:
                        for monitor in monitor_set:
                            if monitor.swigobj is not None:
                                monitor.remove()
            finally:
                if is_cylindrical:
                    self.sim.change_m(original_m)
                if original_k_point:
                    self.sim.change_k_point(original_k_point)

        # update optimizer's state
        self.current_state = "ADJ"

    def calculate_source_gradient(self):
        """Gather the adjoint field over each differentiable source's support.

        Returns a dict keyed by source name (or index), each holding a dict
        keyed by the parameter names the source declared in `differentiable`.
        """
        if not self.differentiable_sources:
            return {}

        out = {}
        for ar in range(len(self.objective_functions)):
            for si, src in enumerate(self.differentiable_sources):
                monitors = self.adjoint_source_monitors[ar][si]
                # the adjoint field carries the normalization the objective
                # quantity applied when it placed the adjoint source, so the
                # transpose has to use that same quantity's phase
                scale = source_gradient.source_grad_scale(
                    self.sim,
                    src,
                    self.frequencies,
                    self.objective_arguments[0]._adj_src_phase(),
                )
                grads = self._source_gradient_for(src, monitors, scale)
                key = source_gradient.source_key(src, si)
                if len(self.objective_functions) == 1:
                    out[key] = grads
                else:
                    out.setdefault(key, []).append(grads)
        return out

    def _source_gradient_for(self, src, monitors, scale):
        """Gradients for one differentiable source, keyed by parameter name."""
        beam_names = [name for name in src.differentiable if name.startswith("beam_")]
        other_names = [
            name for name in src.differentiable if not name.startswith("beam_")
        ]

        grads = {}
        plain = [n for n in other_names if n != "amp_data"]
        if plain:
            currents = monitors[0].gather(self.sim, scale)
            everything = source_gradient.contract(src, currents)
            grads.update({name: everything[name] for name in plain})
        if "amp_data" in other_names:
            # the transpose of the trilinear interpolation amp_file_func does
            currents = monitors[0].gather(self.sim, scale)
            grads["amp_data"] = np.squeeze(
                np.array(
                    [
                        source_gradient.amp_data_gradient(
                            self.sim,
                            src,
                            monitors[0].positions(self.sim),
                            np.asarray(currents)[f_index],
                        )
                        for f_index in range(len(np.atleast_1d(self.frequencies)))
                    ]
                )
            )

        if beam_names:
            # Contract the per-point cotangents onto the beam's own parameters,
            # one frequency at a time so the rows stay separable.
            cotangents_by_frequency = []
            for f_index in range(len(np.atleast_1d(self.frequencies))):
                cotangents = {}
                for monitor in monitors:
                    values = monitor.gather(self.sim, scale)
                    cotangents[monitor.component] = (
                        monitor.positions(self.sim),
                        np.asarray(values)[f_index].ravel(),
                    )
                cotangents_by_frequency.append(
                    source_gradient.beam_parameter_gradients(
                        self.sim,
                        src,
                        beam_names,
                        cotangents,
                        source_gradient.normal_index(src),
                        center=monitors[0].volume.center,
                    )
                )
            for name in beam_names:
                grads[name] = np.squeeze(
                    np.array([row[name] for row in cotangents_by_frequency])
                )
        return grads

    def calculate_geometry_gradient(self):
        """Shape derivatives for each flagged geometric object."""
        if not self.differentiable_objects:
            return {}
        out = {}
        for ar in range(len(self.objective_functions)):
            for oi, (object_index, obj) in enumerate(self.differentiable_objects):
                grads = geometry_gradient.gradient(
                    self.sim,
                    obj,
                    object_index,
                    self.forward_geometry_monitors[oi],
                    self.adjoint_geometry_monitors[ar][oi],
                    self.frequencies,
                )
                key = geometry_gradient.object_key(obj, object_index)
                if len(self.objective_functions) == 1:
                    out[key] = grads
                else:
                    out.setdefault(key, []).append(grads)
        return out

    def calculate_gradient(self):
        self.source_gradient = self.calculate_source_gradient()
        self.geometry_gradient = self.calculate_geometry_gradient()

        # Iterate through all design regions and calculate gradient
        if self.pade_samples:
            self.gradient = self._streamed_gradients
        else:
            self.gradient = [
                [
                    dr.get_gradient(
                        self.sim,
                        self.adjoint_design_region_monitors[ar][dri],
                        self.forward_design_region_monitors[dri],
                        self.frequencies,
                        self.finite_difference_step,
                    )
                    for dri, dr in enumerate(self.design_regions)
                ]
                for ar in range(len(self.objective_functions))
            ]

        for dri in range(self.num_design_regions):
            for i in range(3):
                # note that dft_fields::remove calls delete on its chunks, and the
                # destructor ~dft_chunk automatically removes it from the fields object
                self.forward_design_region_monitors[dri][i].remove()

        # Cleanup list of lists
        if len(self.gradient) == 1:
            self.gradient = self.gradient[0]  # only one objective function
            if len(self.gradient) == 1:
                self.gradient = self.gradient[
                    0
                ]  # only one objective function and one design region
        elif len(self.gradient[0]) == 1:
            self.gradient = [
                g[0] for g in self.gradient
            ]  # multiple objective functions but one design region
        # Return optimizer's state to initialization
        self.current_state = "INIT"

    def calculate_fd_gradient(
        self,
        num_gradients: int = 1,
        db: float = 1e-4,
        design_variables_idx: int = 0,
        filter: Callable = None,
    ) -> List[float]:
        """
        Estimate central difference gradients.

        Parameters
        ----------
        num_gradients ... : scalar
            number of gradients to estimate. Randomly sampled from parameters.
        db ... : scalar
            finite difference step size
        design_variables_idx ... : scalar
            which design region to pull design variables from

        Returns
        -----------
        fd_gradient ... : lists
            [number of objective functions][number of gradients]

        """
        if filter is None:
            filter = lambda x: x
        if num_gradients < self.num_design_params[design_variables_idx]:
            # randomly choose indices to loop estimate
            fd_gradient_idx = np.random.choice(
                self.num_design_params[design_variables_idx],
                num_gradients,
                replace=False,
            )
        elif num_gradients == self.num_design_params[design_variables_idx]:
            fd_gradient_idx = range(self.num_design_params[design_variables_idx])
        else:
            raise ValueError(
                "The requested number of gradients must be less than or equal to the total number of design parameters."
            )

        assert db < 0.2, "The step size of finite difference is too large."

        # cleanup simulation object
        self.sim.reset_meep()
        self.sim.change_sources(self.forward_sources)

        # preallocate result vector
        fd_gradient = []

        for k in fd_gradient_idx:

            b0 = np.ones((self.num_design_params[design_variables_idx],))
            b0[:] = self.design_regions[design_variables_idx].design_parameters.weights
            # -------------------------------------------- #
            # left function evaluation
            # -------------------------------------------- #
            self.sim.reset_meep()

            # assign new design vector
            in_interior = True  # b0[k] is not too close to the boundaries 0 and 1
            if b0[k] < db or b0[k] + db > 1:
                in_interior = False  # b0[k] is too close to 0 or 1

            if b0[k] >= db:
                b0[k] -= db
            self.design_regions[design_variables_idx].update_design_parameters(b0)

            # initialize design monitors
            self.forward_monitors = [
                (
                    m.register_monitors(
                        self.frequencies, pade_samples=self.pade_samples
                    )
                    if self.pade_samples
                    else m.register_monitors(self.frequencies)
                )
                for m in self.objective_arguments
            ]

            if any(isinstance(m, LDOS) for m in self.objective_arguments):
                self.sim.run(
                    mp.dft_ldos(self.frequencies),
                    *self.step_funcs,
                    until_after_sources=mp.stop_when_energy_decayed(
                        dt=1, decay_by=1e-11
                    ),
                )
            else:
                self.sim.run(
                    *self.step_funcs,
                    until_after_sources=mp.stop_when_dft_decayed(
                        self.decay_by, self.minimum_run_time, self.maximum_run_time
                    ),
                )

            # record final objective function value
            results_list = [m() for m in self.objective_arguments]
            fm = [fi(*results_list) for fi in self.objective_functions]

            # -------------------------------------------- #
            # right function evaluation
            # -------------------------------------------- #
            self.sim.reset_meep()

            # assign new design vector
            b0[k] += 2 * db if in_interior else db
            self.design_regions[design_variables_idx].update_design_parameters(b0)

            # initialize design monitors
            self.forward_monitors = [
                (
                    m.register_monitors(
                        self.frequencies, pade_samples=self.pade_samples
                    )
                    if self.pade_samples
                    else m.register_monitors(self.frequencies)
                )
                for m in self.objective_arguments
            ]

            # add monitor used to track dft convergence
            if any(isinstance(m, LDOS) for m in self.objective_arguments):
                self.sim.run(
                    mp.dft_ldos(self.frequencies),
                    *self.step_funcs,
                    until_after_sources=mp.stop_when_energy_decayed(
                        dt=1, decay_by=1e-11
                    ),
                )
            else:
                self.sim.run(
                    *self.step_funcs,
                    until_after_sources=mp.stop_when_dft_decayed(
                        self.decay_by, self.minimum_run_time, self.maximum_run_time
                    ),
                )

            # record final objective function value
            results_list = [m() for m in self.objective_arguments]
            fp = [fi(*results_list) for fi in self.objective_functions]

            # -------------------------------------------- #
            # estimate derivative
            # -------------------------------------------- #
            fd_gradient.append(
                [
                    np.squeeze((fp[fi] - fm[fi]) / db / (2 if in_interior else 1))
                    for fi in range(len(self.objective_functions))
                ]
            )

        # Cleanup singleton dimensions
        if len(fd_gradient) == 1:
            fd_gradient = fd_gradient[0]

        return fd_gradient, fd_gradient_idx

    def update_design(self, rho_vector: List[float], beta: float = None) -> None:
        """Update the design permittivity function.

        rho_vector ....... a list of numpy arrays that maps to each design region
        """
        for bi, b in enumerate(self.design_regions):
            if np.array(rho_vector[bi]).ndim > 1:
                raise ValueError(
                    "Each vector of design variables must contain only one dimension."
                )
            b.update_design_parameters(rho_vector[bi])
            if beta:
                b.update_beta(beta)

        self.sim.reset_meep()
        self.current_state = "INIT"

    def get_objective_arguments(self) -> List[float]:
        """Return list of evaluated objective arguments."""
        return [m.get_evaluation() for m in self.objective_arguments]

    def plot2D(self, init_opt=False, **kwargs) -> None:
        """Produce a graphical visualization of the geometry and/or fields,
        as appropriately autodetermined based on the current state of
        progress.
        """

        if init_opt:
            self.prepare_forward_run()

        self.sim.plot2D(**kwargs)


def atleast_3d(*arys):
    from numpy import array, asanyarray, newaxis

    """
    Modified version of numpy's `atleast_3d`

    Keeps one dimensional array data in first dimension, as
    opposed to moving it to the second dimension as numpy's
    version does. Keeps the meep dimensionality convention.

    View inputs as arrays with at least three dimensions.
    Parameters
    ----------
    arys1, arys2, ... : array_like
        One or more array-like sequences.  Non-array inputs are converted to
        arrays.  Arrays that already have three or more dimensions are
        preserved.
    Returns
    -------
    res1, res2, ... : ndarray
        An array, or list of arrays, each with ``a.ndim >= 3``.  Copies are
        avoided where possible, and views with three or more dimensions are
        returned.  For example, a 1-D array of shape ``(N,)`` becomes a view
        of shape ``(N, 1, 1)``, and a 2-D array of shape ``(M, N)`` becomes a
        view of shape ``(M, N, 1)``.
    """
    res = []
    for ary in arys:
        ary = asanyarray(ary)
        if ary.ndim == 0:
            result = ary.reshape(1, 1, 1)
        elif ary.ndim == 1:
            result = ary[:, newaxis, newaxis]
        elif ary.ndim == 2:
            result = ary[:, :, newaxis]
        else:
            result = ary
        res.append(result)
    return res[0] if len(res) == 1 else res
