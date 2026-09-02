import math
import numbers
from typing import Any, Callable, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as onp
from autograd import make_vjp

import meep as mp

from . import ObjectiveQuantity

# Meep field components used to compute adjoint sensitivities
_ADJOINT_FIELD_COMPONENTS = [mp.Dx, mp.Dy, mp.Dz]
_ADJOINT_FIELD_COMPONENTS_CYL = [mp.Dr, mp.Dp, mp.Dz]

# The frequency axis in the array returned by `mp._get_gradient()`
_GRADIENT_FREQ_AXIS = 1

# The returned axis order from get_dft_array in cylindrical coordinates
_FREQ_AXIS = 0
_RHO_AXIS = 2
_PHI_AXIS = 3
_Z_AXIS = 1

# default finite difference step size when calculating Aᵤ
FD_DEFAULT = 1e-3


def validate_pade_options(pade_samples: int, pade_tolerance) -> int:
    """Validate the shared adjoint Padé controls and return a plain int."""
    if isinstance(pade_samples, (bool, onp.bool_)) or not isinstance(
        pade_samples, numbers.Integral
    ):
        raise TypeError("pade_samples must be an integer")
    pade_samples = int(pade_samples)
    if pade_samples < 0 or 0 < pade_samples < 4:
        raise ValueError("pade_samples must be 0 or at least 4")
    if pade_tolerance is not None:
        if pade_samples < 4:
            raise ValueError("pade_tolerance requires pade_samples >= 4")
        if (
            not onp.isscalar(pade_tolerance)
            or not onp.isfinite(pade_tolerance)
            or pade_tolerance <= 0
        ):
            raise ValueError("pade_tolerance must be a positive finite scalar")
    return pade_samples


def validate_finite_sources(sources: Iterable[mp.Source]) -> None:
    """Reject sources that cannot provide a complete post-source Padé window."""
    for source in sources:
        time_profile = source.src
        if isinstance(time_profile, mp.GaussianSource):
            end_time = (
                time_profile.start_time + 2 * time_profile.width * time_profile.cutoff
            )
        else:
            end_time = getattr(time_profile, "end_time", None)
        try:
            finite_end_time = math.isfinite(end_time)
        except TypeError:
            finite_end_time = False
        if not finite_end_time or end_time >= 1e19:
            raise ValueError("automatic Padé stopping requires finite-duration sources")


class DesignRegion:
    def __init__(
        self,
        design_parameters: Iterable[onp.ndarray],
        volume: mp.Volume = None,
        size: mp.Vector3 = None,
        center: mp.Vector3 = mp.Vector3(),
    ):
        self.volume = volume or mp.Volume(center=center, size=size)
        self.size = self.volume.size
        self.center = self.volume.center
        self.design_parameters = design_parameters
        self.num_design_params = design_parameters.num_params

    def update_design_parameters(self, design_parameters) -> None:
        self.design_parameters.update_weights(design_parameters)

    def update_beta(self, beta: float) -> None:
        self.design_parameters.beta = beta

    def get_gradient(
        self,
        sim: mp.Simulation,
        fields_a: List[mp.DftFields],
        fields_f: List[mp.DftFields],
        frequencies: List[float],
        finite_difference_step: float,
    ) -> onp.ndarray:
        num_freqs = onp.array(frequencies).size
        """We have the option to linearly scale the gradients up front
        using the scalegrad parameter (leftover from MPB API). Not
        currently needed for any existing feature (but available for
        future use)"""
        scalegrad = 1
        grad = onp.zeros((num_freqs, self.num_design_params))  # preallocate
        vol = sim._fit_volume_to_simulation(self.volume)
        # compute the gradient
        mp._get_gradient(
            grad,
            scalegrad,
            fields_a[0].swigobj,
            fields_a[1].swigobj,
            fields_a[2].swigobj,
            fields_f[0].swigobj,
            fields_f[1].swigobj,
            fields_f[2].swigobj,
            sim.gv,
            onp.array(frequencies),
            sim.geps,
            finite_difference_step,
            getattr(self.design_parameters, "grid_id", -1),
        )
        return onp.squeeze(grad).T


def _check_if_cylindrical(sim: mp.Simulation) -> bool:
    return sim.is_cylindrical or (sim.dimensions == mp.CYLINDRICAL)


def _compute_components(sim: mp.Simulation) -> List[int]:
    return (
        _ADJOINT_FIELD_COMPONENTS_CYL
        if _check_if_cylindrical(sim)
        else _ADJOINT_FIELD_COMPONENTS
    )


def calculate_vjps(
    simulation: mp.Simulation,
    design_regions: List[DesignRegion],
    frequencies: List[float],
    fwd_fields: List[List[mp.DftFields]],
    adj_fields: List[List[mp.DftFields]],
    design_variable_shapes: List[Tuple[int, ...]],
    sum_freq_partials: bool = True,
    finite_difference_step: float = FD_DEFAULT,
) -> List[onp.ndarray]:
    """Calculates the VJP for a given set of forward and adjoint fields."""
    vjps = [
        design_region.get_gradient(
            simulation,
            adj_fields[i],
            fwd_fields[i],
            frequencies,
            finite_difference_step,
        )
        for i, design_region in enumerate(design_regions)
    ]
    # `DesignRegion.get_gradient` squeezes the frequency axis away when there is
    # only one frequency, so restore it rather than indexing an axis that may not
    # be there.
    num_freqs = onp.asarray(frequencies).size
    vjps = [vjp.reshape(-1, num_freqs) for vjp in vjps]
    if sum_freq_partials:
        vjps = [
            onp.sum(vjp, axis=_GRADIENT_FREQ_AXIS).reshape(shape)
            for vjp, shape in zip(vjps, design_variable_shapes)
        ]
    else:
        vjps = [
            vjp.reshape(shape + (-1,))
            for vjp, shape in zip(vjps, design_variable_shapes)
        ]
    return vjps


def register_monitors(
    monitors: List[ObjectiveQuantity],
    frequencies: List[float],
    pade_samples: int = 0,
) -> None:
    """Registers a list of monitors."""
    for monitor in monitors:
        if pade_samples:
            monitor.register_monitors(frequencies, pade_samples=pade_samples)
        else:
            monitor.register_monitors(frequencies)


def install_design_region_monitors(
    simulation: mp.Simulation,
    design_regions: List[DesignRegion],
    frequencies: List[float],
    decimation_factor: int = 0,
    pade_samples: int = 0,
) -> List[List[mp.DftFields]]:
    """Installs DFT field monitors at the design regions of the simulation."""
    return [
        [
            simulation.add_dft_fields(
                [comp],
                frequencies,
                where=design_region.volume,
                yee_grid=True,
                decimation_factor=decimation_factor,
                persist=True,
                pade_samples=pade_samples,
            )
            for comp in _compute_components(simulation)
        ]
        for design_region in design_regions
    ]


def gather_monitor_values(
    monitors: Sequence[ObjectiveQuantity],
) -> Union[onp.ndarray, Tuple[onp.ndarray, ...]]:
    """Gathers the values of a list of objective quantities.

    Args:
      monitors: the objective quantities.

    Returns:
      If every monitor yields one value per frequency -- as `EigenmodeCoefficient`
      and `LDOS` do -- a rank-2 ndarray whose dimensions are (monitor,
      frequency). Otherwise, for example when a `FourierFields` monitor
      contributes a whole plane of values, a tuple with one array per monitor.
    """
    values = [onp.atleast_1d(onp.asarray(monitor())) for monitor in monitors]
    shapes = {value.shape for value in values}
    is_stackable = bool(values) and len(shapes) == 1 and values[0].ndim == 1
    return onp.stack(values) if is_stackable else tuple(values)


def validate_and_update_design(
    design_regions: List[DesignRegion], design_variables: Iterable[onp.ndarray]
) -> None:
    """Validate the design regions and variables.

    In particular the design variable should be 1,2,3-D and the design region
    shape should match the design variable shape after dimension expansion.
    The arguments are modified in place.

    Args:
      design_regions: List of mpa.DesignRegion,
      design_variables: Iterable with numpy arrays representing design variables.

    Raises:
      ValueError if the validation of dimensions fails.
    """
    for i, (design_region, design_variable) in enumerate(
        zip(design_regions, design_variables)
    ):
        if design_variable.ndim not in [1, 2, 3]:
            raise ValueError(
                f"Design variables should be 1D, 2D, or 3D, but the design variable at index {i} had a shape of {design_variable.shape}."
            )

        design_region_shape = tuple(
            int(x) for x in design_region.design_parameters.grid_size
        )
        design_variable_padded_shape = design_variable.shape + (1,) * (
            3 - design_variable.ndim
        )
        if design_variable_padded_shape != design_region_shape:
            raise ValueError(
                f"The design variable at index {i} with a shape of {design_variable.shape} is incompatible with the associated design region, which has a shape of {design_region_shape}."
            )

        design_variable = onp.asarray(design_variable, dtype=onp.float64)
        # Update the design variable in Meep
        design_region.update_design_parameters(design_variable.flatten())


_VJP_BACKENDS: List[Tuple[Callable[[Any], bool], Callable]] = []


def register_vjp_backend(
    matches: Callable[[Any], bool],
    vjp: Callable[[Callable, Tuple, Any], Sequence[onp.ndarray]],
) -> None:
    """Registers a way to differentiate objective functions of some flavor.

    This is how support for an autodifferentiation framework other than autograd
    is added without the rest of `meep.adjoint` importing it, and without the
    user having to wrap or annotate their objective function. `wrapper` calls
    this for JAX when it is imported, which happens only if JAX is installed.

    Args:
      matches: given the value an objective function returned, whether this
        backend should differentiate that function. Recognizing the framework's
        own array type is the intended test.
      vjp: called as `vjp(objective_function, args, cotangent)` and returning one
        cotangent per element of `args`, as plain NumPy arrays.
    """
    _VJP_BACKENDS.append((matches, vjp))


def _select_vjp_backend(value) -> Optional[Callable]:
    """Returns the backend that claims `value`, or None for autograd."""
    if value is None:
        return None
    for matches, vjp in _VJP_BACKENDS:
        if matches(value):
            return vjp
    return None


def objective_vjp(
    objective_function: Callable,
    args: Sequence[onp.ndarray],
    cotangent,
    value=None,
) -> Tuple[onp.ndarray, ...]:
    """Pulls a cotangent back through an objective function.

    This is a *single* reverse pass that yields the cotangent with respect to
    every element of `args` at once.

    Objective functions are differentiated with autograd by default. Two things
    override that, in order:

    1. a `vjp(cotangent, *args)` method on the objective function itself, for a
       hand-written or analytic pullback;
    2. a backend registered with `register_vjp_backend` whose predicate accepts
       `value`, which is how an objective function written with another
       autodifferentiation framework is differentiated by that framework without
       the user having to declare anything. `wrapper` registers one for JAX.

    Args:
      objective_function: the objective function, called as
        `objective_function(*args)`.
      args: the values of the objective quantities, in registration order.
      cotangent: the seed of the reverse pass, in the vector space of the
        objective function's output.
      value: what `objective_function(*args)` returned, used to select a backend.
        If omitted, autograd is used.

    Returns:
      One cotangent per element of `args`, each with that element's shape.
    """
    args = tuple(args)
    user_vjp = getattr(objective_function, "vjp", None)
    if callable(user_vjp):
        cotangents = tuple(user_vjp(cotangent, *args))
    else:
        backend = _select_vjp_backend(value)
        if backend is not None:
            cotangents = tuple(backend(objective_function, args, cotangent))
        else:
            vjp, _ = make_vjp(lambda packed: objective_function(*packed))(args)
            cotangents = tuple(vjp(cotangent))
    if len(cotangents) != len(args):
        raise ValueError(
            f"The objective function's vjp returned {len(cotangents)} "
            f"cotangents for {len(args)} objective arguments."
        )
    return cotangents


def create_adjoint_sources(
    monitors: Sequence[ObjectiveQuantity],
    monitor_values_grad: Sequence[onp.ndarray],
) -> List[mp.Source]:
    """Places the adjoint sources for one adjoint simulation.

    Args:
      monitors: the objective quantities registered in the forward run.
      monitor_values_grad: one cotangent per monitor, each with the shape of
        that monitor's own value.

    Returns:
      The list of adjoint sources.
    """
    # Everything downstream -- `_adj_src_scale`, `FilteredSource`, and the
    # `std::complex<double>` buffers the C++ source constructors read -- works in
    # double precision, so upcast here rather than round-tripping through the
    # field precision.
    cotangents = [onp.asarray(dj, dtype=onp.complex128) for dj in monitor_values_grad]
    if len(cotangents) != len(monitors):
        raise ValueError(
            f"Got {len(cotangents)} cotangents for {len(monitors)} monitors."
        )
    if not any(onp.any(dj) for dj in cotangents):
        raise RuntimeError(
            "The gradient of all monitor values is zero, which "
            "means that no adjoint sources can be placed to set "
            "up an adjoint simulation in Meep. Possible causes "
            "could be:\n\n"
            " * the forward simulation was not run for long enough "
            "to allow the input pulse(s) to reach the monitors"
            " * the monitor values are disconnected from the "
            "objective function output."
        )
    adjoint_sources = []
    for monitor, dj in zip(monitors, cotangents):
        if onp.any(dj):
            adjoint_sources += monitor.place_adjoint_source(dj)

    # An empty list here is not an error under MPI.  `fourier_sourcedata`
    # returns one entry per *local* chunk that intersects the monitor, and
    # `IndexedSource` is inherently per-chunk, so a rank whose chunks miss the
    # monitor legitimately places nothing -- which is every rank but one for a
    # point monitor.  The condition worth checking is global and is already
    # checked above: the cotangents are identical on every rank, so if the
    # objective really is disconnected from the monitors, the RuntimeError
    # fires everywhere at once with a message that says so.
    return adjoint_sources
