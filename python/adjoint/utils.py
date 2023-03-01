from typing import Callable, Iterable, List, Tuple

import meep as mp
import numpy as onp
from autograd import tensor_jacobian_product

from . import ObjectiveQuantity, OptimizationProblem

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


def _make_at_least_nd(x: onp.ndarray, dims: int = 3) -> onp.ndarray:
    """Makes an array have at least the specified number of dimensions."""
    return onp.reshape(x, x.shape + onp.maximum(dims - x.ndim, 0) * (1,))


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
) -> None:
    """Registers a list of monitors."""
    for monitor in monitors:
        monitor.register_monitors(frequencies)


def install_design_region_monitors(
    simulation: mp.Simulation,
    design_regions: List[DesignRegion],
    frequencies: List[float],
    decimation_factor: int = 0,
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
            )
            for comp in _compute_components(simulation)
        ]
        for design_region in design_regions
    ]


def gather_monitor_values(monitors: List[ObjectiveQuantity]) -> onp.ndarray:
    """Gathers the mode monitor overlap values as a rank 2 ndarray.

    Args:
      monitors: the mode monitors.

    Returns:
      a rank-2 ndarray, where the dimensions are (monitor, frequency), of dtype
      complex128.  Note that these values refer to the mode as oriented (i.e. they
      are unidirectional).
    """
    monitor_values = [monitor() for monitor in monitors]
    monitor_values = onp.array(monitor_values)
    assert monitor_values.ndim in [1, 2]
    monitor_values = _make_at_least_nd(monitor_values, 2)
    return monitor_values


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


def create_adjoint_sources(
    monitors: Iterable[ObjectiveQuantity], monitor_values_grad: onp.ndarray
) -> List[mp.Source]:
    monitor_values_grad = onp.asarray(
        monitor_values_grad,
        dtype=onp.complex64 if mp.is_single_precision() else onp.complex128,
    )
    if not onp.any(monitor_values_grad):
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
    for monitor_idx, monitor in enumerate(monitors):
        # `dj` for each monitor will have a shape of (num frequencies,)
        dj = onp.asarray(
            monitor_values_grad[monitor_idx],
            dtype=onp.complex64 if mp.is_single_precision() else onp.complex128,
        )
        if onp.any(dj):
            adjoint_sources += monitor.place_adjoint_source(dj)
    assert adjoint_sources
    return adjoint_sources


def get_epigraph_nlopt(
    mpa_opt: OptimizationProblem,
    mapping: Callable,
    callback: Callable | None = None,
) -> tuple[Callable, Callable]:
    """Helper function to create the NLopt objective and vector constraints for an
    epigraph formulation of a meep adjoint optimization problem.

    Parameters
    ----------
    mpa_opt : OptimizationProblem
        Instance of a meep adjoint optimization problem.
    mapping : function
        A Python function that maps the optimization variables to a design
        distribution. Expects the following signature::

            mapping(x: np.ndarray) -> np.ndarray

        If `mapping` takes additional arguments, you can wrap it in a lambda before
        passing it to this function.
    callback : function, optional
        A function that is called after every evaluation of the constraints, useful for
        logging during optimization.
        Expects the following signature::

            callback(t: float, v: np.ndarray, f0: np.ndarray, grad: np.ndarray) -> None

        with the epigraph dummy objective `t`, design weights `v`, the current meep
        objective values `f0` and the gradients `grad`.

    Returns
    -------
    _objective : function
        The epigraph dummy objective.
    _constraints: function
        The epigraph constraints.
    """

    def _objective(x: onp.ndarray, gd: onp.ndarray) -> float:
        if gd.size > 0:
            gd[0] = 1
            gd[1:] = 0
        return x[0]

    def _constraints(result: onp.ndarray, x: onp.ndarray, gd: onp.ndarray) -> None:
        t, v = x[0], x[1:]

        f0, grad = mpa_opt([mapping(v)])

        # f0 -> (obj_funs, wavelengths), grad -> (obj_funs, dofs, wavelengths)
        f0 = onp.atleast_2d(f0)
        grad = onp.reshape(grad, (f0.shape[0], -1, f0.shape[-1]))

        f0 = onp.concatenate(f0, axis=0)
        grad = onp.concatenate(grad, axis=-1)

        for k in range(grad.shape[-1]):
            grad[:, k] = tensor_jacobian_product(mapping, 0)(v, grad[:, k])

        if gd.size > 0:
            gd[:, 0] = -1
            gd[:, 1:] = grad.T

        result[:] = onp.real(f0) - t

        if callback:
            callback(t, v, f0, grad)

    return _objective, _constraints
