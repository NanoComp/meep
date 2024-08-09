"""Shape optimization of a multilayer stack over a broad bandwidth.

The 1D stack consists of alternating materials of index N_LAYER_1 and N_LAYER_2
in the arrangement: N_LAYER_2, N_LAYER_1, N_LAYER_2, N_LAYER_1, ..., N_LAYER_2.

The design parameters are the N layer thicknesses: [t_1, t_2, ..., t_N]. N must
be odd.

The design objective involves minimizing the largest integral of the DFT fields
in the multilayer stack over two wavelengths (worst-case optimization).

Tutorial Reference:
https://meep.readthedocs.io/en/latest/Python_Tutorials/Adjoint_Solver/#shape-optimization-of-a-multilayer-stack
"""

import copy
from typing import Callable, List, Tuple

from autograd.extend import primitive, defvjp
from autograd import numpy as npa
from autograd import tensor_jacobian_product
import meep as mp
import meep.adjoint as mpa
import nlopt
import numpy as np


RESOLUTION_UM = 800
AIR_UM = 1.0
PML_UM = 1.0
NUM_LAYERS = 9
N_LAYER = (1.0, 1.3)
LAYER_PERTURBATION_UM = 1.0 / RESOLUTION_UM
DESIGN_WAVELENGTHS_UM = (0.95, 1.05)
MAX_LAYER_UM = max(DESIGN_WAVELENGTHS_UM) / (4 * min(N_LAYER))
DESIGN_REGION_RESOLUTION_UM = 10 * RESOLUTION_UM
DESIGN_REGION_UM = mp.Vector3(0, 0, NUM_LAYERS * MAX_LAYER_UM)
NZ_DESIGN_GRID = int(DESIGN_REGION_UM.z * DESIGN_REGION_RESOLUTION_UM) + 1
NZ_SIM_GRID = int(DESIGN_REGION_UM.z * RESOLUTION_UM) + 1
MAX_OPT_ITERATIONS = 30
NUM_OPT_REPEAT = 10

num_wavelengths = len(DESIGN_WAVELENGTHS_UM)
frequencies = [1 / wavelength_um for wavelength_um in DESIGN_WAVELENGTHS_UM]


def str_from_list(list_: List[float]) -> str:
    return "[" + ", ".join(f"{val:.4f}" for val in list_) + "]"


def design_region_to_grid(nz: int) -> np.ndarray:
    """Returns the coordinates of the 1D grid for the design region.

    Args:
      nz: number of grid points.

    Returns:
      The 1D coordinates of the grid points.
    """
    z_grid = np.linspace(-0.5 * DESIGN_REGION_UM.z, 0.5 * DESIGN_REGION_UM.z, nz)

    return z_grid


@primitive
def levelset_and_smoothing(layer_thickness_um: np.ndarray) -> np.ndarray:
    """Returns the density weights for a multilayer stack as a levelset.

    Args:
      layer_thickness_um: thickness of each layer in the stack.

    Returns:
      The density weights as a flattened (1D) array.
    """
    air_padding_um = 0.5 * (DESIGN_REGION_UM.z - np.sum(layer_thickness_um))

    weights = np.zeros(NZ_DESIGN_GRID)

    # Air padding at left edge
    z_start = 0
    z_end = int(air_padding_um * DESIGN_REGION_RESOLUTION_UM)
    weights[z_start:z_end] = 0

    z_start = z_end
    for j in range(NUM_LAYERS):
        z_end = z_start + int(layer_thickness_um[j] * DESIGN_REGION_RESOLUTION_UM)
        weights[z_start:z_end] = 1 if (j % 2 == 0) else 0
        z_start = z_end

    # Air padding at right edge
    z_end = z_start + int(air_padding_um * DESIGN_REGION_RESOLUTION_UM)
    weights[z_start:z_end] = 0

    # Smooth the design weights by downsampling from the design grid
    # to the simulation grid using bilinear interpolation.
    z_sim_grid = design_region_to_grid(NZ_SIM_GRID)
    z_design_grid = design_region_to_grid(NZ_DESIGN_GRID)
    smoothed_weights = np.interp(z_sim_grid, z_design_grid, weights)

    return smoothed_weights.flatten()


def levelset_and_smoothing_vjp(
    ans: np.ndarray, layer_thickness_um: np.ndarray
) -> Callable[[np.ndarray], np.ndarray]:
    """Returns a function for computing the vector-Jacobian product."""

    total_layer_thickness_um = np.sum(layer_thickness_um)
    air_padding_um = 0.5 * (DESIGN_REGION_UM.z - total_layer_thickness_um)

    jacobian = np.zeros((NZ_SIM_GRID, NUM_LAYERS))

    z_design_grid = design_region_to_grid(NZ_DESIGN_GRID)
    z_sim_grid = design_region_to_grid(NZ_SIM_GRID)

    for i in range(NUM_LAYERS):
        weights = np.zeros(NZ_DESIGN_GRID)

        # Air padding at left edge
        z_start = 0
        z_end = int(air_padding_um * DESIGN_REGION_RESOLUTION_UM)
        weights[z_start:z_end] = 0

        z_start = z_end
        for j in range(NUM_LAYERS):
            layer_um = layer_thickness_um[j]
            if j == i:
                layer_um += LAYER_PERTURBATION_UM
            z_end = z_start + int(layer_um * DESIGN_REGION_RESOLUTION_UM)
            weights[z_start:z_end] = 1 if (j % 2 == 0) else 0
            z_start = z_end

        # Air padding at right edge
        z_end = z_start + int(air_padding_um * DESIGN_REGION_RESOLUTION_UM)
        weights[z_start:z_end] = 0

        # Smooth the design weights by downsampling from the design grid
        # to the simulation grid using bilinear interpolation.
        smoothed_weights = np.interp(z_sim_grid, z_design_grid, weights)

        jacobian[:, i] = (smoothed_weights - ans) / LAYER_PERTURBATION_UM

    return lambda g: np.tensordot(g, jacobian, axes=1)


def multilayer_stack() -> mpa.OptimizationProblem:
    """Sets up the adjoint optimization of a multilayer stack.

    Returns:
      A `meep.adjoint.Optimization` callable object.
    """
    pml_layers = [mp.PML(thickness=PML_UM)]

    size_z_um = PML_UM + AIR_UM + DESIGN_REGION_UM.z + AIR_UM + PML_UM
    cell_size = mp.Vector3(0, 0, size_z_um)

    frequency_center = np.mean(frequencies)

    # Set source bandwidth to be larger than the range of design wavelengths.
    frequency_width = 1.2 * (np.max(frequencies) - np.min(frequencies))

    src_cmpt = mp.Ex
    src_pt = mp.Vector3(0, 0, -0.5 * size_z_um + PML_UM)
    sources = [
        mp.Source(
            mp.GaussianSource(frequency_center, fwidth=frequency_width),
            component=src_cmpt,
            center=src_pt,
        )
    ]

    mat_1 = mp.Medium(index=N_LAYER[0])
    mat_2 = mp.Medium(index=N_LAYER[1])

    matgrid = mp.MaterialGrid(
        mp.Vector3(0, 0, NZ_SIM_GRID),
        mat_1,
        mat_2,
        weights=np.ones(NZ_SIM_GRID),
        do_averaging=False,
    )

    matgrid_region = mpa.DesignRegion(
        matgrid,
        volume=mp.Volume(center=mp.Vector3(), size=DESIGN_REGION_UM),
    )

    geometry = [
        mp.Block(
            material=matgrid, size=matgrid_region.size, center=matgrid_region.center
        )
    ]

    sim = mp.Simulation(
        resolution=RESOLUTION_UM,
        cell_size=cell_size,
        dimensions=1,
        boundary_layers=pml_layers,
        sources=sources,
        geometry=geometry,
    )

    obj_args = [
        mpa.FourierFields(
            sim, volume=matgrid_region.volume, component=mp.Ex, yee_grid=True
        )
    ]

    def obj_func(dft_ex: np.ndarray) -> Tuple[float]:
        """Objective function for the optimization.

        Args:
          dft_ex: the discrete Fourier-transformed Ex fields as a 2D array of
            dimension (num_wavelengths, NX_DESIGN_GRID * NY_DESIGN_GRID).

        Returns:
          A tuple of the log of the integrals of Ex in the stack for each
          wavelength.
        """
        return npa.log(npa.sum(npa.absolute(dft_ex) ** 2, axis=1))

    opt = mpa.OptimizationProblem(
        simulation=sim,
        objective_functions=obj_func,
        objective_arguments=obj_args,
        design_regions=[matgrid_region],
        frequencies=frequencies,
    )

    return opt


def obj_func(epigraph_and_layer_thickness_um: np.ndarray, grad: np.ndarray) -> float:
    """Objective function for the epigraph formulation.

    Args:
      epigraph_and_layer_thickness_um: 1D array containing epigraph variable
        (first element) and design weights (remaining elements).
      grad: the gradient as a flattened (1D) array, modified in place.

    Returns:
      The scalar epigraph variable.
    """
    epigraph = epigraph_and_layer_thickness_um[0]

    if grad.size > 0:
        grad[0] = 1
        grad[1:] = 0

    return epigraph


def epigraph_constraint(
    result: float, epigraph_and_layer_thickness_um: np.ndarray, gradient: np.ndarray
) -> None:
    """Constraint function for the epigraph formulation.

    Args:
      result: evaluation of the constraint function, modified in place.
      epigraph_and_layer_thickness_um: 1D array containing the epigraph variable
        (first element) and the layer thicknesses (remaining elements).
      gradient: the Jacobian matrix with dimensions (num_wavelengths,
        1 + NUM_LAYERS), modified in place.
    """
    epigraph = epigraph_and_layer_thickness_um[0]
    layer_thickness_um = epigraph_and_layer_thickness_um[1:]
    design_weights = levelset_and_smoothing(layer_thickness_um)

    obj_val, grad = opt([design_weights])

    # Backpropagate the gradients through the levelset construction function.
    grad_backpropagate = np.zeros((NUM_LAYERS, num_wavelengths))
    defvjp(levelset_and_smoothing, levelset_and_smoothing_vjp)
    for k in range(num_wavelengths):
        grad_backpropagate[:, k] = tensor_jacobian_product(levelset_and_smoothing, 0)(
            layer_thickness_um, grad[:, k]
        )

    # TODO (oskooi): determine why factor of 0.5 is necessary for correctness.
    grad_backpropagate = 0.5 * grad_backpropagate

    if gradient.size > 0:
        gradient[:, 0] = -1  # gradient w.r.t. epigraph variable
        gradient[:, 1:] = grad_backpropagate.T

    result[:] = obj_val - epigraph

    print(
        f"iteration:, {current_iteration[0]:2d}, {str_from_list(result)}, "
        f"{str_from_list(layer_thickness_um)}"
    )

    obj_func_history.append(obj_val)
    epigraph_variable_history.append(epigraph)
    current_iteration[0] += 1


if __name__ == "__main__":
    opt = multilayer_stack()

    # Specify lower and upper bounds for layer thicknesses based on layer
    # thicknesses for a quarter-wavelength stack at the mean frequency.
    mean_wavelength_um = 1 / np.mean(1 / np.array(DESIGN_WAVELENGTHS_UM))
    mean_layer_thickness_um = np.zeros(2)
    mean_layer_thickness_um[0] = mean_wavelength_um / (4 * N_LAYER[0])
    mean_layer_thickness_um[1] = mean_wavelength_um / (4 * N_LAYER[1])
    layer_thickness_um_lower_bound = np.zeros(NUM_LAYERS)
    layer_thickness_um_upper_bound = np.zeros(NUM_LAYERS)
    fraction_thickness = 0.05
    layer_thickness_um_lower_bound[0::2] = (
        1 - fraction_thickness
    ) * mean_layer_thickness_um[1]
    layer_thickness_um_upper_bound[0::2] = (
        1 + fraction_thickness
    ) * mean_layer_thickness_um[1]
    layer_thickness_um_lower_bound[1::2] = (
        1 - fraction_thickness
    ) * mean_layer_thickness_um[0]
    layer_thickness_um_upper_bound[1::2] = (
        1 + fraction_thickness
    ) * mean_layer_thickness_um[0]

    # Insert bounds for the epigraph variable.
    epigraph_and_layer_thickness_um_lower_bound = np.insert(
        layer_thickness_um_lower_bound, 0, 0
    )
    epigraph_and_layer_thickness_um_upper_bound = np.insert(
        layer_thickness_um_upper_bound, 0, np.inf
    )
    epigraph_tolerance = [0.01] * num_wavelengths

    min_obj_val = np.inf
    min_epigraph_and_layer_thickness_um = np.zeros(NUM_LAYERS + 1)
    obj_func_history = []
    epigraph_variable_history = []

    for j in range(NUM_OPT_REPEAT):
        rng = np.random.RandomState()
        layer_thickness_um = layer_thickness_um_lower_bound + rng.rand(NUM_LAYERS) * (
            layer_thickness_um_upper_bound - layer_thickness_um_lower_bound
        )

        # Execute a single forward run before the start of the optimization and
        # set the initial epigraph variable to slightly larger than the
        # largest value of the objective function over the wavelengths.
        design_weights = levelset_and_smoothing(layer_thickness_um)
        epigraph, _ = opt([design_weights], need_gradient=False)
        fraction_max_epigraph = 0.2
        epigraph_initial = (1 + fraction_max_epigraph) * np.amax(epigraph)
        epigraph_and_layer_thickness_um_upper_bound[0] = epigraph_initial

        epigraph_and_layer_thickness_um = np.concatenate(
            ([epigraph_initial], layer_thickness_um)
        )

        solver = nlopt.opt(nlopt.LD_CCSAQ, NUM_LAYERS + 1)
        solver.set_lower_bounds(epigraph_and_layer_thickness_um_lower_bound)
        solver.set_upper_bounds(epigraph_and_layer_thickness_um_upper_bound)
        solver.set_min_objective(obj_func)
        solver.set_maxeval(MAX_OPT_ITERATIONS)
        solver.add_inequality_mconstraint(epigraph_constraint, epigraph_tolerance)
        solver.set_ftol_rel(0.05)

        obj_func_history[:] = []
        epigraph_variable_history[:] = []
        current_iteration = [0]

        epigraph_and_layer_thickness_um = solver.optimize(
            epigraph_and_layer_thickness_um
        )
        optimal_obj_val = solver.last_optimum_value()
        return_code = solver.last_optimize_result()

        print(f"optimal_obj_val:, {j:2d}, {optimal_obj_val}, {return_code}")
        print(
            f"optimal_layer_thickness_um:, {j:2d}, "
            f"{str_from_list(epigraph_and_layer_thickness_um[1:])}"
        )

        if optimal_obj_val < min_obj_val:
            min_obj_val = optimal_obj_val
            min_epigraph_and_layer_thickness_um = epigraph_and_layer_thickness_um
            min_obj_func_history = copy.deepcopy(obj_func_history)
            min_epigraph_variable_history = copy.deepcopy(epigraph_variable_history)

    # Save important optimization parameters and output for post processing.
    np.savez(
        "optimal_design.npz",
        RESOLUTION_UM=RESOLUTION_UM,
        AIR_UM=AIR_UM,
        PML_UM=PML_UM,
        NUM_LAYERS=NUM_LAYERS,
        N_LAYER=N_LAYER,
        DESIGN_WAVELENGTHS_UM=DESIGN_WAVELENGTHS_UM,
        MAX_LAYER_UM=MAX_LAYER_UM,
        DESIGN_REGION_UM=DESIGN_REGION_UM,
        DESIGN_REGION_RESOLUTION_UM=DESIGN_REGION_RESOLUTION_UM,
        NZ_DESIGN_GRID=NZ_DESIGN_GRID,
        NZ_SIM_GRID=NZ_SIM_GRID,
        MAX_OPT_ITERATIONS=MAX_OPT_ITERATIONS,
        epigraph_and_layer_thickness_um_lower_bound=epigraph_and_layer_thickness_um_lower_bound,
        epigraph_and_layer_thickness_um_upper_bound=epigraph_and_layer_thickness_um_upper_bound,
        obj_func_history=min_obj_func_history,
        epigraph_variable_history=min_epigraph_variable_history,
        epigraph_variable=min_epigraph_and_layer_thickness_um[0],
        layer_thickness_um=min_epigraph_and_layer_thickness_um[1:],
        epigraph_tolerance=epigraph_tolerance,
        optimal_obj_val=min_obj_val,
    )
