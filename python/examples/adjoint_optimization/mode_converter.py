"""Topology optimization of a waveguide mode converter.

The worst-case optimization is based on minimizing the maximum
of {R,1-T} where R (reflectance) is $|S_{11}|^2$ for mode 1
and T (transmittance) is $|S_{21}|^2$ for mode 2 across six
different wavelengths. The optimization uses the method of moving
asymptotes (MMA) algorithm from NLopt. The minimum linewidth
constraint is based on A.M. Hammond et al., Optics Express,
Vol. 29, pp. 23916-23938, (2021). doi.org/10.1364/OE.431188
"""

from typing import List, NamedTuple, Tuple

from autograd import numpy as npa, tensor_jacobian_product, grad
import matplotlib.pyplot as plt
import meep as mp
import meep.adjoint as mpa
import nlopt
import numpy as np


RESOLUTION_UM = 50
WAVELENGTH_MIN_UM = 1.26
WAVELENGTH_MAX_UM = 1.30
WAVEGUIDE_UM = mp.Vector3(0.4, 3.0, 0)
PADDING_UM = 0.6
PML_UM = 1.0
DESIGN_WAVELENGTHS_UM = (1.265, 1.270, 1.275, 1.285, 1.290, 1.295)
DESIGN_REGION_UM = mp.Vector3(1.6, 1.6, mp.inf)
DESIGN_REGION_RESOLUTION_UM = int(2 * RESOLUTION_UM)
NX_DESIGN_GRID = int(DESIGN_REGION_UM.x * DESIGN_REGION_RESOLUTION_UM) + 1
NY_DESIGN_GRID = int(DESIGN_REGION_UM.y * DESIGN_REGION_RESOLUTION_UM) + 1
MIN_LENGTH_UM = 0.15
SIGMOID_THRESHOLD_INTRINSIC = 0.5
SIGMOID_THRESHOLD_EROSION = 0.75
SIGMOID_THRESHOLD_DILATION = 1 - SIGMOID_THRESHOLD_EROSION
MODE_SYMMETRY = mp.ODD_Z
SILICON = mp.Medium(index=3.5)
SILICON_DIOXIDE = mp.Medium(index=1.5)

cell_um = mp.Vector3(
    PML_UM + WAVEGUIDE_UM.x + DESIGN_REGION_UM.x + WAVEGUIDE_UM.x + PML_UM,
    PML_UM + PADDING_UM + DESIGN_REGION_UM.y + PADDING_UM + PML_UM,
    0,
)
filter_radius_um = mpa.get_conic_radius_from_eta_e(
    MIN_LENGTH_UM, SIGMOID_THRESHOLD_EROSION
)
frequency_min = 1 / WAVELENGTH_MAX_UM
frequency_max = 1 / WAVELENGTH_MIN_UM
frequency_center = 0.5 * (frequency_min + frequency_max)
frequency_width = frequency_max - frequency_min
pml_layers = [mp.PML(thickness=PML_UM)]
frequencies = [1 / wavelength_um for wavelength_um in DESIGN_WAVELENGTHS_UM]
src_pt = mp.Vector3(-0.5 * cell_um.x + PML_UM, 0, 0)
refl_pt = mp.Vector3(-0.5 * cell_um.x + PML_UM + 0.5 * WAVEGUIDE_UM.x)
tran_pt = mp.Vector3(0.5 * cell_um.x - PML_UM - 0.5 * WAVEGUIDE_UM.x)
stop_cond = mp.stop_when_fields_decayed(50, mp.Ez, refl_pt, 1e-6)


def str_from_list(list_: List[float]) -> str:
    return "[" + ", ".join(f"{val:.4f}" for val in list_) + "]"


def border_masks() -> Tuple[np.ndarray, np.ndarray]:
    """Return border masks for the design region.

    The masks are used to prevent violations on constraints on the
    minimum feature size at the boundaries of the design region.

    Returns:
      A 2-tuple of 2D arrays for border masks for Si and SiO2.
    """
    x_grid = np.linspace(
        -DESIGN_REGION_UM.x / 2,
        DESIGN_REGION_UM.x / 2,
        NX_DESIGN_GRID,
    )
    y_grid = np.linspace(
        -DESIGN_REGION_UM.y / 2,
        DESIGN_REGION_UM.y / 2,
        NY_DESIGN_GRID,
    )
    xy_grid_x, xy_grid_y = np.meshgrid(
        x_grid,
        y_grid,
        sparse=True,
        indexing="ij",
    )

    left_waveguide_port = (xy_grid_x <= -DESIGN_REGION_UM.x / 2 + filter_radius_um) & (
        np.abs(xy_grid_y) <= WAVEGUIDE_UM.y / 2
    )
    right_waveguide_port = (xy_grid_x >= DESIGN_REGION_UM.x / 2 - filter_radius_um) & (
        np.abs(xy_grid_y) <= WAVEGUIDE_UM.y / 2
    )
    silicon_mask = left_waveguide_port | right_waveguide_port

    border_mask = (
        (xy_grid_x <= -DESIGN_REGION_UM.x / 2 + filter_radius_um)
        | (xy_grid_x >= DESIGN_REGION_UM.x / 2 - filter_radius_um)
        | (xy_grid_y <= -DESIGN_REGION_UM.y / 2 + filter_radius_um)
        | (xy_grid_y >= DESIGN_REGION_UM.y / 2 - filter_radius_um)
    )
    silicon_dioxide_mask = border_mask.copy()
    silicon_dioxide_mask[silicon_mask] = False

    return silicon_mask, silicon_dioxide_mask


def filter_projection(
    weights: np.ndarray, sigmoid_threshold: float, sigmoid_bias: float
) -> np.ndarray:
    """A differentiable function to filter and project the design weights.

    Args:
      weights: design weights as a flattened (1D) array.
      sigmoid_threshold: erosion/dilation parameter for the projection.
      sigmoid_bias: bias parameter for the projection. 0 is no projection.

    Returns:
      The mapped design weights as a 1D array.
    """
    silicon_mask, silicon_dioxide_mask = border_masks()

    weights_masked = npa.where(
        silicon_mask.flatten(),
        1,
        npa.where(
            silicon_dioxide_mask.flatten(),
            0,
            weights,
        ),
    )

    weights_filtered = mpa.conic_filter(
        weights_masked,
        filter_radius_um,
        DESIGN_REGION_UM.x,
        DESIGN_REGION_UM.y,
        DESIGN_REGION_RESOLUTION_UM,
    )

    if sigmoid_bias == 0:
        return weights_filtered.flatten()
    else:
        weights_projected = mpa.tanh_projection(
            weights_filtered,
            sigmoid_bias,
            sigmoid_threshold,
        )
        return weights_projected.flatten()


def obj_func(epigraph_and_weights: np.ndarray, grad: np.ndarray) -> float:
    """Objective function for the epigraph formulation.

    Args:
      epigraph_and_weights: 1D array containing epigraph variable (first
        element) and design weights (remaining elements).
      grad: the gradient as a flattened (1D) array, modified in place.

    Returns:
      The scalar epigraph variable.
    """
    epigraph = epigraph_and_weights[0]

    if grad.size > 0:
        grad[0] = 1
        grad[1:] = 0

    return epigraph


def epigraph_constraint(
    result: np.ndarray,
    epigraph_and_weights: np.ndarray,
    gradient: np.ndarray,
    sigmoid_threshold: float,
    sigmoid_bias: float,
    use_epsavg: bool,
) -> None:
    """Constraint function for the epigraph formulation.

    Args:
      result: the result of the function evaluation, modified in place.
      epigraph_and_weights: 1D array containing the epigraph variable (first
        element) and design weights (remaining elements).
      gradient: the Jacobian matrix with dimensions (1 + NX_DESIGN_GRID *
        NY_DESIGN_GRID, 2 * num. wavelengths), modified in place.
      sigmoid_threshold: erosion/dilation parameter for projection.
      sigmoid_bias: bias parameter for projection.
      use_epsavg: whether to use subpixel smoothing.
    """
    epigraph = epigraph_and_weights[0]
    weights = epigraph_and_weights[1:]

    obj_val, grad = opt(
        [
            filter_projection(
                weights, sigmoid_threshold, 0 if use_epsavg else sigmoid_bias
            )
        ]
    )

    reflectance = obj_val[0]
    transmittance = obj_val[1]
    obj_val_merged = np.concatenate((reflectance, transmittance))
    obj_val_merged_str = str_from_list(obj_val_merged)

    grad_reflectance = grad[0]
    grad_transmittance = grad[1]
    nfrq = len(frequencies)
    grad = np.zeros((NX_DESIGN_GRID * NY_DESIGN_GRID, 2 * nfrq))
    grad[:, :nfrq] = grad_reflectance
    grad[:, nfrq:] = grad_transmittance

    # backpropagate the gradients through filter and projection function
    for k in range(2 * nfrq):
        grad[:, k] = tensor_jacobian_product(filter_projection, 0)(
            weights,
            sigmoid_threshold,
            sigmoid_bias,
            grad[:, k],
        )

    if gradient.size > 0:
        gradient[:, 0] = -1  # gradient w.r.t. epigraph variable
        gradient[:, 1:] = grad.T  # gradient w.r.t. each frequency objective

    result[:] = np.real(obj_val_merged) - epigraph

    objfunc_history.append(np.real(obj_val_merged))
    epivar_history.append(epigraph)

    print(
        f"iteration:, {cur_iter[0]:3d}, sigmoid_bias: {sigmoid_bias:2d}, "
        f"epigraph: {epigraph:.5f}, obj. func.: {obj_val_merged_str}"
    )

    cur_iter[0] = cur_iter[0] + 1


def line_width_and_spacing_constraint(
    result: np.ndarray,
    epigraph_and_weights: np.ndarray,
    gradient: np.ndarray,
    sigmoid_bias: float,
) -> float:
    """Constraint function for the minimum linewidth.

    Args:
      result: the result of the function evaluation modified in place.
      epigraph_and_weights: 1D array containing epigraph variable (first
        element) and design weights (remaining elements).
      gradient: the Jacobian matrix, modified in place.
      sigmoid_bias: bias parameter for projection.

    Returns:
      The value of the constraint function (a scalar).
    """
    epigraph = epigraph_and_weights[0]
    weights = epigraph_and_weights[1:]
    a1 = 1e-3  # hyper parameter (primary)
    b1 = 0  # hyper parameter (secondary)
    gradient[:, 0] = -a1

    filter_func = lambda a: mpa.conic_filter(
        a.reshape(NX_DESIGN_GRID, NY_DESIGN_GRID),
        filter_radius_um,
        DESIGN_REGION_UM.x,
        DESIGN_REGION_UM.y,
        DESIGN_REGION_RESOLUTION_UM,
    )

    threshold_func = lambda a: mpa.tanh_projection(
        a, sigmoid_bias, SIGMOID_THRESHOLD_INTRINSIC
    )

    # hyper parameter (constant factor and exponent)
    c0 = 1e7 * (filter_radius_um * 1 / RESOLUTION_UM) ** 4

    M1 = lambda a: mpa.constraint_solid(
        a, c0, SIGMOID_THRESHOLD_EROSION, filter_func, threshold_func, 1
    )
    M2 = lambda a: mpa.constraint_void(
        a, c0, SIGMOID_THRESHOLD_DILATION, filter_func, threshold_func, 1
    )

    g1 = grad(M1)(weights)
    g2 = grad(M2)(weights)

    result[0] = M1(weights) - a1 * epigraph - b1
    result[1] = M2(weights) - a1 * epigraph - b1

    gradient[0, 1:] = g1.flatten()
    gradient[1, 1:] = g2.flatten()

    t1 = (M1(weights) - b1) / a1
    t2 = (M2(weights) - b1) / a1

    print(f"line_width_and_spacing_constraint:, {result[0]}, {result[1]}, {t1}, {t2}")

    return max(t1, t2)


def straight_waveguide() -> (np.ndarray, NamedTuple):
    """Computes the DFT fields from the mode source in a straight waveguide.

    The DFT fields are used as normalization of the reflectance measurement
    during the optimization.

    Returns:
      A 2-tuple consisting of a 1d array of DFT fields and DFT fields object
      returned by `meep.get_flux_data`.
    """
    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(frequency_center, fwidth=frequency_width),
            size=mp.Vector3(0, cell_um.y, 0),
            center=src_pt,
            eig_band=1,
            eig_parity=MODE_SYMMETRY,
        )
    ]

    geometry = [
        mp.Block(
            size=mp.Vector3(mp.inf, WAVEGUIDE_UM.y, mp.inf),
            center=mp.Vector3(),
            material=SILICON,
        )
    ]

    sim = mp.Simulation(
        resolution=RESOLUTION_UM,
        default_material=SILICON_DIOXIDE,
        cell_size=cell_um,
        sources=sources,
        geometry=geometry,
        boundary_layers=pml_layers,
        k_point=mp.Vector3(),
    )

    refl_mon = sim.add_mode_monitor(
        frequencies,
        mp.ModeRegion(center=refl_pt, size=mp.Vector3(0, cell_um.y, 0)),
        yee_grid=True,
    )

    sim.run(until_after_sources=stop_cond)

    res = sim.get_eigenmode_coefficients(
        refl_mon,
        [1],
        eig_parity=MODE_SYMMETRY,
    )

    coeffs = res.alpha
    input_flux = np.abs(coeffs[0, :, 0]) ** 2
    input_flux_data = sim.get_flux_data(refl_mon)

    return input_flux, input_flux_data


def mode_converter_optimization(
    input_flux: np.ndarray,
    input_flux_data: NamedTuple,
    use_damping: bool,
    use_epsavg: bool,
    sigmoid_bias: float,
) -> mpa.OptimizationProblem:
    """Sets up the adjoint optimization of the waveguide mode converter.

    Args:
      input_flux: 1D array of DFT fields from the normalization run.
      input_flux_data: DFT fields object returned by `meep.get_flux_data`.
      use_damping: whether to use the damping feature of `MaterialGrid`.
      use_epsavg: whether to use subpixel smoothing in `MaterialGrid`.

    Returns:
      A `meep.adjoint.OptimizationProblem` class object.
    """
    matgrid = mp.MaterialGrid(
        mp.Vector3(NX_DESIGN_GRID, NY_DESIGN_GRID, 0),
        SILICON_DIOXIDE,
        SILICON,
        weights=np.ones((NX_DESIGN_GRID, NY_DESIGN_GRID)),
        beta=sigmoid_bias if use_epsavg else 0,
        do_averaging=True if use_epsavg else False,
        damping=0.02 * 2 * np.pi * frequency_center if use_damping else 0,
    )

    matgrid_region = mpa.DesignRegion(
        matgrid,
        volume=mp.Volume(center=mp.Vector3(), size=DESIGN_REGION_UM),
    )

    matgrid_geometry = [
        mp.Block(
            center=matgrid_region.center,
            size=matgrid_region.size,
            material=matgrid,
        )
    ]

    geometry = [
        mp.Block(
            center=mp.Vector3(),
            size=mp.Vector3(mp.inf, WAVEGUIDE_UM.y, mp.inf),
            material=SILICON,
        )
    ]

    geometry += matgrid_geometry

    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(frequency_center, fwidth=frequency_width),
            size=mp.Vector3(0, cell_um.y, 0),
            center=src_pt,
            eig_band=1,
            eig_parity=MODE_SYMMETRY,
        ),
    ]

    sim = mp.Simulation(
        resolution=RESOLUTION_UM,
        default_material=SILICON_DIOXIDE,
        cell_size=cell_um,
        sources=sources,
        geometry=geometry,
        boundary_layers=pml_layers,
        k_point=mp.Vector3(),
    )

    obj_list = [
        mpa.EigenmodeCoefficient(
            sim,
            mp.Volume(
                center=refl_pt,
                size=mp.Vector3(0, cell_um.y, 0),
            ),
            1,
            forward=False,
            eig_parity=MODE_SYMMETRY,
            subtracted_dft_fields=input_flux_data,
        ),
        mpa.EigenmodeCoefficient(
            sim,
            mp.Volume(
                center=tran_pt,
                size=mp.Vector3(0, cell_um.y, 0),
            ),
            2,
            eig_parity=MODE_SYMMETRY,
        ),
    ]

    def J1(refl_mon, tran_mon):
        return npa.power(npa.abs(refl_mon), 2) / input_flux

    def J2(refl_mon, tran_mon):
        return 1 - npa.power(npa.abs(tran_mon), 2) / input_flux

    opt = mpa.OptimizationProblem(
        simulation=sim,
        objective_functions=[J1, J2],
        objective_arguments=obj_list,
        design_regions=[matgrid_region],
        frequencies=frequencies,
    )

    return opt


if __name__ == "__main__":
    input_flux, input_flux_data = straight_waveguide()

    silicon_mask, silicon_dioxide_mask = border_masks()

    num_weights = NX_DESIGN_GRID * NY_DESIGN_GRID

    # Initial design weights.
    epigraph_and_weights = np.ones((num_weights,)) * 0.5
    epigraph_and_weights[silicon_mask.flatten()] = 1.0
    epigraph_and_weights[silicon_dioxide_mask.flatten()] = 0.0

    # Lower and upper bounds for design weights.
    weights_lower_bound = np.zeros((num_weights,))
    weights_lower_bound[silicon_mask.flatten()] = 1.0
    weights_upper_bound = np.ones((num_weights,))
    weights_upper_bound[silicon_dioxide_mask.flatten()] = 0.0

    # Insert epigraph variable with arbitrary initial value and bounds into the
    # design array. The actual value is determined by the objective and
    # constraint functions below.
    epigraph_and_weights = np.insert(epigraph_and_weights, 0, 1.2)
    weights_lower_bound = np.insert(weights_lower_bound, 0, -np.inf)
    weights_upper_bound = np.insert(weights_upper_bound, 0, +np.inf)

    objfunc_history = []
    epivar_history = []
    cur_iter = [0]

    sigmoid_bias_threshold = 64  # threshold beta above which to use subpixel smoothing
    sigmoid_biases = [8, 16, 32, 64, 128, 256]
    max_evals = [80, 80, 100, 120, 120, 100]
    tolerance_epigraph = np.array([1e-4] * 2 * len(frequencies))  # R, 1-T
    tolerance_width_and_spacing = np.array([1e-8] * 2)  # line width, line spacing

    for sigmoid_bias, max_eval in zip(sigmoid_biases, max_evals):
        solver = nlopt.opt(nlopt.LD_MMA, num_weights + 1)
        solver.set_lower_bounds(weights_lower_bound)
        solver.set_upper_bounds(weights_upper_bound)
        solver.set_min_objective(obj_func)
        solver.set_maxeval(max_eval)
        solver.set_param("dual_ftol_rel", 1e-7)
        solver.add_inequality_mconstraint(
            lambda result_, epigraph_and_weights_, grad_: epigraph_constraint(
                result_,
                epigraph_and_weights_,
                grad_,
                SIGMOID_THRESHOLD_INTRINSIC,
                sigmoid_bias,
                False if sigmoid_bias < sigmoid_bias_threshold else True,
            ),
            tolerance_epigraph,
        )
        solver.set_param("verbosity", 1)

        if sigmoid_bias < sigmoid_bias_threshold:
            use_epsavg = False
        else:
            use_epsavg = True

        opt = mode_converter_optimization(
            input_flux,
            input_flux_data,
            True,  # use_damping
            use_epsavg,
            sigmoid_bias,
        )

        # Apply the minimum linewidth constraint
        # only in the final epoch to an initial
        # binary design from the previous epoch.
        if sigmoid_bias == sigmoid_biases[-1]:
            line_width_and_spacing = np.zeros(2)
            grad_line_width_and_spacing = np.zeros((2, num_weights + 1))
            linewidth_constraint_val = line_width_and_spacing_constraint(
                line_width_and_spacing,
                epigraph_and_weights,
                grad_line_width_and_spacing,
                sigmoid_bias,
            )
            solver.add_inequality_mconstraint(
                lambda result_, epigraph_and_weights_, grad_: line_width_and_spacing_constraint(
                    result_,
                    epigraph_and_weights_,
                    grad_,
                    sigmoid_bias,
                ),
                tolerance_width_and_spacing,
            )

        # Execute a single forward run before the start of each
        # epoch and manually set the initial epigraph variable to
        # slightly larger than the largest value of the objective
        # function over the six wavelengths and the lengthscale
        # constraint (final epoch only).
        epigraph_initial = opt(
            [
                filter_projection(
                    epigraph_and_weights[1:],
                    SIGMOID_THRESHOLD_INTRINSIC,
                    sigmoid_bias if sigmoid_bias < sigmoid_bias_threshold else 0,
                ),
            ],
            need_gradient=False,
        )
        epigraph_initial = np.concatenate(
            (epigraph_initial[0][0], epigraph_initial[0][1])
        )
        epigraph_initial_str = str_from_list(epigraph_initial)

        epigraph_and_weights[0] = np.amax(epigraph_initial)
        if sigmoid_bias == sigmoid_biases[-1]:
            epigraph_and_weights[0] = 1.05 * max(
                epigraph_and_weights[0], linewidth_constraint_val
            )
        print(
            f"epigraph-calibration:, {sigmoid_bias}, {epigraph_initial_str}, "
            f"{epigraph_and_weights[0]}"
        )

        epigraph_and_weights[:] = solver.optimize(epigraph_and_weights)

        optimal_design_weights = filter_projection(
            epigraph_and_weights[1:],
            SIGMOID_THRESHOLD_INTRINSIC,
            sigmoid_bias,
        ).reshape(NX_DESIGN_GRID, NY_DESIGN_GRID)

        # Save the unmapped weights and a bitmap image
        # of the design weights at the end of each epoch.
        fig, ax = plt.subplots()
        ax.imshow(
            optimal_design_weights,
            cmap="binary",
            interpolation="none",
        )
        ax.set_axis_off()
        if mp.am_master():
            fig.savefig(
                f"optimal_design_beta{sigmoid_bias}.png",
                dpi=150,
                bbox_inches="tight",
            )
            # Save the final (unmapped) design as a 2D array in CSV format
            np.savetxt(
                f"unmapped_design_weights_beta{sigmoid_bias}.csv",
                epigraph_and_weights[1:].reshape(NX_DESIGN_GRID, NY_DESIGN_GRID),
                fmt="%4.2f",
                delimiter=",",
            )

    # Save important optimization parameters and output as
    # separate fields for post processing.
    with open("optimal_design.npz", "wb") as datafile:
        np.savez(
            datafile,
            RESOLUTION_UM=RESOLUTION_UM,
            DESIGN_WAVELENGTHS_UM=DESIGN_WAVELENGTHS_UM,
            WAVEGUIDE_UM=WAVEGUIDE_UM,
            PADDING_UM=PADDING_UM,
            PML_UM=PML_UM,
            DESIGN_REGION_UM=DESIGN_REGION_UM,
            DESIGN_REGION_RESOLUTION_UM=DESIGN_REGION_RESOLUTION_UM,
            NX_DESIGN_GRID=NX_DESIGN_GRID,
            NY_DESIGN_GRID=NY_DESIGN_GRID,
            MIN_LENGTH_UM=MIN_LENGTH_UM,
            sigmoid_biases=sigmoid_biases,
            max_eval=max_eval,
            objfunc_history=objfunc_history,
            epivar_history=epivar_history,
            epigraph_variable=epigraph_and_weights[0],
            unmapped_design_weights=epigraph_and_weights[1:],
            optimal_design_weights=optimal_design_weights,
        )
