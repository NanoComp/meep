"""waveguide_crossing.py - Using meep's adjoint solver, designs a waveguide
crossing that maximizes transmission at a single frequency. 

Two approaches are demonstrated: (1) start with the nominal crossing shape and
perform shape optimization via meep's smoothed projection feature; (2) run a
full topology optimization starting from an initial grayscale. Evolve the design
until β=∞.

Importantly, this particular example highlights some of the ways one can use the
novel smoothed projection function to perform both shape and topology
optimization.
"""

from typing import Callable, List, Optional, Tuple

import meep.adjoint as mpa
import nlopt
import numpy as np
from autograd import grad
from autograd import numpy as npa
from autograd import tensor_jacobian_product
from matplotlib import pyplot as plt

import meep as mp

mp.quiet()

DEFAULT_MIN_LENGTH = 0.09
DEFAULT_DESIGN_REGION_WIDTH = 3.0
DEFAULT_DESIGN_REGION_HEIGHT = 3.0
DEFAULT_WAVEGUIDE_WIDTH = 0.5
DEFAULT_ETA = 0.5
DEFAULT_ETA_E = 0.75
DEFAULT_MAX_EVAL = 30


def build_optimization_problem(
    resolution: float,
    beta: float,
    use_smoothed_projection: bool,
    min_length: float = DEFAULT_MIN_LENGTH,
    dx: float = DEFAULT_DESIGN_REGION_WIDTH,
    dy: float = DEFAULT_DESIGN_REGION_HEIGHT,
    waveguide_width: float = DEFAULT_WAVEGUIDE_WIDTH,
    eta: float = DEFAULT_ETA,
    eta_e: float = DEFAULT_ETA_E,
    damping_factor: float = 0.0,
) -> Tuple[mpa.OptimizationProblem, Callable]:
    """Build the waveguide-crossing optimization problem.

    The waveguide crossing is a cononical inverse design problem with both
    shape- and topology-optimization implementations. The idea is to find the
    optimal structure that maximizes transmission from one side to the other. It
    exhibits C4 symmetry, and generally resembles the following structure:

         |  |
         |  |
    -----    ------
    -----    ------
         |  |
         |  |

    Args:
        resolution: Simulation resolution in pixels/micron.
        beta: Tanh function projection strength parameter, ranging from [0,∞].
        use_smoothed_projection: Whether or not to use the smoothed projection.
        min_length: Minimum length scale in microns.
        dx: Design region width in microns.
        dy: Design region height in microns.
        waveguide_width: Waveguide width in microns.
        eta: Projection function threshold parameter.
        eta_e: Projection function eroded threshold parameter.
        damping_factor: The material grid damping scalar factor.

    Returns:
        The corresponding optimization problem object and the mapping function
        that applies the linear and nonlinear transformations.
    """
    # Map the design region resolution to the yee grid, which is twice the standard resolution.
    design_region_resolution = int(2 * resolution)

    # pml thickness
    dpml = 1.0

    filter_radius = mpa.get_conic_radius_from_eta_e(min_length, eta_e)

    sxy = dx + 1 + 2 * dpml

    silicon = mp.Medium(epsilon=12)
    cell_size = mp.Vector3(sxy, sxy, 0)
    boundary_layers = [mp.PML(thickness=dpml)]

    eig_parity = mp.EVEN_Y + mp.ODD_Z

    design_region_size = mp.Vector3(dx, dy)
    Nx = int(design_region_resolution * design_region_size.x) + 1
    Ny = int(design_region_resolution * design_region_size.y) + 1

    waveguide_geometry = [
        mp.Block(material=silicon, size=mp.Vector3(mp.inf, waveguide_width, mp.inf)),
        mp.Block(material=silicon, size=mp.Vector3(waveguide_width, mp.inf, mp.inf)),
    ]

    # Source centered in optical c-band
    fcen = 1 / 1.55
    df = 0.23 * fcen
    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(-0.5 * sxy + dpml + 0.1, 0),
            size=mp.Vector3(0, sxy - 2 * dpml),
            eig_band=1,
            eig_parity=eig_parity,
        )
    ]

    damping = damping_factor * fcen
    matgrid = mp.MaterialGrid(
        mp.Vector3(Nx, Ny),
        mp.air,
        silicon,
        weights=np.ones((Nx, Ny)),
        beta=0,  # disable meep's internal smoothing
        do_averaging=False,  # disable meep's internal mg smoothing
        damping=damping,
    )

    matgrid_region = mpa.DesignRegion(
        matgrid,
        volume=mp.Volume(
            center=mp.Vector3(),
            size=mp.Vector3(design_region_size.x, design_region_size.y, 0),
        ),
    )

    matgrid_geometry = [
        mp.Block(
            center=matgrid_region.center, size=matgrid_region.size, material=matgrid
        )
    ]

    geometry = waveguide_geometry + matgrid_geometry

    sim = mp.Simulation(
        resolution=resolution,
        cell_size=cell_size,
        boundary_layers=boundary_layers,
        sources=sources,
        geometry=geometry,
    )

    frequencies = [fcen]

    obj_list = [
        mpa.EigenmodeCoefficient(
            sim,
            mp.Volume(
                center=mp.Vector3(-0.5 * sxy + dpml + 0.2),
                size=mp.Vector3(0, sxy - 2 * dpml, 0),
            ),
            1,
            eig_parity=eig_parity,
        ),
        mpa.EigenmodeCoefficient(
            sim,
            mp.Volume(
                center=mp.Vector3(0.5 * sxy - dpml - 0.2),
                size=mp.Vector3(0, sxy - 2 * dpml, 0),
            ),
            1,
            eig_parity=eig_parity,
        ),
    ]

    def J(input, output):
        """Simple objective function to minimize loss."""
        return 1 - npa.power(npa.abs(output / input), 2)

    opt = mpa.OptimizationProblem(
        simulation=sim,
        maximum_run_time=500,
        objective_functions=J,
        objective_arguments=obj_list,
        design_regions=[matgrid_region],
        frequencies=frequencies,
    )

    def mapping(x: npa.ndarray):
        """Applies the smoothing and projection."""
        x = x.reshape(Nx, Ny)
        x = mpa.conic_filter(
            x,
            filter_radius,
            design_region_size.x,
            design_region_size.y,
            design_region_resolution,
        )

        # Enforce the square symmetry
        x = (x + npa.rot90(x) + npa.rot90(x, 2) + npa.rot90(x, 3)) / 4

        # Only used the smoothed projection if prompted.
        if use_smoothed_projection:
            x = mpa.smoothed_projection(
                x, beta=beta, eta=eta, resolution=design_region_resolution
            )
        else:
            x = mpa.tanh_projection(x, beta=beta, eta=eta)

        return x.flatten()

    return opt, mapping


def nlopt_fom(
    x: np.ndarray,
    gradient: np.ndarray,
    opt: mpa.OptimizationProblem,
    mapping: Callable,
    data: List,
    results: List,
):
    """Wrapper for NLopt FOM.

    Args:
        x: Degrees of freedom array.
        gradient: Gradient of FOM.
        opt: Optimization problem object.
        mapping: Mapping function.
        data: Structure to store the simulated design each iteration.
        results: Structure to store the simulated FOM each iteration.

    Returns:
        The FOM value at the current iteration.
    """
    grid_size = opt.design_regions[0].design_parameters.grid_size
    Nx, Ny = int(grid_size.x), int(grid_size.y)

    f0, dJ_du = opt([mapping(x)])
    backprop_gradient = tensor_jacobian_product(mapping, 0)(x, dJ_du)
    if gradient.size > 0:
        gradient[:] = backprop_gradient

    data.append(
        [
            np.squeeze(x.copy().reshape(Nx, Ny)),
            np.squeeze(mapping(x).copy().reshape(Nx, Ny)),
            np.squeeze(backprop_gradient.reshape(Nx, Ny)),
        ]
    )

    print(
        f"FOM: {np.real(f0)} | x NaNs:{np.sum(np.isnan(backprop_gradient))} | grad NaNs: {np.sum(np.isnan(backprop_gradient))}"
    )
    results.append(np.real(f0))

    return float(np.real(f0))


def _plot_optimization_results(
    data: List,
    results: List,
    num_samples: int = 4,
):
    samples = np.linspace(0, len(results), num_samples, dtype=int, endpoint=False)
    plt.figure(figsize=(5.25, 4.0))
    plt.subplot(2, 1, 1)
    plt.plot(results, "o-")
    plt.xlabel("Optimization Iteration")
    plt.ylabel("FOM")
    for k in range(len(samples)):
        plt.subplot(2, 4, 5 + k)
        plt.imshow(data[samples[k]], cmap="binary", vmin=0.0, vmax=1.0)
        plt.axis("off")
        plt.title(f"It. {samples[k]+1}")
    plt.tight_layout()
    plt.show()


def run_shape_optimization(
    resolution: float,
    beta: float,
    maxeval: int,
    use_smoothed_projection: bool = True,
    damping_factor: float = 0.0,
    dx: float = DEFAULT_DESIGN_REGION_WIDTH,
    dy: float = DEFAULT_DESIGN_REGION_HEIGHT,
    waveguide_width: float = DEFAULT_WAVEGUIDE_WIDTH,
    output_filename_prefix: Optional[str] = None,
    plot_results: bool = True,
):
    """Run shape optimization using a cross as a starting guess.

    Args:
        resolution: Simulation resolution in pixels/micron.
        beta: Tanh function projection strength parameter, ranging from [0,∞].
        maxeval: Maximum number of optimization iterations to run.
        use_smoothed_projection: Whether or not to use the smoothed projection.
        dx: Design region width in microns.
        dy: Design region height in microns.
        waveguide_width: Waveguide width in microns.
        output_filename_prefix: The filename prefix that will store the
            optimization results. If `None`, no file is saved.
        plot_results: Whether or not to plot results.

    Returns:
        The design and FOM result arrays, along with the optimization problem
        and mapping function.

    """
    # Initialize the optimization problem and mapping functions
    opt, mapping = build_optimization_problem(
        resolution=resolution,
        beta=beta,
        dx=dx,
        dy=dy,
        waveguide_width=waveguide_width,
        use_smoothed_projection=use_smoothed_projection,
        damping_factor=damping_factor,
    )

    # pull number of parameters from the design region
    n = opt.design_regions[0].design_parameters.weights.size
    grid_size = opt.design_regions[0].design_parameters.grid_size
    Nx, Ny = int(grid_size.x), int(grid_size.y)

    # Set up the optimizer
    algorithm = nlopt.LD_CCSAQ
    solver = nlopt.opt(algorithm, n)
    solver.set_lower_bounds(0)
    solver.set_upper_bounds(1)
    solver.set_maxeval(maxeval)

    # initial guess, which is just a simple waveguide crossing
    x = np.linspace(-dx / 2, dy / 2, Nx)
    y = np.linspace(-dx / 2, dy / 2, Ny)
    X, Y = np.meshgrid(x, y)
    mask = (np.abs(Y) <= waveguide_width / 2.0) + (np.abs(X) <= waveguide_width / 2.0)
    x0 = np.zeros((Nx, Ny))
    x0[mask.T] = 1
    x0 = mapping(x0)

    # Create empty datastructures we can use to log the results
    data = []
    results = []

    # prepare the optimizer objective function. We need to wrap the above fom in
    # a format that nlopt expects.
    nlopt_fom_simple = lambda x, gradient: nlopt_fom(
        x, gradient, opt=opt, mapping=mapping, data=data, results=results
    )
    solver.set_min_objective(nlopt_fom_simple)

    # Run the optimization
    x_final = solver.optimize(x0)

    # Log the final results
    opt.update_design([x_final.copy(), mapping(x_final)])
    f0, final_grad = opt(need_gradient=False)
    final_backprop_gradient = tensor_jacobian_product(mapping, 0)(x_final, final_grad)
    results.append(np.real(f0))
    data.append(
        [
            np.squeeze(x_final.copy().reshape(Nx, Ny)),
            np.squeeze(mapping(x_final).copy().reshape(Nx, Ny)),
            np.squeeze(final_backprop_gradient.reshape(Nx, Ny)),
        ]
    )

    if plot_results:
        _plot_optimization_results(data, results)

    # Save to disk
    if mp.am_really_master() and (output_filename_prefix is not None):
        np.savez(output_filename_prefix + "_data.npz", data=data, results=results)

    return data, results, opt, mapping


def run_topology_optimization(
    resolution: float,
    beta_evolution: List[float],
    maxeval: int,
    use_smoothed_projection: bool = True,
    dx: float = DEFAULT_DESIGN_REGION_WIDTH,
    dy: float = DEFAULT_DESIGN_REGION_HEIGHT,
    waveguide_width: float = DEFAULT_WAVEGUIDE_WIDTH,
    output_filename_prefix: Optional[str] = None,
    damping_factor: float = 0.0,
):
    """Run shape optimization using a cross as a starting guess.

    Args:
        resolution: Simulation resolution in pixels/micron.
        beta_evolution: List of Tanh function projection strength parameter,
            ranging from [0,∞], for each optimization epoch.
        maxeval: Maximum number of optimization iterations to run.
        use_smoothed_projection: Whether or not to use the smoothed projection.
        dx: Design region width in microns.
        dy: Design region height in microns.
        waveguide_width: Waveguide width in microns.
        output_filename_prefix: The filename prefix that will store the
            optimization results. If `None`, no file is saved.

    Returns:
        The design and FOM result arrays, along with the optimization problem
        and mapping function.

    """
    # Initialize the optimization problem and mapping functions so that we can
    # set up the problem.
    opt, mapping = build_optimization_problem(
        resolution=resolution,
        beta=beta_evolution[0],
        dx=dx,
        dy=dy,
        waveguide_width=waveguide_width,
        use_smoothed_projection=use_smoothed_projection,
        damping_factor=damping_factor,
    )

    # pull number of parameters from the design region
    n = opt.design_regions[0].design_parameters.weights.size
    grid_size = opt.design_regions[0].design_parameters.grid_size
    Nx, Ny = int(grid_size.x), int(grid_size.y)

    # Set up the optimizer
    algorithm = nlopt.LD_CCSAQ
    solver = nlopt.opt(algorithm, n)
    solver.set_lower_bounds(0)
    solver.set_upper_bounds(1)
    solver.set_maxeval(maxeval)

    # initial guess, which is just a uniform gray region
    x0 = 0.5 * np.ones((Nx, Ny))
    x0 = mapping(x0)

    # Create empty datastructures we can use to log the results
    data = []
    results = []

    for beta in beta_evolution:
        # Re-initialize the optimization problem and mapping functions
        opt, mapping = build_optimization_problem(
            resolution=resolution,
            beta=beta,
            dx=dx,
            dy=dy,
            waveguide_width=waveguide_width,
            use_smoothed_projection=use_smoothed_projection,
            damping_factor=damping_factor,
        )

        # prepare the optimizer objective function. We need to wrap the above fom in
        # a format that nlopt expects.
        nlopt_fom_simple = lambda x, gradient: nlopt_fom(
            x, gradient, opt=opt, mapping=mapping, data=data, results=results
        )
        solver.set_min_objective(nlopt_fom_simple)

        # Run the optimization
        x0 = solver.optimize(x0)

    x_final = x0

    # Log the final results
    opt.update_design([mapping(x_final)])
    f0, final_grad = opt(need_gradient=False)
    final_backprop_gradient = tensor_jacobian_product(mapping, 0)(x_final, final_grad)
    results.append(np.real(f0))
    data.append(
        [
            np.squeeze(x_final.copy().reshape(Nx, Ny)),
            np.squeeze(mapping(x_final).copy().reshape(Nx, Ny)),
            np.squeeze(final_backprop_gradient.reshape(Nx, Ny)),
        ]
    )

    _plot_optimization_results(data, results)

    # Save to disk
    if mp.am_really_master() and (output_filename_prefix is not None):
        np.savez(output_filename_prefix + "_data.npz", data=data, results=results)

    return data, results, opt, mapping


def analyze_gradient_convergence(
    resolution: float,
    beta_range: List[float],
    dx: float = DEFAULT_DESIGN_REGION_WIDTH,
    dy: float = DEFAULT_DESIGN_REGION_HEIGHT,
    waveguide_width: float = DEFAULT_WAVEGUIDE_WIDTH,
) -> Tuple[List[float], List[float]]:
    """Analyze the norm of the gradient vs beta.

    This experiment plots the norm of the gradient vector as a function of beta
    for both the smoothed projection and standard projection functions. As seen
    here, the smoothed projection function has a well-defined.

    Args:
        resolution: Simulation resolution in pixels/micron.
        beta_range: List of projection threshold parameters to test.
        dx: Design region width in microns.
        dy: Design region height in microns.
        waveguide_width: Waveguide width in microns.

    Returns:
        The norm of the gradient array when using smoothing and when not using
        smoothing as a function of beta.
    """

    # Initialize the optimization problem and mapping functions
    opt, mapping = build_optimization_problem(
        resolution=resolution,
        beta=beta_range[0],
        dx=dx,
        dy=dy,
        waveguide_width=waveguide_width,
        use_smoothed_projection=True,
    )

    # pull number of parameters from the design region
    n = opt.design_regions[0].design_parameters.weights.size
    grid_size = opt.design_regions[0].design_parameters.grid_size
    Nx, Ny = int(grid_size.x), int(grid_size.y)

    # initial guess, which is just a simple waveguide crossing
    x = np.linspace(-dy / 2, dy / 2, Nx)
    y = np.linspace(-dx / 2, dy / 2, Ny)
    X, Y = np.meshgrid(x, y)
    mask = (np.abs(Y) <= waveguide_width / 2.0) + (np.abs(X) <= waveguide_width / 2.0)
    x0 = np.zeros((Nx, Ny))
    x0[mask.T] = 1
    x0 = mapping(x0)

    # Initialize storage structures
    with_smoothing = []
    without_smoothing = []

    for beta in beta_range:
        # Compute the norm with smoothing
        opt, mapping = build_optimization_problem(
            resolution=resolution,
            beta=beta,
            dx=dx,
            dy=dy,
            waveguide_width=waveguide_width,
            use_smoothed_projection=True,
        )
        _, dJ_du = opt([mapping(x0)])
        backprop_gradient = tensor_jacobian_product(mapping, 0)(x0, dJ_du)
        norm_gradient_smoothing = np.linalg.norm(backprop_gradient.flatten())
        with_smoothing.append(norm_gradient_smoothing)

        # Compute the norm without smoothing
        opt, mapping = build_optimization_problem(
            resolution=resolution,
            beta=beta,
            dx=dx,
            dy=dy,
            waveguide_width=waveguide_width,
            use_smoothed_projection=False,
        )
        _, dJ_du = opt([mapping(x0)])
        backprop_gradient = tensor_jacobian_product(mapping, 0)(x0, dJ_du)
        norm_gradient_no_smoothing = np.linalg.norm(backprop_gradient.flatten())
        without_smoothing.append(norm_gradient_no_smoothing)

    # plot results
    plt.figure(figsize=(5.2, 2.5), constrained_layout=True)
    plt.loglog(beta_range, without_smoothing, "o-", label="W/o smoothing")
    plt.loglog(beta_range, with_smoothing, "o-", label="W/ smoothing")
    plt.legend()
    plt.xlabel("β")
    plt.ylabel("|df/dx|")
    plt.show()

    return with_smoothing, without_smoothing


def analyze_FOM_convergence(
    resolution: float, beta: float, maxeval: int, damping_factor: float = 0.0
) -> Tuple[List[float], List[float]]:
    """Analyze the convergence of the new projection method.

    Because the smoothed projection has a well-defined gradient for all values
    of beta, it should exhibit better convergence properties than the standard
    projection step, which grows increasingly more ill-conditioned as
    beta->infinity. This experiment runs the same shape optimization problem for
    both projection functions and plots the results.

    Args:
        resolution: Simulation resolution in pixels/micron.
        beta: Tanh function projection strength parameter, ranging from [0,∞].
        maxeval: Maximum number of optimization iterations to run.

    Returns:
        The FOM evolution for both the smoothed and not smoothed case as a
        function of optimization iteration.

    """
    print("Running shape optimization WITHOUT smoothing...")
    _, results, _, _ = run_shape_optimization(
        beta=beta,
        resolution=resolution,
        maxeval=maxeval,
        use_smoothed_projection=False,
        plot_results=False,
        output_filename_prefix=f"without_smoothing_grad_{beta}",
        damping_factor=damping_factor,
    )
    print("Running shape optimization WITH smoothing...")
    (
        _,
        results_smoothed,
        _,
        _,
    ) = run_shape_optimization(
        beta=beta,
        resolution=resolution,
        maxeval=maxeval,
        use_smoothed_projection=True,
        plot_results=False,
        output_filename_prefix=f"with_smoothing_grad_{beta}",
        damping_factor=damping_factor,
    )

    plt.figure(figsize=(5.2, 2.0), constrained_layout=True)
    plt.loglog(results, "o-", label="W/o smoothing")
    plt.loglog(results_smoothed, "o-", label="W/ smoothing")
    plt.legend()
    plt.xlabel("Optimization iteration")
    plt.ylabel("FOM")
    plt.show()

    return results, results_smoothed


if __name__ == "__main__":

    analyze_gradient_convergence(
        beta_range=np.logspace(1, 3, base=10, num=10), resolution=20
    )
