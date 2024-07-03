"""Shape optimization of a multilayer stack.

The 1D stack consists of alternating materials of index N_LAYER_1 and N_LAYER_2
in the arrangement: N_LAYER_1, N_LAYER_2, N_LAYER_1, N_LAYER_2, ..., N_LAYER_1.

The design parameters are the N layer thicknesses: [t_1, t_2, ..., t_N].
"""

from typing import Callable

from autograd.extend import primitive, defvjp
from autograd import numpy as npa
from autograd import tensor_jacobian_product
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import meep as mp
import meep.adjoint as mpa
import numpy as np


RESOLUTION_UM = 200
WAVELENGTH_UM = 1.0
AIR_UM = 1.0
PML_UM = 1.0
NUM_LAYERS = 11
MIN_LAYER_UM = 0.1
MAX_LAYER_UM = 0.9
N_LAYER_1 = 3.5
N_LAYER_2 = 1.0  # TODO (oskooi): support arbitrary use case (not just air).
LAYER_PERTURBATION_UM = 1.0 / RESOLUTION_UM
DESIGN_REGION_RESOLUTION_UM = 10 * RESOLUTION_UM

DEBUG_OUTPUT = True

design_region_size = mp.Vector3(0, 0, NUM_LAYERS * MAX_LAYER_UM)
nz_design_grid = int(design_region_size.z * DESIGN_REGION_RESOLUTION_UM) + 1
nz_sim_grid = int(design_region_size.z * RESOLUTION_UM) + 1
num_unit_cells = int(0.5 * (NUM_LAYERS - 1))


def design_region_to_grid(nz: int) -> np.ndarray:
    """Returns the 1D coordinates of the grid for the design region.

    Args:
      nz: number of grid points.

    Returns:
      The coordinates of the grid points.
    """

    z_grid = np.linspace(-0.5 * design_region_size.z, +0.5 * design_region_size.z, nz)

    return z_grid


@primitive
def levelset_and_smoothing(layer_thickness_um: np.ndarray) -> np.ndarray:
    """Returns the density weights for a multilayer stack as a levelset.

    Args:
      layer_thickness_um: thickness of each layer in the stack.

    Returns:
      The density weights as a flattened (1D) array.
    """

    air_padding_um = 0.5 * (design_region_size.z - np.sum(layer_thickness_um))

    weights = np.zeros(nz_design_grid)
    z_grid = design_region_to_grid(nz_design_grid)

    # Air padding at left edge.
    z_start = 0
    z_end = int(air_padding_um * DESIGN_REGION_RESOLUTION_UM)
    weights[z_start:z_end] = 0

    z_start = z_end
    for j in range(NUM_LAYERS):
        z_end = z_start + int(layer_thickness_um[j] * DESIGN_REGION_RESOLUTION_UM)
        weights[z_start:z_end] = 1 if (j % 2 == 0) else 0
        z_start = z_end

    # Air padding at right edge.
    weights[z_start:] = 0

    # Smooth the design weights by downsampling from the design grid
    # to the simulation grid using bilinear interpolation.
    z_sim_grid = design_region_to_grid(nz_sim_grid)
    smoothed_weights = np.interp(z_sim_grid, z_grid, weights)

    if DEBUG_OUTPUT:
        fig, ax = plt.subplots()
        ax.plot(z_sim_grid, smoothed_weights, "b-")
        ax.plot(z_grid, weights, "r-")
        ax.set_xlabel(r"z ($\mu$m)")
        ax.set_ylabel("design weights (smoothed)")
        if mp.am_master():
            fig.savefig("multilayer_stack_levelset.png", dpi=150, bbox_inches="tight")

    return smoothed_weights.flatten()


def levelset_and_smoothing_vjp(
    ans: np.ndarray, layer_thickness_um: np.ndarray
) -> Callable[[np.ndarray], np.ndarray]:
    """Returns a function for computing the vector-Jacobian product."""

    air_padding_um = 0.5 * (design_region_size.z - np.sum(layer_thickness_um))

    jacobian = np.zeros((nz_sim_grid, NUM_LAYERS))

    z_grid = design_region_to_grid(nz_design_grid)
    z_sim_grid = design_region_to_grid(nz_sim_grid)

    for i in range(NUM_LAYERS):
        weights = np.zeros(nz_design_grid)

        # Air padding at left edge.
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

        # Air padding at right edge.
        weights[z_start:] = 0

        # Smooth the design weights by downsampling from the design grid
        # to the simulation grid using bilinear interpolation.
        smoothed_weights = np.interp(z_sim_grid, z_grid, weights)

        jacobian[:, i] = (smoothed_weights - ans) / LAYER_PERTURBATION_UM

    if DEBUG_OUTPUT:
        fig, ax = plt.subplots()
        im = ax.imshow(
            np.transpose(jacobian),
            cmap="inferno",
            interpolation="none",
            aspect="equal",
        )
        ax.set_title(r"$\partial \rho_{smooth-levelset} / \partial t$")
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)
        if mp.am_master():
            fig.savefig("multilayer_stack_gradient.png", dpi=150, bbox_inches="tight")

    return lambda g: np.tensordot(g, jacobian, axes=1)


def input_flux() -> float:
    """Returns the flux generated by the source."""

    frequency = 1 / WAVELENGTH_UM
    pml_layers = [mp.PML(direction=mp.Z, thickness=PML_UM)]

    size_z_um = PML_UM + AIR_UM + design_region_size.z + AIR_UM + PML_UM
    cell_size = mp.Vector3(0, 0, size_z_um)

    src_cmpt = mp.Ex
    src_pt = mp.Vector3(0, 0, -0.5 * size_z_um + PML_UM)
    sources = [
        mp.Source(
            mp.GaussianSource(frequency, fwidth=0.1 * frequency),
            component=src_cmpt,
            center=src_pt,
        )
    ]

    sim = mp.Simulation(
        resolution=RESOLUTION_UM,
        cell_size=cell_size,
        dimensions=1,
        boundary_layers=pml_layers,
        sources=sources,
    )

    tran_pt = mp.Vector3(0, 0, 0.5 * size_z_um - PML_UM)

    flux_mon = sim.add_flux(
        frequency, 0, 1, mp.FluxRegion(center=tran_pt, size=mp.Vector3())
    )

    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(25.0, src_cmpt, tran_pt, 1e-6)
    )

    flux = mp.get_fluxes(flux_mon)[0]

    return flux


def multilayer_stack(norm_flux: float) -> mpa.OptimizationProblem:
    """Sets up the adjoint optimization of a multilayer stack.

    Args:
      norm_flux: input flux to use for normalizing the objective function.

    Returns:
      A meep.adjoint.Optimization object for the simulation.
    """

    frequency = 1 / WAVELENGTH_UM
    pml_layers = [mp.PML(direction=mp.Z, thickness=PML_UM)]

    size_z_um = PML_UM + AIR_UM + design_region_size.z + AIR_UM + PML_UM
    cell_size = mp.Vector3(0, 0, size_z_um)

    src_cmpt = mp.Ex
    src_pt = mp.Vector3(0, 0, -0.5 * size_z_um + PML_UM)
    sources = [
        mp.Source(
            mp.GaussianSource(frequency, fwidth=0.1 * frequency),
            component=src_cmpt,
            center=src_pt,
        )
    ]

    mat_1 = mp.Medium(index=N_LAYER_1)
    mat_2 = mp.Medium(index=N_LAYER_2)

    matgrid = mp.MaterialGrid(
        mp.Vector3(0, 0, nz_sim_grid),
        mat_1,
        mat_2,
        weights=np.ones(nz_sim_grid),
        do_averaging=False,
    )

    matgrid_region = mpa.DesignRegion(
        matgrid,
        volume=mp.Volume(center=mp.Vector3(), size=design_region_size),
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

    tran_pt = mp.Vector3(0, 0, 0.5 * size_z_um - PML_UM)

    obj_args = [
        mpa.FourierFields(sim, mp.Volume(center=tran_pt), mp.Ex),
        mpa.FourierFields(sim, mp.Volume(center=tran_pt), mp.Hy),
    ]

    def obj_func(dft_ex, dft_hy):
        return npa.real(npa.conj(dft_ex) * dft_hy) / norm_flux

    opt = mpa.OptimizationProblem(
        simulation=sim,
        objective_functions=obj_func,
        objective_arguments=obj_args,
        design_regions=[matgrid_region],
        frequencies=[frequency],
    )

    return opt


if __name__ == "__main__":
    layer_thickness_um = np.array(
        [0.3, 0.1, 0.2, 0.1, 0.6, 0.4, 0.2, 0.1, 0.4, 0.3, 0.2]
    )

    smoothed_design_weights = levelset_and_smoothing(layer_thickness_um)

    norm_flux = input_flux()

    opt = multilayer_stack(norm_flux)

    obj_val_unperturbed, grad_unperturbed = opt(
        [smoothed_design_weights], need_gradient=True
    )
    print(f"obj_val_unperturbed[0] = {obj_val_unperturbed[0]}")

    defvjp(levelset_and_smoothing, levelset_and_smoothing_vjp)
    grad_backprop = tensor_jacobian_product(levelset_and_smoothing, 0)(
        layer_thickness_um, grad_unperturbed
    )

    if DEBUG_OUTPUT:
        fig, ax = plt.subplots()
        ax.plot(grad_unperturbed)
        ax.set_title(r"$\partial F / \partial \rho_{smooth-levelset}$")
        if mp.am_master():
            fig.savefig(
                "gradient_wrt_smoothed_design_weights.png", dpi=150, bbox_inches="tight"
            )

    perturbation_um = 1e-3 * np.ones(NUM_LAYERS)
    perturbed_design_weights = levelset_and_smoothing(
        layer_thickness_um + perturbation_um
    )
    perturbed_design_weights = perturbed_design_weights.flatten()

    obj_val_perturbed, _ = opt([perturbed_design_weights], need_gradient=False)
    print(f"obj_val_perturbed[0] = {obj_val_perturbed[0]}")

    adj_directional_deriv = (perturbation_um[None, :] @ grad_backprop)[0]
    fnd_directional_deriv = obj_val_perturbed[0] - obj_val_unperturbed[0]
    print(f"adj_directional_deriv:, {adj_directional_deriv}")
    print(f"fnd_directional_deriv:, {fnd_directional_deriv}")
    rel_err = abs(
        (fnd_directional_deriv - adj_directional_deriv) / fnd_directional_deriv
    )
    print(
        f"dir-deriv:, {fnd_directional_deriv:.8f} (finite difference), "
        f"{adj_directional_deriv:.8f} (adjoint), {rel_err:.6f} (error)"
    )
