"""Computes the adjoint gradient of a level set.

This is a 2D example for computing the gradient of the diffraction efficiency
of the first transmitted order of a 1D grating with respect to the grating
height. The adjoint gradient is validated using the brute-force finite
difference via the directional derivative. The grating structure is
represented as a level set.
"""

from enum import Enum
from typing import Callable, Tuple

from autograd.extend import primitive, defvjp
from autograd import numpy as npa
from autograd import tensor_jacobian_product
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import meep as mp
import meep.adjoint as mpa
import numpy as np
from scipy.interpolate import RegularGridInterpolator


RESOLUTION_UM = 100
WAVELENGTH_UM = 0.53
GRATING_PERIOD_UM = 1.28
GRATING_HEIGHT_UM = 0.52
GRATING_DUTY_CYCLE = 0.65
GRATING_HEIGHT_PERTURBATION_UM = 2.0 / RESOLUTION_UM
DIFFRACTION_ORDER = 1
SUBSTRATE_UM = 1.5
PML_UM = 1.0
AIR_UM = 1.0
N_GLASS = 1.5
DESIGN_REGION_RESOLUTION_UM = 10 * RESOLUTION_UM

Polarization = Enum("Polarization", "S P")

design_region_size = mp.Vector3(
    GRATING_HEIGHT_UM + GRATING_HEIGHT_PERTURBATION_UM, GRATING_PERIOD_UM, 0
)
nx_design_grid = int(design_region_size.x * DESIGN_REGION_RESOLUTION_UM) + 1
ny_design_grid = int(design_region_size.y * DESIGN_REGION_RESOLUTION_UM) + 1
nx_sim_grid = int(design_region_size.x * RESOLUTION_UM) + 1
ny_sim_grid = int(design_region_size.y * RESOLUTION_UM) + 1

DEBUG_OUTPUT = False


def design_region_to_meshgrid(nx: int, ny: int) -> Tuple[np.ndarray, np.ndarray]:
    """Returns the 2D coordinates of the meshgrid for the design region.

    Args:
      nx, ny: the number of grid points in the x and y directions, respectively.

    Returns:
      The coordinates of the x and y grid points (2D arrays) as a 2-tuple.
    """

    xcoord = np.linspace(-0.5 * design_region_size.x, +0.5 * design_region_size.x, nx)
    ycoord = np.linspace(-0.5 * design_region_size.y, +0.5 * design_region_size.y, ny)
    xv, yv = np.meshgrid(xcoord, ycoord, indexing="ij")

    return xv, yv


@primitive
def levelset_and_smoothing(grating_height_um: float) -> np.ndarray:
    """Returns the density weights for a binary grating as a levelset.

    Args:
      grating_height_um: the height of the grating.

    Returns:
      The density weights as a flattened (1D) array.
    """

    xv, yv = design_region_to_meshgrid(nx_design_grid, ny_design_grid)
    weights = np.where(
        np.abs(yv) <= 0.5 * GRATING_DUTY_CYCLE * GRATING_PERIOD_UM,
        np.where(xv <= xv[0][0] + grating_height_um, 1, 0),
        0,
    )
    xcoord = np.linspace(
        -0.5 * design_region_size.x, +0.5 * design_region_size.x, nx_design_grid
    )
    ycoord = np.linspace(
        -0.5 * design_region_size.y, +0.5 * design_region_size.y, ny_design_grid
    )
    interp = RegularGridInterpolator((xcoord, ycoord), weights)

    # Smooth the design weights by downsampling from the design grid
    # to the simulation grid using bilinear interpolation.
    sim_grid_xv, sim_grid_yv = design_region_to_meshgrid(nx_sim_grid, ny_sim_grid)
    smoothed_weights = interp((sim_grid_xv, sim_grid_yv))

    if DEBUG_OUTPUT:
        fig, ax = plt.subplots()
        im = ax.imshow(
            np.transpose(smoothed_weights),
            cmap="binary",
            interpolation="none",
            aspect="equal",
        )
        ax.set_title(r"$\rho_{smooth-levelset}(h)$")
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)
        if mp.am_master():
            fig.savefig("smoothed_levelset.png", dpi=150, bbox_inches="tight")

    return smoothed_weights.flatten()


def levelset_and_smoothing_vjp(
    ans: np.ndarray,
    grating_height_um: float,
) -> Callable[[np.ndarray], np.ndarray]:
    """Returns a function for computing the vector-Jacobian product.

       The Jacobian is computed manually using a finite difference.

    Args:
      ans: the design weights for the smoothed grating (no perturbation).
      grating_height_um: the height of the grating.

    Returns:
      An anonymous function for computing the vector-Jacobian product.
    """

    xv, yv = design_region_to_meshgrid(nx_design_grid, ny_design_grid)
    weights = np.where(
        np.abs(yv) <= 0.5 * GRATING_DUTY_CYCLE * GRATING_PERIOD_UM,
        np.where(
            xv <= xv[0][0] + grating_height_um + GRATING_HEIGHT_PERTURBATION_UM, 1, 0
        ),
        0,
    )
    xcoord = np.linspace(
        -0.5 * design_region_size.x, +0.5 * design_region_size.x, nx_design_grid
    )
    ycoord = np.linspace(
        -0.5 * design_region_size.y, +0.5 * design_region_size.y, ny_design_grid
    )
    interp = RegularGridInterpolator((xcoord, ycoord), weights)

    # Smooth the design weights by downsampling from the design grid
    # to the simulation grid using bilinear interpolation.
    xv_sim_grid, yv_sim_grid = design_region_to_meshgrid(nx_sim_grid, ny_sim_grid)
    smoothed_weights = interp((xv_sim_grid, yv_sim_grid))
    smoothed_weights = smoothed_weights.flatten()

    jacobian = (smoothed_weights - ans) / GRATING_HEIGHT_PERTURBATION_UM

    if DEBUG_OUTPUT:
        fig, ax = plt.subplots()
        im = ax.imshow(
            np.transpose(jacobian.reshape(nx_sim_grid, ny_sim_grid)),
            cmap="inferno",
            interpolation="none",
            aspect="equal",
        )
        ax.set_title(r"$\partial \rho_{smooth-levelset} / \partial h$")
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)
        if mp.am_master():
            fig.savefig("gradient_smooth_levelset.png", dpi=150, bbox_inches="tight")

    return lambda g: np.tensordot(g, jacobian, axes=1)


def grating_1d(pol: Polarization) -> mpa.OptimizationProblem:
    """Sets up the adjoint optimization of a 1D grating.

    Args:
      pol: the polarization state (S or P).

    Returns:
      A meep.adjoint.OptimizationProblem object for the simulation.
    """

    frequency = 1 / WAVELENGTH_UM
    pml_layers = [mp.PML(direction=mp.X, thickness=PML_UM)]
    glass = mp.Medium(index=N_GLASS)
    size_x_um = (
        PML_UM
        + SUBSTRATE_UM
        + GRATING_HEIGHT_UM
        + GRATING_HEIGHT_PERTURBATION_UM
        + AIR_UM
        + PML_UM
    )
    size_y_um = GRATING_PERIOD_UM
    cell_size = mp.Vector3(size_x_um, size_y_um, 0)
    k_point = mp.Vector3()

    if pol.name == "S":
        eig_parity = mp.ODD_Z
        src_cmpt = mp.Ez
    else:
        eig_parity = mp.EVEN_Z
        src_cmpt = mp.Hz

    src_pt = mp.Vector3(-0.5 * size_x_um + PML_UM, 0, 0)
    sources = [
        mp.Source(
            mp.GaussianSource(frequency, fwidth=0.1 * frequency),
            component=src_cmpt,
            center=src_pt,
            size=mp.Vector3(0, size_y_um, 0),
        )
    ]

    matgrid = mp.MaterialGrid(
        mp.Vector3(nx_sim_grid, ny_sim_grid),
        mp.air,
        glass,
        weights=np.ones((nx_sim_grid, ny_sim_grid)),
        do_averaging=False,
    )

    matgrid_region = mpa.DesignRegion(
        matgrid,
        volume=mp.Volume(
            center=mp.Vector3(
                (-0.5 * size_x_um + PML_UM + SUBSTRATE_UM + 0.5 * design_region_size.x),
                0,
                0,
            ),
            size=design_region_size,
        ),
    )

    geometry = [
        mp.Block(
            material=glass,
            size=mp.Vector3(PML_UM + SUBSTRATE_UM, mp.inf, mp.inf),
            center=mp.Vector3(-0.5 * size_x_um + 0.5 * (PML_UM + SUBSTRATE_UM), 0, 0),
        ),
        mp.Block(
            material=matgrid,
            size=matgrid_region.size,
            center=matgrid_region.center,
        ),
    ]

    sim = mp.Simulation(
        resolution=RESOLUTION_UM,
        cell_size=cell_size,
        boundary_layers=pml_layers,
        k_point=k_point,
        sources=sources,
        geometry=geometry,
    )

    tran_pt = mp.Vector3(0.5 * size_x_um - PML_UM, 0, 0)

    kdiff = mp.Vector3(
        (frequency**2 - (DIFFRACTION_ORDER / GRATING_PERIOD_UM) ** 2) ** 0.5,
        DIFFRACTION_ORDER / GRATING_PERIOD_UM,
        0,
    )

    obj_args = [
        mpa.EigenmodeCoefficient(
            sim,
            mp.Volume(
                center=tran_pt,
                size=mp.Vector3(0, size_y_um, 0),
            ),
            mode=1,
            kpoint_func=lambda *not_used: kdiff,
            eig_parity=eig_parity,
            eig_vol=mp.Volume(center=tran_pt, size=mp.Vector3(0, 1 / RESOLUTION_UM, 0)),
        ),
    ]

    def obj_func(mode_coeff):
        return npa.abs(mode_coeff) ** 2

    opt = mpa.OptimizationProblem(
        simulation=sim,
        objective_functions=obj_func,
        objective_arguments=obj_args,
        design_regions=[matgrid_region],
        frequencies=[frequency],
    )

    return opt


if __name__ == "__main__":
    # Only the S polarization is supported (for now).
    opt = grating_1d(Polarization.S)

    smoothed_design_weights = levelset_and_smoothing(GRATING_HEIGHT_UM)

    obj_val_unperturbed, grad_unperturbed = opt(
        [smoothed_design_weights], need_gradient=True
    )

    if DEBUG_OUTPUT:
        fig, ax = plt.subplots()
        im = ax.imshow(
            np.transpose(grad_unperturbed.reshape(nx_sim_grid, ny_sim_grid)),
            cmap="inferno",
            interpolation="none",
            aspect="equal",
        )
        ax.set_title(r"$\partial F / \partial \rho_{smooth-levelset}$")
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)
        if mp.am_master():
            fig.savefig(
                "gradient_wrt_smoothed_design_weights.png", dpi=150, bbox_inches="tight"
            )

        fig, ax = plt.subplots()
        opt.plot2D(init_opt=False, ax=ax)
        if mp.am_master():
            fig.savefig("sim_layout.png", dpi=150, bbox_inches="tight")

    defvjp(levelset_and_smoothing, levelset_and_smoothing_vjp)
    grad_backprop = tensor_jacobian_product(levelset_and_smoothing, 0)(
        GRATING_HEIGHT_UM,
        grad_unperturbed,
    )

    perturbed_design_weights = levelset_and_smoothing(
        GRATING_HEIGHT_UM + GRATING_HEIGHT_PERTURBATION_UM
    )
    perturbed_design_weights = perturbed_design_weights.flatten()

    obj_val_perturbed, _ = opt([perturbed_design_weights], need_gradient=False)

    adj_directional_deriv = GRATING_HEIGHT_PERTURBATION_UM * grad_backprop
    fnd_directional_deriv = obj_val_perturbed[0] - obj_val_unperturbed[0]
    rel_err = abs(
        (fnd_directional_deriv - adj_directional_deriv) / fnd_directional_deriv
    )
    print(
        f"dir-deriv:, {fnd_directional_deriv:.8f} (finite difference), "
        f"{adj_directional_deriv:.8f} (adjoint), {rel_err:.6f} (error)"
    )
