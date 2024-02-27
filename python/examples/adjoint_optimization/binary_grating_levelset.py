"""Validates the adjoint gradient of a level set.

This is a 2D test which uses subpixel smoothing of the MaterialGrid for
computing the gradient of the diffraction efficiency of the first transmitted
order of a 1D grating. The adjoint gradient is validated using the brute-force
finite difference via the directional derivative. The grating structure is
represented as a level set.
"""

from enum import Enum
from typing import Tuple

from autograd.extend import primitive, defvjp
from autograd import numpy as anp
from autograd import tensor_jacobian_product
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import meep as mp
import meep.adjoint as mpa
from scipy.interpolate import interp1d
import skfmm


RESOLUTION_UM = 200
WAVELENGTH_UM = 0.53
GRATING_PERIOD_UM = 1.2
PADDING_UM = 1.5

Polarization = Enum("Polarization", "S P")


def design_region_to_meshgrid(
    design_region_size: mp.Vector3, nx: int, ny: int
) -> Tuple[anp.ndarray, anp.ndarray]:
    """Returns the 2D coordinates of the meshgrid for the design region."""

    xcoord = anp.linspace(-0.5 * design_region_size.x, +0.5 * design_region_size.x, nx)
    ycoord = anp.linspace(-0.5 * design_region_size.y, +0.5 * design_region_size.y, ny)
    xv, yv = anp.meshgrid(xcoord, ycoord, indexing="ij")

    return xv, yv


def signed_distance(weights: anp.ndarray) -> anp.ndarray:
    """Maps the 2D weights using a signed-distance function."""

    # Create signed distance function.
    sd = skfmm.distance(weights - 0.5)

    # Linear interpolation of zero-levelset onto 0.5-levelset.
    xp = [anp.min(sd.flatten()), 0, anp.max(sd.flatten())]
    yp = [0, 0.5001, 1]

    return anp.interp(sd.flatten(), xp, yp)


@primitive
def levelset_and_smoothing(
    grating_height_um: float,
    grating_duty_cycle: float,
    design_region_size: mp.Vector3,
    nx: int,
    ny: int,
) -> anp.ndarray:
    """Generates the grating as a levelset with signed-distance smoothing."""

    xv, yv = design_region_to_meshgrid(design_region_size, nx, ny)

    weights = anp.where(
        anp.abs(yv) <= 0.5 * grating_duty_cycle * GRATING_PERIOD_UM,
        anp.where(xv <= xv[0][0] + grating_height_um, 1, 0),
        0,
    )

    smoothed_weights = signed_distance(weights)

    return smoothed_weights


def levelset_and_smoothing_vjp(
    ans: anp.ndarray,
    grating_height_um: float,
    grating_duty_cycle: float,
    design_region_size: mp.Vector3,
    nx: int,
    ny: int,
) -> anp.ndarray:
    """Computes the vector-Jacobian product."""

    xv, yv = design_region_to_meshgrid(design_region_size, nx, ny)

    # pixel dimensions
    delta_x = design_region_size.x / (nx - 1)

    weights = anp.where(
        anp.abs(yv) <= 0.5 * grating_duty_cycle * GRATING_PERIOD_UM,
        anp.where(xv <= xv[0][0] + grating_height_um + delta_x, 1, 0),
        0,
    )

    jacobian = (signed_distance(weights) - ans) / delta_x
    print(f"jacobian:, {anp.linalg.norm(jacobian)}")

    extent = [xv[0][0], xv[0][-1], yv[0][0], yv[0][-1]]
    fig, ax = plt.subplots()
    im = ax.imshow(
        anp.transpose(jacobian.reshape(nx, ny)),
        extent=extent,
        cmap="hot_r",
        interpolation="none",
        aspect="equal",
    )
    ax.set_title("gradient")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im, cax=cax)
    if mp.am_master():
        fig.savefig("grating_gradient.png", dpi=150, bbox_inches="tight")

    return lambda g: anp.tensordot(g, jacobian, axes=1)


def grating_1d(
    pol: Polarization,
    grating_height_um: float,
    design_region_size: mp.Vector3,
    nx: int,
    ny: int,
) -> mpa.OptimizationProblem:
    """Sets up the adjoint optimization of a 1D grating."""

    frequency = 1 / WAVELENGTH_UM
    substrate_um = 2.0
    pml_um = 1.0
    pml_layers = [mp.PML(direction=mp.X, thickness=pml_um)]

    n_glass = 1.5
    glass = mp.Medium(index=n_glass)

    size_x_um = pml_um + substrate_um + grating_height_um + PADDING_UM + pml_um
    size_y_um = GRATING_PERIOD_UM
    cell_size = mp.Vector3(size_x_um, size_y_um, 0)

    k_point = mp.Vector3()

    if pol.name == "S":
        eig_parity = mp.ODD_Z
        src_cmpt = mp.Ez
    else:
        eig_parity = mp.EVEN_Z
        src_cmpt = mp.Hz

    src_pt = mp.Vector3(-0.5 * size_x_um + pml_um, 0, 0)
    sources = [
        mp.Source(
            mp.GaussianSource(frequency, fwidth=0.1 * frequency),
            component=src_cmpt,
            center=src_pt,
            size=mp.Vector3(0, size_y_um, 0),
        )
    ]

    matgrid = mp.MaterialGrid(
        mp.Vector3(nx, ny),
        mp.air,
        glass,
        weights=anp.ones((nx, ny)),
        do_averaging=True,
        beta=anp.inf,
    )

    matgrid_region = mpa.DesignRegion(
        matgrid,
        volume=mp.Volume(
            center=mp.Vector3(
                (-0.5 * size_x_um + pml_um + substrate_um + 0.5 * design_region_size.x),
                0,
                0,
            ),
            size=design_region_size,
        ),
    )

    geometry = [
        mp.Block(
            material=glass,
            size=mp.Vector3(pml_um + substrate_um, mp.inf, mp.inf),
            center=mp.Vector3(-0.5 * size_x_um + 0.5 * (pml_um + substrate_um), 0, 0),
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

    tran_pt = mp.Vector3(0.5 * size_x_um - pml_um, 0, 0)

    order_y = 1
    kdiff = mp.Vector3(
        (frequency**2 - (order_y / GRATING_PERIOD_UM) ** 2) ** 0.5,
        order_y / GRATING_PERIOD_UM,
        0,
    )
    print(f"kdiff = ({kdiff.x:.5f}, {kdiff.y:.5f}, {kdiff.z:.5f})")

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

    def obj_func(tran_mon):
        return anp.abs(tran_mon) ** 2

    opt = mpa.OptimizationProblem(
        simulation=sim,
        objective_functions=obj_func,
        objective_arguments=obj_args,
        design_regions=[matgrid_region],
        frequencies=[frequency],
    )

    return opt


if __name__ == "__main__":
    grating_height_um = 0.52
    grating_duty_cycle = 0.65
    height_perturbation_um = 2.0 / RESOLUTION_UM

    design_region_size = mp.Vector3(
        grating_height_um + PADDING_UM, GRATING_PERIOD_UM, 0
    )
    design_region_resolution_um = int(2 * RESOLUTION_UM)
    nx = int(design_region_size.x * design_region_resolution_um) + 1
    ny = int(design_region_size.y * design_region_resolution_um) + 1

    opt = grating_1d(
        Polarization.P,
        grating_height_um,
        design_region_size,
        nx,
        ny,
    )

    smoothed_design_weights = levelset_and_smoothing(
        grating_height_um,
        grating_duty_cycle,
        design_region_size,
        nx,
        ny,
    )

    obj_value_unperturbed, grad_unperturbed = opt(
        [smoothed_design_weights],
        need_gradient=True,
    )

    fig, ax = plt.subplots()
    opt.plot2D(init_opt=False, ax=ax)
    if mp.am_master():
        fig.savefig("grating_1D_plot2D.png", dpi=150, bbox_inches="tight")

    defvjp(levelset_and_smoothing, levelset_and_smoothing_vjp)
    grad_backprop = tensor_jacobian_product(levelset_and_smoothing, 0)(
        grating_height_um,
        grating_duty_cycle,
        design_region_size,
        nx,
        ny,
        grad_unperturbed,
    )

    perturbed_design_weights = levelset_and_smoothing(
        grating_height_um + height_perturbation_um,
        grating_duty_cycle,
        design_region_size,
        nx,
        ny,
    )
    perturbed_design_weights = perturbed_design_weights.flatten()

    obj_value_perturbed, _ = opt([perturbed_design_weights], need_gradient=False)

    adj_directional_deriv = height_perturbation_um * grad_backprop
    fnd_directional_deriv = obj_value_perturbed[0] - obj_value_unperturbed[0]
    rel_err = abs(
        (fnd_directional_deriv - adj_directional_deriv) / fnd_directional_deriv
    )
    print(
        f"directional-deriv:, {fnd_directional_deriv:.8f} (finite difference), "
        f"{adj_directional_deriv:.8f} (adjoint), {rel_err:.6f} (error)"
    )
