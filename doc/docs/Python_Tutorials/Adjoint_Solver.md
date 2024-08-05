---
# Adjoint Solver
---

Meep contains an adjoint-solver module for efficiently computing the
gradient of an arbitrary function of the mode coefficients
($S$-parameters), DFT fields, local density of states (LDOS), and
"far" fields with respect to $\varepsilon$ on a discrete spatial grid
(a [`MaterialGrid`](../Python_User_Interface.md#materialgrid) class
object) at multiple frequencies over a broad bandwidth. Regardless of
the number of degrees of freedom for the grid points, just **two**
distinct timestepping runs are required. The first run is the
"forward" calculation to compute the objective function and the DFT
fields of the design region. The second run is the "adjoint"
calculation to compute the gradient of the objective function with
respect to the design variables. The adjoint run involves a special
type of current source distribution used to compute the DFT fields of
the design region. The gradient is computed in post processing using
the DFT fields from the forward and adjoint runs. The gradient
calculation is fully automated. The theoretical and computational
details of the adjoint-solver module are described in this
publication:

- A. M. Hammond, A. Oskooi, M. Chen, Z. Lin, S. G. Johnson, and
  S. E. Ralph, “[High-performance hybrid time/frequency-domain
  topology optimization for large-scale photonics inverse
  design](https://doi.org/10.1364/OE.442074),” *Optics Express*,
  vol. 30, no. 3, pp. 4467–4491 (2022).

Much of the functionality of the adjoint solver is implemented in
Python using [autograd](https://github.com/HIPS/autograd) as well as
[JAX](https://github.com/google/jax).

The adjoint solver supports inverse design and [topology
optimization](https://en.wikipedia.org/wiki/Topology_optimization) by
providing the functionality to wrap an optimization library around the
gradient computation. The adjoint solver also supports enforcing
minimum feature sizes on the optimal design in 1D (line width and
spacing) and 2D (arbitrary-shaped holes and islands). This is
demonstrated in several tutorials below.

[TOC]

Broadband Waveguide Mode Converter with Minimum Feature Size
------------------------------------------------------------

This example demonstrates some of the advanced functionality of the
adjoint solver including worst-case (minimax) optimization across
multiple wavelengths, multiple objective functions, and design
constraints on the minimum line width and line spacing. The design
problem involves a broadband waveguide mode converter in 2D with
minimum feature size. This is based on [M.F. Schubert et al., ACS
Photonics, Vol. 9, pp. 2327-36,
(2022)](https://doi.org/10.1021/acsphotonics.2c00313).

The mode converter must satisfy two separate design objectives: (1)
minimize the reflectance $R$ into the fundamental mode of the input
port ($S_{11}$) and (2) maximize the transmittance $T$ into the
second-order mode of the output port ($S_{21}$). There are different
ways to define this multi-objective and multi-wavelength
optimization. The approach taken here is *worst-case* optimization
whereby we minimize the maximum (the worst case) of $R$ and $1-T$
across six wavelengths in the $O$-band for telecommunications. This is
known as *minimax* optimization. The fundamental mode launched from
the input port has $E_z$ polarization given a 2D cell in the $xy$
plane.

The challenge with minimax optimization is that the $\max$ objective
function is not everywhere differentiable. This property would seem to
preclude the use of gradient-based optimization algorithms for this
problem which involves twelve independent functions ($R$ and $1-T$ for
each of six wavelengths). Fortunately, there is a workaround: the
problem can be reformulated as a differentiable problem using a
so-called "epigraph" formulation: introducing a new "dummy"
optimization variable $t$ and adding each independent function
$f_k(x)$ as a new nonlinear constraint $t \ge f_k(x)$. See the [NLopt
documentation](https://nlopt.readthedocs.io/en/latest/NLopt_Introduction/#equivalent-formulations-of-optimization-problems)
for an overview of this approach. (Note: this tutorial example requires NLopt [version 2.7.0](https://github.com/stevengj/nlopt/releases/tag/v2.7.0) or higher.) The minimax/epigraph approach is
also covered in the [near-to-far field
tutorial](https://nbviewer.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/06-Near2Far-Epigraph.ipynb).

In this example, we use a minimum feature size of 150 nm for the
linewidth *and* linespacing. The implementation of these constraints
in the adjoint-solver module is based on [A.M. Hammond et al., Optics
Express, Vol. 29, pp. 23916-38
(2021)](https://doi.org/10.1364/OE.431188).

There are six important items to highlight in the set up of this
optimization problem:

1. The lengthscale constraint is activated only in the final epoch. It
   is often helpful to binarize the design using a large $\beta$ value
   for the projection operator (hyperbolic tangent) before this final
   epoch. This is because the lengthscale constraint forces
   binarization which could induce large changes in an initial
   grayscale design and thus irrevocably spoil the performance of the
   final design. Note that regardless of the value of $\beta$,
   projecting the design weights $u$ will produce grayscale values
   between 0 and 1 whenever $u \approx \eta \pm \frac{1}{\beta}$.

2. The initial value of the epigraph variable of the final epoch (in
  which the minimum feature size constraint is imposed) should take
  into account the value of the constraint itself. This ensures a
  feasible starting point for the [method of moving asymptotes (MMA)
  optimization
  algorithm](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#mma-method-of-moving-asymptotes-and-ccsa),
  (which is based on the conservative convex separable approximation
  (CCSA) algorithm).

3. The edge of the design region is padded by a filter radius (rather
  than e.g., a single pixel) to produce measured minimum feature sizes
  of the final design that are consistent with the imposed constraint.

4. The hyperparameters of the feature-size constraint function
  (`a1`, `b1`, and `c0` in the `glc` function of the script below), need to
  be chosen carefully to produce final designs which do not
  significantly degrade the performance of the unconstrained designs
  at the start of the final epoch.

5. Damping of the design weights is used for the early epochs in which
  the $\beta$ parameter of the projection function is small (< ~50) and
  the design is mostly grayscale in order to induce binarization.

6. The subpixel-smoothing feature of the [`MaterialGrid`](../Python_User_Interface.md#materialgrid)
  is necessary whenever the $\beta$ parameter of the projection function
  is large (> ~50) and thus the design is binary (or nearly so). Without
  subpixel smoothing, then the gradients are nearly zero
 except for $u \approx \eta \pm 1/\beta$ where the derivatives
 and second derivatives blow up, causing optimization algorithms to break down. When subpixel
  smoothing is enabled (`do_averaging=True`), the weights are projected
  *internally* using the `beta` parameter. For this reason, any
  preprocessing (i.e., mapping) of the weights outside of the `MaterialGrid`
  should apply only a filter to the weights but must not perform any
  projection.

A schematic of the final design and the simulation layout is shown
below. The minimum feature size of the final design, measured using a
[ruler](https://github.com/NanoComp/photonics-opt-testbed/tree/main/ruler),
is 165 nm. This value is consistent with the imposed constraint since
it is approximately within one design pixel (10 nm).

![](../images/mode_converter_sim_layout.png#center)

A plot of the reflectance and transmittance spectrum in linear and log
(dB) scales of the final design is shown below. The worst-case
reflectance is -17.7 dB at a wavelength of 1.295 μm. The worst-case
transmittance is -2.1 dB at a wavelength of 1.265 μm.

![](../images/mode_converter_refl_tran_spectra.png#center)

A plot of the objective-function history (worst case or maximum value
across the six wavelengths) for this design is also shown. The
"spikes" present in the plot are a normal feature of
nonlinear-optimization algorithms. The algorithm may take too large a
step which turns out to make the objective function worse. This means
the algorithm then has to "backtrack" and take a smaller step. This
occurs in the CCSA algorithm by increasing a penalty term. Also, even
though the worst-case objective function is constant during most of
the first epoch which may indicate the optimizer is making no
progress, in fact the optimizer is working to improve the objective
function at the other (non worst-case) wavelengths.

![](../images/mode_converter_objfunc_hist.png#center)

Finally, here are additional designs generated using constraints on
the minimum feature size of 50 nm, 70 nm, 90 nm, and 225 nm. The
designs with smaller minimum feature sizes are clearly distinguishable
from the designs with the larger ones.

![](../images/mode_converter_designs.png#center)

The table below shows a comparison of the imposed constraint on the
minimum feature size of the optimizer versus the measured minimum
linewidth and linespacing for the five designs presented in this
tutorial. There is fairly consistent agreement in the constraint and
measured values except for the design with the largest minimum feature
size (225 nm vs. 277 nm). For cases in which the measured minimum
feature size is significantly *larger* than the constraint used in the
optimization, this could be an indication that the final design can be
improved by a better choice of the hyperparameters in the feature-size
constraint function. Generally, one should expect the constraint and
measured values to agree within a length of about one to two
design-region pixels (10 nm, in this example).

| constraint (nm) | measured linewidth (nm) | measured linespacing (nm) |
|:---------------:|:-----------------------:|:-------------------------:|
|        50       |            85           |             60            |
|        70       |           103           |             85            |
|        90       |           109           |            134            |
|       150       |           221           |            165            |
|       225       |           277           |            681            |

Finally, a plot of the worst-case reflectance and transmittance versus
the imposed constraint on the minimum feature size is shown below for
the five designs. The general trend of decreasing performance (i.e.,
increasing reflectance and decreasing transmittance) with increasing
minimum feature size is evident.

![](../images/mode_converter_worst_case_refl_tran.png#center)

The script is
[python/examples/adjoint_optimization/mode_converter.py](https://github.com/NanoComp/meep/tree/master/python/examples/adjoint_optimization/mode_converter.py). The runtime of this script using three cores of an Intel Xeon at 2.0 GHz is approximately 14 hours.

```py
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
num_wavelengths = len(DESIGN_WAVELENGTHS_UM)
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


def filter_and_project(
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
      result: evaluation of this constraint function, modified in place.
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
            filter_and_project(
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
    grad = np.zeros((NX_DESIGN_GRID * NY_DESIGN_GRID, 2 * num_wavelengths))
    grad[:, :num_wavelengths] = grad_reflectance
    grad[:, num_wavelengths:] = grad_transmittance

    # Backpropagate the gradients through the filter and project function.
    for k in range(2 * num_wavelengths):
        grad[:, k] = tensor_jacobian_product(filter_and_project, 0)(
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
        f"epigraph: {epigraph:.5f}, obj. func.: {obj_val_merged_str}, "
        f"epigraph constraint: {str_from_list(result)}"
    )

    cur_iter[0] = cur_iter[0] + 1


def line_width_and_spacing_constraint(
    result: np.ndarray,
    epigraph_and_weights: np.ndarray,
    gradient: np.ndarray,
    sigmoid_bias: float,
) -> float:
    """Constraint function for the minimum line width and spacing.

    Args:
      result: evaluation of this constraint function, modified in place.
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
      A 2-tuple consisting of (1) a 1D array of DFT fields and (2) the DFT
      fields object returned by `meep.get_flux_data`.
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

    # Initial design weights (arbitrary constant value).
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

    # Threshold beta above which to use subpixel smoothing.
    sigmoid_bias_threshold = 64

    sigmoid_biases = [8, 16, 32, 64, 128, 256]
    max_evals = [80, 80, 100, 120, 120, 100]
    epigraph_tolerance = np.array([1e-4] * 2 * num_wavelengths)  # R, 1-T
    tolerance_width_and_spacing = np.array([1e-8] * 2)  # line width, line spacing

    for sigmoid_bias, max_eval in zip(sigmoid_biases, max_evals):
        solver = nlopt.opt(nlopt.LD_CCSAQ, num_weights + 1)
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
            epigraph_tolerance,
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

        # Apply the constraint for the minimum line width and spacing only in
        # the final epoch to an initial binary design from the previous epoch.
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

        # Execute a single forward run before the start of each epoch and
        # manually set the initial epigraph variable to slightly larger than
        # the largest value of the objective function over the six wavelengths
        # and the lengthscale constraint (final epoch only).
        epigraph_initial = opt(
            [
                filter_and_project(
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

        epigraph_and_weights[0] = np.max(epigraph_initial)
        fraction_max_epigraph = 0.05
        if sigmoid_bias == sigmoid_biases[-1]:
            epigraph_and_weights[0] = (1 + fraction_max_epigraph) * max(
                epigraph_and_weights[0], linewidth_constraint_val
            )
        print(
            f"epigraph-calibration:, {sigmoid_bias}, "
            f"{str_from_list(epigraph_initial)}, {epigraph_and_weights[0]}"
        )

        epigraph_and_weights[:] = solver.optimize(epigraph_and_weights)

        optimal_design_weights = filter_and_project(
            epigraph_and_weights[1:],
            SIGMOID_THRESHOLD_INTRINSIC,
            sigmoid_bias,
        ).reshape(NX_DESIGN_GRID, NY_DESIGN_GRID)

        # Save the unmapped weights and a bitmap image of the design weights
        # at the end of each epoch.
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

    # Save important optimization parameters and output for post processing.
    np.savez(
        "optimal_design.npz",
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
```


Derivatives with Respect to Shape Parameters
--------------------------------------------

It is also possible to compute the derivative of the Meep outputs with respect to a geometric parameter via a [level-set](https://en.wikipedia.org/wiki/Level_set) formulation (an [implicit-function](https://en.wikipedia.org/wiki/Implicit_function) representation of a material discontinuity) using the density-based adjoint solver. This is useful for applications involving [shape optimization](https://en.wikipedia.org/wiki/Shape_optimization) of explicitly parameterized geometries.

As a demonstration, we will compute the derivative of the diffraction efficiency of the first transmitted order with $\mathcal{S}$ polarization (electric field out of the plane, i.e. $E_z$) of a 1D binary grating with respect to the grating height. The accuracy of the adjoint derivative is validated using a brute-force finite-difference approximation.

The calculation of the diffraction efficiency involves [mode decomposition](Python_User_Interface/#mode-decomposition) of planewaves which is described in [Tutorial/Transmissive Diffraction Spectrum for Planewave at Normal Incidence](Mode_Decomposition.md#transmissive-diffraction-spectrum-for-planewave-at-normal-incidence). An important aspect of setting up this simulation is specifying a diffraction order (planewave) in 2D using the `meep.adjoint.EigenmodeCoefficient` object. This involves three components: (1) the wavevector via the `kpoint_func` parameter, (2) the polarization via the `eig_parity` parameter, and (3) the `eig_vol` parameter must be set to have a *length of one pixel in the periodic direction*. The latter is necessary for MPB to correctly interpret the wavevector in the Cartesian basis rather than the reciprocal-lattice basis defined by the grating period.

Calculating the derivative of the diffraction efficiency ($F$) with respect to the grating height ($h$) involves using the chain rule to obtain a product of two derivatives (Jacobians): (1) $\partial F / \partial \rho$ and (2) $\partial \rho / \partial h$, where $\rho$ are the density weights (2D grid of $N_x \times N_y = N$ points) used to represent the binary grating as a level set. (1) is the adjoint gradient computed by Meep. (2) requires the level set to be *smoothed* to ensure that the derivative $\partial \rho / \partial h$ exists. (2) is implemented in Autograd using a custom vector–Jacobian product (vJP) which is used as part of reverse-mode differentiation (backpropagation) to compute $\partial F / \partial h$. An overview of this calculation is shown below.

![](../images/levelset_gradient_backpropagation.png#center)

The derivative $\partial \rho / \partial h$ can be approximated by a finite difference using a one-pixel perturbation applied to the grating height; this is computationally convenient because it greatly simplifies the construction of $\rho(h)$ as explained below. (A finite difference involves two function evaluations of $\rho(h)$, but the cost for this is negligible since it involves no Meep simulations.) The construction of $\rho(h)$ involves two steps: first, constructing a simple binary image $b(h)$ at a high resolution; and second, smoothing $b(h)$ into a continuous level-set function $\rho(h)$. This smoothing of the image can be performed using a number of different methods including a [signed-distance function](https://en.wikipedia.org/wiki/Signed_distance_function) or convolution filter. In this example, the smoothing is based simply on downsampling the image from a high-resolution grid (10$\times$ the resolution of the simulation grid) to the lower-resolution simulation grid using bilinear interpolation, which leads to "gray" pixels at the boundaries between materials in a way that changes continuously with $h$. Only these boundary pixels have nonzero derivatives in the Jacobian $\partial\rho/\partial h$ in this case. This calculation is summarized in the schematic below.

![](../images/levelset_jacobian_matrix.png#center)

A schematic of the simulation layout showing the design region containing the grating is shown below. The adjoint gradient $\partial F / \partial \rho$ computed by Meep at a resolution of 400 pixels/$\mu$m is shown next to this schematic.

![](../images/levelset_adjoint_gradient.png#center)

The simulation script is in [python/examples/adjoint_optimization/binary_grating_levelset.py](https://github.com/NanoComp/meep/tree/master/python/examples/adjoint_optimization/binary_grating_levelset.py). Running this script at five different resolutions produces the output below for the derivative computed using the finite difference and the adjoint derivative as well as the relative error from these two results.

```
RESOLUTION_UM = 50
deriv:, -0.02535355 (finite difference), -0.01469133 (adjoint), 0.420542 (error)

RESOLUTION_UM = 100
deriv:, -0.00604817 (finite difference), -0.00473221 (adjoint), 0.217580 (error)

RESOLUTION_UM = 200
deriv:, -0.00284452 (finite difference), -0.00252470 (adjoint), 0.112432 (error)

RESOLUTION_UM = 400
deriv:, -0.00140221 (finite difference), -0.00132065 (adjoint), 0.058165 (error)

RESOLUTION_UM = 800
deriv:, -0.00069606 (finite difference), -0.00067547 (adjoint), 0.029583 (error)
```

A logarithmic plot of the relative error vs. grid resolution based on these results demonstrates linear convergence. This is expected for fields right on the boundary of a discontinuous interface.

Currently, we recommend using this procedure only for the $\mathcal{S}$ ($E_z$) polarization, since for the $\mathcal{P}$ polarization (electric field in the $xy$ plane), there appear to be large discretization errors in the adjoint gradient which [we are currently investigating](https://github.com/NanoComp/meep/pull/2792#issuecomment-2171942477). This is why this demonstration only involved the $\mathcal{S}$ polarization.

![](../images/levelset_adjoint_gradient_error_vs_resolution.png#center)

```py
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
GRATING_HEIGHT_PERTURBATION_UM = 1.0 / RESOLUTION_UM
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
        f"deriv:, {fnd_directional_deriv:.8f} (finite difference), "
        f"{adj_directional_deriv:.8f} (adjoint), {rel_err:.6f} (error)"
    )
```

Compact Notebook Tutorials of Basic Features
--------------------------------------------

As an alternative to the first tutorial which combined multiple
features into a single demonstration, there are six notebook tutorials
that demonstrate various basic features of the adjoint solver.

- [Introduction](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/01-Introduction.ipynb)

- [Waveguide Bend Optimization](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/02-Waveguide_Bend.ipynb)

- [Filtering and Thresholding Design Parameters and Broadband Objective Functions](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/03-Filtered_Waveguide_Bend.ipynb)

- [Design of a Symmetric Broadband Splitter](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/04-Splitter.ipynb)

- [Metalens Optimization with Near2Far](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/05-Near2Far.ipynb)

- [Near2Far Optimization with Epigraph Formulation](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/06-Near2Far-Epigraph.ipynb)

- [Connectivity Constraint](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/07-Connectivity-Constraint.ipynb)
