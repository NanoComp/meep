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
[python/examples/adjoint_optimization/mode_converter.py](https://github.com/NanoComp/meep/tree/master/python/examples/adjoint_optimization/mode_converter.py). The
runtime of this script using three Intel Xeon 2.0 GHz processors is
approximately 14 hours.

```py
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from autograd import numpy as npa, tensor_jacobian_product, grad
import nlopt
import meep as mp
import meep.adjoint as mpa
from typing import NamedTuple

resolution = 50  # pixels/μm

w = 0.4  # waveguide width
l = 3.0  # waveguide length (on each side of design region)
dpad = 0.6  # padding length above/below design region
dpml = 1.0  # PML thickness
dx = 1.6  # length of design region
dy = 1.6  # width of design region

sx = dpml + l + dx + l + dpml
sy = dpml + dpad + dy + dpad + dpml
cell_size = mp.Vector3(sx, sy, 0)

pml_layers = [mp.PML(thickness=dpml)]

# wavelengths for minimax optimization
wvls = (1.265, 1.270, 1.275, 1.285, 1.290, 1.295)
frqs = [1 / wvl for wvl in wvls]

minimum_length = 0.15  # minimum length scale (μm)
eta_i = 0.5  # blueprint design field thresholding point (between 0 and 1)
eta_e = 0.75  # erosion design field thresholding point (between 0 and 1)
eta_d = 1 - eta_e  # dilation design field thresholding point (between 0 and 1)
filter_radius = mpa.get_conic_radius_from_eta_e(minimum_length, eta_e)
print(f"filter_radius:, {filter_radius:.6f}")

# pulsed source center frequency and bandwidth
wvl_min = 1.26
wvl_max = 1.30
frq_min = 1 / wvl_max
frq_max = 1 / wvl_min
fcen = 0.5 * (frq_min + frq_max)
df = frq_max - frq_min

eig_parity = mp.ODD_Z
src_pt = mp.Vector3(-0.5 * sx + dpml, 0, 0)

nSiO2 = 1.5
SiO2 = mp.Medium(index=nSiO2)
nSi = 3.5
Si = mp.Medium(index=nSi)

design_region_size = mp.Vector3(dx, dy, 0)
design_region_resolution = int(2 * resolution)
Nx = int(design_region_size.x * design_region_resolution) + 1
Ny = int(design_region_size.y * design_region_resolution) + 1

# impose a bit "mask" of thickness equal to the filter radius
# around the edges of the design region in order to prevent
# violations of the minimum feature size constraint.

x_g = np.linspace(
    -design_region_size.x / 2,
    design_region_size.x / 2,
    Nx,
)
y_g = np.linspace(
    -design_region_size.y / 2,
    design_region_size.y / 2,
    Ny,
)
X_g, Y_g = np.meshgrid(
    x_g,
    y_g,
    sparse=True,
    indexing="ij",
)

left_wg_mask = (X_g <= -design_region_size.x / 2 + filter_radius) & (
    np.abs(Y_g) <= w / 2
)
right_wg_mask = (X_g >= design_region_size.x / 2 - filter_radius) & (
    np.abs(Y_g) <= w / 2
)
Si_mask = left_wg_mask | right_wg_mask

border_mask = (
    (X_g <= -design_region_size.x / 2 + filter_radius)
    | (X_g >= design_region_size.x / 2 - filter_radius)
    | (Y_g <= -design_region_size.y / 2 + filter_radius)
    | (Y_g >= design_region_size.y / 2 - filter_radius)
)
SiO2_mask = border_mask.copy()
SiO2_mask[Si_mask] = False

refl_pt = mp.Vector3(-0.5 * sx + dpml + 0.5 * l)
tran_pt = mp.Vector3(0.5 * sx - dpml - 0.5 * l)

stop_cond = mp.stop_when_fields_decayed(50, mp.Ez, refl_pt, 1e-8)


def mapping(x: np.ndarray, eta: float, beta: float) -> np.ndarray:
    """A differentiable mapping function which applies, in order,
       the following sequence of transformations to the design weights:
       (1) a bit mask for the boundary pixels, (2) convolution with a
       conic filter, and (3) projection via a hyperbolic tangent (if
       necessary).

    Args:
      x: design weights as a 1d array of size Nx*Ny.
      eta: erosion/dilation parameter for the projection.
      beta: bias parameter for the projection. A value of 0 is no projection.

    Returns:
      The mapped design weights as a 1d array.
    """
    x = npa.where(
        Si_mask.flatten(),
        1,
        npa.where(
            SiO2_mask.flatten(),
            0,
            x,
        ),
    )

    filtered_field = mpa.conic_filter(
        x,
        filter_radius,
        design_region_size.x,
        design_region_size.y,
        design_region_resolution,
    )

    if beta == 0:
        return filtered_field.flatten()

    else:
        projected_field = mpa.tanh_projection(
            filtered_field,
            beta,
            eta,
        )

        return projected_field.flatten()


def f(x: np.ndarray, grad: np.ndarray) -> float:
    """Objective function for the epigraph formulation.

    Args:
      x: 1d array of size 1+Nx*Ny containing epigraph variable (first element)
         and design weights (remaining Nx*Ny elements).
      grad: the gradient as a 1d array of size 1+Nx*Ny modified in place.

    Returns:
      The epigraph variable (a scalar).
    """
    t = x[0]  # epigraph variable
    v = x[1:]  # design weights
    if grad.size > 0:
        grad[0] = 1
        grad[1:] = 0
    return t


def c(result: np.ndarray, x: np.ndarray, gradient: np.ndarray, eta: float,
      beta: float, use_epsavg: bool):
    """Constraint function for the epigraph formulation.

    Args:
      result: the result of the function evaluation modified in place.
      x: 1d array of size 1+Nx*Ny containing epigraph variable (first
         element) and design weights (remaining Nx*Ny elements).
      gradient: the Jacobian matrix with dimensions (1+Nx*Ny,
                2*num. wavelengths) modified in place.
      eta: erosion/dilation parameter for projection.
      beta: bias parameter for projection.
      use_epsavg: whether to use subpixel smoothing.
    """
    t = x[0]  # epigraph variable
    v = x[1:]  # design weights

    f0, dJ_du = opt([mapping(v, eta, 0 if use_epsavg else beta)])

    f0_reflection = f0[0]
    f0_transmission = f0[1]
    f0_merged = np.concatenate((f0_reflection, f0_transmission))
    f0_merged_str = '[' + ','.join(str(ff) for ff in f0_merged) + ']'

    dJ_du_reflection = dJ_du[0]
    dJ_du_transmission = dJ_du[1]
    nfrq = len(frqs)
    my_grad = np.zeros((Nx * Ny, 2 * nfrq))
    my_grad[:, :nfrq] = dJ_du_reflection
    my_grad[:, nfrq:] = dJ_du_transmission

    # backpropagate the gradients through mapping function
    for k in range(2 * nfrq):
        my_grad[:, k] = tensor_jacobian_product(mapping, 0)(
            v,
            eta,
            beta,
            my_grad[:, k],
        )

    if gradient.size > 0:
        gradient[:, 0] = -1  # gradient w.r.t. epigraph variable ("t")
        gradient[:, 1:] = my_grad.T  # gradient w.r.t. each frequency objective

    result[:] = np.real(f0_merged) - t

    objfunc_history.append(np.real(f0_merged))
    epivar_history.append(t)

    print(
        f"iteration:, {cur_iter[0]:3d}, eta: {eta}, beta: {beta:2d}, "
        f"t: {t:.5f}, obj. func.: {f0_merged_str}"
    )

    cur_iter[0] = cur_iter[0] + 1


def glc(result: np.ndarray, x: np.ndarray, gradient: np.ndarray,
        beta: float) -> float:
    """Constraint function for the minimum linewidth.

    Args:
      result: the result of the function evaluation modified in place.
      x: 1d array of size 1+Nx*Ny containing epigraph variable (first
         element) and design weights (remaining elements).
      gradient: the Jacobian matrix with dimensions (1+Nx*Ny,
                num. wavelengths) modified in place.
      beta: bias parameter for projection.

    Returns:
      The value of the constraint function (a scalar).
    """
    t = x[0]  # dummy parameter
    v = x[1:]  # design parameters
    a1 = 1e-3  # hyper parameter (primary)
    b1 = 0  # hyper parameter (secondary)
    gradient[:, 0] = -a1

    filter_f = lambda a: mpa.conic_filter(
        a.reshape(Nx, Ny),
        filter_radius,
        design_region_size.x,
        design_region_size.y,
        design_region_resolution,
    )
    threshold_f = lambda a: mpa.tanh_projection(a, beta, eta_i)

    # hyper parameter (constant factor and exponent)
    c0 = 1e7 * (filter_radius * 1 / resolution) ** 4

    M1 = lambda a: mpa.constraint_solid(a, c0, eta_e, filter_f, threshold_f, 1)
    M2 = lambda a: mpa.constraint_void(a, c0, eta_d, filter_f, threshold_f, 1)

    g1 = grad(M1)(v)
    g2 = grad(M2)(v)

    result[0] = M1(v) - a1 * t - b1
    result[1] = M2(v) - a1 * t - b1

    gradient[0, 1:] = g1.flatten()
    gradient[1, 1:] = g2.flatten()

    t1 = (M1(v) - b1) / a1
    t2 = (M2(v) - b1) / a1

    print(f"glc:, {result[0]}, {result[1]}, {t1}, {t2}")

    return max(t1, t2)


def straight_waveguide() -> (np.ndarray, NamedTuple):
    """Computes the DFT fields from the mode source in a straight waveguide
       for use as normalization of the reflectance measurement during the
       optimization.

    Returns:
      A 2-tuple consisting of a 1d array of DFT fields and DFT fields object
      returned by `meep.get_flux_data`.
    """
    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            size=mp.Vector3(0, sy, 0),
            center=src_pt,
            eig_band=1,
            eig_parity=eig_parity,
        )
    ]

    geometry = [
        mp.Block(
            size=mp.Vector3(mp.inf, w, mp.inf),
            center=mp.Vector3(),
            material=Si,
        )
    ]

    sim = mp.Simulation(
        resolution=resolution,
        default_material=SiO2,
        cell_size=cell_size,
        sources=sources,
        geometry=geometry,
        boundary_layers=pml_layers,
        k_point=mp.Vector3(),
    )

    refl_mon = sim.add_mode_monitor(
        frqs,
        mp.ModeRegion(center=refl_pt, size=mp.Vector3(0, sy, 0)),
        yee_grid=True,
    )

    sim.run(until_after_sources=stop_cond)

    res = sim.get_eigenmode_coefficients(
        refl_mon,
        [1],
        eig_parity=eig_parity,
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
        beta: float) -> mpa.OptimizationProblem:
    """Sets up the adjoint optimization of the waveguide mode converter.

    Args:
      input_flux: 1d array of DFT fields from normalization run.
      input_flux_data: DFT fields object returned by `meep.get_flux_data`.
      use_damping: whether to use the damping feature of `MaterialGrid`.
      use_epsavg: whether to use subpixel smoothing in `MaterialGrid`.

    Returns:
      A `meep.adjoint.OptimizationProblem` class object.
    """
    matgrid = mp.MaterialGrid(
        mp.Vector3(Nx, Ny, 0),
        SiO2,
        Si,
        weights=np.ones((Nx, Ny)),
        beta=beta if use_epsavg else 0,
        do_averaging=True if use_epsavg else False,
        damping=0.02 * 2 * np.pi * fcen if use_damping else 0,
    )

    matgrid_region = mpa.DesignRegion(
        matgrid,
        volume=mp.Volume(
            center=mp.Vector3(),
            size=mp.Vector3(design_region_size.x, design_region_size.y, mp.inf),
        ),
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
            size=mp.Vector3(mp.inf, w, mp.inf),
            material=Si,
        )
    ]

    geometry += matgrid_geometry

    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            size=mp.Vector3(0, sy, 0),
            center=src_pt,
            eig_band=1,
            eig_parity=eig_parity,
        ),
    ]

    sim = mp.Simulation(
        resolution=resolution,
        default_material=SiO2,
        cell_size=cell_size,
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
                size=mp.Vector3(0, sy, 0),
            ),
            1,
            forward=False,
            eig_parity=eig_parity,
            subtracted_dft_fields=input_flux_data,
        ),
        mpa.EigenmodeCoefficient(
            sim,
            mp.Volume(
                center=tran_pt,
                size=mp.Vector3(0, sy, 0),
            ),
            2,
            eig_parity=eig_parity,
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
        frequencies=frqs,
    )

    return opt


if __name__ == "__main__":
    input_flux, input_flux_data = straight_waveguide()

    algorithm = nlopt.LD_MMA

    # number of design parameters
    n = Nx * Ny

    # initial design parameters
    x = np.ones((n,)) * 0.5
    x[Si_mask.flatten()] = 1.0  # set the edges of waveguides to silicon
    x[SiO2_mask.flatten()] = 0.0  # set the other edges to SiO2

    # lower and upper bounds for design weights
    lb = np.zeros((n,))
    lb[Si_mask.flatten()] = 1.0
    ub = np.ones((n,))
    ub[SiO2_mask.flatten()] = 0.0

    # insert epigraph variable initial value (arbitrary) and bounds into the
    # design array. the actual value is determined by the objective and
    # constraint functions below.
    x = np.insert(x, 0, 1.2)
    lb = np.insert(lb, 0, -np.inf)
    ub = np.insert(ub, 0, +np.inf)

    objfunc_history = []
    epivar_history = []
    cur_iter = [0]

    beta_thresh = 64 # threshold beta above which to use subpixel smoothing
    betas = [8, 16, 32, 64, 128, 256]
    max_evals = [80, 80, 100, 120, 120, 100]
    tol_epi = np.array([1e-4] * 2 * len(frqs))  # R, 1-T
    tol_lw = np.array([1e-8] * 2)  # line width, line spacing

    for beta, max_eval in zip(betas, max_evals):
        solver = nlopt.opt(algorithm, n + 1)
        solver.set_lower_bounds(lb)
        solver.set_upper_bounds(ub)
        solver.set_min_objective(f)
        solver.set_maxeval(max_eval)
        solver.set_param("dual_ftol_rel", 1e-7)
        solver.add_inequality_mconstraint(
            lambda rr, xx, gg: c(
                rr,
                xx,
                gg,
                eta_i,
                beta,
                False if beta < beta_thresh else True,
            ),
            tol_epi,
        )
        solver.set_param("verbosity", 1)

        opt = mode_converter_optimization(
            input_flux,
            input_flux_data,
            True,  # use_damping
            False if beta < beta_thresh else True,  # use_epsavg
            beta,
        )

        # apply the minimum linewidth constraint
        # only in the final epoch to an initial
        # binary design from the previous epoch.
        if beta == betas[-1]:
            res = np.zeros(2)
            grd = np.zeros((2, n + 1))
            t = glc(res, x, grd, beta)
            solver.add_inequality_mconstraint(
                lambda rr, xx, gg: glc(
                    rr,
                    xx,
                    gg,
                    beta,
                ),
                tol_lw,
            )

        # execute a single forward run before the start of each
        # epoch and manually set the initial epigraph variable to
        # slightly larger than the largest value of the objective
        # function over the six wavelengths and the lengthscale
        # constraint (final epoch only).
        t0 = opt(
            [
                mapping(
                    x[1:],
                    eta_i,
                    beta if beta < beta_thresh else 0,
                ),
            ],
            need_gradient=False,
        )
        t0 = np.concatenate((t0[0][0], t0[0][1]))
        t0_str = '[' + ','.join(str(tt) for tt in t0) + ']'
        x[0] = np.amax(t0)
        x[0] = 1.05 * (max(x[0], t) if beta == betas[-1] else x[0])
        print(f"data:, {beta}, {t0_str}, {x[0]}")

        x[:] = solver.optimize(x)

        optimal_design_weights = mapping(
            x[1:],
            eta_i,
            beta,
        ).reshape(Nx, Ny)

        # save the unmapped weights and a bitmap image
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
                f"optimal_design_beta{beta}.png",
                dpi=150,
                bbox_inches="tight",
            )
            # save the final (unmapped) design as a 2d array in CSV format
            np.savetxt(
                f"unmapped_design_weights_beta{beta}.csv",
                x[1:].reshape(Nx, Ny),
                fmt="%4.2f",
                delimiter=",",
            )

    # save all the important optimization parameters and output
    # as separate arrays in a single file for post processing.
    with open("optimal_design.npz", "wb") as fl:
        np.savez(
            fl,
            Nx=Nx,
            Ny=Ny,
            design_region_size=(dx, dy),
            design_region_resolution=design_region_resolution,
            betas=betas,
            max_eval=max_eval,
            objfunc_history=objfunc_history,
            epivar_history=epivar_history,
            t=x[0],
            unmapped_design_weights=x[1:],
            minimum_length=minimum_length,
            optimal_design_weights=optimal_design_weights,
        )
```

Binary Grating as a Level Set
-----------------------------

This is a demonstration of computing the gradient of a level set using the adjoint solver.

The simulation script is in [python/examples/adjoint_optimization/binary_grating_levelset.py](https://github.com/NanoComp/meep/tree/master/python/examples/adjoint_optimization/binary_grating_levelset.py).

![](../images/levelset_gradient_backpropagation.png#center)

![](../images/levelset_jacobian_matrix.png#center)

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
