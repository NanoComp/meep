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

Shape Optimization of a Multilayer Stack
----------------------------------------

We extend the demonstration of the shape derivative from the previous tutorial to perform shape optimization of a multilayer stack over a broad bandwidth. The 1D design problem is shown in the schematic below and involves finding the layer thicknesses for a fixed number of layers (9) which minimize the integrated field intensity $\int |E_x|^2$, which demonstrates a `FourierFields` objective function.  (This is equivalent to minimizing absorption if $\varepsilon$ had a small imaginary part, and is related but not precisely equivalent to minimizing transmission or maximizing reflection.)  In particular, we minimize the worst case (largest) of the intensities at two wavelengths: $\lambda_1$ = 0.95 μm and $\lambda_2$ = 1.05 μm, to roughly emulate a broadband problem.

The stack consists of two materials of alternating refractive index $n_A$ = 1.3 and $n_B$ = 1.0. The layers are arranged as $n_A$, $n_B$, $n_A$, $n_B$, ..., $n_A$. The semi-infinite regions to the left and right of the stack are vacuum.

![](../images/multilayer_stack_shape_opt.png#center)

A reference design to compare our results against is a [quarter-wavelength stack](https://en.wikipedia.org/wiki/Distributed_Bragg_reflector). The mean wavelength of $\lambda_1$ and $\lambda_2$ is $\lambda = 1.0 μm for which the quarter-wavelength layer thicknesses are $\lambda / (4 n_A)$ = 0.19 μm and $\lambda / (4 n_B)$ = 0.25 μm. These values are used to specify upper and lower bounds for the layer thicknesses. This is important given the non-convex nature of this particular design problem.

(For several layers and a single wavelength, the optimizer can easily find a non-quarter-wavelength structure with very low transmission, and once it gets to that point the gradients will be almost zero so probably it will converge very slowly towards a quarter-wavelength solution because the optimizer is generally very slow to try to reduce transmittance from 1e-3 to 1e-6, for example.)

The worst-case optimization is implemented using the [epigraph formulation](https://nlopt.readthedocs.io/en/latest/NLopt_Introduction/#equivalent-formulations-of-optimization-problems) and the Conservative Convex Separable Approximation (CCSA) algorithm. The optimization is run ten times with random initial designs from which the local optima with the smallest objective value is chosen. For the run involving the best design, a plot of the objective function vs. iteration number for the two wavelengths is shown below. Also included in this plot is the epigraph variable. These results demonstrate that the stack is being optimized for $\lambda$ = 0.95 μm and *not* $\lambda$ = 1.05 μm. We observed this trend for various local optima. The optimizer converges to the local optima using 13 iterations. The nine layer thicknesses (μm) of this design are: 0.1822, 0.2369, 0.1822, 0.2369, 0.1822, 0.2369, 0.1910, 0.2504, 0.1931. Since these are 1D simulations, the runtime for each design iteration is generally fast (about ten seconds using one core of an Intel Xeon i7-7700K CPU @ 4.20GHz).

In some cases (not shown), the optimizer may take some steps that make things *worse* during the early iterations. This is a common phenomenon in many algorithms — initially, the optimizer takes steps that are too large and then has to backtrack. After a few iterations, the algorithm has more information about a "good" step size.

![](../images/multilayer_opt_obj_func_history.png#center)

Finally, a plot of $|E_x|^2$ in the stack (design region) normalized by the fields from the quarter-wavelength stack is shown below. The plot shows that while the DFT fields for each of the two wavelengths are decaying through the stack, the transmittance values are different: T($\lambda_1$ = 0.95 μm) = 0.2574 vs. T($\lambda_2$ = 1.05 μm) = 0.3734. The transmittance for the quarter-wavelength stack at $\lambda \approx$ 1.0 μm is 0.2522. Note that the magnitude of $|E_x|^2$ at the right edge of the stack is larger for $\lambda_2$ = 1.05 μm than $\lambda_1$ = 0.95 μm. This is consistent with the transmittance being larger for $\lambda_2$ = 1.05 than for $\lambda_1$ = 0.95 μm.

![](../images/broadband_field_decay_in_multilayer_stack.png#center)

The script is [python/examples/adjoint_optimization/multilayer_opt.py](https://github.com/NanoComp/meep/tree/master/python/examples/adjoint_optimization/multilayer_opt.py).

```py
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
    z_grid = np.linspace(
        -0.5 * DESIGN_REGION_UM.z,
        0.5 * DESIGN_REGION_UM.z,
        nz
    )

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
        ans: np.ndarray,
        layer_thickness_um: np.ndarray
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
        do_averaging=False
    )

    matgrid_region = mpa.DesignRegion(
        matgrid,
        volume=mp.Volume(
            center=mp.Vector3(),
            size=DESIGN_REGION_UM
        ),
    )

    geometry = [
        mp.Block(
            material=matgrid,
            size=matgrid_region.size,
            center=matgrid_region.center
        )
    ]

    sim = mp.Simulation(
        resolution=RESOLUTION_UM,
        cell_size=cell_size,
        dimensions=1,
        boundary_layers=pml_layers,
        sources=sources,
        geometry=geometry
    )

    obj_args = [
        mpa.FourierFields(
            sim,
            volume=matgrid_region.volume,
            component=mp.Ex,
            yee_grid=True
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
        return npa.log(npa.sum(npa.absolute(dft_ex)**2, axis=1))

    opt = mpa.OptimizationProblem(
        simulation=sim,
        objective_functions=obj_func,
        objective_arguments=obj_args,
        design_regions=[matgrid_region],
        frequencies=frequencies
    )

    return opt


def obj_func(
        epigraph_and_layer_thickness_um: np.ndarray,
        grad: np.ndarray
) -> float:
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
        result: float,
        epigraph_and_layer_thickness_um: np.ndarray,
        gradient: np.ndarray
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
            layer_thickness_um,
            grad[:, k]
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
    fraction_thickness = 0.08
    layer_thickness_um_lower_bound[0::2] = ((1 - fraction_thickness) *
                                            mean_layer_thickness_um[1])
    layer_thickness_um_upper_bound[0::2] = ((1 + fraction_thickness) *
                                            mean_layer_thickness_um[1])
    layer_thickness_um_lower_bound[1::2] = ((1 - fraction_thickness) *
                                            mean_layer_thickness_um[0])
    layer_thickness_um_upper_bound[1::2] = ((1 + fraction_thickness) *
                                            mean_layer_thickness_um[0])

    # Insert bounds for the epigraph variable.
    epigraph_and_layer_thickness_um_lower_bound = np.insert(
        layer_thickness_um_lower_bound, 0, 0
    )
    epigraph_and_layer_thickness_um_upper_bound = np.insert(
        layer_thickness_um_upper_bound, 0, np.inf
    )
    epigraph_tolerance = [0.1] * num_wavelengths

    min_obj_val = np.inf
    min_epigraph_and_layer_thickness_um = np.zeros(NUM_LAYERS + 1)
    obj_func_history = []
    epigraph_variable_history = []

    for j in range(NUM_OPT_REPEAT):
        rng = np.random.RandomState()
        layer_thickness_um = (
            layer_thickness_um_lower_bound + rng.rand(NUM_LAYERS) *
            (layer_thickness_um_upper_bound - layer_thickness_um_lower_bound)
        )

        # Execute a single forward run before the start of the optimization and
        # set the initial epigraph variable to slightly larger than the
        # largest value of the objective function over the wavelengths.
        design_weights = levelset_and_smoothing(layer_thickness_um)
        epigraph, _ = opt(
            [design_weights],
            need_gradient=False
        )
        fraction_max_epigraph = 0.2
        epigraph_initial = (1 + fraction_max_epigraph) * np.amax(epigraph)
        epigraph_and_layer_thickness_um_upper_bound[0] = epigraph_initial

        epigraph_and_layer_thickness_um = np.concatenate(
            (
                [epigraph_initial],
                layer_thickness_um
            )
        )

        solver = nlopt.opt(nlopt.LD_CCSAQ, NUM_LAYERS + 1)
        solver.set_lower_bounds(epigraph_and_layer_thickness_um_lower_bound)
        solver.set_upper_bounds(epigraph_and_layer_thickness_um_upper_bound)
        solver.set_min_objective(obj_func)
        solver.set_maxeval(MAX_OPT_ITERATIONS)
        solver.add_inequality_mconstraint(
            epigraph_constraint,
            epigraph_tolerance
        )
        solver.set_ftol_rel(0.02)

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
            print(f"new_minima:, {j:2d}, {min_obj_val}")


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
        optimal_obj_val=min_obj_val
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
