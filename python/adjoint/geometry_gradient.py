"""Adjoint gradients with respect to a geometric object's centre and size.

The design gradient finite-differences `eff_chi1inv_row` -- the subpixel
smoothed inverse permittivity over a voxel -- with respect to a design weight.
That function is a function of the *geometry*, so perturbing an object's centre
or size instead gives the shape derivative from the same machinery, at no cost
in additional timestepping. It is the same pattern as the design weights and as
the Gaussian beam parameters: finite-difference a cheap analytic map, contract
against the adjoint field.

Two things differ from the design gradient.

A design weight is local, so it can be perturbed inside the point loop. A centre
or a size is global, so the perturbation hoists out: move the geometry once,
visit every point, restore. Twelve passes for six parameters, not twelve
perturbations per point.

And the step has units. `utils.FD_DEFAULT` is a dimensionless perturbation of a
weight; here it is a length, and the scale that makes sense is a fraction of a
pixel. `DEFAULT_STEP_PIXELS` sets it relative to the grid rather than absolutely.

Accuracy: a shape derivative is one order worse than the fields
------------------------------------------------------------------
Subpixel smoothing makes the permittivity -- and so the objective -- a
second-order accurate, *continuous* function of an object's position. Measured
on a block of index 2.5 at resolution 20, turning smoothing off makes the
objective piecewise constant in position, swinging 77% across a single pixel;
turning it on makes it continuous, swinging 14%.

But second-order accuracy in the value is only first-order accuracy in the
derivative. Writing the discretization error as a function of the sub-pixel
phase,

    J(p) = J_exact(p) + Delta^2 E(p / Delta)

with E an O(1) oscillatory function -- which is why the residual has a period of
exactly one pixel -- differentiating gives

    dJ/dp = dJ_exact/dp + Delta^2 (1/Delta) E'(p / Delta)
          = dJ_exact/dp + O(Delta)

Measured: the per-pixel swing falls 9.9% -> 4.0% -> 1.9% at resolutions
20 -> 40 -> 80, halving with each doubling rather than quartering. This is
inherent to differentiating a smoothed staircase and cannot be recovered by
improving the adjoint, because the error is in the function being
differentiated, not in how it is differentiated.

Two consequences for testing. The adjoint should be checked against a finite
difference of the *discrete* objective, since both differentiate the same
smoothed function and should agree closely; that is a statement about the
implementation. Agreement with the physically intended derivative is limited to
O(Delta) and is a statement about the discretization.
"""

import warnings
from typing import List, Optional

import numpy as np

import meep as mp

# Perturbation for the geometric finite difference, as a fraction of a pixel.
# An absolute step is the wrong idea: subpixel smoothing varies over a pixel, so
# the step has to be defined against the grid.
#
# Measured on a block whose faces sit mid-voxel, against a finite difference of
# the objective: 0.2 px is 10% out, 0.05 px is 6.5%, 0.02 px is 3.4%, 0.005 px
# is 1.3% and 0.001 px is 0.7%, still improving. The error is truncation, not
# roundoff, over that whole range, so the step wants to be small.
DEFAULT_STEP_PIXELS = 0.002

# How far past an object the gradient monitor reaches, in pixels. The support of
# d(epsilon)/d(parameter) is the shell of voxels the boundary sweeps through, so
# it extends outside the object itself.
MONITOR_PAD_PIXELS = 2.0

# Indices `geometry_addgradient` uses, matching the geom_param enumeration in
# meepgeom.cpp.
PARAMETER_INDEX = {
    "center": (0, 1, 2),
    "size": (3, 4, 5),
}

DIFFERENTIABLE_PARAMS = tuple(PARAMETER_INDEX)


def differentiable_objects(sim: mp.Simulation) -> List:
    """The geometric objects flagged for differentiation, with their indices."""
    return [
        (index, obj)
        for index, obj in enumerate(sim.geometry)
        if getattr(obj, "differentiable", ())
    ]


def object_key(obj, index: int):
    """The key an object's gradient appears under: its name, or its position."""
    name = getattr(obj, "name", None)
    return name if name is not None else index


def _check_supported(obj) -> None:
    """Refuse the cases the diagonal-only contraction would get wrong."""
    if not isinstance(obj, mp.Block):
        raise NotImplementedError(
            f"Geometry gradients are implemented for mp.Block, not "
            f"{type(obj).__name__}. 'center' is meaningful for any object, but "
            "'size' is not, and the contraction currently assumes axis-aligned "
            "faces."
        )
    axes = (obj.e1, obj.e2, obj.e3)
    expected = (mp.Vector3(1, 0, 0), mp.Vector3(0, 1, 0), mp.Vector3(0, 0, 1))
    for axis, unit in zip(axes, expected):
        if abs(axis.x - unit.x) + abs(axis.y - unit.y) + abs(axis.z - unit.z) > 1e-12:
            raise NotImplementedError(
                "Geometry gradients assume an axis-aligned block, so that "
                "'size' means what the three components suggest. The "
                "off-diagonal terms subpixel smoothing produces at a tilted "
                "boundary *are* contracted; it is the parameterization that is "
                "restricted, not the physics."
            )


def check_faces_off_voxel_edges(sim: mp.Simulation, obj) -> None:
    """Warn when a face sits exactly on a voxel edge, where no derivative exists.

    A voxel's filling fraction is piecewise linear in the position of the
    boundary crossing it, with a kink each time the boundary reaches a voxel
    edge. The smoothed permittivity -- and so the objective -- is therefore C0
    but not C1 in an object's position, and *on* an edge the left and right
    derivatives differ.

    A central difference straddles that kink and returns a weighted mixture of
    the two, with weights that depend on the step. Measured on such a face, the
    reported gradient swept monotonically from -0.19 to +0.12 as the step went
    from 0.2 px to 0.001 px, passing through the true value without settling.
    Moving the same face half a pixel made the sweep converge to 0.7%.

    This is easy to hit by accident, because round sizes at round resolutions
    land on edges: a 0.8-wide block at resolution 20 has faces exactly 8 pixels
    from its centre.
    """
    dx = 1.0 / sim.resolution
    size = (obj.size.x, obj.size.y, obj.size.z)
    center = (obj.center.x, obj.center.y, obj.center.z)
    for axis, letter in enumerate("xyz"):
        if size[axis] == 0:
            continue
        for face in (center[axis] - size[axis] / 2, center[axis] + size[axis] / 2):
            # Pixel centres sit at integer multiples of dx, so pixel *edges*
            # are at half-integers. A face at a pixel centre is straddled by
            # that pixel and smooths normally; a face at a pixel edge is
            # straddled by nothing, every neighbouring pixel is full or empty,
            # and the derivative there is one-sided.
            offset = abs((face / dx) - round(face / dx))
            if abs(offset - 0.5) < 1e-6:
                warnings.warn(
                    f"The {letter} face of "
                    f"{getattr(obj, 'name', None) or 'this object'} at "
                    f"{face:.6g} lies on a pixel edge, where the smoothed "
                    "permittivity has a kink and the derivative with respect "
                    "to position is one-sided. The reported gradient will "
                    "depend on the finite-difference step rather than "
                    f"converging. Shift it by about {dx / 2:.6g} along "
                    + letter
                    + " so the face lands on a pixel centre, or change the "
                    "resolution.",
                    RuntimeWarning,
                    stacklevel=3,
                )


def check_smoothing(sim: mp.Simulation) -> None:
    """Refuse to differentiate a staircase.

    Without subpixel smoothing the effective permittivity is a step function of
    an object's position: shifting a boundary changes nothing until it crosses a
    voxel edge, then changes by the full material contrast. A finite difference
    of that is zero almost everywhere and a spike in a few voxels, and the sum
    is not a derivative. It would look like a plausible number rather than an
    error, so refuse instead.
    """
    if not getattr(sim, "eps_averaging", True):
        raise ValueError(
            "Geometry gradients need subpixel smoothing, but this simulation "
            "has eps_averaging=False. Without it the permittivity is a step "
            "function of an object's position and the derivative does not "
            "exist at the grid scale."
        )


def install_geometry_monitors(
    sim: mp.Simulation,
    objects,
    frequencies,
    decimation_factor: int = 0,
    pade_samples: int = 0,
) -> List[List[mp.DftFields]]:
    """DFT monitors over each differentiable object's support.

    Unlike a design region, an ordinary object has no monitor of its own, so the
    fields it needs have to be recorded in both the forward and adjoint runs.
    """
    from . import utils

    monitors = []
    for _, obj in objects:
        # Pad beyond the object. Moving a boundary changes the smoothed
        # permittivity in every voxel the boundary passes through, and that
        # includes voxels whose centres lie *outside* the original extent. A
        # monitor covering only the object misses them, which costs little
        # where the field is continuous across the boundary and a great deal
        # where it is not -- the normal E field jumps by the index contrast.
        pad = MONITOR_PAD_PIXELS / sim.resolution

        # `mp.inf` is how one writes "spans the cell", and is the natural way to
        # describe a layer -- but it is 1e20, not a flag, and
        # `_fit_volume_to_simulation` passes it straight through. A DFT volume
        # of that extent fails inside meep with "impossible(?) looping
        # boundaries", so clamp to the cell before fitting.
        cell = sim.cell_size
        padded = mp.Vector3(
            *(
                min(s + 2 * pad, c) if s else 0.0
                for s, c in zip(
                    (obj.size.x, obj.size.y, obj.size.z), (cell.x, cell.y, cell.z)
                )
            )
        )
        volume = sim._fit_volume_to_simulation(
            mp.Volume(center=obj.center, size=padded)
        )
        monitors.append(
            [
                sim.add_dft_fields(
                    [component],
                    frequencies,
                    where=volume,
                    yee_grid=True,
                    decimation_factor=decimation_factor,
                    persist=True,
                    pade_samples=pade_samples,
                )
                for component in utils._compute_components(sim)
            ]
        )
    return monitors


def gradient(
    sim: mp.Simulation,
    obj,
    object_index: int,
    forward_fields: List[mp.DftFields],
    adjoint_fields: List[mp.DftFields],
    frequencies,
    step: Optional[float] = None,
) -> dict:
    """dJ/d(parameter) for one object, keyed by the names it declared."""
    _check_supported(obj)
    check_smoothing(sim)
    check_faces_off_voxel_edges(sim, obj)

    names = [n for n in obj.differentiable]
    indices = []
    for name in names:
        indices.extend(PARAMETER_INDEX[name])
    indices = np.asarray(indices, dtype=np.intc)

    frequencies = np.asarray(frequencies, dtype=np.float64)
    out = np.zeros((frequencies.size, indices.size), dtype=np.float64)

    if step is None:
        step = DEFAULT_STEP_PIXELS / sim.resolution

    mp._get_geometry_gradient(
        out,
        1.0,
        adjoint_fields[0].swigobj,
        adjoint_fields[1].swigobj,
        adjoint_fields[2].swigobj,
        forward_fields[0].swigobj,
        forward_fields[1].swigobj,
        forward_fields[2].swigobj,
        sim.gv,
        frequencies,
        sim.geps,
        object_index,
        indices,
        step,
    )

    grads = {}
    for position, name in enumerate(names):
        columns = slice(3 * position, 3 * position + 3)
        grads[name] = np.squeeze(out[:, columns])
    return grads
