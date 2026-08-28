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
"""

from typing import List, Optional

import numpy as np

import meep as mp

# Perturbation for the geometric finite difference, as a fraction of a pixel.
# An absolute step is the wrong idea: subpixel smoothing varies over a pixel, so
# the step has to be defined against the grid. Too small and the smoothed
# permittivity has not changed measurably; too large and the difference samples
# curvature rather than a derivative.
DEFAULT_STEP_PIXELS = 0.02

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
        padded = mp.Vector3(
            obj.size.x + 2 * pad if obj.size.x else 0.0,
            obj.size.y + 2 * pad if obj.size.y else 0.0,
            obj.size.z + 2 * pad if obj.size.z else 0.0,
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
