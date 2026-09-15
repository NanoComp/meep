"""Adjoint gradients with respect to source amplitudes.

Meep's design gradient contracts the forward and adjoint fields inside a design
region.  The source gradient is simpler: writing Maxwell's equations as
`A(rho) E = -i omega J`, we have `dE/dJ = -i omega A^-1`, so

    dJ_obj/dJ = (dJ_obj/dE) (dE/dJ) = -i omega lambda

which is just the adjoint field sampled where the source sits.  No additional
simulation is needed -- only a DFT monitor over the source's support during the
adjoint run that is already being performed.

The gathering is done by `dft_fields::fourier_sourcegradient`, the transpose of
the `fourier_sourcedata` scatter that places adjoint sources.  Writing the two
as an exact transpose pair is what makes the array ordering automatic: the
cotangent comes back indexed exactly like the array that would be handed in.

Meep's contract stops at the discrete current amplitudes.  Turning `dJ_obj/dA`
into a derivative with respect to a beam waist or a density is the caller's
job -- either JAX's, via `MeepJaxWrapper`, or the user's.
"""

from typing import List, Optional

import numpy as np

import meep as mp

# Parameters this module can currently evaluate.  `source.differentiable` is
# validated against a wider list at construction time; see
# `meep.source._validate_differentiable`.
IMPLEMENTED_PARAMS = (
    "currents",
    "amplitude",
    "amp_data",
    "beam_x0",
    "beam_kdir",
    "beam_w0",
    "beam_E0",
)


def differentiable_sources(sim: mp.Simulation) -> List:
    """Return the sources in `sim` that were flagged for differentiation."""
    return [s for s in sim.sources if getattr(s, "differentiable", ())]


def source_key(source, index: int):
    """The key a source's gradient appears under.

    Uses the source's `name` when it has one, and its position in
    `Simulation.sources` otherwise.  Names are preferred because they survive
    the source list being reordered, and because they are meaningful across
    processes.
    """
    name = getattr(source, "name", None)
    return name if name is not None else index


class SourceGradientMonitor:
    """A DFT monitor over one differentiable source's support.

    Installed during the adjoint run.  The adjoint field it records *is* the
    gradient with respect to the source's currents, up to the per-frequency
    scale factor supplied by the caller.
    """

    def __init__(
        self,
        sim: mp.Simulation,
        source,
        frequencies: np.ndarray,
        decimation_factor: Optional[int] = 0,
        component: Optional[int] = None,
        pade_samples: int = 0,
    ):
        self.source = source
        # A Gaussian beam has no single component: it places four, so each gets
        # its own monitor and the caller says which one this is.
        self.component = source.component if component is None else component
        self._frequencies = np.asarray(frequencies)

        if getattr(source, "amp_func", None) is not None or getattr(
            source, "amp_func_file", ""
        ):
            # An amp_func is evaluated inside Meep at each grid point, so the
            # discrete amplitudes it produces are not visible here.  The
            # currents cotangent would be correct but the user has no array to
            # apply it to, and `amplitude` would silently ignore the profile.
            raise NotImplementedError(
                "Differentiating a source defined by `amp_func` or "
                "`amp_func_file` is not supported: the profile is evaluated "
                "inside Meep, so there is no array for the caller to apply a "
                "cotangent to. Use `amp_data` instead."
            )

        self.volume = sim._fit_volume_to_simulation(
            mp.Volume(center=source.center, size=source.size)
        )
        self._check_not_in_pml(sim)
        self.decimation_factor = decimation_factor
        self.pade_samples = pade_samples
        self._monitor = None

    def _check_not_in_pml(self, sim: mp.Simulation) -> None:
        """Refuse a source whose adjoint field would be absorbed.

        Inside a PML the adjoint field is damped, so the gradient would come
        back finite, smooth, and wrong.  That is worth an error rather than a
        surprise.
        """
        if not sim.boundary_layers:
            return
        thickness = max(
            (getattr(bl, "thickness", 0.0) for bl in sim.boundary_layers), default=0.0
        )
        if thickness <= 0:
            return

        half_cell = np.array([sim.cell_size.x, sim.cell_size.y, sim.cell_size.z]) / 2
        center = np.array(
            [self.volume.center.x, self.volume.center.y, self.volume.center.z]
        )
        half_size = (
            np.array([self.volume.size.x, self.volume.size.y, self.volume.size.z]) / 2
        )

        for i, dim in enumerate("xyz"):
            if half_cell[i] == 0:
                continue  # not a simulated direction
            overshoot = (np.abs(center[i]) + half_size[i]) - (half_cell[i] - thickness)
            if overshoot > 0:
                raise ValueError(
                    f"Differentiable source extends {overshoot:.4g} into the PML "
                    f"along {dim}. The adjoint field is absorbed there, so the "
                    "gradient would be wrong without reporting an error. Move "
                    "the source away from the boundary."
                )

    def register(self, sim: mp.Simulation) -> None:
        """Install the DFT monitor. Called once per adjoint run."""
        # yee_grid=True so the monitor samples exactly the points a volume
        # source drives, rather than voxel centers.
        self._monitor = sim.add_dft_fields(
            [self.component],
            self._frequencies,
            where=self.volume,
            yee_grid=True,
            decimation_factor=self.decimation_factor,
            pade_samples=self.pade_samples,
        )

    def get_pade_error(self):
        """Delegate to the underlying DFT monitor.

        The adjoint stopping criterion is handed every monitor whose
        convergence the gradient depends on, including this one -- otherwise a
        run could stop before the source cotangent has settled. This class wraps
        a `DftFields` rather than being one, so forward the query.
        """
        if self._monitor is None:
            raise RuntimeError(
                "the monitor has not been registered yet; call register() first"
            )
        return self._monitor.get_pade_error()

    def num_points(self, sim: mp.Simulation) -> int:
        """How many grid points this source actually drives.

        The *unreduced* count. A zero-width source volume straddles two Yee
        planes of a magnetic component, each carrying half the current, and
        collapsing that direction would merge them -- preserving the total
        weight but losing how the amplitude varies between the two, which is
        exactly what distinguishes a beam from a uniform sheet.
        """
        dims = sim.fields.dft_monitor_full_size(
            self._monitor.swigobj, self.volume.swigobj, self.component
        )
        return int(np.prod(dims))

    def shape(self, sim: mp.Simulation):
        """(nfreq, num points)."""
        return (len(self._frequencies), self.num_points(sim))

    def gather(self, sim: mp.Simulation, scale: np.ndarray) -> np.ndarray:
        """Return dJ_obj/d(currents), summed over all processes.

        `scale` is the per-frequency factor relating a current amplitude to the
        adjoint field, broadcast over the spatial axes.
        """
        num_points = self.num_points(sim)
        grad = np.zeros(num_points * len(self._frequencies), dtype=np.complex128)

        self._monitor.swigobj.volume_source_gradient(
            self.volume.swigobj, self.component, sim.fields, grad
        )

        grad = grad.reshape(len(self._frequencies), num_points)
        # volume_source_gradient negates nothing, while the calibration that
        # fixed `scale` was done against an electric source. Electric and
        # magnetic components therefore differ by a sign here -- getting this
        # wrong leaves each component individually plausible but silently
        # subtracts one contribution from the other where a source drives
        # several, as a Gaussian beam does.
        sign = -1.0 if mp.is_electric(self.component) else 1.0
        # No half-timestep phase correction here, despite step.cpp evaluating
        # magnetic sources at t and electric ones at t + dt/2. Applying one
        # breaks the single-component checks by exactly 1/cos(w dt/2), so
        # whatever absorbs that offset is already accounted for in the
        # calibrated adjoint source phase.
        grad *= sign * np.asarray(scale).reshape(-1, 1)
        return grad.reshape(self.shape(sim))

    def positions(self, sim: mp.Simulation) -> np.ndarray:
        """(num points, 3) positions matching `gather`'s ordering.

        Emitted from the same C++ loop as the cotangent rather than
        reconstructed here, so the two cannot disagree about which point is
        which -- which is exactly the kind of half-pixel mismatch that produces
        a smooth, plausible, wrong gradient.
        """
        num_points = self.num_points(sim)
        out = np.zeros(3 * num_points, dtype=np.float64)
        self._monitor.swigobj.monitor_positions(
            self.volume.swigobj, self.component, sim.fields, out
        )
        return out.reshape(num_points, 3)


def install_source_gradient_monitors(
    sim: mp.Simulation,
    sources: List,
    frequencies: np.ndarray,
    decimation_factor: Optional[int] = 0,
    pade_samples: int = 0,
) -> List[SourceGradientMonitor]:
    """Install a DFT monitor over each differentiable source's support."""
    monitors = []
    for source in sources:
        # one monitor per component the source actually drives; a Gaussian beam
        # places four
        group = [
            SourceGradientMonitor(
                sim,
                source,
                frequencies,
                decimation_factor,
                component=component,
                pade_samples=pade_samples,
            )
            for component in source_components(source, sim.dimensions)
        ]
        for m in group:
            m.register(sim)
        monitors.append(group)
    return monitors


def is_gaussian_beam(source) -> bool:
    return hasattr(source, "beam_w0") and hasattr(source, "beam_kdir")


def normal_index(source) -> int:
    """Which axis the source plane is normal to, as 0, 1 or 2."""
    size = [source.size.x, source.size.y, source.size.z]
    zeros = [i for i, s in enumerate(size) if s == 0]
    if not zeros:
        raise ValueError(
            "A Gaussian beam source must be a plane (a line in 2D), so one of "
            "its size components has to be zero."
        )
    return zeros[0]


def beam_places(source, component, dimensions: int) -> bool:
    """Whether `add_volume_source_check` actually places this component.

    It declines several, and contracting a sensitivity for a source that was
    never placed is silently wrong rather than an error:

      - components along the plane normal;
      - in 2D, whichever parity the beam does not excite.

    Mirrors sources.cpp:495.
    """
    normal = normal_index(source)
    axis = {mp.X: 0, mp.Y: 1, mp.Z: 2}[mp.component_direction(component)]
    if axis == normal:
        return False
    if dimensions == 2:
        e0 = source.beam_E0
        has_tm = abs(complex(e0.z)) > 0
        has_te = abs(complex(e0.x)) > 0 or abs(complex(e0.y)) > 0
        tm = component in (mp.Ez, mp.Hx, mp.Hy)
        if has_te and tm:
            return False
        if has_tm and not tm:
            return False
    return True


def source_components(source, dimensions: int = 3):
    """Which field components a source actually drives."""
    if is_gaussian_beam(source):
        return [
            component
            for component, _, _ in beam_component_map(normal_index(source))
            if beam_places(source, component, dimensions)
        ]
    return [source.component]


def time_profile_dtft(
    sim: mp.Simulation, src_time, frequencies: np.ndarray
) -> np.ndarray:
    """DTFT of a source's time profile, in the convention `_adj_src_scale` uses.

    Meep relates a source's *amplitude* to the field it produces through the
    Fourier transform of its time profile, so converting the adjoint field into
    a derivative with respect to that amplitude requires this factor.
    """
    dt = sim.fields.dt
    T = sim.meep_time()
    y = np.array([src_time.swigobj.current(t, dt) for t in np.arange(0, T, dt)])
    freqs = np.asarray(frequencies)
    phase = np.exp(1j * 2 * np.pi * freqs[:, np.newaxis] * np.arange(y.size) * dt)
    return (phase @ y) * dt / np.sqrt(2 * np.pi)


def source_grad_scale(
    sim: mp.Simulation,
    source,
    frequencies: np.ndarray,
    adj_src_phase: np.ndarray,
) -> np.ndarray:
    """Per-frequency factor relating the adjoint field to dJ_obj/d(amplitude).

    This is the transpose of `ObjectiveQuantity._adj_src_scale`, which maps a
    field cotangent onto an adjoint source amplitude.  Going the other way, the
    adjoint field at the source is multiplied by the *forward* source's own
    transform, divided by the volume element and the discrete-time `i*omega`,
    and multiplied by the same unit-modulus phase the adjoint source divided
    out.

    The returned gradient follows the convention JAX and autograd use for a
    real function of a complex input, namely
    `g = dJ/d(Re a) - i dJ/d(Im a)`, so that it chains without adjustment.
    """
    frequencies = np.asarray(frequencies)
    dt = sim.fields.dt

    # discrete-time derivative, matching _adj_src_scale
    iomega = (1.0 - np.exp(-1j * (2 * np.pi * frequencies) * dt)) * (1.0 / dt)

    # src_vol_chunkloop multiplies the amplitude by gv.a once per *zero-width*
    # direction, to keep the integrated current fixed as a delta function is
    # resolved (`data.amp *= gv.a` in sources.cpp). That is a factor of the
    # resolution per delta direction, not per dimension: a point source in 2D
    # carries a^2 and a line source a^1. Using 1/dV here instead would be right
    # for the point and wrong by a factor of the resolution for the line.
    num_dims = sim._infer_dimensions(sim.k_point)
    size = [source.size.x, source.size.y, source.size.z][:num_dims]
    num_delta = sum(1 for extent in size if extent == 0)

    fwd_dtft = time_profile_dtft(sim, source.src, frequencies)

    # The component-dependent sign is applied in `SourceGradientMonitor.gather`,
    # since it differs between electric and magnetic components and this scale
    # is shared across all the components one source drives.
    scale = np.asarray(adj_src_phase) * fwd_dtft * sim.resolution**num_delta / iomega

    if sim.using_real_fields():
        # real fields keep only Re[J], halving the amplitude at +omega
        scale *= 2
    return scale


# The four component sources `fields::add_volume_source(src, where, beam)`
# places, as (source component, amplitude sign, which beam field is evaluated).
# With n the index of the plane normal and np1/np2 the two tangential axes:
#   K = n x H  goes on the electric components, N = -n x E on the magnetic ones.
_E_COMPONENTS = (mp.Ex, mp.Ey, mp.Ez)
_H_COMPONENTS = (mp.Hx, mp.Hy, mp.Hz)


def beam_component_map(normal_index: int):
    """Which sources a Gaussian beam places, mirroring sources.cpp:526."""
    np1 = (normal_index + 1) % 3
    np2 = (normal_index + 2) % 3
    return (
        (_E_COMPONENTS[np2], +1.0, _H_COMPONENTS[np1]),
        (_E_COMPONENTS[np1], -1.0, _H_COMPONENTS[np2]),
        (_H_COMPONENTS[np2], -1.0, _E_COMPONENTS[np1]),
        (_H_COMPONENTS[np1], +1.0, _E_COMPONENTS[np2]),
    )


def _beam_at(sim, source, positions, overrides=None, center=None):
    """Evaluate a Gaussian beam's six field components at each position.

    Rebuilt from the source's own parameters so that a perturbed copy can be
    evaluated without touching the simulation, which is what makes the
    parameter derivatives cost no FDTD runs at all.
    """
    values = {
        "beam_x0": source.beam_x0,
        "beam_kdir": source.beam_kdir,
        "beam_w0": source.beam_w0,
        "beam_E0": source.beam_E0,
    }
    values.update(overrides or {})

    dims, cyl = sim.dimensions, sim.is_cylindrical
    # add_volume_source measures positions from the centre of the volume Meep
    # actually used, which grid snapping can move off source.center by up to
    # half a pixel. Using the wrong one shifts every sensitivity.
    origin = source.center if center is None else center
    beam = mp.gaussianbeam(
        mp.py_v3_to_vec(dims, values["beam_x0"], cyl),
        mp.py_v3_to_vec(dims, values["beam_kdir"], cyl),
        float(values["beam_w0"]),
        source.src.swigobj.frequency().real,
        sim.fields.get_eps(mp.py_v3_to_vec(dims, source.center, cyl)).real,
        sim.fields.get_mu(mp.py_v3_to_vec(dims, source.center, cyl)).real,
        np.array(
            [values["beam_E0"].x, values["beam_E0"].y, values["beam_E0"].z],
            dtype=np.complex128,
        ),
    )

    out = np.zeros((len(positions), 6), dtype=np.complex128)
    buffer = np.zeros(6, dtype=np.complex128)
    for i, position in enumerate(positions):
        # gaussianbeam_ampfunc is handed the position relative to the source
        # volume's center, so the same offset has to be applied here.
        relative = mp.Vector3(*position) - origin
        beam.get_fields(buffer, mp.py_v3_to_vec(dims, relative, cyl))
        out[i] = buffer
    return out


def _perturbations(name, value, step):
    """The plus and minus variations of one beam parameter.

    Vector parameters are varied one component at a time, so the derivative
    comes back with the shape of the parameter.
    """
    if name == "beam_w0":
        yield (), value + step, value - step
        return
    for axis, letter in enumerate("xyz"):
        delta = mp.Vector3(**{letter: step})
        yield (axis,), value + delta, value - delta


def beam_parameter_gradients(
    sim: mp.Simulation,
    source,
    names,
    cotangents: dict,
    normal_index: int,
    step: float = 1e-6,
    center=None,
) -> dict:
    """Contract per-point cotangents onto a Gaussian beam's own parameters.

    The map from a beam's parameters to the currents it places is an analytic
    function Meep evaluates itself, so a central difference over it is both
    appropriate and cheap: it costs no FDTD runs, only re-evaluating the beam.
    Meep already takes the same approach one level up, finite-differencing the
    material grid inside `material_grids_addgradient`.

    Only `J^T lambda` is ever needed, so each directional derivative is
    contracted as soon as it is formed and no dense Jacobian is built. For a
    source plane with many points that matters.

    Args:
        cotangents: maps a source component to (positions, cotangent), where
            cotangent is dJ/d(amplitude) at each of those positions.
        normal_index: 0, 1 or 2 for a plane normal to x, y or z.
    """
    mapping = beam_component_map(normal_index)
    field_index = {c: i for i, c in enumerate(_E_COMPONENTS + _H_COMPONENTS)}

    out = {}
    for name in names:
        value = {
            "beam_x0": source.beam_x0,
            "beam_kdir": source.beam_kdir,
            "beam_w0": source.beam_w0,
            "beam_E0": source.beam_E0,
        }[name]

        entries = {}
        for key, plus, minus in _perturbations(name, value, step):
            total = 0.0 + 0.0j
            for component, sign, evaluated in mapping:
                # skipped when add_volume_source_check declined to place it
                if component not in cotangents:
                    continue
                positions, cotangent = cotangents[component]
                if not len(positions):
                    continue
                column = field_index[evaluated]
                high = _beam_at(sim, source, positions, {name: plus}, center)[:, column]
                low = _beam_at(sim, source, positions, {name: minus}, center)[:, column]
                sensitivity = sign * (high - low) / (2 * step)
                total += np.sum(np.asarray(cotangent).ravel() * sensitivity)
            entries[key] = total

        if name == "beam_w0":
            out[name] = entries[()]
        else:
            gradient = np.array([entries[(axis,)] for axis in range(3)])
            if name == "beam_kdir":
                # Only the direction of beam_kdir is meaningful -- its length is
                # ignored -- so the component along it is not a derivative of
                # anything. Projecting it out leaves the part that is.
                axis = np.array([value.x, value.y, value.z], dtype=float)
                norm = np.linalg.norm(axis)
                if norm > 0:
                    axis = axis / norm
                    gradient = gradient - axis * np.dot(axis, gradient)
            out[name] = gradient
    return out


def _mirrorindex(i: np.ndarray, n: int) -> np.ndarray:
    """`mirrorindex` from fields.cpp:755."""
    return np.where(i >= n, 2 * n - 1 - i, np.where(i < 0, -1 - i, i))


def _map_coordinates(r: np.ndarray, n: int):
    """`map_coordinates` from fields.cpp:759, vectorized over points.

    Returns the two bracketing indices and the weight of the second, matching
    the `do_fabs` behaviour `linear_interpolate` relies on.
    """
    if n == 1:
        zero = np.zeros(r.shape, dtype=int)
        return zero, zero, np.zeros(r.shape)
    r = np.where(r < 0.0, -r, np.where(r > 1.0, 1.0 - r, r))
    i1 = _mirrorindex(np.floor(r * n).astype(int), n)
    d = r * n - i1 - 0.5
    i2 = _mirrorindex(np.where(d >= 0.0, i1 + 1, i1 - 1), n)
    return i1, i2, np.abs(d)


def amp_data_gradient(
    sim: mp.Simulation,
    source,
    positions: np.ndarray,
    cotangent: np.ndarray,
) -> np.ndarray:
    """Contract a per-point cotangent onto a source's `amp_data` array.

    `Source(amp_data=...)` reaches the grid through `amp_file_func`, which
    trilinearly interpolates the user's array at each grid point. This is the
    transpose of that interpolation: each point's cotangent is scattered back
    onto the (at most) eight array entries that fed it, with the same weights.

    It is a plain scatter rather than anything new in C++, because the hard
    part -- the cotangent with respect to the amplitude Meep actually applied
    at each grid point, together with where that point is -- is already done by
    `volume_source_gradient`.
    """
    data = np.asarray(source.amp_data)
    shape = tuple(data.shape) + (1,) * (3 - data.ndim)
    nx, ny, nz = shape

    positions = np.asarray(positions, dtype=float)
    relative = positions - np.array(
        [source.center.x, source.center.y, source.center.z], dtype=float
    )
    extent = np.array([source.size.x, source.size.y, source.size.z], dtype=float)
    # `amp_file_func` maps a position onto [0, 1] across the source volume, and
    # pins the coordinate to the centre in any direction of zero extent.
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(
            extent > 0, 0.5 + relative / np.where(extent > 0, extent, 1.0), 0.0
        )

    ix1, ix2, wx = _map_coordinates(r[:, 0], nx)
    iy1, iy2, wy = _map_coordinates(r[:, 1], ny)
    iz1, iz2, wz = _map_coordinates(r[:, 2], nz)

    grad = np.zeros(nx * ny * nz, dtype=np.complex128)
    values = np.asarray(cotangent).ravel()
    for ix, fx in ((ix1, 1.0 - wx), (ix2, wx)):
        for iy, fy in ((iy1, 1.0 - wy), (iy2, wy)):
            for iz, fz in ((iz1, 1.0 - wz), (iz2, wz)):
                flat = (ix * ny + iy) * nz + iz
                np.add.at(grad, flat, values * fx * fy * fz)
    return grad.reshape(data.shape)


def contract(source, currents_grad: np.ndarray, source_amplitudes=None) -> dict:
    """Contract the currents cotangent onto the source's declared parameters.

    Every requested parameter is a linear function of the currents, so no
    finite difference over Meep's source construction is involved here.  The
    nonlinear built-in parameters (`beam_w0` and friends) are rejected at
    construction time; see `meep.source._validate_differentiable`.
    """
    out = {}
    for name in source.differentiable:
        if name == "currents":
            out[name] = currents_grad
        elif name == "amplitude":
            # `amplitude` scales every point identically, so its derivative is
            # the sum of the currents cotangent against the unit profile. Exact,
            # with no finite difference.
            if source_amplitudes is None:
                profile = np.ones(currents_grad.shape[1:])
            else:
                profile = np.asarray(source_amplitudes) / source.amplitude
            out[name] = np.sum(
                currents_grad * np.conj(profile)[np.newaxis, ...],
                axis=tuple(range(1, currents_grad.ndim)),
            )
        else:  # pragma: no cover - blocked by _validate_differentiable
            raise NotImplementedError(
                f"'{name}' is not implemented in this version of the source "
                "adjoint; see meep.source._validate_differentiable."
            )
    return out
