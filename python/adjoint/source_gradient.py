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
IMPLEMENTED_PARAMS = ("currents", "amplitude")


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
    ):
        self.source = source
        self.component = source.component
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
                "`amp_func_file` is not supported; supply the amplitudes as an "
                "array (`amp_data`) or drive the source from JAX instead."
            )

        self.volume = sim._fit_volume_to_simulation(
            mp.Volume(center=source.center, size=source.size)
        )
        self._check_not_in_pml(sim)
        self.decimation_factor = decimation_factor
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
        center = np.array([self.volume.center.x, self.volume.center.y, self.volume.center.z])
        half_size = np.array([self.volume.size.x, self.volume.size.y, self.volume.size.z]) / 2

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
        )

    def shape(self, sim: mp.Simulation):
        """(nfreq,) + the monitor's spatial shape, with trailing 1s dropped."""
        dims = sim.fields.dft_monitor_size(
            self._monitor.swigobj, self.volume.swigobj, self.component
        )
        dims = [d for d in dims if d > 1] or [1]
        return (len(self._frequencies), *dims)

    def gather(self, sim: mp.Simulation, scale: np.ndarray) -> np.ndarray:
        """Return dJ_obj/d(currents), summed over all processes.

        `scale` is the per-frequency factor relating a current amplitude to the
        adjoint field, broadcast over the spatial axes.
        """
        dims = sim.fields.dft_monitor_size(
            self._monitor.swigobj, self.volume.swigobj, self.component
        )
        num_points = int(np.prod(dims))
        grad = np.zeros(num_points * len(self._frequencies), dtype=np.complex128)

        self._monitor.swigobj.fourier_sourcegradient(
            self.volume.swigobj, self.component, sim.fields, grad
        )

        grad = grad.reshape(len(self._frequencies), num_points)
        grad *= np.asarray(scale).reshape(-1, 1)
        return grad.reshape(self.shape(sim))


def install_source_gradient_monitors(
    sim: mp.Simulation,
    sources: List,
    frequencies: np.ndarray,
    decimation_factor: Optional[int] = 0,
) -> List[SourceGradientMonitor]:
    """Install a DFT monitor over each differentiable source's support."""
    monitors = [
        SourceGradientMonitor(sim, s, frequencies, decimation_factor) for s in sources
    ]
    for m in monitors:
        m.register(sim)
    return monitors


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

    num_dims = sim._infer_dimensions(sim.k_point)
    dV = 1 / sim.resolution**num_dims

    fwd_dtft = time_profile_dtft(sim, source.src, frequencies)

    scale = np.asarray(adj_src_phase) * fwd_dtft / (dV * iomega)

    if sim.using_real_fields():
        # real fields keep only Re[J], halving the amplitude at +omega
        scale *= 2
    return scale


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
