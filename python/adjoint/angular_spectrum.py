"""Differentiable angular-spectrum propagation through stratified media.

Meep's near-to-far transformation requires the near surface to sit in a
homogeneous medium -- `dft_near2far` aborts otherwise -- so a structure whose
radiation crosses a material interface, such as a grating coupler radiating up
through a cladding into air, has to keep that interface inside the FDTD cell.

This module propagates the tangential DFT fields on a planar monitor through an
arbitrary layer stack analytically instead, in JAX. The layers above the monitor
leave the simulation, which both shrinks the cell and makes the stack itself a
differentiable design parameter: thicknesses, indices, and the position and angle
of whatever the field is finally coupled into are all ordinary JAX values.

The method is the plane-wave decomposition of the measured fields. On a plane in
a homogeneous layer, the tangential `E` and `H` together determine the up- and
down-going plane-wave amplitudes at every transverse wavevector; each amplitude
is then carried through the stack by a scattering-matrix recursion and
recombined. Recording both `E` and `H` is what makes the up/down split possible,
and is why this is strictly better posed than an open near-to-far surface, which
has no way to reject radiation heading the wrong way.

Two entry points, neither requiring the other:

    # Post-processing an ordinary forward run. NumPy in, NumPy out.
    monitor = sim.add_dft_fields([mp.Ez, mp.Hx], frequencies, where=plane)
    sim.run(...)
    asp = mpa.AngularSpectrum.from_monitor(sim, monitor, stack)
    efficiency = asp.overlap_monitor(sim, monitor, fiber_mode, distance=300.0)

    # Inside an objective function, differentiated through the adjoint.
    fields = asp.take(args)
    efficiency = asp.overlap(fields, fiber_mode, distance=300.0)

Only planar monitors normal to a coordinate axis are supported, and the monitor
must lie in a homogeneous region. Cylindrical coordinates are not supported.
"""

import math
from typing import Callable, Dict, NamedTuple, Optional, Sequence, Tuple, Union

import jax
import jax.numpy as jnp
import numpy as onp

import meep as mp

# Every layer is given an infinitesimal loss before its longitudinal wavevector
# is computed. `kz = sqrt(n^2 k0^2 - kt^2)` vanishes on the light line, where its
# derivative diverges: `jnp.sqrt(0.0)` evaluates to zero happily but its VJP is
# infinite, which becomes NaN and destroys the whole gradient while leaving the
# objective value looking perfect. A grid point lands exactly on the light line
# whenever the padded monitor width is an integer number of wavelengths in the
# medium, which is not a rare coincidence when cell sizes and wavelengths are
# chosen as round numbers.
DEFAULT_LOSS_REGULARIZATION = 1e-12

# Polarization indices used throughout. In 3D these are the s and p components
# relative to the plane of incidence; in 2D with an out-of-plane E field only S
# is populated, and with an in-plane E field only P.
S_POLARIZATION = 0
P_POLARIZATION = 1


class Layer(NamedTuple):
    """One layer of a stack.

    Attributes:
        index: the refractive index. A scalar, an array with one entry per
            frequency, or a callable mapping frequency to index.
        thickness: the thickness in Meep units, or None for the semi-infinite
            layer that terminates the stack.
    """

    index: Union[complex, float, onp.ndarray, Callable]
    thickness: Optional[float] = None

    @classmethod
    def from_medium(cls, medium: mp.Medium, thickness: Optional[float] = None):
        """Builds a layer from a `meep.Medium`, evaluating its index per frequency."""
        return cls(
            index=lambda frequency: onp.sqrt(complex(medium.epsilon(frequency)[0][0])),
            thickness=thickness,
        )

    def indices(self, frequencies: onp.ndarray) -> jnp.ndarray:
        """Returns the index at each frequency, as a (num frequencies,) array."""
        if callable(self.index):
            values = jnp.asarray([self.index(f) for f in frequencies])
        else:
            values = jnp.asarray(self.index)
        return jnp.broadcast_to(values, (len(frequencies),))


class Stack(NamedTuple):
    """The layers above (or below) the monitor, ordered away from it.

    The first layer is the one containing the monitor, and the last must be
    semi-infinite, i.e. have `thickness=None`.
    """

    layers: Sequence[Layer]
    loss_regularization: float = DEFAULT_LOSS_REGULARIZATION

    def validate(self) -> None:
        if len(self.layers) < 1:
            raise ValueError("A stack needs at least one layer.")
        if self.layers[-1].thickness is not None:
            raise ValueError(
                "The last layer of a stack terminates it and must be "
                "semi-infinite, i.e. constructed with thickness=None; got "
                f"{self.layers[-1].thickness}."
            )
        for i, layer in enumerate(self.layers[:-1]):
            if layer.thickness is None:
                raise ValueError(
                    f"Layer {i} of {len(self.layers)} has thickness=None, but "
                    "only the final layer may be semi-infinite."
                )


class TangentialFields(NamedTuple):
    """Tangential fields on the monitor plane.

    Attributes:
        E: maps a Meep field component to a (num frequencies,) + spatial array.
        H: likewise for the magnetic components.
        normal: the coordinate direction normal to the plane, `mp.X`, `mp.Y`, or
            `mp.Z`.
        sign: +1 when the radiation of interest travels along +normal, -1 for
            -normal. A monitor below the structure, propagating into the
            substrate, uses -1.
    """

    E: Dict[int, jnp.ndarray]
    H: Dict[int, jnp.ndarray]
    normal: int
    sign: int = 1


def _as_concrete(value) -> Optional[float]:
    """Returns `value` as a float, or None if it is a JAX tracer.

    Distances and thicknesses may be traced, so checks on them have to be
    skipped rather than forced.
    """
    try:
        return float(value)
    except (TypeError, jax.errors.ConcretizationTypeError):
        return None


def _safe_sqrt(argument: jnp.ndarray) -> jnp.ndarray:
    """Square root with the branch chosen so that the imaginary part is >= 0.

    Evanescent orders must decay away from the monitor rather than grow, which
    fixes the branch. The caller is expected to have regularized the argument so
    that it never lands exactly on zero; see DEFAULT_LOSS_REGULARIZATION.
    """
    root = jnp.sqrt(jnp.asarray(argument, dtype=jnp.complex128))
    return jnp.where(jnp.imag(root) < 0, -root, root)


def _longitudinal_wavevector(
    index: jnp.ndarray, k0: jnp.ndarray, kt: jnp.ndarray
) -> jnp.ndarray:
    """kz for each (frequency, transverse wavevector), shaped (nfreq, nkt)."""
    kt = jnp.asarray(kt)
    transverse = jnp.sum(jnp.square(kt.reshape(kt.shape[0], -1)), axis=-1)
    return _safe_sqrt(
        jnp.square(index)[:, None] * jnp.square(k0)[:, None] - transverse[None, :]
    )


def _admittances(
    index: jnp.ndarray, kz: jnp.ndarray, k0: jnp.ndarray
) -> Tuple[jnp.ndarray, jnp.ndarray]:
    """Wave admittances for s and p polarization.

    In Meep units the vacuum impedance is 1, so with omega = k0 these reduce to
    `Y_s = kz / k0` and `Y_p = n^2 k0 / kz`.
    """
    admittance_s = kz / k0[:, None]
    admittance_p = jnp.square(index)[:, None] * k0[:, None] / kz
    return admittance_s, admittance_p


def _interface_coefficients(
    admittance_in: jnp.ndarray, admittance_out: jnp.ndarray
) -> Tuple[jnp.ndarray, jnp.ndarray]:
    """Fresnel reflection and transmission at one interface, in admittances."""
    total = admittance_in + admittance_out
    reflection = (admittance_in - admittance_out) / total
    transmission = 2 * admittance_in / total
    return reflection, transmission


def _stack_transmission(
    admittances: Sequence[jnp.ndarray],
    wavevectors: Sequence[jnp.ndarray],
    thicknesses: Sequence[float],
) -> Tuple[jnp.ndarray, jnp.ndarray]:
    """Transmission and reflection of a stack, by scattering-matrix recursion.

    Accumulating a transfer matrix instead would overflow as soon as a layer is
    thick or an order is evanescent, since that formulation contains a growing
    exponential; the scattering recursion below only ever forms `exp(i kz d)`
    with non-negative imaginary `kz`, which decays.

    Args:
        admittances: the wave admittance in each layer, one (nfreq, nkt) array
            per layer, for a single polarization.
        wavevectors: the longitudinal wavevector in each layer, likewise.
        thicknesses: the thickness of each layer except the last.

    Returns:
        The amplitude transmission and reflection coefficients of the whole
        stack, referenced to the front face of the first layer.
    """
    # Start from the terminating interface and recurse toward the monitor. At
    # each step `reflection` is the reflection looking into the remaining stack
    # and `transmission` the accumulated transmission through it.
    reflection, transmission = _interface_coefficients(admittances[-2], admittances[-1])
    for j in range(len(admittances) - 3, -1, -1):
        phase = jnp.exp(1j * wavevectors[j + 1] * thicknesses[j + 1])
        interface_r, interface_t = _interface_coefficients(
            admittances[j], admittances[j + 1]
        )
        # Redheffer star product of the interface with the layer already
        # accumulated behind it.
        backward_r = -interface_r  # reflection of the interface from the far side
        denominator = 1 - backward_r * reflection * phase**2
        transmission = interface_t * phase * transmission / denominator
        reflection = interface_r + (
            interface_t * (2 - interface_t) * reflection * phase**2 / denominator
        )
    return transmission, reflection


def _single_interface_limit(
    admittances: Sequence[jnp.ndarray],
) -> Tuple[jnp.ndarray, jnp.ndarray]:
    """Transmission and reflection when the stack is a single interface."""
    reflection, transmission = _interface_coefficients(admittances[0], admittances[1])
    return transmission, reflection


class PropagationResult(NamedTuple):
    """The plane-wave amplitudes produced by a propagation.

    Attributes:
        kt: the transverse wavevectors, (num kt,) in 2D or (num kt, 2) in 3D.
        amplitudes: (num frequencies, num kt, 2) electric field amplitudes in
            the s and p basis, at the target plane.
        kz: (num frequencies, num kt) longitudinal wavevector in the output
            layer, used for flux and overlap integrals.
        index: (num frequencies,) index of the output layer.
    """

    kt: jnp.ndarray
    amplitudes: jnp.ndarray
    kz: jnp.ndarray
    index: jnp.ndarray


# Sign and normalization conventions, stated once because they are the usual
# source of silent errors.
#
# Meep uses exp(-i omega t), so Faraday's law reads curl E = i omega mu H. In
# Meep units mu = 1 and omega = k0. For a two-dimensional simulation in the x-y
# plane with a monitor line normal to y and an out-of-plane electric field, a
# plane wave exp(i(kx x + ky y)) therefore satisfies
#
#     H_x = ky E_z / k0
#
# so an up-going wave (ky = +kz) has H_x = +Y_s E_z and a down-going one has
# H_x = -Y_s E_z, with Y_s = kz / k0. That is what separates the two:
#
#     E_up = (E_z + H_x / Y_s) / 2      E_down = (E_z - H_x / Y_s) / 2
#
# Reflection and transmission coefficients are those of the *tangential* field
# components, which is the admittance convention. For p polarization this is the
# negative of the more familiar form written in terms of the full electric
# vector; the two agree on |r|, on the Brewster angle, and on energy
# conservation, and differ only in the sign of r_p.


# The tangential axes of a plane, ordered right-handed so that u_hat x v_hat =
# n_hat. That ordering is what makes the decomposition below sign-free:
# n_hat x H_t has components (-H_v, +H_u) in the same basis.
_PLANE_AXES = {
    mp.X: ((mp.Ey, mp.Ez), (mp.Hy, mp.Hz)),
    mp.Y: ((mp.Ez, mp.Ex), (mp.Hz, mp.Hx)),
    mp.Z: ((mp.Ex, mp.Ey), (mp.Hx, mp.Hy)),
}


def _tangential_components(normal: int, dimensions: int) -> Tuple[Tuple[int, ...], ...]:
    """The E and H components tangential to a plane with the given normal."""
    if normal not in _PLANE_AXES:
        raise ValueError(f"Unsupported monitor normal {normal}.")
    if dimensions == 2 and normal == mp.Z:
        raise ValueError(
            "In a 2D simulation the monitor plane must be normal to x or y, "
            f"but got normal={normal}."
        )
    if dimensions not in (2, 3):
        raise NotImplementedError(
            "Angular-spectrum propagation supports 2D and 3D Cartesian "
            f"simulations, not {dimensions}D. Cylindrical coordinates are not "
            "supported."
        )
    return _PLANE_AXES[normal]


class Mode(NamedTuple):
    """A target field to project onto, defined by its angular spectrum.

    Attributes:
        spectrum: called as `spectrum(propagator, kt, k0, index)` and returning
            a (num frequencies, num kt, 2) array of tangential electric field
            amplitudes in the s and p directions. Carrying both is what lets a
            linearly polarized beam be represented in 3D, where its s and p
            content varies with azimuth.
    """

    spectrum: Callable


def gaussian_mode(
    waist: float,
    tilt_deg: float = 0.0,
    offset: float = 0.0,
    polarization: int = S_POLARIZATION,
    tilt_azimuth_deg: float = 0.0,
) -> Mode:
    """A tilted, laterally offset, linearly polarized Gaussian, e.g. a fiber mode.

    The spectrum is written in closed form rather than sampled, so the tilt and
    the offset are exact continuous parameters and differentiable, instead of
    being quantized by the monitor grid.

    Args:
        waist: the 1/e field radius at the target plane, in Meep units. For a
            fiber quoted by mode-field diameter, this is MFD / 2.
        tilt_deg: the angle from the plane normal, in degrees.
        offset: the lateral displacement of the beam center at the target plane,
            along the tilt direction.
        polarization: the linear polarization direction. `S_POLARIZATION` means
            along the second tangential axis, `P_POLARIZATION` the first; in 2D
            these are the out-of-plane and in-plane electric fields
            respectively. In 3D the beam is linearly polarized along that axis
            and its s and p content follows the azimuth of each wavevector.
        tilt_azimuth_deg: which way the beam is tilted, measured in the
            tangential plane. Ignored in 2D, where there is only one direction
            to tilt in.

    Returns:
        A `Mode`.
    """

    def spectrum(propagator, kt, k0, index):
        kt = jnp.asarray(kt)
        transverse = (
            jnp.asarray(index).real
            * jnp.asarray(k0)
            * jnp.sin(jnp.deg2rad(jnp.asarray(tilt_deg)))
        )
        azimuth = jnp.deg2rad(jnp.asarray(tilt_azimuth_deg))
        direction = (
            jnp.array([1.0])
            if kt.shape[-1] == 1
            else jnp.stack([jnp.cos(azimuth), jnp.sin(azimuth)])
        )
        center = transverse[:, None] * direction
        detuning = kt[None, :, :] - center[:, None, :]
        envelope = jnp.exp(
            -jnp.sum(jnp.square(detuning * waist), axis=-1) / 4.0
        ) * jnp.exp(-1j * jnp.sum(detuning * direction * offset, axis=-1))

        rotation = propagator._azimuth()
        if rotation is None:
            weights = (
                jnp.array([1.0, 0.0])
                if polarization == S_POLARIZATION
                else jnp.array([0.0, 1.0])
            )
            components = jnp.broadcast_to(weights, envelope.shape + (2,))
        else:
            cosine, sine = rotation
            if polarization == S_POLARIZATION:
                # polarized along v_hat: v_hat . s_hat = cos, v_hat . p_hat = sin
                components = jnp.stack([cosine, sine], axis=-1)
            else:
                components = jnp.stack([-sine, cosine], axis=-1)
            components = jnp.broadcast_to(components[None, :, :], envelope.shape + (2,))
        return envelope[..., None] * components

    return Mode(spectrum=spectrum)


class AngularSpectrum:
    """Propagates monitor fields through a layer stack.

    Attributes:
        stack: the layers above the monitor, ordered away from it.
        frequencies: the monitor frequencies, in Meep units.
        pitch: the spacing of the monitor samples.
        num_points: the number of monitor samples.
        normal: the coordinate direction normal to the monitor plane.
        sign: +1 if the radiation of interest travels along +normal, else -1.
    """

    def __init__(
        self,
        stack: Stack,
        frequencies: Sequence[float],
        pitch: float,
        num_points: int,
        normal: int = mp.Y,
        sign: int = 1,
        pad_factor: int = 4,
        kt: Optional[onp.ndarray] = None,
    ):
        """Initializes a propagator.

        Args:
            stack: the layers above the monitor.
            frequencies: the monitor frequencies.
            pitch: the monitor sample spacing, normally 1 / resolution.
            num_points: the number of monitor samples.
            normal: the direction normal to the monitor plane.
            sign: the direction of the radiation of interest along that normal.
            pad_factor: how much to zero-pad before transforming. The transform
                is periodic, so a field that has not decayed by the edges of the
                padded window wraps around; padding pushes the replicas apart.
                Ignored when `kt` is given.
            kt: transverse wavevectors to evaluate, instead of the uniform grid
                a padded transform would produce. Evaluating a chosen set costs
                O(num_points * num kt) and sidesteps padding entirely, which is
                worthwhile when only the wavevectors within some numerical
                aperture matter.
        """
        stack.validate()
        if sign not in (1, -1):
            raise ValueError(f"sign must be +1 or -1, got {sign}.")

        # A monitor plane is a line in a 2D simulation and a rectangle in a 3D
        # one, so the transverse space is one- or two-dimensional. Everything
        # downstream of the transform works on a flat list of transverse
        # wavevectors and is indifferent to which it was.
        self.num_points = (
            (int(num_points),)
            if onp.ndim(num_points) == 0
            else tuple(int(n) for n in num_points)
        )
        self.pitch = (
            (float(pitch),) * len(self.num_points)
            if onp.ndim(pitch) == 0
            else tuple(float(p) for p in pitch)
        )
        if len(self.pitch) != len(self.num_points):
            raise ValueError(
                f"pitch has {len(self.pitch)} entries but num_points has "
                f"{len(self.num_points)}."
            )
        self.transverse_dimensions = len(self.num_points)
        _tangential_components(normal, self.transverse_dimensions + 1)

        self.stack = stack
        self.frequencies = onp.asarray(frequencies, dtype=float)
        self.normal = normal
        self.sign = int(sign)
        self.pad_factor = int(pad_factor)

        self._k0 = 2 * onp.pi * self.frequencies
        if kt is None:
            self._padded = tuple(n * self.pad_factor for n in self.num_points)
            axes = [
                2 * onp.pi * onp.fft.fftfreq(n, d=d)
                for n, d in zip(self._padded, self.pitch)
            ]
            grids = onp.meshgrid(*axes, indexing="ij")
            self._kt = jnp.asarray(onp.stack([g.ravel() for g in grids], axis=-1))
            self._uniform = True
        else:
            kt = onp.asarray(kt, dtype=float)
            self._kt = jnp.asarray(kt.reshape(kt.shape[0], -1))
            self._uniform = False
            self._padded = None

        # Regularizing the index keeps kz off the light line, where its
        # derivative is unbounded; see DEFAULT_LOSS_REGULARIZATION.
        regularizer = 1 + 1j * self.stack.loss_regularization
        self._indices = [
            layer.indices(self.frequencies) * regularizer for layer in stack.layers
        ]
        self._wavevectors = [
            _longitudinal_wavevector(index, self._k0, self._kt)
            for index in self._indices
        ]
        self._admittances = [
            _admittances(index, kz, self._k0)
            for index, kz in zip(self._indices, self._wavevectors)
        ]
        # Deliberately not coerced to float: a thickness may be a JAX tracer, so
        # that an antireflection coating or a cladding depth can be optimized
        # alongside the design.
        self._thicknesses = [
            0.0 if layer.thickness is None else layer.thickness
            for layer in stack.layers
        ]

    @property
    def kt(self) -> jnp.ndarray:
        """The transverse wavevectors at which the spectrum is evaluated."""
        return self._kt

    @property
    def stack_thickness(self):
        """The total thickness of the stack, excluding the semi-infinite layer.

        Not necessarily a float: layer thicknesses may be JAX values.
        """
        total = 0.0
        for thickness in self._thicknesses[:-1]:
            total = total + thickness
        return total

    def coordinates(self):
        """The transverse sample coordinates, centered on zero.

        A single array for a line monitor, a list of two for a rectangular one.
        """
        axes = [
            (onp.arange(n) - (n - 1) / 2) * d
            for n, d in zip(self.num_points, self.pitch)
        ]
        return axes[0] if self.transverse_dimensions == 1 else axes

    def _coordinate_axes(self):
        """`coordinates`, always as a list, one array per transverse axis."""
        axes = self.coordinates()
        return [axes] if self.transverse_dimensions == 1 else axes

    def _sample_positions(self) -> onp.ndarray:
        """Every sample position, flattened to (num samples, num transverse)."""
        grids = onp.meshgrid(*self._coordinate_axes(), indexing="ij")
        return onp.stack([g.ravel() for g in grids], axis=-1)

    def _transform(self, values: jnp.ndarray) -> jnp.ndarray:
        """Transforms monitor samples to (num frequencies, num kt)."""
        values = jnp.asarray(values)
        transverse = tuple(range(-self.transverse_dimensions, 0))
        measure = float(onp.prod(self.pitch))
        if not self._uniform:
            flat = values.reshape(values.shape[: transverse[0]] + (-1,))
            phase = jnp.exp(
                -1j * jnp.einsum("kd,xd->kx", self._kt, self._sample_positions())
            )
            return jnp.einsum("kx,...x->...k", phase, flat) * measure
        padded = jnp.zeros(
            values.shape[: transverse[0]] + self._padded, dtype=jnp.complex128
        )
        padded = padded.at[
            (Ellipsis,) + tuple(slice(0, n) for n in self.num_points)
        ].set(values)
        spectrum = jnp.fft.fftn(padded, axes=transverse) * measure
        spectrum = spectrum.reshape(spectrum.shape[: transverse[0]] + (-1,))
        # The samples are centered on zero, so undo the phase ramp implied by
        # having placed them at indices starting from zero.
        origins = jnp.asarray([axis[0] for axis in self._coordinate_axes()])
        return spectrum * jnp.exp(-1j * (self._kt @ origins))

    def _tangential_spectra(self, fields: TangentialFields):
        """Transforms the tangential fields, returning them on the (u, v) axes.

        Returns `(E_u, E_v, H_u, H_v)`, each (num frequencies, num kt), with a
        missing component treated as zero. In 2D exactly one of the two
        polarizations is populated and the other stays zero, which costs nothing
        and spares the caller declaring which one their source excites.
        """
        if fields.normal != self.normal:
            raise ValueError(
                f"The fields are on a plane normal to {fields.normal} but this "
                f"propagator was built for {self.normal}."
            )
        (e_u, e_v), (h_u, h_v) = _PLANE_AXES[self.normal]
        # E_u pairs with H_v, not H_u: the decomposition contracts E_t against
        # n_hat x H_t, which swaps the two tangential axes. For a y-normal plane
        # that means Ez goes with Hx and Ex with Hz.
        pairs = ((e_u, h_v), (e_v, h_u))
        present = [(fields.E.get(e), fields.H.get(h)) for e, h in pairs]
        if all(e is None and h is None for e, h in present):
            raise ValueError(
                "No tangential field components were supplied; nothing to " "propagate."
            )
        for (electric, magnetic), (e, h) in zip(present, pairs):
            if (electric is None) != (magnetic is None):
                raise ValueError(
                    "Separating up-going from down-going radiation needs both "
                    f"tangential fields, but only one of "
                    f"{mp.component_name(e)} and {mp.component_name(h)} was "
                    "supplied."
                )
        zero = None
        for electric, magnetic in present:
            if electric is not None:
                zero = jnp.zeros_like(self._transform(jnp.asarray(electric)))
                break
        # Returned as (E_u, H_v, E_v, H_u), matching how they pair up.
        return tuple(
            zero if value is None else self._transform(jnp.asarray(value))
            for pair in present
            for value in pair
        )

    def _azimuth(self):
        """cos and sin of the angle from u_hat to the transverse wavevector.

        At normal incidence the plane of incidence is undefined and any
        orthogonal pair will do, so the azimuth is pinned to zero there rather
        than left to produce a nan from atan2(0, 0).
        """
        if self.transverse_dimensions == 1:
            # The transverse wavevector lies along a single axis, so there is no
            # azimuth to speak of and the (u, v) axes are already the s and p
            # directions.
            return None
        magnitude = jnp.linalg.norm(self._kt, axis=-1)
        safe = jnp.where(magnitude > 0, magnitude, 1.0)
        cosine = jnp.where(magnitude > 0, self._kt[:, 0] / safe, 1.0)
        sine = jnp.where(magnitude > 0, self._kt[:, 1] / safe, 0.0)
        return cosine, sine

    def decompose(self, fields: TangentialFields):
        """Splits the monitor fields into up- and down-going spectra.

        The tangential fields are rotated into the s and p directions of each
        transverse wavevector, where the admittance is a scalar, and separated
        using

            E_up_s = (E_s - H_p / Y_s) / 2      E_up_p = (E_p + H_s / Y_p) / 2

        which is `E_up = (E_t - (n_hat x H_t) / Y) / 2` written out in that
        basis. In 2D the wavevector lies along one axis, so the rotation is the
        identity and (u, v) are already (s, p) up to the ordering below.

        Returns:
            `(up, down)`, each a dict mapping polarization to a
            (num frequencies, num kt) array of tangential electric field
            amplitudes at the monitor plane.
        """
        electric_u, magnetic_v, electric_v, magnetic_u = self._tangential_spectra(
            fields
        )
        azimuth = self._azimuth()
        admittance_s, admittance_p = self._admittances[0]
        # `fields.sign` flips which branch counts as outgoing, for a monitor
        # facing down into a substrate.
        outgoing = fields.sign

        if azimuth is None:
            # One transverse direction, so there is no azimuth and the
            # right-handed axes are already a valid s/p pair: u_hat is
            # perpendicular to the wavevector, v_hat lies along it. Using
            # (n_hat x H)_u = -H_v and (n_hat x H)_v = +H_u directly,
            #
            #     E_up_u = (E_u + H_v / Y_s) / 2   E_up_v = (E_v - H_u / Y_p) / 2
            #
            # Note the signs are *not* those of the rotated case below: the 3D
            # convention puts s_hat along -u_hat here, since the wavevector runs
            # along v_hat, and the two bases therefore differ by a sign.
            electric_s, magnetic_s = electric_u, magnetic_u
            electric_p, magnetic_p = electric_v, magnetic_v
            cross_s = -outgoing * magnetic_p / admittance_s
            cross_p = outgoing * magnetic_s / admittance_p
        else:
            cosine, sine = azimuth
            electric_s = -electric_u * sine + electric_v * cosine
            electric_p = electric_u * cosine + electric_v * sine
            magnetic_s = -magnetic_u * sine + magnetic_v * cosine
            magnetic_p = magnetic_u * cosine + magnetic_v * sine
            cross_s = outgoing * magnetic_p / admittance_s
            cross_p = -outgoing * magnetic_s / admittance_p

        up = {
            S_POLARIZATION: 0.5 * (electric_s - cross_s),
            P_POLARIZATION: 0.5 * (electric_p - cross_p),
        }
        down = {
            S_POLARIZATION: 0.5 * (electric_s + cross_s),
            P_POLARIZATION: 0.5 * (electric_p + cross_p),
        }
        return up, down

    def _transmission(self, polarization: int):
        """Transmission and reflection of the stack for one polarization."""
        admittances = [pair[polarization] for pair in self._admittances]
        if len(admittances) == 2:
            return _single_interface_limit(admittances)
        return _stack_transmission(admittances, self._wavevectors, self._thicknesses)

    def spectrum(self, fields: TangentialFields, distance: float):
        """The up-going spectrum carried to a plane `distance` from the monitor.

        Args:
            fields: the monitor fields.
            distance: how far along the outgoing normal to propagate, measured
                from the monitor plane. It must reach at least to the far side of
                the stack.

        Returns:
            A `PropagationResult`.
        """
        remaining = distance - self.stack_thickness
        concrete = _as_concrete(remaining)
        if concrete is not None and concrete < 0:
            raise ValueError(
                f"distance={distance} does not clear the stack, which is "
                f"{_as_concrete(self.stack_thickness)} thick. The target plane "
                "has to lie in the semi-infinite layer."
            )
        up, _ = self.decompose(fields)
        amplitudes = {}
        for polarization, amplitude in up.items():
            transmission, _ = self._transmission(polarization)
            # Monitor to the first interface, through the stack, then onward in
            # the terminating layer.
            to_interface = jnp.exp(1j * self._wavevectors[0] * self._thicknesses[0])
            beyond = jnp.exp(1j * self._wavevectors[-1] * remaining)
            amplitudes[polarization] = amplitude * to_interface * transmission * beyond
        stacked = jnp.stack(
            [
                amplitudes.get(
                    polarization,
                    jnp.zeros_like(next(iter(amplitudes.values()))),
                )
                for polarization in (S_POLARIZATION, P_POLARIZATION)
            ],
            axis=-1,
        )
        return PropagationResult(
            kt=self._kt,
            amplitudes=stacked,
            kz=self._wavevectors[-1],
            index=self._indices[-1],
        )

    def _weights(self, result: PropagationResult) -> jnp.ndarray:
        """Per-wavevector power weights, (num frequencies, num kt, 2).

        Power carried along the normal by a plane-wave component is
        `Re(Y) |E|^2 / 2`, so this is the measure that makes overlaps unitary.
        Evanescent components carry none and drop out.
        """
        admittance_s, admittance_p = _admittances(
            self._indices[-1], result.kz, self._k0
        )
        measure = self._spectral_measure()
        return (
            0.5
            * measure
            * jnp.stack([jnp.real(admittance_s), jnp.real(admittance_p)], axis=-1)
        )

    def _spectral_measure(self) -> float:
        """The dk / 2pi factor that turns a spectral sum into a real-space integral."""
        if self._uniform:
            measure = 1.0
            for padded, pitch in zip(self._padded, self.pitch):
                measure *= (2 * onp.pi / (padded * pitch)) / (2 * onp.pi)
            return measure
        # A chosen set of wavevectors is assumed uniform along each axis;
        # trapezoid weights would be more careful but the set is normally a grid.
        measure = 1.0
        for axis in range(self._kt.shape[-1]):
            values = onp.unique(onp.asarray(self._kt[:, axis]))
            spacing = onp.mean(onp.diff(values)) if values.size > 1 else 1.0
            measure *= spacing / (2 * onp.pi)
        return measure

    def power(self, fields: TangentialFields, distance: Optional[float] = None):
        """Outgoing power through the target plane, one value per frequency."""
        result = self.spectrum(
            fields, self.stack_thickness if distance is None else distance
        )
        weights = self._weights(result)
        return jnp.sum(weights * jnp.abs(result.amplitudes) ** 2, axis=(1, 2))

    def overlap(
        self,
        fields: TangentialFields,
        mode: Mode,
        distance: float,
        incident_power: Optional[jnp.ndarray] = None,
    ):
        """Fraction of the incident power coupled into `mode`.

        The propagated field is projected onto the normalized mode using the
        power inner product, so the result is bounded by one and equals one when
        the field is a pure multiple of the mode and all of the incident power
        reaches the target plane.

        Args:
            fields: the monitor fields.
            mode: the target, e.g. `gaussian_mode(...)`.
            distance: distance from the monitor to the target plane.
            incident_power: the power to normalize against, one value per
                frequency. Defaults to the power crossing the target plane, in
                which case the result is the modal purity of the radiation
                rather than an end-to-end efficiency.

        Returns:
            A (num frequencies,) array.
        """
        result = self.spectrum(fields, distance)
        weights = self._weights(result)
        field_amplitude = result.amplitudes
        mode_amplitude = mode.spectrum(self, result.kt, self._k0, result.index)

        axes = (1, 2)
        cross = jnp.sum(weights * field_amplitude * jnp.conj(mode_amplitude), axis=axes)
        mode_norm = jnp.sum(weights * jnp.abs(mode_amplitude) ** 2, axis=axes)
        coupled = jnp.abs(cross) ** 2 / mode_norm
        if incident_power is None:
            incident_power = jnp.sum(weights * jnp.abs(field_amplitude) ** 2, axis=axes)
        return coupled / incident_power

    def report(self, fields: TangentialFields) -> Dict[str, jnp.ndarray]:
        """Diagnostics for the four ways a monitor is usually placed wrongly.

        Returns a dict with, per frequency:

        `downgoing_fraction`
            power heading back toward the structure. Large means something above
            the monitor is scattering, or the monitor sits inside the near field.
        `evanescent_fraction`
            power above the light line, which does not propagate. Large means
            the monitor is too close to the structure.
        `edge_amplitude`
            field magnitude at the ends of the monitor relative to its peak. A
            transform is periodic, so a field that has not decayed by the edges
            wraps around; large means the monitor is too narrow.
        """
        up, down = self.decompose(fields)
        admittance = self._admittances[0]
        measure = self._spectral_measure()
        propagating = jnp.abs(jnp.imag(self._wavevectors[0])) < 1e-8 * jnp.abs(
            jnp.real(self._wavevectors[0]) + 1e-30
        )

        def spectral_power(spectra, mask=None):
            total = 0.0
            for polarization, amplitude in spectra.items():
                weight = 0.5 * measure * jnp.real(admittance[polarization])
                weighted = weight * jnp.abs(amplitude) ** 2
                if mask is not None:
                    weighted = jnp.where(mask, weighted, 0.0)
                total = total + jnp.sum(weighted, axis=-1)
            return total

        up_power = spectral_power(up)
        down_power = spectral_power(down)
        up_propagating = spectral_power(up, propagating)

        edges = []
        for values in list(fields.E.values()):
            magnitude = jnp.abs(jnp.asarray(values))
            axes = tuple(range(-self.transverse_dimensions, 0))
            peak = jnp.max(magnitude, axis=axes)
            border = jnp.zeros_like(peak)
            for axis in axes:
                border = jnp.maximum(
                    border,
                    jnp.max(
                        jnp.take(magnitude, jnp.array([0, -1]), axis=axis),
                        axis=axes,
                    ),
                )
            edges.append(border / jnp.where(peak > 0, peak, 1.0))

        total = up_power + down_power
        return {
            "downgoing_fraction": down_power / jnp.where(total > 0, total, 1.0),
            "evanescent_fraction": 1.0
            - up_propagating / jnp.where(up_power > 0, up_power, 1.0),
            "edge_amplitude": jnp.max(jnp.stack(edges, axis=0), axis=0),
        }

    def propagate(self, fields: TangentialFields, distance: float, coordinates=None):
        """The tangential electric field at the target plane, in real space.

        Args:
            fields: the monitor fields.
            distance: distance from the monitor to the target plane.
            coordinates: where to evaluate, as one array per transverse
                dimension. Defaults to the monitor's own coordinates; pass a
                wider range to see a beam that has spread.

        Returns:
            A dict mapping each tangential electric component to an array of
            shape (num frequencies,) + the coordinate shape.
        """
        result = self.spectrum(fields, distance)
        if coordinates is None:
            axes = self._coordinate_axes()
        else:
            axes = (
                [onp.asarray(coordinates)]
                if self.transverse_dimensions == 1
                else [onp.asarray(a) for a in coordinates]
            )
        grids = onp.meshgrid(*axes, indexing="ij")
        positions = jnp.asarray(onp.stack([g.ravel() for g in grids], axis=-1))

        rotation = self._azimuth()
        amplitude_s = result.amplitudes[..., S_POLARIZATION]
        amplitude_p = result.amplitudes[..., P_POLARIZATION]
        if rotation is None:
            amplitude_u, amplitude_v = amplitude_s, amplitude_p
        else:
            cosine, sine = rotation
            # The rotation into (s, p) is a reflection, hence its own inverse.
            amplitude_u = -amplitude_s * sine + amplitude_p * cosine
            amplitude_v = amplitude_s * cosine + amplitude_p * sine

        phase = jnp.exp(1j * jnp.einsum("xd,kd->xk", positions, result.kt))
        measure = self._spectral_measure()
        (e_u, e_v), _ = _PLANE_AXES[self.normal]
        shape = tuple(len(a) for a in axes)
        outputs = {}
        for component, amplitude in ((e_u, amplitude_u), (e_v, amplitude_v)):
            if not jnp.any(amplitude):
                continue
            values = jnp.einsum("xk,fk->fx", phase, amplitude) * measure
            outputs[component] = values.reshape(values.shape[:1] + shape)
        return outputs

    # ---------------------------------------------------------------- Meep glue
    #
    # Two ways in. `from_monitor` post-processes an ordinary forward run and
    # deals in NumPy, so a user who only wants a far field never meets JAX.
    # `for_design` builds the `FourierFields` an objective function needs, for
    # use inside `OptimizationProblem` or a `MeepJaxWrapper` loss.

    @staticmethod
    def _plane_geometry(simulation: mp.Simulation, volume: mp.Volume):
        """Infers the normal, sample pitch, and sample counts of a planar volume."""
        size = [volume.size.x, volume.size.y, volume.size.z]
        dimensions = simulation.dimensions
        if dimensions not in (2, 3):
            raise NotImplementedError(
                "Angular-spectrum propagation supports 2D and 3D Cartesian "
                f"simulations, not {dimensions}D."
            )
        in_plane = [0, 1] if dimensions == 2 else [0, 1, 2]
        flat = [axis for axis in in_plane if size[axis] == 0]
        if len(flat) != 1:
            raise ValueError(
                "The monitor must be a plane normal to one coordinate axis, "
                f"i.e. a Volume with exactly one zero size; got size={size} in "
                f"{dimensions}D."
            )
        normal = (mp.X, mp.Y, mp.Z)[flat[0]]
        tangential = [axis for axis in in_plane if axis != flat[0]]
        metadata = simulation.get_array_metadata(vol=volume)

        pitches, counts = [], []
        for axis in tangential:
            samples = onp.asarray(metadata[axis])
            if samples.size < 2:
                raise ValueError(
                    "The monitor needs at least two sample points along every "
                    f"tangential axis, but has {samples.size} along axis {axis}."
                )
            pitches.append(float(onp.mean(onp.diff(samples))))
            counts.append(int(samples.size))

        # Meep orders the array axes x, y, z; the propagator wants them in the
        # right-handed (u, v) order of the plane, which for a y-normal plane is
        # (z, x) rather than (x, z).
        if normal == mp.Y and dimensions == 3:
            pitches, counts = pitches[::-1], counts[::-1]
        return normal, tuple(pitches), tuple(counts)

    @staticmethod
    def _assert_homogeneous(
        simulation: mp.Simulation, volume: mp.Volume, tolerance: float = 1e-6
    ):
        """Rejects a monitor that does not lie in a homogeneous region.

        The plane-wave decomposition assumes a single index at the monitor. This
        mirrors the check `dft_near2far` makes for the same reason, and catches a
        monitor accidentally clipping a waveguide or a PML.
        """
        epsilon = onp.asarray(simulation.get_array(vol=volume, component=mp.Dielectric))
        spread = float(onp.max(epsilon) - onp.min(epsilon))
        if spread > tolerance * max(1.0, float(onp.max(epsilon))):
            raise ValueError(
                "The monitor plane does not lie in a homogeneous medium: "
                f"epsilon varies from {onp.min(epsilon):.6g} to "
                f"{onp.max(epsilon):.6g} along it. Angular-spectrum propagation "
                "needs a single index at the monitor, so move the plane clear of "
                "any structure."
            )
        return float(onp.mean(epsilon))

    @classmethod
    def from_monitor(
        cls,
        simulation: mp.Simulation,
        monitor,
        stack: Stack,
        volume: mp.Volume,
        sign: int = 1,
        **kwargs,
    ):
        """Builds a propagator matching an existing `dft_fields` monitor.

        Args:
            simulation: the simulation the monitor belongs to.
            monitor: the object returned by `add_dft_fields`.
            stack: the layers above the monitor.
            volume: the same volume the monitor was registered with.
            sign: +1 if the radiation of interest travels along +normal.
            **kwargs: forwarded to the constructor, e.g. `pad_factor`.
        """
        normal, pitch, num_points = cls._plane_geometry(simulation, volume)
        cls._assert_homogeneous(simulation, volume)
        frequencies = onp.asarray(monitor.freq)
        return cls(
            stack,
            frequencies,
            pitch,
            num_points,
            normal=normal,
            sign=sign,
            **kwargs,
        )

    def fields_from_monitor(
        self, simulation: mp.Simulation, monitor
    ) -> TangentialFields:
        """Reads the tangential DFT fields off a monitor into a `TangentialFields`."""
        electric, magnetic = {}, {}
        # `DftFields` keeps only the component count, but `DftObj` retains the
        # arguments `add_dft_fields` was called with, the first of which is the
        # component list. Fall back to trying every tangential component if that
        # ever stops being true.
        registered = getattr(monitor, "args", None)
        components = (
            set(registered[0])
            if registered
            else {
                component for group in _PLANE_AXES[self.normal] for component in group
            }
        )
        (e_u, e_v), (h_u, h_v) = _PLANE_AXES[self.normal]
        for e_component, h_component in ((e_u, h_u), (e_v, h_v)):
            for component, target in (
                (e_component, electric),
                (h_component, magnetic),
            ):
                if component not in components:
                    continue
                values = onp.array(
                    [
                        simulation.get_dft_array(monitor, component, i)
                        for i in range(len(self.frequencies))
                    ]
                )
                if values.shape[1:] != self.num_points:
                    raise ValueError(
                        f"The monitor returned {values.shape[1:]} samples for "
                        f"{mp.component_name(component)} but the propagator was "
                        f"built for {self.num_points}. The volume Meep actually "
                        "used may have been snapped to the grid; build the "
                        "propagator from the same volume that was registered."
                    )
                target[component] = values
        return TangentialFields(
            E=electric, H=magnetic, normal=self.normal, sign=self.sign
        )

    def propagate_monitor(self, simulation, monitor, distance, coordinates=None):
        """`propagate`, reading from a monitor and returning NumPy arrays."""
        fields = self.fields_from_monitor(simulation, monitor)
        return {
            component: onp.asarray(values)
            for component, values in self.propagate(
                fields, distance, coordinates
            ).items()
        }

    def overlap_monitor(self, simulation, monitor, mode, distance, incident_power=None):
        """`overlap`, reading from a monitor and returning a NumPy array."""
        fields = self.fields_from_monitor(simulation, monitor)
        return onp.asarray(self.overlap(fields, mode, distance, incident_power))

    def power_monitor(self, simulation, monitor, distance=None):
        """`power`, reading from a monitor and returning a NumPy array."""
        fields = self.fields_from_monitor(simulation, monitor)
        return onp.asarray(self.power(fields, distance))

    def report_monitor(self, simulation, monitor) -> Dict[str, onp.ndarray]:
        """`report`, reading from a monitor and returning NumPy arrays."""
        fields = self.fields_from_monitor(simulation, monitor)
        return {key: onp.asarray(value) for key, value in self.report(fields).items()}

    def objective_arguments(self, simulation, volume, **kwargs):
        """The `FourierFields` an objective function needs, for the adjoint path.

        Both tangential components of both polarizations are registered. The
        unused ones cost two extra DFT line monitors in the forward run and
        nothing in the adjoint, since a monitor whose cotangent is identically
        zero places no adjoint source, and registering them removes the need for
        the user to work out which polarization their source excites.
        """
        from . import FourierFields

        self._objective_components = [
            component for group in _PLANE_AXES[self.normal] for component in group
        ]
        return [
            FourierFields(simulation, volume, component, yee_grid=False, **kwargs)
            for component in self._objective_components
        ]

    def take(self, args) -> TangentialFields:
        """Repacks the leading objective-function arguments into `TangentialFields`.

        `OptimizationProblem` passes objective arguments positionally, so this
        consumes the ones `objective_arguments` produced and leaves the rest for
        the caller's own monitors.
        """
        if not hasattr(self, "_objective_components"):
            raise RuntimeError(
                "Call objective_arguments() before take(), so that the "
                "component order is known."
            )
        electric, magnetic = {}, {}
        for component, values in zip(self._objective_components, args):
            target = magnetic if mp.is_magnetic(component) else electric
            target[component] = values
        return TangentialFields(
            E=electric, H=magnetic, normal=self.normal, sign=self.sign
        )

    def __len__(self) -> int:
        """How many leading objective arguments `take` consumes."""
        return 4
