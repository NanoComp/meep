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
            index=lambda frequency: onp.sqrt(
                complex(medium.epsilon(frequency)[0][0])
            ),
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
    reflection, transmission = _interface_coefficients(
        admittances[-2], admittances[-1]
    )
    for j in range(len(admittances) - 3, -1, -1):
        phase = jnp.exp(1j * wavevectors[j + 1] * thicknesses[j + 1])
        interface_r, interface_t = _interface_coefficients(
            admittances[j], admittances[j + 1]
        )
        # Redheffer star product of the interface with the layer already
        # accumulated behind it.
        backward_r = -interface_r  # reflection of the interface from the far side
        denominator = 1 - backward_r * reflection * phase**2
        transmission = (
            interface_t * phase * transmission / denominator
        )
        reflection = interface_r + (
            interface_t * (2 - interface_t) * reflection * phase**2 / denominator
        )
    return transmission, reflection


def _single_interface_limit(
    admittances: Sequence[jnp.ndarray],
) -> Tuple[jnp.ndarray, jnp.ndarray]:
    """Transmission and reflection when the stack is a single interface."""
    reflection, transmission = _interface_coefficients(
        admittances[0], admittances[1]
    )
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


def _tangential_components(normal: int, dimensions: int) -> Tuple[Tuple[int, ...], ...]:
    """The E and H components tangential to a plane with the given normal."""
    if dimensions == 2:
        if normal == mp.Y:
            return (mp.Ex, mp.Ez), (mp.Hx, mp.Hz)
        if normal == mp.X:
            return (mp.Ey, mp.Ez), (mp.Hy, mp.Hz)
        raise ValueError(
            "In a 2D simulation the monitor plane must be normal to x or y, "
            f"but got normal={normal}."
        )
    raise NotImplementedError(
        "Angular-spectrum propagation currently supports 2D simulations only. "
        "The three-dimensional case additionally needs the s/p rotation by "
        "azimuth, with its removable singularity at normal incidence, and is "
        "not implemented here."
    )


# For an up-going plane wave the tangential fields satisfy H_t = Y (n_hat x E_t)
# in both polarizations, the cross product supplying the orientation, so
#
#     E_up = (E_t - (n_hat x H_t) / Y) / 2
#
# uniformly. Written out in a tangential basis this becomes one (E component,
# partnering H component, sign, polarization) tuple per polarization, with
# E_up = (E + sign * H / Y) / 2.
_DECOMPOSITION = {
    mp.Y: (
        (mp.Ez, mp.Hx, +1, S_POLARIZATION),
        (mp.Ex, mp.Hz, -1, P_POLARIZATION),
    ),
    mp.X: (
        (mp.Ez, mp.Hy, -1, S_POLARIZATION),
        (mp.Ey, mp.Hz, +1, P_POLARIZATION),
    ),
}


class Mode(NamedTuple):
    """A target field to project onto, defined by its angular spectrum.

    Attributes:
        spectrum: called as `spectrum(kt, k0, index)` and returning a
            (num frequencies, num kt) array of tangential electric field
            amplitudes.
        polarization: which polarization the amplitudes belong to.
    """

    spectrum: Callable
    polarization: int = S_POLARIZATION


def gaussian_mode(
    waist: float,
    tilt_deg: float = 0.0,
    offset: float = 0.0,
    polarization: int = S_POLARIZATION,
) -> Mode:
    """A tilted, laterally offset Gaussian, e.g. a fiber mode.

    The spectrum is written in closed form rather than sampled, so the tilt and
    the offset are exact continuous parameters and differentiable, instead of
    being quantized by the monitor grid.

    Args:
        waist: the 1/e field radius at the target plane, in Meep units. For a
            fiber quoted by mode-field diameter, this is MFD / 2.
        tilt_deg: the angle from the plane normal, in degrees.
        offset: the lateral displacement of the beam center at the target plane.
        polarization: S_POLARIZATION or P_POLARIZATION.

    Returns:
        A `Mode`.
    """

    def spectrum(kt, k0, index):
        kt = jnp.asarray(kt)
        center = (
            jnp.asarray(index).real[:, None]
            * jnp.asarray(k0)[:, None]
            * jnp.sin(jnp.deg2rad(jnp.asarray(tilt_deg)))
        )
        detuning = kt[None, :] - center
        return jnp.exp(-jnp.square(detuning * waist) / 4.0) * jnp.exp(
            -1j * kt[None, :] * offset
        )

    return Mode(spectrum=spectrum, polarization=polarization)


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
        if normal not in _DECOMPOSITION:
            raise ValueError(f"Unsupported monitor normal {normal}.")
        if sign not in (1, -1):
            raise ValueError(f"sign must be +1 or -1, got {sign}.")

        self.stack = stack
        self.frequencies = onp.asarray(frequencies, dtype=float)
        self.pitch = float(pitch)
        self.num_points = int(num_points)
        self.normal = normal
        self.sign = int(sign)
        self.pad_factor = int(pad_factor)

        self._k0 = 2 * onp.pi * self.frequencies
        if kt is None:
            padded = self.num_points * self.pad_factor
            self._kt = jnp.asarray(
                2 * onp.pi * onp.fft.fftfreq(padded, d=self.pitch)
            )
            self._uniform = True
            self._padded = padded
        else:
            self._kt = jnp.asarray(kt, dtype=float)
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
        self._thicknesses = [
            0.0 if layer.thickness is None else float(layer.thickness)
            for layer in stack.layers
        ]

    @property
    def kt(self) -> jnp.ndarray:
        """The transverse wavevectors at which the spectrum is evaluated."""
        return self._kt

    @property
    def stack_thickness(self) -> float:
        """The total thickness of the stack, excluding the semi-infinite layer."""
        return float(sum(self._thicknesses[:-1]))

    def coordinates(self) -> onp.ndarray:
        """The transverse coordinates of the monitor samples, centered on zero."""
        return (onp.arange(self.num_points) - (self.num_points - 1) / 2) * self.pitch

    def _transform(self, values: jnp.ndarray) -> jnp.ndarray:
        """Transforms (num frequencies, num points) samples to (num freq, num kt)."""
        values = jnp.asarray(values)
        if self._uniform:
            padded = jnp.zeros(
                values.shape[:-1] + (self._padded,), dtype=jnp.complex128
            )
            padded = padded.at[..., : self.num_points].set(values)
            spectrum = jnp.fft.fft(padded, axis=-1) * self.pitch
            # The samples are centered on zero, so undo the phase ramp implied by
            # having placed them at indices 0..num_points-1.
            origin = self.coordinates()[0]
            return spectrum * jnp.exp(-1j * self._kt * origin)
        phase = jnp.exp(-1j * self._kt[:, None] * self.coordinates()[None, :])
        return jnp.einsum("kx,...x->...k", phase, values) * self.pitch

    def _polarization_terms(self, fields: TangentialFields):
        """Yields (polarization, E samples, H samples, sign) for what is present."""
        if fields.normal != self.normal:
            raise ValueError(
                f"The fields are on a plane normal to {fields.normal} but this "
                f"propagator was built for {self.normal}."
            )
        for e_component, h_component, sign, polarization in _DECOMPOSITION[
            self.normal
        ]:
            e_values = fields.E.get(e_component)
            h_values = fields.H.get(h_component)
            if e_values is None and h_values is None:
                continue
            if e_values is None or h_values is None:
                raise ValueError(
                    "Separating up-going from down-going radiation needs both "
                    f"tangential fields, but only one of {mp.component_name(e_component)}"
                    f" and {mp.component_name(h_component)} was supplied."
                )
            yield polarization, jnp.asarray(e_values), jnp.asarray(h_values), sign

    def decompose(self, fields: TangentialFields):
        """Splits the monitor fields into up- and down-going spectra.

        Returns:
            `(up, down)`, each a dict mapping polarization to a
            (num frequencies, num kt) array of tangential electric field
            amplitudes at the monitor plane.
        """
        up, down = {}, {}
        admittance = self._admittances[0]
        for polarization, e_values, h_values, sign in self._polarization_terms(
            fields
        ):
            e_spectrum = self._transform(e_values)
            h_spectrum = self._transform(h_values)
            # `sign` already carries the orientation of n_hat x H_t; `fields.sign`
            # flips which branch counts as outgoing for a downward-facing monitor.
            scaled = fields.sign * sign * h_spectrum / admittance[polarization]
            up[polarization] = 0.5 * (e_spectrum + scaled)
            down[polarization] = 0.5 * (e_spectrum - scaled)
        if not up:
            raise ValueError(
                "No tangential field components were supplied; nothing to "
                "propagate."
            )
        return up, down

    def _transmission(self, polarization: int):
        """Transmission and reflection of the stack for one polarization."""
        admittances = [pair[polarization] for pair in self._admittances]
        if len(admittances) == 2:
            return _single_interface_limit(admittances)
        return _stack_transmission(
            admittances, self._wavevectors, self._thicknesses
        )

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
        if onp.any(onp.asarray(remaining) < 0):
            raise ValueError(
                f"distance={distance} does not clear the stack, which is "
                f"{self.stack_thickness} thick. The target plane has to lie in "
                "the semi-infinite layer."
            )
        up, _ = self.decompose(fields)
        amplitudes = {}
        for polarization, amplitude in up.items():
            transmission, _ = self._transmission(polarization)
            # Monitor to the first interface, through the stack, then onward in
            # the terminating layer.
            to_interface = jnp.exp(
                1j * self._wavevectors[0] * self._thicknesses[0]
            )
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
            return (2 * onp.pi / (self._padded * self.pitch)) / (2 * onp.pi)
        spacing = jnp.diff(self._kt)
        # Trapezoid weights would be more careful, but a chosen kt set is
        # normally uniform; require that rather than silently mis-weighting.
        return jnp.mean(spacing) / (2 * onp.pi)

    def power(self, fields: TangentialFields, distance: Optional[float] = None):
        """Outgoing power through the target plane, one value per frequency."""
        result = self.spectrum(fields, self.stack_thickness if distance is None else distance)
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
        weights = self._weights(result)[..., mode.polarization]
        field_amplitude = result.amplitudes[..., mode.polarization]
        mode_amplitude = mode.spectrum(result.kt, self._k0, result.index)

        cross = jnp.sum(weights * field_amplitude * jnp.conj(mode_amplitude), axis=-1)
        mode_norm = jnp.sum(weights * jnp.abs(mode_amplitude) ** 2, axis=-1)
        coupled = jnp.abs(cross) ** 2 / mode_norm
        if incident_power is None:
            incident_power = jnp.sum(
                weights * jnp.abs(field_amplitude) ** 2, axis=-1
            )
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
        for _, e_values, _, _ in self._polarization_terms(fields):
            magnitude = jnp.abs(e_values)
            peak = jnp.max(magnitude, axis=-1)
            edge = jnp.maximum(magnitude[..., 0], magnitude[..., -1])
            edges.append(edge / jnp.where(peak > 0, peak, 1.0))

        total = up_power + down_power
        return {
            "downgoing_fraction": down_power / jnp.where(total > 0, total, 1.0),
            "evanescent_fraction": 1.0
            - up_propagating / jnp.where(up_power > 0, up_power, 1.0),
            "edge_amplitude": jnp.max(jnp.stack(edges, axis=0), axis=0),
        }

    def propagate(
        self, fields: TangentialFields, distance: float, coordinates=None
    ) -> Dict[int, jnp.ndarray]:
        """The tangential electric field at the target plane, in real space.

        Args:
            fields: the monitor fields.
            distance: distance from the monitor to the target plane.
            coordinates: where to evaluate. Defaults to the monitor's own
                coordinates; pass a wider range to see a beam that has spread.

        Returns:
            A dict mapping each tangential electric component to a
            (num frequencies, num coordinates) array.
        """
        result = self.spectrum(fields, distance)
        if coordinates is None:
            coordinates = self.coordinates()
        coordinates = jnp.asarray(coordinates)
        measure = self._spectral_measure()
        phase = jnp.exp(1j * result.kt[None, :] * coordinates[:, None])
        outputs = {}
        for e_component, _, _, polarization in _DECOMPOSITION[self.normal]:
            if not jnp.any(result.amplitudes[..., polarization]):
                continue
            outputs[e_component] = (
                jnp.einsum(
                    "xk,fk->fx", phase, result.amplitudes[..., polarization]
                )
                * measure
            )
        return outputs
