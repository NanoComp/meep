"""Tests for angular-spectrum propagation through stratified media.

Everything except `TestAgainstNearToFar` runs without an FDTD simulation, using
analytic oracles: Fresnel coefficients, the Brewster angle, a quarter-wave
antireflection coating, an independently written transfer matrix, energy
conservation, and the spreading of a Gaussian beam.
"""

import math
import unittest

import numpy as onp

import meep as mp

try:
    import meep.adjoint as mpa
except ImportError:
    import adjoint as mpa

from utils import ApproxComparisonTestCase

try:
    import jax

    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp
    from meep.adjoint import angular_spectrum as asm
except ImportError:
    jax = None

mp.verbosity(0)

WAVELENGTH = 1.55
N_OXIDE, N_AIR, N_HIGH = 1.444, 1.0, 1.9


def _wavevector(index, k0, kt):
    return asm._longitudinal_wavevector(jnp.array([index + 0j]), k0, kt)


def _admittance(index, kz, k0, polarization):
    return asm._admittances(jnp.array([index + 0j]), kz, k0)[polarization]


@unittest.skipIf(jax is None, "jax is not installed")
class TestStackSolver(ApproxComparisonTestCase):
    """The scattering recursion, against closed-form references."""

    def setUp(self):
        self.k0 = onp.array([2 * onp.pi / WAVELENGTH])

    def _coefficients(self, indices, thicknesses, angle_deg, polarization):
        kt = jnp.array([indices[0] * self.k0[0] * math.sin(math.radians(angle_deg))])
        wavevectors = [_wavevector(n, self.k0, kt) for n in indices]
        admittances = [
            _admittance(n, kz, self.k0, polarization)
            for n, kz in zip(indices, wavevectors)
        ]
        if len(indices) == 2:
            return asm._single_interface_limit(admittances)
        return asm._stack_transmission(admittances, wavevectors, thicknesses)

    def test_single_interface_matches_fresnel(self):
        """Reflection at one interface, swept to near grazing.

        The p-polarization reference carries a sign: these are the coefficients
        of the *tangential* field, which is the admittance convention, and it is
        the negative of the form written for the full electric vector.
        """
        for polarization, name in (
            (asm.S_POLARIZATION, "s"),
            (asm.P_POLARIZATION, "p"),
        ):
            for angle in (0, 10, 30, 45, 60, 80, 89):
                theta = math.radians(angle)
                _, reflection = self._coefficients(
                    [N_OXIDE, N_AIR], [0.0], angle, polarization
                )
                cos_in = math.cos(theta)
                cos_out = onp.sqrt(
                    complex(1 - (N_OXIDE * math.sin(theta) / N_AIR) ** 2)
                )
                if polarization == asm.S_POLARIZATION:
                    expected = (N_OXIDE * cos_in - N_AIR * cos_out) / (
                        N_OXIDE * cos_in + N_AIR * cos_out
                    )
                else:
                    expected = -(N_AIR * cos_in - N_OXIDE * cos_out) / (
                        N_AIR * cos_in + N_OXIDE * cos_out
                    )
                self.assertAlmostEqual(
                    complex(reflection[0, 0]), expected, places=12,
                    msg=f"{name} polarization at {angle} deg",
                )

    def test_brewster_angle(self):
        """p-polarized reflection vanishes at atan(n_out / n_in)."""
        angle = math.degrees(math.atan(N_AIR / N_OXIDE))
        _, reflection = self._coefficients(
            [N_OXIDE, N_AIR], [0.0], angle, asm.P_POLARIZATION
        )
        self.assertLess(abs(complex(reflection[0, 0])), 1e-12)

    def test_quarter_wave_antireflection_coating(self):
        """An index-matched quarter wave nulls the reflection at normal incidence."""
        index = math.sqrt(N_OXIDE * N_AIR)
        thickness = WAVELENGTH / (4 * index)
        _, reflection = self._coefficients(
            [N_OXIDE, index, N_AIR], [0.0, thickness], 0.0, asm.S_POLARIZATION
        )
        self.assertLess(abs(complex(reflection[0, 0])), 1e-12)

    def test_matches_independent_transfer_matrix(self):
        """Against a transfer-matrix implementation written from scratch here.

        The propagator deliberately does not use transfer matrices, since they
        overflow on a thick layer or an evanescent order, so this is a genuinely
        independent route to the same numbers in a regime where both are stable.
        """
        indices, thicknesses = [N_OXIDE, N_HIGH, N_AIR], [0.0, 0.37]

        def transfer_matrix(angle_deg, polarization):
            kt = indices[0] * self.k0[0] * math.sin(math.radians(angle_deg))
            wavevectors = [
                complex(onp.sqrt(complex((n * self.k0[0]) ** 2 - kt**2)))
                for n in indices
            ]
            if polarization == asm.S_POLARIZATION:
                admittances = [kz / self.k0[0] for kz in wavevectors]
            else:
                admittances = [
                    n**2 * self.k0[0] / kz for n, kz in zip(indices, wavevectors)
                ]
            matrix = onp.eye(2, dtype=complex)
            for j in range(len(indices) - 1):
                r = (admittances[j] - admittances[j + 1]) / (
                    admittances[j] + admittances[j + 1]
                )
                t = 2 * admittances[j] / (admittances[j] + admittances[j + 1])
                matrix = matrix @ (onp.array([[1, r], [r, 1]], dtype=complex) / t)
                if j + 1 < len(indices) - 1:
                    phase = onp.exp(1j * wavevectors[j + 1] * thicknesses[j + 1])
                    matrix = matrix @ onp.array(
                        [[1 / phase, 0], [0, phase]], dtype=complex
                    )
            return 1 / matrix[0, 0], matrix[1, 0] / matrix[0, 0]

        for polarization in (asm.S_POLARIZATION, asm.P_POLARIZATION):
            for angle in (0, 25, 55):
                transmission, reflection = self._coefficients(
                    indices, thicknesses, angle, polarization
                )
                expected_t, expected_r = transfer_matrix(angle, polarization)
                self.assertAlmostEqual(
                    complex(transmission[0, 0]), expected_t, places=12
                )
                self.assertAlmostEqual(
                    complex(reflection[0, 0]), expected_r, places=12
                )

    def test_energy_is_conserved(self):
        """Reflected plus transmitted power equals the incident, losslessly."""
        indices, thicknesses = [N_OXIDE, N_HIGH, N_AIR], [0.0, 0.37]
        for polarization in (asm.S_POLARIZATION, asm.P_POLARIZATION):
            for angle in (0, 20, 40):
                transmission, reflection = self._coefficients(
                    indices, thicknesses, angle, polarization
                )
                kt = jnp.array(
                    [indices[0] * self.k0[0] * math.sin(math.radians(angle))]
                )
                y_in = _admittance(
                    indices[0], _wavevector(indices[0], self.k0, kt), self.k0,
                    polarization,
                )
                y_out = _admittance(
                    indices[-1], _wavevector(indices[-1], self.k0, kt), self.k0,
                    polarization,
                )
                reflected = abs(complex(reflection[0, 0])) ** 2
                transmitted = abs(complex(transmission[0, 0])) ** 2 * float(
                    onp.real(complex(y_out[0, 0])) / onp.real(complex(y_in[0, 0]))
                )
                self.assertAlmostEqual(reflected + transmitted, 1.0, places=12)


def _uniform_propagator(index=N_OXIDE, num_points=512, pitch=0.05, pad_factor=8):
    """A propagator with no interface, i.e. plain homogeneous propagation."""
    stack = mpa.Stack([mpa.Layer(index, 0.0), mpa.Layer(index)])
    return mpa.AngularSpectrum(
        stack, [1 / WAVELENGTH], pitch, num_points, normal=mp.Y,
        pad_factor=pad_factor,
    )


def _up_going(propagator, values):
    """Builds the magnetic partner that makes `values` purely up-going.

    The relation is per wavevector, so the partner has to be constructed in the
    spectral domain; scaling by the admittance of one launch angle in real space
    leaves spurious down-going content behind.
    """
    spectrum = propagator._transform(jnp.asarray(values))
    admittance = propagator._admittances[0][asm.S_POLARIZATION]
    coordinates = propagator.coordinates()
    magnetic = (
        jnp.einsum(
            "xk,fk->fx",
            jnp.exp(1j * propagator.kt[None, :] * coordinates[:, None]),
            spectrum * admittance,
        )
        * propagator._spectral_measure()
    )
    return mpa.TangentialFields(
        E={mp.Ez: jnp.asarray(values)}, H={mp.Hx: magnetic}, normal=mp.Y
    )


@unittest.skipIf(jax is None, "jax is not installed")
class TestDecomposition(ApproxComparisonTestCase):
    """Separating up- from down-going radiation."""

    def test_residual_tracks_edge_truncation(self):
        """A purely up-going field yields no down-going part, up to truncation.

        The transform is periodic, so whatever amplitude survives at the ends of
        the monitor wraps around and appears as spurious down-going content. The
        two should therefore track each other, which is what makes the
        `edge_amplitude` diagnostic meaningful.
        """
        propagator = _uniform_propagator()
        x = propagator.coordinates()
        for width in (5.0, 3.0):
            values = (onp.exp(-((x / width) ** 2)))[None, :]
            up, down = propagator.decompose(_up_going(propagator, values))
            residual = float(
                jnp.max(jnp.abs(down[asm.S_POLARIZATION]))
                / jnp.max(jnp.abs(up[asm.S_POLARIZATION]))
            )
            edge = float(onp.exp(-((x.max() / width) ** 2)))
            self.assertLess(residual, max(10 * edge, 1e-9), f"width {width}")

    def test_down_going_field_is_recognized(self):
        """Flipping the magnetic field flips which branch the power lands in."""
        propagator = _uniform_propagator()
        x = propagator.coordinates()
        values = (onp.exp(-((x / 3.0) ** 2)))[None, :]
        fields = _up_going(propagator, values)
        flipped = mpa.TangentialFields(
            E=fields.E, H={mp.Hx: -fields.H[mp.Hx]}, normal=mp.Y
        )
        up, down = propagator.decompose(flipped)
        self.assertLess(
            float(
                jnp.max(jnp.abs(up[asm.S_POLARIZATION]))
                / jnp.max(jnp.abs(down[asm.S_POLARIZATION]))
            ),
            1e-6,
        )


@unittest.skipIf(jax is None, "jax is not installed")
class TestPropagation(ApproxComparisonTestCase):
    """Propagation against analytic beam physics."""

    def test_round_trip_at_zero_distance(self):
        propagator = _uniform_propagator()
        x = propagator.coordinates()
        values = (onp.exp(-((x / 3.0) ** 2)))[None, :]
        output = propagator.propagate(_up_going(propagator, values), 0.0)[mp.Ez]
        self.assertLess(float(jnp.max(jnp.abs(output - values))), 1e-8)

    def test_gaussian_beam_spreading(self):
        """The 1/e width follows w0 sqrt(1 + (z/zR)^2) exactly.

        The padded window has to be wide enough for the spread beam: the
        transform is periodic, so a beam wider than the window wraps onto itself
        and the answer is nonsense rather than merely inaccurate.
        """
        waist = 3.0
        rayleigh = math.pi * waist**2 * N_OXIDE / WAVELENGTH
        for distance, pad_factor in ((50.0, 16), (200.0, 32)):
            propagator = _uniform_propagator(pad_factor=pad_factor)
            x = propagator.coordinates()
            values = (onp.exp(-((x / waist) ** 2)))[None, :]
            expected = waist * math.sqrt(1 + (distance / rayleigh) ** 2)
            samples = onp.linspace(-4 * expected, 4 * expected, 3001)
            profile = onp.abs(
                onp.asarray(
                    propagator.propagate(
                        _up_going(propagator, values), distance,
                        coordinates=samples,
                    )[mp.Ez][0]
                )
            )
            above = samples[profile >= profile.max() / onp.e]
            measured = (above.max() - above.min()) / 2
            self.assertAlmostEqual(measured / expected, 1.0, places=6)

    def test_refuses_to_stop_short_of_the_stack(self):
        stack = mpa.Stack([mpa.Layer(N_OXIDE, 5.0), mpa.Layer(N_AIR)])
        propagator = mpa.AngularSpectrum(
            stack, [1 / WAVELENGTH], 0.05, 128, normal=mp.Y
        )
        x = propagator.coordinates()
        fields = _up_going(propagator, (onp.exp(-((x / 1.0) ** 2)))[None, :])
        with self.assertRaisesRegex(ValueError, "does not clear the stack"):
            propagator.spectrum(fields, 2.0)


@unittest.skipIf(jax is None, "jax is not installed")
class TestOverlap(ApproxComparisonTestCase):
    """Projection onto a target mode."""

    def test_self_overlap_is_unity(self):
        propagator = _uniform_propagator()
        x = propagator.coordinates()
        fields = _up_going(propagator, (onp.exp(-((x / 3.0) ** 2)))[None, :])
        result = propagator.spectrum(fields, 100.0)
        itself = mpa.Mode(
            spectrum=lambda kt, k0, index: result.amplitudes[..., asm.S_POLARIZATION]
        )
        self.assertAlmostEqual(
            float(propagator.overlap(fields, itself, 100.0)[0]), 1.0, places=10
        )

    def test_recovers_the_launch_angle(self):
        """A beam launched at an angle matches a mode tilted to that angle.

        The mode spectrum is written in closed form, so the tilt is a continuous
        parameter rather than something quantized by the monitor grid.
        """
        propagator = _uniform_propagator()
        x = propagator.coordinates()
        k0 = 2 * onp.pi / WAVELENGTH
        for angle in (0.0, 8.0, 20.0):
            transverse = N_OXIDE * k0 * math.sin(math.radians(angle))
            values = (onp.exp(1j * transverse * x) * onp.exp(-((x / 20.0) ** 2)))[
                None, :
            ]
            fields = _up_going(propagator, values)
            candidates = onp.linspace(angle - 4, angle + 4, 33)
            best = max(
                candidates,
                key=lambda t: float(
                    propagator.overlap(
                        fields, mpa.gaussian_mode(20.0, tilt_deg=float(t)), 1e-9
                    )[0]
                ),
            )
            self.assertAlmostEqual(float(best), angle, places=6)


@unittest.skipIf(jax is None, "jax is not installed")
class TestGradients(ApproxComparisonTestCase):
    """Differentiability, including the failure that motivates the regularizer."""

    def test_gradient_is_finite_on_the_light_line(self):
        """A wavevector landing exactly on the light line must not poison the VJP.

        kz vanishes there and its derivative is unbounded, so an unregularized
        sqrt evaluates to zero happily and returns an infinite cotangent. The
        objective value looks perfect while every gradient is NaN. The window
        below is chosen so that a grid point lands exactly on n k0.
        """
        num_points, pad_factor = 64, 1
        # Place a sample exactly on the light line: kt = 2 pi m / (N dx) equals
        # n k0 = 2 pi n / lambda when N dx = m lambda / n.
        pitch = 8 * WAVELENGTH / N_OXIDE / num_points
        propagator = _uniform_propagator(
            num_points=num_points, pitch=pitch, pad_factor=pad_factor
        )
        on_light_line = onp.min(
            onp.abs(onp.asarray(propagator.kt) - N_OXIDE * 2 * onp.pi / WAVELENGTH)
        )
        self.assertLess(on_light_line, 1e-9, "the test window missed the light line")

        x = propagator.coordinates()
        values = jnp.asarray((onp.exp(-((x / 0.5) ** 2)))[None, :])

        def objective(scale):
            fields = _up_going(propagator, scale * values)
            return jnp.real(
                jnp.sum(propagator.spectrum(fields, 1.0).amplitudes)
            )

        gradient = jax.grad(objective)(1.0)
        self.assertTrue(onp.isfinite(float(gradient)), "gradient is not finite")

    def test_differentiable_in_stack_and_mode_parameters(self):
        """Layer thickness and fiber tilt carry gradients, not just the fields."""
        x = onp.linspace(-6.4, 6.4, 256)
        values = (onp.exp(-((x / 2.0) ** 2)))[None, :]

        def objective(thickness, tilt):
            stack = mpa.Stack([mpa.Layer(N_OXIDE, thickness), mpa.Layer(N_AIR)])
            propagator = mpa.AngularSpectrum(
                stack, [1 / WAVELENGTH], 0.05, 256, normal=mp.Y, pad_factor=8
            )
            fields = _up_going(propagator, values)
            return propagator.overlap(
                fields, mpa.gaussian_mode(4.0, tilt_deg=tilt), thickness + 10.0
            )[0]

        gradients = jax.grad(objective, argnums=(0, 1))(2.0, 5.0)
        for name, value in zip(("thickness", "tilt"), gradients):
            self.assertTrue(onp.isfinite(float(value)), name)
            self.assertNotEqual(float(value), 0.0, name)


@unittest.skipIf(jax is None, "jax is not installed")
class TestValidation(unittest.TestCase):
    """The constructor rejects stacks and monitors that cannot work."""

    def test_stack_must_terminate(self):
        with self.assertRaisesRegex(ValueError, "semi-infinite"):
            mpa.Stack([mpa.Layer(N_OXIDE, 1.0), mpa.Layer(N_AIR, 1.0)]).validate()

    def test_only_the_last_layer_may_be_semi_infinite(self):
        with self.assertRaisesRegex(ValueError, "only the final layer"):
            mpa.Stack(
                [mpa.Layer(N_OXIDE), mpa.Layer(N_HIGH, 1.0), mpa.Layer(N_AIR)]
            ).validate()

    def test_both_tangential_fields_are_required(self):
        """One field alone cannot distinguish up-going from down-going."""
        propagator = _uniform_propagator(num_points=64)
        fields = mpa.TangentialFields(
            E={mp.Ez: onp.ones((1, 64))}, H={}, normal=mp.Y
        )
        with self.assertRaisesRegex(ValueError, "both"):
            propagator.decompose(fields)

    def test_three_dimensions_is_rejected_clearly(self):
        with self.assertRaisesRegex(NotImplementedError, "2D"):
            asm._tangential_components(mp.Z, 3)


@unittest.skipIf(jax is None, "jax is not installed")
class TestAgainstNearToFar(ApproxComparisonTestCase):
    """Against Meep's own near-to-far transformation, in a homogeneous medium.

    This is the end-to-end check on the Meep adapter: sample pitch, coordinate
    origin and the phase ramp that comes with it. What matters is that the
    disagreement falls with resolution -- a mistake in any of those produces a
    fixed error instead.
    """

    def _run(self, resolution):
        cell_x, cell_y, pml = 40.0, 8.0, 1.0
        frequency = 1 / WAVELENGTH
        source = [
            mp.Source(
                mp.GaussianSource(frequency, fwidth=0.1 * frequency),
                component=mp.Ez, center=mp.Vector3(0, -1.0),
                size=mp.Vector3(12.0, 0),
                # Narrow and apodized, so the field has decayed by the ends of
                # the monitor. A bare dipole never does, and near2far tolerates
                # that while a periodic transform does not.
                amp_func=lambda p: onp.exp(-((p.x / 1.5) ** 2)),
            )
        ]
        simulation = mp.Simulation(
            resolution=resolution, cell_size=mp.Vector3(cell_x, cell_y),
            boundary_layers=[mp.PML(pml)], sources=source,
            default_material=mp.Medium(index=N_OXIDE),
        )
        height, width = 1.0, cell_x - 2 * pml - 0.4
        volume = mp.Volume(
            center=mp.Vector3(0, height), size=mp.Vector3(width, 0)
        )
        monitor = simulation.add_dft_fields(
            [mp.Ez, mp.Hx], [frequency], where=volume
        )
        near2far = simulation.add_near2far(
            [frequency],
            mp.Near2FarRegion(
                center=mp.Vector3(0, height), size=mp.Vector3(width, 0), weight=+1
            ),
        )
        simulation.run(
            until_after_sources=mp.stop_when_dft_decayed(1e-10, minimum_run_time=30)
        )

        stack = mpa.Stack([mpa.Layer(N_OXIDE, 0.0), mpa.Layer(N_OXIDE)])
        propagator = mpa.AngularSpectrum.from_monitor(
            simulation, monitor, stack, volume, pad_factor=16
        )
        distance = 25.0
        samples = onp.linspace(-8, 8, 33)
        got = propagator.propagate_monitor(
            simulation, monitor, distance, coordinates=samples
        )[mp.Ez][0]
        want = onp.array(
            [
                simulation.get_farfield(near2far, mp.Vector3(x, height + distance))[2]
                for x in samples
            ]
        )
        report = propagator.report_monitor(simulation, monitor)
        return (
            onp.abs(got - want).max() / onp.abs(want).max(),
            float(report["edge_amplitude"][0]),
        )

    def test_agrees_and_converges(self):
        coarse, coarse_edge = self._run(20)
        fine, _ = self._run(40)
        # If the field had not decayed at the monitor edges the comparison would
        # be meaningless, so assert that first.
        self.assertLess(coarse_edge, 1e-6)
        self.assertLess(coarse, 3e-2)
        self.assertLess(fine, coarse / 2, "did not converge with resolution")


if __name__ == "__main__":
    unittest.main()
