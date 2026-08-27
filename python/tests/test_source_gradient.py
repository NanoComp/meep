"""Tests for adjoint gradients with respect to source amplitudes.

The finite-difference comparisons all pin the run length rather than using
`stop_when_dft_decayed`. That is deliberate: an adaptive stop makes a perturbed
run end at a different time, and the difference scales with the perturbation,
so the ratio of the adjoint gradient to the finite difference settles at a
fixed wrong value instead of converging. It looks exactly like a wrong
gradient. The runs here are also long enough for the *adjoint* DFT, which rings
considerably longer than the forward one in a lossless background.
"""
import unittest

import numpy as np
from autograd import numpy as npa

import meep as mp
import meep.adjoint as mpa

RES = 20
FCEN = 1.0
RUN = 200.0
SRC_C = mp.Vector3(-1.0, 0)
MON_C = mp.Vector3(1.2, 0)
CELL = mp.Vector3(6, 4)
RHO = 0.5 * np.ones(64)
FD_STEP = 1e-4


def _design_region(sim):
    """An inert design region; OptimizationProblem requires one."""
    mg = mp.MaterialGrid(
        mp.Vector3(8, 8), mp.air, mp.Medium(index=1.5), grid_type="U_MEAN"
    )
    dr = mpa.DesignRegion(
        mg, volume=mp.Volume(center=mp.Vector3(0.2, 0), size=mp.Vector3(0.4, 0.4))
    )
    sim.geometry = [mp.Block(center=dr.center, size=dr.size, material=mg)]
    return dr


def _problem(source, fcen=FCEN, res=RES, comp=mp.Ez):
    sim = mp.Simulation(
        cell_size=CELL,
        resolution=res,
        boundary_layers=[mp.PML(1.0)],
        sources=[source],
        force_complex_fields=True,
    )
    dr = _design_region(sim)
    mon = mpa.FourierFields(sim, mp.Volume(center=MON_C, size=mp.Vector3(0, 0)), comp)
    opt = mpa.OptimizationProblem(
        simulation=sim,
        objective_functions=[lambda f: npa.sum(npa.abs(f) ** 2)],
        objective_arguments=[mon],
        design_regions=[dr],
        frequencies=[fcen],
        minimum_run_time=RUN,
        maximum_run_time=RUN,
    )
    return sim, opt


def _num_source_points(center, size, comp=mp.Ez, res=RES):
    """How many grid points a source region covers."""
    sim = mp.Simulation(
        cell_size=CELL,
        resolution=res,
        boundary_layers=[mp.PML(1.0)],
        force_complex_fields=True,
        sources=[
            mp.Source(
                mp.GaussianSource(FCEN, fwidth=0.2), component=comp, center=center
            )
        ],
    )
    sim.init_sim()
    vol = sim._fit_volume_to_simulation(mp.Volume(center=center, size=size))
    mon = sim.add_dft_fields([comp], [FCEN], where=vol, yee_grid=True)
    sim._evaluate_dft_objects()
    return int(np.prod(sim.fields.dft_monitor_size(mon.swigobj, vol.swigobj, comp)))


def _complex_fd(evaluate, perturb, step=FD_STEP):
    """Central difference in the convention JAX and autograd use.

    For a real objective of a complex parameter both frameworks return
    `dJ/d(Re a) - i dJ/d(Im a)`, so that is what the adjoint result is compared
    against.
    """
    out = {}
    for tag, delta in (("re", step), ("im", 1j * step)):
        out[tag] = (evaluate(perturb(delta)) - evaluate(perturb(-delta))) / (2 * step)
    return out["re"] - 1j * out["im"]


class TestDifferentiableFlag(unittest.TestCase):
    """Validation of `differentiable=`; none of these run a simulation."""

    def _src(self, **kwargs):
        return mp.Source(
            mp.GaussianSource(FCEN, fwidth=0.2),
            component=mp.Ez,
            center=SRC_C,
            **kwargs,
        )

    def test_accepts_universal_parameters(self):
        src = self._src(differentiable=["currents", "amplitude"], name="drive")
        self.assertEqual(src.differentiable, ("currents", "amplitude"))
        self.assertEqual(src.name, "drive")

    def test_default_is_not_differentiable(self):
        self.assertEqual(self._src().differentiable, ())

    def test_rejects_unknown_parameter(self):
        with self.assertRaisesRegex(ValueError, "not a differentiable parameter"):
            self._src(differentiable=["nonsense"])

    def test_rejects_bare_string(self):
        with self.assertRaisesRegex(ValueError, "list of parameter names"):
            self._src(differentiable="currents")

    def test_rejects_duplicates(self):
        with self.assertRaisesRegex(ValueError, "duplicate"):
            self._src(differentiable=["currents", "currents"])

    def test_support_moving_parameters_get_their_own_message(self):
        # These are outside the formulation rather than merely unimplemented,
        # and the error should distinguish the two.
        for name in ("center", "size"):
            with self.assertRaisesRegex(ValueError, "grid points the source occupies"):
                self._src(differentiable=[name])

    def test_unimplemented_parameter_is_not_reported_as_unknown(self):
        beam = dict(
            src=mp.GaussianSource(FCEN, fwidth=0.2),
            center=SRC_C,
            size=mp.Vector3(0, 2),
            beam_kdir=mp.Vector3(1, 0, 0),
            beam_w0=1.0,
            beam_E0=mp.Vector3(0, 0, 1),
        )
        with self.assertRaises(NotImplementedError):
            mp.GaussianBeam3DSource(differentiable=["beam_w0"], **beam)
        # ... but it is still rejected as unknown on a class that lacks it
        with self.assertRaises(ValueError):
            self._src(differentiable=["beam_w0"])


class TestTranspose(unittest.TestCase):
    """`fourier_sourcegradient` must be the exact transpose of the scatter."""

    # Kept clear of the PML boundary so the monitor lands in a single chunk.
    # `test_dot_product_identity` needs that: it pairs the scatter's per-chunk
    # amplitudes with `get_dft_array`, and the local-to-global index mapping
    # that would let it do so across chunks is not reachable from Python. The
    # cross-chunk reduction is covered instead by
    # `test_gather_matches_dft_for_an_aligned_monitor`, which compares two
    # globally assembled arrays, and by the finite-difference tests.
    MON_SIZE = mp.Vector3(0, 1.0)

    def _monitor(self, mon_x=1.0, symmetries=(), comp=mp.Ez, freqs=(1.0, 1.15)):
        sim = mp.Simulation(
            cell_size=CELL,
            resolution=RES,
            boundary_layers=[mp.PML(1.0)],
            sources=[
                mp.Source(
                    mp.GaussianSource(FCEN, fwidth=0.2), component=mp.Ez, center=SRC_C
                )
            ],
            symmetries=list(symmetries),
            force_complex_fields=True,
        )
        vol = mp.Volume(center=mp.Vector3(mon_x, 0), size=self.MON_SIZE)
        mon = sim.add_dft_fields([comp], np.asarray(freqs), where=vol, yee_grid=True)
        sim.run(until=25)
        npts = int(np.prod(sim.fields.dft_monitor_size(mon.swigobj, vol.swigobj, comp)))
        return sim, vol, mon, npts, np.asarray(freqs)

    def test_gather_matches_dft_for_an_aligned_monitor(self):
        # With no symmetry and a grid-aligned electric monitor every weight is
        # 1 except the -1 the scatter applies to electric components, so the
        # gathered gradient must be exactly -dft. This pins ordering and sign.
        sim, vol, mon, npts, freqs = self._monitor()
        grad = np.zeros(npts * len(freqs), dtype=np.complex128)
        mon.swigobj.fourier_sourcegradient(vol.swigobj, mp.Ez, sim.fields, grad)
        grad = grad.reshape(len(freqs), npts)
        for i in range(len(freqs)):
            dft = np.asarray(sim.get_dft_array(mon, mp.Ez, i)).ravel()
            np.testing.assert_allclose(grad[i], -dft, rtol=1e-13, atol=0)

    def test_dot_product_identity(self):
        # <dft, S x> == <S^T dft, x>. If the scatter and the gather disagree
        # about ordering or weights, this fails. No FDTD gradient involved.
        sim, vol, mon, npts, freqs = self._monitor()
        nf = len(freqs)
        rng = np.random.default_rng(0)
        x = (
            rng.standard_normal((nf, npts)) + 1j * rng.standard_normal((nf, npts))
        ).ravel()

        grad = np.zeros(nf * npts, dtype=np.complex128)
        mon.swigobj.fourier_sourcegradient(vol.swigobj, mp.Ez, sim.fields, grad)
        rhs = np.sum(grad * x)

        srcdata = mon.swigobj.fourier_sourcedata(vol.swigobj, mp.Ez, sim.fields, x)
        dft = np.stack(
            [np.asarray(sim.get_dft_array(mon, mp.Ez, i)).ravel() for i in range(nf)]
        )
        lhs = 0j
        for sd in srcdata:
            amp = np.array(sd.amp_arr).reshape(-1, nf)
            self.assertEqual(amp.shape[0], npts)  # single aligned chunk
            for i in range(nf):
                lhs += np.sum(dft[i] * amp[:, i])

        self.assertAlmostEqual(abs(lhs - rhs) / abs(lhs), 0.0, places=12)


class TestArraySource(unittest.TestCase):
    def _field(self, source):
        sim = mp.Simulation(
            cell_size=CELL,
            resolution=RES,
            boundary_layers=[mp.PML(1.0)],
            sources=[source],
            force_complex_fields=True,
        )
        mon = sim.add_dft_fields(
            [mp.Ez],
            [FCEN],
            where=mp.Volume(center=MON_C, size=mp.Vector3(0.2, 0.2)),
        )
        sim.run(until_after_sources=60)
        return np.asarray(sim.get_dft_array(mon, mp.Ez, 0)).ravel()

    def test_one_point_matches_an_ordinary_source(self):
        # Pins both conversions in ArraySource.add_source: the -1 the scatter
        # applies to electric components, and the fact that it places a current
        # density rather than a point amplitude.
        #
        # Both an electric and a magnetic component are checked, because the
        # scatter negates electric components and only those. Undoing that
        # negation for magnetic ones too leaves each component individually
        # plausible but flips the relative sign of the two sheets of an
        # equivalent-current pair, which reverses the direction such a pair
        # radiates in without producing anything that looks like an error.
        t = lambda: mp.GaussianSource(FCEN, fwidth=0.2)
        for component in (mp.Ez, mp.Hz):
            with self.subTest(component=mp.component_name(component)):
                plain = self._field(
                    mp.Source(t(), component=component, center=SRC_C, amplitude=1.0)
                )
                array = self._field(
                    mp.ArraySource(
                        t(),
                        component,
                        amplitudes=np.array([1.0 + 0j]),
                        frequency=FCEN,
                        center=SRC_C,
                        size=mp.Vector3(0, 0),
                    )
                )
                np.testing.assert_allclose(array, plain, rtol=2e-6)

    def test_rejects_wrong_length(self):
        src = mp.ArraySource(
            mp.GaussianSource(FCEN, fwidth=0.2),
            mp.Ez,
            amplitudes=np.ones(3, dtype=complex),
            frequency=FCEN,
            center=SRC_C,
            size=mp.Vector3(0, 0.4),
        )
        sim = mp.Simulation(
            cell_size=CELL,
            resolution=RES,
            boundary_layers=[mp.PML(1.0)],
            sources=[src],
            force_complex_fields=True,
        )
        with self.assertRaisesRegex(ValueError, "grid points"):
            sim.init_sim()


class TestSourceGradient(unittest.TestCase):
    def test_no_flagged_source_leaves_the_return_shape_alone(self):
        src = mp.Source(
            mp.GaussianSource(FCEN, fwidth=0.2), component=mp.Ez, center=SRC_C
        )
        _, opt = _problem(src)
        _, grad = opt([RHO])
        self.assertNotIsInstance(grad, dict)

    def test_rejects_a_source_inside_the_pml(self):
        # The adjoint field is absorbed there, so the gradient would come back
        # finite, smooth, and wrong.
        src = mp.Source(
            mp.GaussianSource(FCEN, fwidth=0.2),
            component=mp.Ez,
            center=mp.Vector3(-2.7, 0),
            differentiable=["amplitude"],
            name="drive",
        )
        _, opt = _problem(src)
        with self.assertRaisesRegex(ValueError, "PML"):
            opt([RHO])

    def _amplitude_case(self, fcen=FCEN, res=RES, amp=1.0, comp=mp.Ez):
        def make(a):
            return mp.Source(
                mp.GaussianSource(fcen, fwidth=0.2 * fcen),
                component=comp,
                center=SRC_C,
                amplitude=a,
                differentiable=["amplitude"],
                name="drive",
            )

        _, opt = _problem(make(amp), fcen, res, comp)
        _, grad = opt([RHO])
        adjoint = np.ravel(grad["drive"]["amplitude"])[0]

        def evaluate(a):
            _, o = _problem(make(a), fcen, res, comp)
            return np.asarray(o([RHO], need_gradient=False)[0]).item()

        reference = _complex_fd(evaluate, lambda d: amp + d)
        return adjoint, reference

    def test_amplitude_gradient(self):
        for label, kwargs in [
            ("baseline", {}),
            ("higher resolution", dict(res=30)),
            ("lower frequency", dict(fcen=0.8)),
            ("higher frequency", dict(fcen=1.3)),
            ("complex amplitude", dict(amp=0.7 - 0.4j)),
            ("magnetic source", dict(comp=mp.Hz)),
        ]:
            with self.subTest(label):
                adjoint, reference = self._amplitude_case(**kwargs)
                self.assertLess(
                    abs(adjoint - reference) / abs(reference),
                    1e-5,
                    f"{label}: adjoint {adjoint} vs finite difference {reference}",
                )

    def test_currents_gradient(self):
        size = mp.Vector3(0, 0.4)
        npts = _num_source_points(SRC_C, size)
        rng = np.random.default_rng(7)
        amps = rng.standard_normal(npts) + 1j * rng.standard_normal(npts)

        def make(a):
            return mp.ArraySource(
                mp.GaussianSource(FCEN, fwidth=0.2),
                mp.Ez,
                amplitudes=a,
                frequency=FCEN,
                center=SRC_C,
                size=size,
                differentiable=["currents"],
                name="sheet",
            )

        _, opt = _problem(make(amps))
        _, grad = opt([RHO])
        adjoint = np.ravel(grad["sheet"]["currents"])
        self.assertEqual(adjoint.size, npts)

        def evaluate(a):
            _, o = _problem(make(a))
            return np.asarray(o([RHO], need_gradient=False)[0]).item()

        for j in (0, npts // 2, npts - 1):
            with self.subTest(point=j):

                def perturb(delta, j=j):
                    out = amps.copy()
                    out[j] += delta
                    return out

                reference = _complex_fd(evaluate, perturb)
                self.assertLess(
                    abs(adjoint[j] - reference) / abs(reference),
                    1e-5,
                    f"point {j}: adjoint {adjoint[j]} vs {reference}",
                )


try:
    import jax

    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp

    _HAVE_JAX = True
except ImportError:
    _HAVE_JAX = False


@unittest.skipUnless(_HAVE_JAX, "JAX is an optional dependency")
class TestJaxRoundTrip(unittest.TestCase):
    """`jax.grad` through an `ArraySource` must match a finite difference.

    This is what pins the cotangent convention. Meep returns
    `dJ/d(Re a) - i dJ/d(Im a)`, which is what JAX and autograd both produce
    for a real function of a complex input, so it chains with no adjustment.
    A convention error here would not raise, it would just be wrong.
    """

    def test_gradient_chains_through_jax(self):
        size = mp.Vector3(0, 0.4)
        npts = _num_source_points(SRC_C, size)

        def wrapper():
            src = mp.ArraySource(
                mp.GaussianSource(FCEN, fwidth=0.2),
                mp.Ez,
                amplitudes=np.ones(npts, dtype=complex),
                frequency=FCEN,
                center=SRC_C,
                size=size,
                name="sheet",
            )
            sim = mp.Simulation(
                cell_size=CELL,
                resolution=RES,
                boundary_layers=[mp.PML(1.0)],
                force_complex_fields=True,
            )
            dr = _design_region(sim)
            mon = mpa.FourierFields(
                sim, mp.Volume(center=MON_C, size=mp.Vector3(0, 0)), mp.Ez
            )
            return mpa.MeepJaxWrapper(
                sim,
                [src],
                [mon],
                [dr],
                [FCEN],
                minimum_run_time=RUN,
                maximum_run_time=RUN,
            )

        rng = np.random.default_rng(11)
        amps0 = jnp.asarray(rng.standard_normal(npts) + 1j * rng.standard_normal(npts))
        rho = 0.5 * jnp.ones((8, 8))

        def loss(rho, amps):
            (dft,) = wrapper()([rho], [amps])
            return jnp.sum(jnp.abs(dft) ** 2)

        _, (_, adjoint) = jax.value_and_grad(loss, argnums=(0, 1))(rho, amps0)
        self.assertEqual(adjoint.shape, (npts,))

        for j in (0, npts // 2, npts - 1):
            with self.subTest(point=j):
                reference = _complex_fd(
                    lambda a: loss(rho, a), lambda d, j=j: amps0.at[j].add(d)
                )
                self.assertLess(
                    abs(adjoint[j] - reference) / abs(reference),
                    1e-5,
                    f"point {j}: jax.grad {adjoint[j]} vs {reference}",
                )


if __name__ == "__main__":
    unittest.main()
