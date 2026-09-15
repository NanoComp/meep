"""Tests for adjoint gradients with respect to a geometric object's shape.

Two conventions run through all of these, and both cost real debugging time to
establish.

Object faces are placed at *quarter*-pixel offsets. Yee components sit half a
pixel apart, so a face at a pixel centre for one component lies exactly on a
pixel edge for another. On an edge no pixel straddles the face: the filling
fraction is 0 or 1 on both sides, meep applies no smoothing there, and the
derivative is one-sided. A quarter-pixel offset avoids both integer and
half-integer multiples of the pixel, so every component straddles every face
and the derivative is two-sided everywhere.

Run lengths are pinned rather than left to `stop_when_dft_decayed`, for the
same reason as the source-gradient tests: an adaptive stop makes a perturbed
run end at a different time, and the difference scales with the perturbation.
"""

import unittest
import warnings

import numpy as np
from autograd import numpy as npa

import meep as mp
import meep.adjoint as mpa

RES = 20
FCEN = 1.0
RUN = 200.0
CELL = mp.Vector3(8, 6)
SRC_C = mp.Vector3(-2.0, 0)
MON_C = mp.Vector3(2.0, 0)
RHO = 0.5 * np.ones(64)
FD_STEP = 1e-3


def quarter_pixel(resolution: float) -> float:
    """The offset that keeps a face clear of every component's pixel edges."""
    return 0.25 / resolution


def _design_region(sim):
    """An inert design region; OptimizationProblem requires one."""
    grid = mp.MaterialGrid(
        mp.Vector3(8, 8), mp.air, mp.Medium(index=1.5), grid_type="U_MEAN"
    )
    region = mpa.DesignRegion(
        grid, volume=mp.Volume(center=mp.Vector3(0, -2.0), size=mp.Vector3(0.4, 0.4))
    )
    sim.geometry = list(sim.geometry) + [
        mp.Block(center=region.center, size=region.size, material=grid)
    ]
    return region


def _problem(block, res=RES, component=mp.Hz, eps_averaging=True):
    sim = mp.Simulation(
        cell_size=CELL,
        resolution=res,
        boundary_layers=[mp.PML(1.0)],
        sources=[
            mp.Source(
                mp.GaussianSource(FCEN, fwidth=0.2), component=component, center=SRC_C
            )
        ],
        geometry=[block],
        eps_averaging=eps_averaging,
        force_complex_fields=True,
    )
    region = _design_region(sim)
    monitor = mpa.FourierFields(
        sim, mp.Volume(center=MON_C, size=mp.Vector3(0, 0)), component
    )
    return mpa.OptimizationProblem(
        simulation=sim,
        objective_functions=[lambda f: npa.sum(npa.abs(f) ** 2)],
        objective_arguments=[monitor],
        design_regions=[region],
        frequencies=[FCEN],
        minimum_run_time=RUN,
        maximum_run_time=RUN,
    )


def _block(center, size, index=2.5, differentiable=None, res=RES):
    """A block whose faces avoid every component's pixel edges."""
    offset = quarter_pixel(res)
    return mp.Block(
        center=mp.Vector3(center.x, center.y + offset),
        size=mp.Vector3(size.x + 2 * offset, size.y),
        material=mp.Medium(index=index),
        differentiable=differentiable,
        name="scatterer",
    )


class TestDifferentiableFlag(unittest.TestCase):
    """Validation of `differentiable=` on a geometric object."""

    def test_accepts_center_and_size(self):
        block = mp.Block(
            center=mp.Vector3(),
            size=mp.Vector3(1, 1),
            differentiable=["center", "size"],
            name="b",
        )
        self.assertEqual(block.differentiable, ("center", "size"))
        self.assertEqual(block.name, "b")

    def test_default_is_not_differentiable(self):
        self.assertEqual(mp.Block(size=mp.Vector3(1, 1)).differentiable, ())

    def test_rejects_unknown_parameter(self):
        with self.assertRaisesRegex(ValueError, "not a differentiable parameter"):
            mp.Block(size=mp.Vector3(1, 1), differentiable=["waist"])

    def test_rejects_bare_string(self):
        with self.assertRaisesRegex(ValueError, "list of parameter names"):
            mp.Block(size=mp.Vector3(1, 1), differentiable="center")

    def test_center_and_size_are_allowed_here_unlike_for_sources(self):
        # A source's cotangent is gathered over a fixed set of grid points, so
        # moving it changes which points it occupies rather than the amplitudes
        # applied to them, and 'center' is rejected there. Epsilon is a
        # function of position, so for geometry the same name is exactly what a
        # derivative with respect to position means.
        with self.assertRaises(ValueError):
            mp.Source(
                mp.GaussianSource(FCEN, fwidth=0.2),
                component=mp.Ez,
                center=SRC_C,
                differentiable=["center"],
            )
        mp.Block(size=mp.Vector3(1, 1), differentiable=["center"])


class TestGuards(unittest.TestCase):
    def test_refuses_without_subpixel_smoothing(self):
        # Without smoothing the permittivity is a step function of position:
        # nothing changes until a boundary crosses a pixel edge, then it
        # changes by the full contrast. A finite difference of that is not a
        # derivative, and returning a number would hide it.
        block = _block(mp.Vector3(), mp.Vector3(1.0, 0.8), differentiable=["center"])
        opt = _problem(block, eps_averaging=False)
        with self.assertRaisesRegex(ValueError, "subpixel smoothing"):
            opt([RHO])

    def test_warns_when_a_face_lies_on_a_pixel_edge(self):
        # Pixel centres are at integer multiples of dx, so edges are at
        # half-integers. A face there is straddled by nothing.
        half = 0.5 / RES
        block = mp.Block(
            center=mp.Vector3(0, half),
            size=mp.Vector3(1.0, 0.8),
            material=mp.Medium(index=2.5),
            differentiable=["center"],
            name="s",
        )
        opt = _problem(block)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            opt([RHO])
        self.assertTrue(
            any("pixel edge" in str(w.message) for w in caught),
            "expected a warning about a face on a pixel edge",
        )


class TestShapeDerivative(unittest.TestCase):
    SIZE = mp.Vector3(1.0, 0.8)

    def _value(self, center, size, res=RES, index=2.5):
        opt = _problem(_block(center, size, index, res=res), res=res)
        return float(np.asarray(opt([RHO], need_gradient=False)[0]).item())

    def _adjoint(self, names, res=RES, index=2.5):
        opt = _problem(
            _block(mp.Vector3(), self.SIZE, index, differentiable=names, res=res),
            res=res,
        )
        _, grad = opt([RHO])
        return grad["scatterer"]

    def test_center_gradient(self):
        for label, res, index in (
            ("baseline", RES, 2.5),
            ("higher resolution", 40, 2.5),
            ("lower contrast", RES, 1.5),
        ):
            with self.subTest(label):
                gradient = self._adjoint(["center"], res, index)["center"]
                adjoint = float(np.real(np.atleast_1d(gradient)[1]))
                d = mp.Vector3(y=FD_STEP)
                reference = (
                    self._value(d, self.SIZE, res, index)
                    - self._value(mp.Vector3() - d, self.SIZE, res, index)
                ) / (2 * FD_STEP)
                self.assertLess(
                    abs(adjoint - reference) / abs(reference),
                    2e-2,
                    f"{label}: adjoint {adjoint} vs finite difference {reference}",
                )

    def test_size_gradient(self):
        adjoint = float(np.real(np.atleast_1d(self._adjoint(["size"])["size"])[1]))
        d = mp.Vector3(y=FD_STEP)
        reference = (
            self._value(mp.Vector3(), self.SIZE + d)
            - self._value(mp.Vector3(), self.SIZE - d)
        ) / (2 * FD_STEP)
        self.assertLess(abs(adjoint - reference) / abs(reference), 2e-2)

    def test_translation_through_a_uniform_medium_is_zero(self):
        # A block of the background material is not there at all, so moving it
        # cannot change anything. Exact, and independent of finite-difference
        # error, which makes it sharper than the comparisons above.
        gradient = self._adjoint(["center"], index=1.0)["center"]
        self.assertLess(float(np.max(np.abs(np.atleast_1d(gradient)))), 1e-12)

    def test_returns_one_entry_per_named_parameter(self):
        gradient = self._adjoint(["center", "size"])
        self.assertEqual(set(gradient), {"center", "size"})
        for name in ("center", "size"):
            self.assertEqual(np.shape(np.atleast_1d(gradient[name])), (3,))

    def test_dispersive_object(self):
        # A metal reflector's position is the motivating case.  Before the
        # dispersive branch was routed, `geometry_addgradient` contracted only
        # the real non-dispersive tensor, so a Drude block's poles contributed
        # nothing and the gradient was of a structure that was not simulated.
        drude = mp.Medium(
            epsilon=1.0,
            E_susceptibilities=[
                mp.DrudeSusceptibility(frequency=1.0, gamma=0.2, sigma=0.6)
            ],
        )
        size = mp.Vector3(1.0, 0.3)
        offset = quarter_pixel(RES)

        def block(center, differentiable=None):
            return mp.Block(
                center=mp.Vector3(center.x, center.y + offset),
                size=mp.Vector3(size.x + 2 * offset, size.y),
                material=drude,
                differentiable=differentiable,
                name="scatterer",
            )

        opt = _problem(block(mp.Vector3(), differentiable=["center"]))
        _, grad = opt([RHO])
        adjoint = float(np.real(np.atleast_1d(grad["scatterer"]["center"])[1]))

        def value(dy):
            return float(
                np.asarray(
                    _problem(block(mp.Vector3(0, dy)))([RHO], need_gradient=False)[0]
                ).real.item()
            )

        reference = (value(FD_STEP) - value(-FD_STEP)) / (2 * FD_STEP)
        self.assertGreater(abs(reference), 1e-6, "objective must respond to the move")
        self.assertLess(
            abs(adjoint - reference) / abs(reference),
            2e-2,
            f"adjoint {adjoint} vs finite difference {reference}",
        )

    def test_object_with_an_infinite_extent(self):
        # `mp.inf` is the idiomatic way to write a layer that spans the cell,
        # and it is 1e20 rather than a flag. A monitor built around one without
        # clamping fails inside meep with "impossible(?) looping boundaries",
        # which is what a stratified stack -- a metal reflector under a grating,
        # say -- runs into immediately.
        offset = quarter_pixel(RES)

        def slab(dy, differentiable=None):
            return mp.Block(
                center=mp.Vector3(0, dy + offset),
                size=mp.Vector3(mp.inf, 0.4),
                material=mp.Medium(index=2.5),
                differentiable=differentiable,
                name="slab",
            )

        opt = _problem(slab(0.0, differentiable=["center"]))
        _, grad = opt([RHO])
        adjoint = float(np.real(np.atleast_1d(grad["slab"]["center"])[1]))

        def value(dy):
            return float(
                np.asarray(_problem(slab(dy))([RHO], need_gradient=False)[0]).item()
            )

        reference = (value(FD_STEP) - value(-FD_STEP)) / (2 * FD_STEP)
        self.assertGreater(abs(reference), 1e-6, "objective must respond to the move")
        self.assertLess(abs(adjoint - reference) / abs(reference), 2e-2)

    def test_no_flagged_object_leaves_the_return_shape_alone(self):
        opt = _problem(_block(mp.Vector3(), self.SIZE))
        _, gradient = opt([RHO])
        self.assertNotIsInstance(gradient, dict)


if __name__ == "__main__":
    unittest.main()
