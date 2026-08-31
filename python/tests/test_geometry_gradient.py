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
        geometry=list(block) if isinstance(block, (list, tuple)) else [block],
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
        # Ex is smoothed over voxels whose y edges lie at half-integer
        # multiples of dx, so a face there is straddled by none of them.
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
            any("voxel edge" in str(w.message) for w in caught),
            "expected a warning about a face on a pixel edge",
        )


class TestShapeDerivative(unittest.TestCase):
    SIZE = mp.Vector3(1.0, 0.8)

    def _value(self, center, size, res=RES, index=2.5, component=mp.Hz):
        opt = _problem(
            _block(center, size, index, res=res), res=res, component=component
        )
        return float(np.asarray(opt([RHO], need_gradient=False)[0]).item())

    def _adjoint(self, names, res=RES, index=2.5, component=mp.Hz):
        opt = _problem(
            _block(mp.Vector3(), self.SIZE, index, differentiable=names, res=res),
            res=res,
            component=component,
        )
        _, grad = opt([RHO])
        return grad["scatterer"]

    def test_both_axes_and_both_polarizations(self):
        # Hz drives the in-plane Ex and Ey, Ez the out-of-plane one; they see
        # different voxels and different smoothed tensor entries.
        unit = (mp.Vector3(x=FD_STEP), mp.Vector3(y=FD_STEP))
        for component in (mp.Hz, mp.Ez):
            grad = self._adjoint(["center", "size"], component=component)
            center = np.real(np.atleast_1d(grad["center"]))
            size = np.real(np.atleast_1d(grad["size"]))
            scale = np.linalg.norm([center[1], size[0], size[1]])
            # The source and monitor sit symmetrically about x = 0, so by
            # reciprocity the objective is even in the block's x position and
            # its derivative there is exactly zero.
            with self.subTest(component=mp.component_name(component), name="center x"):
                self.assertLess(abs(center[0]), 1e-3 * scale)
            for name, axis in (("center", 1), ("size", 0), ("size", 1)):
                with self.subTest(
                    component=mp.component_name(component), name=name, axis=axis
                ):
                    d = unit[axis]
                    if name == "center":
                        plus, minus = (d, self.SIZE), (mp.Vector3() - d, self.SIZE)
                    else:
                        plus = (mp.Vector3(), self.SIZE + d)
                        minus = (mp.Vector3(), self.SIZE - d)
                    reference = (
                        self._value(*plus, component=component)
                        - self._value(*minus, component=component)
                    ) / (2 * FD_STEP)
                    adjoint = (center if name == "center" else size)[axis]
                    self.assertLess(
                        abs(adjoint - reference) / abs(reference),
                        2e-2,
                        f"adjoint {adjoint} vs finite difference {reference}",
                    )

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

    def test_no_flagged_object_leaves_the_return_shape_alone(self):
        opt = _problem(_block(mp.Vector3(), self.SIZE))
        _, gradient = opt([RHO])
        self.assertNotIsInstance(gradient, dict)


class TestSeveralObjects(unittest.TestCase):
    """Two flagged objects, each with its own gradient."""

    def _blocks(self, shift_a=mp.Vector3(), shift_b=mp.Vector3(), flagged=True):
        q = quarter_pixel(RES)
        flags = ["center"] if flagged else None
        return [
            mp.Block(
                center=mp.Vector3(-0.6, q) + shift_a,
                size=mp.Vector3(0.5 + 2 * q, 0.6),
                material=mp.Medium(index=2.5),
                differentiable=flags,
                name="a",
            ),
            mp.Block(
                center=mp.Vector3(0.7 + q, 0.3 + q) + shift_b,
                size=mp.Vector3(0.4, 0.5),
                material=mp.Medium(index=2.0),
                differentiable=flags,
                name="b",
            ),
        ]

    def _value(self, **shifts):
        opt = _problem(self._blocks(flagged=False, **shifts))
        return float(np.asarray(opt([RHO], need_gradient=False)[0]).item())

    def test_each_object_gets_its_own_gradient(self):
        _, grad = _problem(self._blocks())([RHO])
        self.assertLessEqual({"a", "b"}, set(grad))
        checks = (
            ("a", 1, dict(shift_a=mp.Vector3(y=FD_STEP))),
            ("b", 0, dict(shift_b=mp.Vector3(x=FD_STEP))),
        )
        for name, axis, plus in checks:
            with self.subTest(object=name, axis=axis):
                minus = {k: mp.Vector3() - v for k, v in plus.items()}
                reference = (self._value(**plus) - self._value(**minus)) / (2 * FD_STEP)
                adjoint = float(np.real(np.atleast_1d(grad[name]["center"])[axis]))
                self.assertLess(abs(adjoint - reference) / abs(reference), 2e-2)


class TestThreeDimensions(unittest.TestCase):
    """A block in a small 3D cell."""

    RES3 = 12

    def _problem(self, dz, differentiable=None):
        q = quarter_pixel(self.RES3)
        block = mp.Block(
            center=mp.Vector3(0, q, q + dz),
            size=mp.Vector3(0.8 + 2 * q, 0.6, 0.5),
            material=mp.Medium(index=2.5),
            differentiable=differentiable,
            name="block",
        )
        sim = mp.Simulation(
            cell_size=mp.Vector3(4, 3, 3),
            resolution=self.RES3,
            boundary_layers=[mp.PML(0.6)],
            sources=[
                mp.Source(
                    mp.GaussianSource(FCEN, fwidth=0.2),
                    component=mp.Ex,
                    center=mp.Vector3(-1.0, 0.1, 0.1),
                )
            ],
            geometry=[block],
            force_complex_fields=True,
        )
        grid = mp.MaterialGrid(mp.Vector3(4, 4, 1), mp.air, mp.Medium(index=1.5))
        region = mpa.DesignRegion(
            grid,
            volume=mp.Volume(
                center=mp.Vector3(0, -1.0, 0), size=mp.Vector3(0.3, 0.3, 0.3)
            ),
        )
        sim.geometry = list(sim.geometry) + [
            mp.Block(center=region.center, size=region.size, material=grid)
        ]
        monitor = mpa.FourierFields(
            sim, mp.Volume(center=mp.Vector3(1.0, 0.2, -0.1)), mp.Ex
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

    def test_center_gradient_along_z(self):
        _, grad = self._problem(0.0, ["center"])([0.5 * np.ones(16)])
        adjoint = float(np.real(np.atleast_1d(grad["block"]["center"])[2]))

        def value(dz):
            opt = self._problem(dz)
            return float(
                np.asarray(opt([0.5 * np.ones(16)], need_gradient=False)[0]).item()
            )

        reference = (value(FD_STEP) - value(-FD_STEP)) / (2 * FD_STEP)
        self.assertLess(abs(adjoint - reference) / abs(reference), 2e-2)


class TestEdgeWarningPerComponent(unittest.TestCase):
    """Voxel edges depend on the component; only components with field count."""

    def _warnings(self, component):
        # x faces at integer multiples of dx: a voxel edge for Ex (along its
        # own direction) but a voxel centre for Ey and Ez
        block = mp.Block(
            center=mp.Vector3(0, quarter_pixel(RES)),
            size=mp.Vector3(1.0, 0.8),
            material=mp.Medium(index=2.5),
            differentiable=["center"],
            name="s",
        )
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            _problem(block, component=component)([RHO])
        return [str(w.message) for w in caught if "voxel edge" in str(w.message)]

    def test_in_plane_field_warns(self):
        messages = self._warnings(mp.Hz)  # drives Ex and Ey
        self.assertTrue(any(" ex" in m for m in messages), messages)

    def test_out_of_plane_field_does_not(self):
        self.assertEqual(self._warnings(mp.Ez), [])  # drives Ez only


if __name__ == "__main__":
    unittest.main()
