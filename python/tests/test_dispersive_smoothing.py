"""Subpixel averaging of susceptibility amplitudes at a geometric interface.

The instantaneous permittivity has been subpixel-smoothed since forever while
sigma was point sampled, so a dispersive object's response switched on a whole
pixel at a time as it moved.  `eff_sigma_row` tapers each medium's pole
amplitude by the pixel's filling fraction instead.

Every medium here has `eps_infinity = 1`, matching the background, so the
instantaneous permittivity is uniform and contributes no position dependence at
all: the whole signal in these tests comes from sigma.  `eps_averaging=False`
routes back to the point sample, so the before/after comparison comes out of
one binary.
"""

import unittest

import numpy as np

import meep as mp

RES = 20
FCEN = 1.0
RUN = 120.0
DX = 1.0 / RES

DRUDE = mp.Medium(
    epsilon=1.0,
    E_susceptibilities=[mp.DrudeSusceptibility(frequency=1.0, gamma=0.5, sigma=1.0)],
)
LOSSLESS = mp.Medium(
    epsilon=1.0,
    E_susceptibilities=[
        mp.LorentzianSusceptibility(frequency=1.6, gamma=0.0, sigma=0.6)
    ],
)


def transmission(offset, averaging=True, res=RES, material=DRUDE, background=mp.air):
    sim = mp.Simulation(
        cell_size=mp.Vector3(6, 4),
        resolution=res,
        boundary_layers=[mp.PML(1.0)],
        default_material=background,
        sources=[
            mp.Source(
                mp.GaussianSource(FCEN, fwidth=0.2),
                component=mp.Ez,
                center=mp.Vector3(-1.8, 0),
            )
        ],
        geometry=[
            mp.Block(
                center=mp.Vector3(offset, 0),
                size=mp.Vector3(0.5, mp.inf),
                material=material,
            )
        ],
        eps_averaging=averaging,
        force_complex_fields=True,
    )
    mon = sim.add_dft_fields(
        [mp.Ez], FCEN, 0, 1, center=mp.Vector3(1.8, 0), size=mp.Vector3()
    )
    sim.run(until=RUN)
    return float(np.abs(sim.get_dft_array(mon, mp.Ez, 0)) ** 2)


def step_ratio(curve):
    """Largest single step relative to the average step.

    A staircase puts all of its variation in one place, so this is close to the
    number of samples; a smooth curve sampled this finely keeps it near 1.
    """
    steps = np.abs(np.diff(np.asarray(curve)))
    return float(steps.max() / steps.mean()) if steps.mean() > 0 else np.inf


class TestSubpixelContinuity(unittest.TestCase):
    OFFSETS = np.linspace(0, DX, 11)

    def setUp(self):
        mp.verbosity(0)

    def test_point_sampled_sigma_is_a_staircase(self):
        # Not a claim about desired behaviour -- it pins what the comparison
        # below is measuring against, so a change that made both paths smooth
        # for some unrelated reason could not pass vacuously.
        curve = [transmission(o, averaging=False) for o in self.OFFSETS]
        identical = sum(a == b for a, b in zip(curve, curve[1:]))
        self.assertGreaterEqual(
            identical, len(curve) - 3, f"expected a staircase, got {curve}"
        )

    def test_averaged_sigma_varies_continuously(self):
        curve = [transmission(o, averaging=True) for o in self.OFFSETS]
        self.assertEqual(
            sum(a == b for a, b in zip(curve, curve[1:])),
            0,
            f"adjacent offsets should differ: {curve}",
        )
        # Measured ~1.9 for the averaged curve and ~10 (i.e. one step carries
        # everything) for the point-sampled one.
        self.assertLess(step_ratio(curve), 3.0, f"curve looks stepped: {curve}")

    def test_averaging_shrinks_the_discretization_sensitivity(self):
        # The structure is the same slab throughout; all variation across the
        # sweep is discretization error. Smoothing should reduce it a lot.
        span = lambda c: max(c) - min(c)
        rough = span([transmission(o, averaging=False) for o in self.OFFSETS])
        smooth = span([transmission(o, averaging=True) for o in self.OFFSETS])
        self.assertLess(smooth, rough / 10, f"span {smooth} vs {rough}")

    def test_translating_by_a_whole_pixel_reproduces_the_structure(self):
        # Exact invariant: shifting by one pixel maps the grid onto itself, so
        # any discrepancy is a bug in the fill rather than discretization.
        a, b = transmission(0.0), transmission(DX)
        self.assertAlmostEqual(a / b, 1.0, places=5)


class TestExclusions(unittest.TestCase):
    """Pole kinds whose amplitude a filling fraction cannot meaningfully scale.

    A multilevel atom carries per-pixel level populations with their own
    nonlinear dynamics, so a partially filled pixel is not a weaker gain
    medium. A gyrotropic sigma is antisymmetric, so the argument that a
    non-negative combination of two tensors stays positive semidefinite -- the
    reason the linear mixture is safe -- does not apply. Both stay point
    sampled, and these pin that rather than leaving it to a comment.
    """

    def setUp(self):
        mp.verbosity(0)

    def assertStillStepped(self, susceptibility):
        medium = mp.Medium(epsilon=1.0, E_susceptibilities=[susceptibility])
        curve = [
            transmission(o, averaging=True, material=medium)
            for o in (0.1 * DX, 0.3 * DX, 0.5 * DX)
        ]
        self.assertEqual(
            curve[0], curve[1], f"pole should not have been smoothed: {curve}"
        )
        self.assertEqual(curve[1], curve[2], f"pole should not been smoothed: {curve}")

    def test_multilevel_atom_is_not_smoothed(self):
        self.assertStillStepped(
            mp.MultilevelAtom(
                sigma=1,
                transitions=[
                    mp.Transition(1, 2, pumping_rate=0.005, frequency=FCEN, gamma=1e-5)
                ],
                initial_populations=[1],
            )
        )

    def test_gyrotropic_is_not_smoothed(self):
        self.assertStillStepped(
            mp.GyrotropicLorentzianSusceptibility(
                frequency=1.1, gamma=0.2, sigma=0.5, bias=mp.Vector3(0, 0, 1)
            )
        )


class TestInvariants(unittest.TestCase):
    def setUp(self):
        mp.verbosity(0)

    def test_a_block_of_the_background_medium_is_invisible(self):
        # Exact, and independent of any tolerance: an object made of what is
        # already there cannot change the answer no matter where it sits.
        a = transmission(0.0, material=DRUDE, background=DRUDE)
        b = transmission(0.4 * DX, material=DRUDE, background=DRUDE)
        self.assertAlmostEqual(a / b, 1.0, places=9)

    def test_a_lossless_inclusion_does_not_pump_energy(self):
        # Passivity, end to end.  Mixing two media's sigma tensors is only safe
        # if the result stays positive semidefinite; scaling each by a
        # non-negative fill and adding does, but this is the check that would
        # catch it if the construction were ever changed to something that
        # does not.  A closed cavity with gamma = 0 has nowhere to lose energy,
        # so the field must not grow without bound.
        sim = mp.Simulation(
            cell_size=mp.Vector3(3, 3),
            resolution=RES,
            geometry=[
                mp.Block(
                    center=mp.Vector3(0.4 * DX, 0.4 * DX),
                    size=mp.Vector3(1.0, 1.0),
                    material=LOSSLESS,
                )
            ],
            sources=[
                mp.Source(
                    mp.GaussianSource(FCEN, fwidth=0.4),
                    component=mp.Ez,
                    center=mp.Vector3(-0.7, 0.3),
                )
            ],
            eps_averaging=True,
        )
        sim.run(until=50)
        early = sim.field_energy_in_box(mp.Volume(mp.Vector3(), size=mp.Vector3(3, 3)))
        sim.run(until=400)
        late = sim.field_energy_in_box(mp.Volume(mp.Vector3(), size=mp.Vector3(3, 3)))
        self.assertLess(late, 1.05 * early, f"energy grew: {early} -> {late}")


if __name__ == "__main__":
    unittest.main()
