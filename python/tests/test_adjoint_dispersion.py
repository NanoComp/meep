"""Adjoint gradients of dispersive structures against finite differences.

`get_chi1_tensor_disp` used to build the dispersive permittivity from the
*continuum* lineshape -- `lorentzian_susceptibility::chi1` and
`1 + i*sigma/(2*pi*f)` -- while the FDTD timesteps a discrete recurrence. The
adjoint gradient was then the exact derivative of an operator slightly
different from the one being simulated, and disagreed with a finite difference
of the discrete objective at `O((freq*dt)^2)`.

Measured relative disagreement at `res` 20 / 40 / 80, continuum against
discrete:

    D-conductivity   0.455%  0.107%  0.027%   ->  0.001%  0.000%  0.001%
    Drude pole       0.232%  0.058%  0.015%   ->  0.000%  0.000%  0.000%

The `4x`-per-doubling of the left-hand columns is the second-order signature;
the right-hand columns sit at the same roundoff floor as a non-dispersive
structure. `TOL` below is set from those numbers: loose enough to be immune to
roundoff, tight enough that reverting the fix fails at every resolution tested.

Run lengths are pinned rather than left to `stop_when_dft_decayed`: an adaptive
stop makes a perturbed run end at a different time, and the difference scales
with the perturbation.
"""

import unittest

import numpy as np
from autograd import numpy as npa

import meep as mp
import meep.adjoint as mpa

FCEN = 1.0
RUN = 150.0
CELL = mp.Vector3(6, 4)
NG = 4
IDX = 5
U0 = 0.3
FD_STEP = 1e-4

# Comfortably below the 0.232% that the continuum lineshape produces at the
# coarsest resolution tested here, and ~100x above the roundoff floor the fix
# achieves. See the measurements in the module docstring.
TOL = 5e-4

DRUDE = mp.Medium(
    epsilon=1.0,
    E_susceptibilities=[mp.DrudeSusceptibility(frequency=1.0, gamma=0.1, sigma=0.3)],
)
LORENTZ = mp.Medium(
    epsilon=1.0,
    E_susceptibilities=[
        mp.LorentzianSusceptibility(frequency=1.3, gamma=0.2, sigma=0.5)
    ],
)


def _problem(medium2, res, damping=0.0):
    grid = mp.MaterialGrid(
        mp.Vector3(NG, NG), mp.air, medium2, do_averaging=False, beta=0, damping=damping
    )
    region = mpa.DesignRegion(
        grid, volume=mp.Volume(center=mp.Vector3(), size=mp.Vector3(1.0, 1.0))
    )
    sim = mp.Simulation(
        cell_size=CELL,
        resolution=res,
        boundary_layers=[mp.PML(1.0)],
        sources=[
            mp.Source(
                mp.GaussianSource(FCEN, fwidth=0.2),
                component=mp.Ez,
                center=mp.Vector3(-1.8, 0),
            )
        ],
        geometry=[mp.Block(center=region.center, size=region.size, material=grid)],
        eps_averaging=False,
        force_complex_fields=True,
    )
    monitor = mpa.FourierFields(
        sim, mp.Volume(center=mp.Vector3(1.8, 0), size=mp.Vector3()), mp.Ez
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


class TestDispersiveGradient(unittest.TestCase):
    def assertGradientMatchesFD(self, medium2, res=20, damping=0.0):
        rho = U0 * np.ones(NG * NG)
        _, grad = _problem(medium2, res, damping)([rho])
        adjoint = float(np.real(np.atleast_1d(np.squeeze(grad))[IDX]))

        def value(r):
            return float(
                np.asarray(
                    _problem(medium2, res, damping)([r], need_gradient=False)[0]
                ).real.item()
            )

        hi, lo = rho.copy(), rho.copy()
        hi[IDX] += FD_STEP
        lo[IDX] -= FD_STEP
        reference = (value(hi) - value(lo)) / (2 * FD_STEP)

        # Guards against the comparison passing vacuously: a structure the
        # objective does not respond to would make any adjoint value "agree".
        self.assertGreater(abs(reference), 1e-6, "objective does not respond to rho")
        self.assertLess(
            abs(adjoint - reference) / abs(reference),
            TOL,
            f"adjoint {adjoint} vs finite difference {reference}",
        )
        return adjoint, reference

    def test_drude_pole(self):
        self.assertGradientMatchesFD(DRUDE)

    def test_lorentzian_pole(self):
        # The Drude branch of chi1_discrete() drops the omega_0^2 term, so a
        # resonant pole exercises a path the Drude test does not reach.
        self.assertGradientMatchesFD(LORENTZ)

    def test_d_conductivity(self):
        # `damping` adds sigma = u(1-u)*damping to D_conductivity, which reaches
        # the operator through conductivity_factor() rather than through a pole.
        self.assertGradientMatchesFD(mp.Medium(index=1.5), damping=20.0)

    def test_coarse_resolution(self):
        # The sharpest of these. The error being removed is O((freq*dt)^2), so
        # it is largest where dt is largest; a formula that is merely closer
        # than the continuum one would still converge and still pass at res 20.
        self.assertGradientMatchesFD(DRUDE, res=10)
        self.assertGradientMatchesFD(mp.Medium(index=1.5), res=10, damping=20.0)


class TestReportedEpsilonIsContinuum(unittest.TestCase):
    """`chi1` describes the material; only the adjoint wants the discretization.

    Guards the split: the discrete lineshape lives in `chi1_discrete`, and
    anything reporting material properties to the user must keep using the
    continuum one, or `get_epsilon` would start returning resolution-dependent
    numbers for a resolution-independent material.
    """

    def test_get_epsilon_matches_the_medium_at_two_resolutions(self):
        reference = complex(np.squeeze(DRUDE.epsilon(FCEN))[0, 0])
        for res in (10, 40):
            sim = mp.Simulation(
                cell_size=mp.Vector3(2, 2),
                resolution=res,
                default_material=DRUDE,
                force_complex_fields=True,
            )
            sim.init_sim()
            eps = np.mean(sim.get_epsilon(frequency=FCEN))
            self.assertAlmostEqual(eps.real, reference.real, places=5)
            self.assertAlmostEqual(eps.imag, reference.imag, places=5)


if __name__ == "__main__":
    unittest.main()
