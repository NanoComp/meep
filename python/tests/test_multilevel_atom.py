import math
import unittest

from utils import ApproxComparisonTestCase

import meep as mp


class TestMultiLevelAtom(ApproxComparisonTestCase):
    @unittest.skipIf(
        mp.is_single_precision(), "double-precision floating point specific test"
    )
    def test_multilevel_atom(self):
        resolution = 200
        ncav = 1.5
        Lcav = 1
        dpad = 1
        dpml = 1
        sz = Lcav + dpad + dpml

        cell_size = mp.Vector3(z=sz)
        dimensions = 1
        pml_layers = [mp.PML(dpml, side=mp.High)]

        omega_a = 40
        freq_21 = omega_a / (2 * math.pi)

        gamma_perp = 4
        gamma_21 = (2 * gamma_perp) / (2 * math.pi)

        theta = 1
        sigma_21 = 2 * theta * theta * omega_a

        rate_21 = 0.005
        N0 = 28
        Rp = 0.0051

        t1 = mp.Transition(
            1,
            2,
            pumping_rate=Rp,
            frequency=freq_21,
            gamma=gamma_21,
            sigma_diag=mp.Vector3(sigma_21, sigma_21, sigma_21),
        )
        t2 = mp.Transition(2, 1, transition_rate=rate_21)
        ml_atom = mp.MultilevelAtom(
            sigma=1, transitions=[t1, t2], initial_populations=[N0]
        )
        two_level = mp.Medium(index=ncav, E_susceptibilities=[ml_atom])

        geometry = [
            mp.Block(
                center=mp.Vector3(z=(-0.5 * sz) + (0.5 * Lcav)),
                size=mp.Vector3(mp.inf, mp.inf, Lcav),
                material=two_level,
            )
        ]

        sim = mp.Simulation(
            cell_size=cell_size,
            resolution=resolution,
            boundary_layers=pml_layers,
            geometry=geometry,
            dimensions=dimensions,
        )

        def field_func(p):
            return 1 if p.z == (-0.5 * sz) + (0.5 * Lcav) else 0

        boundary = mp.Vector3(z=(-0.5 * sz) + Lcav + (0.5 * dpad))
        envelope = []

        def record_field(sim):
            envelope.append(abs(sim.get_field_point(mp.Ex, boundary).real))

        sim.init_sim()
        sim.initialize_field(mp.Ex, field_func)
        # Let the laser settle, then measure the steady-state envelope at the
        # cavity boundary, as Opt. Express 20, 474 (2012) does.  Steady-state
        # lasing fixes the amplitude but leaves the phase free, so the envelope
        # is reproducible where an instantaneous field sample is not.
        sim.run(until=2400)
        sim.run(record_field, until=100)
        self.assertClose(max(envelope), 0.103733, epsilon=1e-4)


if __name__ == "__main__":
    unittest.main()
