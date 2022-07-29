import math
import unittest

import meep as mp


class TestMultiLevelAtom(unittest.TestCase):
    @unittest.skipIf(
        mp.is_single_precision(), "double-precision floating point specific test"
    )
    def test_multilevel_atom(self):
        resolution = 40
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

        def check_field(sim):
            fp = sim.get_field_point(
                mp.Ex, mp.Vector3(z=(-0.5 * sz) + Lcav + (0.5 * dpad))
            ).real
            self.assertAlmostEqual(fp, -2.7110969214986387)

        sim.init_sim()
        sim.initialize_field(mp.Ex, field_func)
        sim.run(mp.at_end(check_field), until=7000)


if __name__ == "__main__":
    unittest.main()
