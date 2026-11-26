import unittest

import numpy as np

import meep as mp


class TestSolveCw(unittest.TestCase):
    def setUp(self):
        mp.verbosity(0)

    def tearDown(self):
        mp.verbosity(1)

    @staticmethod
    def slab_simulation(**source_kwargs):
        """A 1D cell with a dielectric slab and a CW line source."""
        return mp.Simulation(
            cell_size=mp.Vector3(z=12),
            dimensions=1,
            resolution=40,
            boundary_layers=[mp.PML(2.0)],
            force_complex_fields=True,
            geometry=[
                mp.Block(
                    center=mp.Vector3(z=2),
                    size=mp.Vector3(mp.inf, mp.inf, 1),
                    material=mp.Medium(epsilon=4),
                )
            ],
            sources=[
                mp.Source(
                    mp.ContinuousSource(1.0, **source_kwargs),
                    component=mp.Ex,
                    center=mp.Vector3(z=-1),
                )
            ],
        )

    def test_matches_time_domain(self):
        """The solve_cw fields equal the steady-state time-domain fields."""
        sim = self.slab_simulation(width=5)
        sim.run(until=300)
        ex_td = sim.get_array(component=mp.Ex)
        hy_td = sim.get_array(component=mp.Hy)
        t0 = sim.meep_time()

        self.assertTrue(sim.solve_cw())
        self.assertEqual(sim.meep_time(), t0)

        ex_cw = sim.get_array(component=mp.Ex)
        hy_cw = sim.get_array(component=mp.Hy)

        # Both fields are at the same time step as the time-domain fields, so
        # an error in the phase of the CW solution (e.g., a solution offset by
        # one time step, |1 - exp(-2i*pi*f*dt)| ~ 0.08) is easily detected.
        tol = 1e-3 if mp.is_single_precision() else 1e-4
        for name, cw, td in (("Ex", ex_cw, ex_td), ("Hy", hy_cw, hy_td)):
            err = np.linalg.norm(cw - td) / np.linalg.norm(td)
            self.assertLess(err, tol, f"{name}: relative error {err:.2e}")

    # With more than one process, meep::abort calls MPI_Abort instead of
    # raising a catchable exception, so there is no state left to check.
    @unittest.skipIf(mp.count_processors() > 1, "single-process test")
    def test_state_restored_after_failure(self):
        """A failed solve_cw leaves the time, sources, and DFTs working."""
        # A source that has not yet turned on makes solve_cw raise an
        # exception after it has already switched the fields into CW mode.
        sim = self.slab_simulation(start_time=1)
        dft = sim.add_dft_fields(
            [mp.Ex], 1.0, 0, 1, center=mp.Vector3(), size=mp.Vector3()
        )
        sim.init_sim()
        t0 = sim.meep_time()

        with self.assertRaises(RuntimeError):
            sim.solve_cw()
        self.assertEqual(sim.meep_time(), t0)

        sim.run(until=5)
        self.assertGreater(np.abs(sim.get_array(component=mp.Ex)).max(), 0)
        self.assertGreater(np.abs(sim.get_dft_array(dft, mp.Ex, 0)).max(), 0)


if __name__ == "__main__":
    unittest.main()
