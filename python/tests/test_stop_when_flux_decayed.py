import unittest

import meep as mp
import numpy as np


class TestStopWhenFluxDecayed(unittest.TestCase):
    def setUp(self):
        self.resolution = 20
        self.dpml = 1.0
        self.sxy = 6.0
        self.cell_size = mp.Vector3(self.sxy + 2 * self.dpml, self.sxy + 2 * self.dpml)
        self.fcen = 1.0
        self.df = 0.2
        self.nfreq = 3
        self.sources = [
            mp.Source(
                src=mp.GaussianSource(self.fcen, fwidth=self.df),
                center=mp.Vector3(),
                component=mp.Ez,
            )
        ]
        self.pml_layers = [mp.PML(self.dpml)]

    def _make_sim(self):
        return mp.Simulation(
            cell_size=self.cell_size,
            resolution=self.resolution,
            sources=self.sources,
            boundary_layers=self.pml_layers,
            symmetries=[
                mp.Mirror(mp.X, phase=+1),
                mp.Mirror(mp.Y, phase=+1),
            ],
        )

    def test_single_frequency(self):
        """Simulation terminates and produces non-zero flux with one frequency."""
        sim = self._make_sim()
        flux_mon = sim.add_flux(
            self.fcen,
            0,
            1,
            mp.FluxRegion(
                center=mp.Vector3(0.5 * self.sxy, 0),
                size=mp.Vector3(0, self.sxy),
            ),
        )

        sim.run(until_after_sources=mp.stop_when_flux_decayed(flux_mon, tol=1e-3))

        fluxes = mp.get_fluxes(flux_mon)
        self.assertEqual(len(fluxes), 1)
        self.assertGreater(abs(fluxes[0]), 0)

    def test_multiple_frequencies(self):
        """All frequencies must converge before termination."""
        sim = self._make_sim()
        flux_mon = sim.add_flux(
            self.fcen,
            self.df,
            self.nfreq,
            mp.FluxRegion(
                center=mp.Vector3(0.5 * self.sxy, 0),
                size=mp.Vector3(0, self.sxy),
            ),
        )

        sim.run(until_after_sources=mp.stop_when_flux_decayed(flux_mon, tol=1e-3))

        fluxes = mp.get_fluxes(flux_mon)
        self.assertEqual(len(fluxes), self.nfreq)
        for f in fluxes:
            self.assertGreater(abs(f), 0)

    def test_maximum_run_time(self):
        """Simulation stops at maximum_run_time even if not converged."""
        sim = self._make_sim()
        flux_mon = sim.add_flux(
            self.fcen,
            0,
            1,
            mp.FluxRegion(
                center=mp.Vector3(0.5 * self.sxy, 0),
                size=mp.Vector3(0, self.sxy),
            ),
        )

        max_time = 50
        sim.run(
            until_after_sources=mp.stop_when_flux_decayed(
                flux_mon, tol=1e-30, maximum_run_time=max_time
            )
        )

        self.assertLessEqual(sim.round_time(), max_time + 1)

    def test_minimum_run_time(self):
        """Simulation runs at least minimum_run_time."""
        sim = self._make_sim()
        flux_mon = sim.add_flux(
            self.fcen,
            0,
            1,
            mp.FluxRegion(
                center=mp.Vector3(0.5 * self.sxy, 0),
                size=mp.Vector3(0, self.sxy),
            ),
        )

        min_time = 50
        sim.run(
            until_after_sources=mp.stop_when_flux_decayed(
                flux_mon, tol=1e-3, minimum_run_time=min_time
            )
        )

        self.assertGreaterEqual(sim.round_time(), min_time)

    def test_agrees_with_long_run(self):
        """Flux values from stop_when_flux_decayed match a long fixed-time run."""
        sim1 = self._make_sim()
        flux_mon1 = sim1.add_flux(
            self.fcen,
            self.df,
            self.nfreq,
            mp.FluxRegion(
                center=mp.Vector3(0.5 * self.sxy, 0),
                size=mp.Vector3(0, self.sxy),
            ),
        )
        sim1.run(until_after_sources=mp.stop_when_flux_decayed(flux_mon1, tol=1e-5))
        fluxes_decayed = np.array(mp.get_fluxes(flux_mon1))

        sim2 = self._make_sim()
        flux_mon2 = sim2.add_flux(
            self.fcen,
            self.df,
            self.nfreq,
            mp.FluxRegion(
                center=mp.Vector3(0.5 * self.sxy, 0),
                size=mp.Vector3(0, self.sxy),
            ),
        )
        sim2.run(until_after_sources=500)
        fluxes_long = np.array(mp.get_fluxes(flux_mon2))

        np.testing.assert_allclose(fluxes_decayed, fluxes_long, rtol=1e-3)

    def test_invalid_flux_argument(self):
        """Raises TypeError for non-DftFlux arguments."""
        with self.assertRaises(TypeError):
            mp.stop_when_flux_decayed("not_a_flux_object")

    def test_none_flux_argument(self):
        """Raises ValueError for None flux argument."""
        with self.assertRaises(ValueError):
            mp.stop_when_flux_decayed(None)


if __name__ == "__main__":
    unittest.main()
