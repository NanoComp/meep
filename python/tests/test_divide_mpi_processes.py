import unittest

import meep as mp


@unittest.skipIf(mp.count_processors() < 2, "MPI specific test")
class TestDivideParallelProcesses(unittest.TestCase):
    def test_divide_parallel_processes(self):
        resolution = 20

        sxy = 4
        dpml = 1
        cell = mp.Vector3(sxy + 2 * dpml, sxy + 2 * dpml)

        pml_layers = [mp.PML(dpml)]

        n = mp.divide_parallel_processes(2)
        fcen = 1.0 / (n + 1)

        sources = [
            mp.Source(
                src=mp.GaussianSource(fcen, fwidth=0.2 * fcen),
                center=mp.Vector3(),
                component=mp.Ez,
            )
        ]

        symmetries = [mp.Mirror(mp.X), mp.Mirror(mp.Y)]

        self.sim = mp.Simulation(
            cell_size=cell,
            resolution=resolution,
            sources=sources,
            symmetries=symmetries,
            boundary_layers=pml_layers,
        )

        flux_box = self.sim.add_flux(
            fcen,
            0,
            1,
            mp.FluxRegion(mp.Vector3(y=0.5 * sxy), size=mp.Vector3(sxy)),
            mp.FluxRegion(mp.Vector3(y=-0.5 * sxy), size=mp.Vector3(sxy), weight=-1),
            mp.FluxRegion(mp.Vector3(0.5 * sxy), size=mp.Vector3(y=sxy)),
            mp.FluxRegion(mp.Vector3(-0.5 * sxy), size=mp.Vector3(y=sxy), weight=-1),
            decimation_factor=1,
        )

        self.sim.run(until_after_sources=30)

        tot_flux = mp.get_fluxes(flux_box)[0]

        tot_fluxes = mp.merge_subgroup_data(tot_flux)
        fcens = mp.merge_subgroup_data(fcen)

        self.assertEqual(fcens[0], 1)
        self.assertEqual(fcens[1], 0.5)
        places = 4 if mp.is_single_precision() else 7
        self.assertAlmostEqual(tot_fluxes[0], 9.8628728533, places=places)
        self.assertAlmostEqual(tot_fluxes[1], 19.6537275387, places=places)


if __name__ == "__main__":
    unittest.main()
