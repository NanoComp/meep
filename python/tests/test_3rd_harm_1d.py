import unittest

from utils import ApproxComparisonTestCase

import meep as mp


class Test3rdHarm1d(ApproxComparisonTestCase):
    def setUp(self):
        self.sz = 100
        fcen = 1 / 3.0
        df = fcen / 20.0
        self.amp = 1.0
        self.k = 1e-2
        self.dpml = 1.0
        dimensions = 1
        cell = mp.Vector3(0, 0, self.sz)

        default_material = mp.Medium(index=1, chi3=self.k)

        pml_layers = mp.PML(self.dpml)

        sources = mp.Source(
            mp.GaussianSource(fcen, fwidth=df),
            component=mp.Ex,
            center=mp.Vector3(0, 0, (-0.5 * self.sz) + self.dpml),
            amplitude=self.amp,
        )

        nfreq = 400
        fmin = fcen / 2.0
        fmax = fcen * 4

        self.sim = mp.Simulation(
            cell_size=cell,
            geometry=[],
            sources=[sources],
            boundary_layers=[pml_layers],
            default_material=default_material,
            resolution=20,
            dimensions=dimensions,
        )

        fr = mp.FluxRegion(mp.Vector3(0, 0, (0.5 * self.sz) - self.dpml - 0.5))
        self.trans = self.sim.add_flux(
            0.5 * (fmin + fmax), fmax - fmin, nfreq, fr, decimation_factor=1
        )
        self.trans1 = self.sim.add_flux(fcen, 0, 1, fr, decimation_factor=1)
        self.trans3 = self.sim.add_flux(3 * fcen, 0, 1, fr, decimation_factor=1)

    def test_3rd_harm_1d(self):

        expected_harmonics = [0.01, 1.0, 221.89548712071553, 1.752960413399477]

        self.sim.run(
            until_after_sources=mp.stop_when_fields_decayed(
                50, mp.Ex, mp.Vector3(0, 0, (0.5 * self.sz) - self.dpml - 0.5), 1e-6
            )
        )

        harmonics = [
            self.k,
            self.amp,
            mp.get_fluxes(self.trans1)[0],
            mp.get_fluxes(self.trans3)[0],
        ]

        tol = 3e-5 if mp.is_single_precision() else 1e-7
        self.assertClose(expected_harmonics, harmonics, epsilon=tol)


if __name__ == "__main__":
    unittest.main()
