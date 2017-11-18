import unittest

import meep as mp


class TestFieldFunctions(unittest.TestCase):

    def setUp(self):
        resolution = 20

        cell = mp.Vector3(10, 10, 0)

        pml_layers = mp.PML(1.0)

        self.fcen = 1.0
        df = 1.0

        sources = mp.Source(src=mp.GaussianSource(self.fcen, fwidth=df), center=mp.Vector3(),
                            component=mp.Ez)

        symmetries = [mp.Mirror(mp.X), mp.Mirror(mp.Y)]

        self.sim = mp.Simulation(resolution=resolution,
                                 cell_size=cell,
                                 boundary_layers=[pml_layers],
                                 sources=[sources],
                                 symmetries=symmetries)

    def test_field_functions(self):
        self.sim.run(until=200)

        def f(r, ex, hz, eps):
            return (r.x * r.norm() + ex) - (eps * hz)

        cs = [mp.Ex, mp.Hz, mp.Dielectric]
        res1 = self.sim.integrate_field_function(cs, f)
        vol = mp.Volume(size=mp.Vector3(1), center=mp.Vector3())
        res2 = self.sim.integrate_field_function(cs, f, vol)
        res3 = self.sim.max_abs_field_function(cs, f, vol)

        self.assertAlmostEqual(res1, complex(-6.938893903907228e-18, 0))
        self.assertAlmostEqual(res2, complex(0, 0))
        self.assertAlmostEqual(res3, 0.27593732304637586)


if __name__ == '__main__':
    unittest.main()
