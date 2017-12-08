import unittest

import meep as mp


def f(r, ex, hz, eps):
    return (r.x * r.norm() + ex) - (eps * hz)


class TestFieldFunctions(unittest.TestCase):

    cs = [mp.Ex, mp.Hz, mp.Dielectric]
    vol = mp.Volume(size=mp.Vector3(1), center=mp.Vector3())

    def setUp(self):
        resolution = 20

        cell = mp.Vector3(10, 10, 0)

        pml_layers = mp.PML(1.0)

        fcen = 1.0
        df = 1.0

        sources = mp.Source(src=mp.GaussianSource(fcen, fwidth=df), center=mp.Vector3(),
                            component=mp.Ez)

        symmetries = [mp.Mirror(mp.X), mp.Mirror(mp.Y)]

        self.sim = mp.Simulation(resolution=resolution,
                                 cell_size=cell,
                                 boundary_layers=[pml_layers],
                                 sources=[sources],
                                 symmetries=symmetries)

    def test_integrate_field_function(self):
        self.sim.run(until=200)

        res1 = self.sim.integrate_field_function(self.cs, f)
        res2 = self.sim.integrate_field_function(self.cs, f, self.vol)

        self.assertAlmostEqual(res1, complex(-6.938893903907228e-18, 0))
        self.assertAlmostEqual(res2, complex(0, 0))

        self.sim.output_field_function("weird-function", self.cs, f)

    def test_max_abs_field_function(self):
        self.sim.run(until=200)

        self.cs = [mp.Ex, mp.Hz, mp.Dielectric]
        res = self.sim.max_abs_field_function(self.cs, f, self.vol)
        self.assertAlmostEqual(res, 0.27593732304637586)

if __name__ == '__main__':
    unittest.main()
