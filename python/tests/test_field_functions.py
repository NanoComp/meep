import unittest

import meep as mp


def f(r, ex, hz, eps):
    return (r.x * r.norm() + ex) - (eps * hz)


def f2(r, ez1, ez2):
    return ez1.conjugate() * ez2


class TestFieldFunctions(unittest.TestCase):

    cs = [mp.Ex, mp.Hz, mp.Dielectric]
    vol = mp.Volume(size=mp.Vector3(1), center=mp.Vector3())

    @classmethod
    def setUpClass(cls):
        cls.temp_dir = mp.make_output_directory()

    @classmethod
    def tearDownClass(cls):
        mp.delete_directory(cls.temp_dir)

    def init(self):
        resolution = 20

        cell = mp.Vector3(10, 10, 0)

        pml_layers = mp.PML(1.0)

        fcen = 1.0
        df = 1.0

        sources = mp.Source(
            src=mp.GaussianSource(fcen, fwidth=df), center=mp.Vector3(), component=mp.Ez
        )

        symmetries = [mp.Mirror(mp.X), mp.Mirror(mp.Y)]

        return mp.Simulation(
            resolution=resolution,
            cell_size=cell,
            boundary_layers=[pml_layers],
            sources=[sources],
            symmetries=symmetries,
        )

    def init2(self):
        n = 3.4
        w = 1
        r = 1
        pad = 4
        dpml = 2
        sxy = 2 * (r + w + pad + dpml)
        cell = mp.Vector3(sxy, sxy)

        geometry = [
            mp.Cylinder(radius=r + w, height=mp.inf, material=mp.Medium(index=n)),
            mp.Cylinder(radius=r, height=mp.inf, material=mp.air),
        ]

        pml_layers = [mp.PML(dpml)]
        resolution = 5
        fcen = 0.118
        df = 0.010

        sources = [
            mp.Source(
                src=mp.GaussianSource(fcen, fwidth=df),
                component=mp.Ez,
                center=mp.Vector3(r + 0.1),
            )
        ]

        symmetries = [mp.Mirror(mp.Y)]

        return mp.Simulation(
            cell_size=cell,
            resolution=resolution,
            geometry=geometry,
            boundary_layers=pml_layers,
            sources=sources,
            symmetries=symmetries,
        )

    def test_integrate_field_function(self):
        sim = self.init()
        sim.use_output_directory(self.temp_dir)
        sim.run(until=200)

        res1 = sim.integrate_field_function(self.cs, f)
        res2 = sim.integrate_field_function(self.cs, f, self.vol)

        self.assertAlmostEqual(res1, complex(-6.938893903907228e-18, 0.0))
        self.assertAlmostEqual(res2, 0.0j)

        sim.output_field_function("weird-function", self.cs, f)

    def test_integrate2_field_function(self):
        sim = self.init2()
        sim.run(until_after_sources=10)
        fields2 = sim.fields
        sim.reset_meep()
        sim.run(until_after_sources=10)

        res1 = sim.integrate2_field_function(fields2, [mp.Ez], [mp.Ez], f2)
        places = 6 if mp.is_single_precision() else 7
        self.assertAlmostEqual(res1, 0.17158099566244897, places=places)

    def test_max_abs_field_function(self):
        sim = self.init()
        sim.run(until=200)

        self.cs = [mp.Ex, mp.Hz, mp.Dielectric]
        res = sim.max_abs_field_function(self.cs, f, self.vol)
        self.assertAlmostEqual(res, 0.27593732304637586)


if __name__ == "__main__":
    unittest.main()
