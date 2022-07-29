import unittest

import meep as mp


# Material function that recreates the ellipsoid-in-cylinder configuration of
# examples/cyl-ellipsoid.py
def my_material_func(p):
    R1X = 0.5
    R1Y = 1.0
    R2 = 3.0

    x = p.x
    y = p.y

    # test for point inside inner ellipsoid
    if (x**2 / (R1X**2) + y**2 / (R1Y**2)) < 1.0:
        nn = 1.0
    elif (x**2 / (R2**2) + y**2 / (R2**2)) < 1.0:
        nn = 3.5
    else:
        nn = 1.0

    return mp.Medium(epsilon=nn**2)


def my_epsilon_func(p):
    R1X = 0.5
    R1Y = 1.0
    R2 = 3.0

    x = p.x
    y = p.y

    if (x**2 / (R1X**2) + y**2 / (R1Y**2)) < 1.0:
        return 1.0
    elif (x**2 / (R2**2) + y**2 / (R2**2)) < 1.0:
        return 3.5
    return 1.0


class TestUserMaterials(unittest.TestCase):
    def setUp(self):
        self.resolution = 10
        self.cell = mp.Vector3(10, 10)
        self.symmetries = [mp.Mirror(mp.X), mp.Mirror(mp.Y)]
        self.boundary_layers = [mp.PML(1.0)]
        self.sources = [
            mp.Source(
                src=mp.GaussianSource(0.2, fwidth=0.1),
                component=mp.Ez,
                center=mp.Vector3(),
            )
        ]

    def test_user_material_func(self):
        sim = mp.Simulation(
            cell_size=self.cell,
            resolution=self.resolution,
            symmetries=self.symmetries,
            boundary_layers=self.boundary_layers,
            sources=self.sources,
            material_function=my_material_func,
        )

        sim.run(until=200)
        fp = sim.get_field_point(mp.Ez, mp.Vector3(x=1))

        places = 3 if mp.is_single_precision() else 7
        self.assertAlmostEqual(fp, 4.816403627871773e-4, places=places)

    def test_epsilon_func(self):
        sim = mp.Simulation(
            cell_size=self.cell,
            resolution=self.resolution,
            symmetries=self.symmetries,
            boundary_layers=self.boundary_layers,
            sources=self.sources,
            epsilon_func=my_epsilon_func,
        )

        sim.run(until=100)
        fp = sim.get_field_point(mp.Ez, mp.Vector3(x=1))

        self.assertAlmostEqual(fp, -7.895783750440999e-4)

    def test_geometric_obj_with_user_material(self):
        geometry = [mp.Cylinder(5, material=my_material_func)]

        sim = mp.Simulation(
            cell_size=self.cell,
            resolution=self.resolution,
            symmetries=self.symmetries,
            geometry=geometry,
            boundary_layers=self.boundary_layers,
            sources=self.sources,
        )
        sim.run(until=200)
        fp = sim.get_field_point(mp.Ez, mp.Vector3(x=1))

        places = 3 if mp.is_single_precision() else 7
        self.assertAlmostEqual(fp, 4.816403627871773e-4, places=places)

    def test_geometric_obj_with_epsilon_func(self):
        geometry = [mp.Cylinder(5, epsilon_func=my_epsilon_func)]

        sim = mp.Simulation(
            cell_size=self.cell,
            resolution=self.resolution,
            symmetries=self.symmetries,
            geometry=geometry,
            boundary_layers=self.boundary_layers,
            sources=self.sources,
        )
        sim.run(until=100)
        fp = sim.get_field_point(mp.Ez, mp.Vector3(x=1))

        self.assertAlmostEqual(fp, -7.895783750440999e-4)


if __name__ == "__main__":
    unittest.main()
