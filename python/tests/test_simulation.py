import unittest

import meep as mp
from meep.geom import Cylinder, Lattice, Vector3
from meep.simulation import Pml, Simulation, no_size


class TestSimulation(unittest.TestCase):

    # def setUp(self):
    #     c = Cylinder(1)
    #     self.sim = Simulation(
    #         cell=Lattice(size=Vector3(16, 16, no_size)),
    #         geometry=[c],
    #         sources=[],
    #         resolution=10,
    #         pml_layers=[Pml(2)]
    #     )

    #     def say_hi():
    #         print("Hi")

    #     self.say_hi = say_hi

    # def test_at_beginning(self):
    #     self.sim.run(self.sim.at_beginning(self.say_hi), until=300)

    # def test_at_every(self):
    #     self.sim.run(self.sim.at_every(0.1, self.say_hi), until=5)

    def test_adaptive_integration(self):
        integration_args = [
            lambda u: u * u,
            [0.0],
            [1.0],
            1e-9,
            1e-4,
            50000,
        ]

        integral, abstol = mp.py_adaptive_integration2(*integration_args)
        print(integral, abstol)
        self.assertEqual(integral, 0.333333)
        # integral, abstol = mp.py_adaptive_integration(mp.py_pml_profile2u, *integration_args)


if __name__ == '__main__':
    unittest.main()
