import unittest

import meep as mp
from meep.geom import Cylinder, Lattice, Medium, Vector3
from meep.simulation import Mirror, Pml, Simulation, no_size, py_v3_to_vec
from meep.source import GaussianSource, Source


class TestSimulation(unittest.TestCase):

    def setUp(self):
        c = Cylinder(1)
        self.sim = Simulation(
            cell=Lattice(size=Vector3(16, 16, no_size)),
            geometry=[c],
            sources=[],
            resolution=10,
            pml_layers=[Pml(2)]
        )

        def say_hi():
            print("Hi")

        self.say_hi = say_hi

    def test_at_beginning(self):
        self.sim.run(self.sim.at_beginning(self.say_hi), until=300)

    def test_at_every(self):
        self.sim.run(self.sim.at_every(0.1, self.say_hi), until=5)


if __name__ == '__main__':
    unittest.main()
