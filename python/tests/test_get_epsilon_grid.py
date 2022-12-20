import unittest

import meep.adjoint as mpa
import numpy as np
import parameterized
from meep.materials import Co, SiN

import meep as mp


class TestGetEpsilonGrid(unittest.TestCase):
    def setUp(self):
        resolution = 60
        self.cell_size = mp.Vector3(1.0, 1.0, 0)

        matgrid_resolution = 200
        matgrid_size = mp.Vector3(1.0, 1.0, mp.inf)
        Nx, Ny = int(matgrid_size.x * matgrid_resolution), int(
            matgrid_size.y * matgrid_resolution
        )
        x = np.linspace(-0.5 * matgrid_size.x, 0.5 * matgrid_size.x, Nx)
        y = np.linspace(-0.5 * matgrid_size.y, 0.5 * matgrid_size.y, Ny)
        xv, yv = np.meshgrid(x, y)
        rad = 0.201943
        w = 0.104283
        weights = np.logical_and(
            np.sqrt(np.square(xv) + np.square(yv)) > rad,
            np.sqrt(np.square(xv) + np.square(yv)) < rad + w,
        )

        matgrid = mp.MaterialGrid(
            mp.Vector3(Nx, Ny),
            mp.air,
            mp.Medium(index=3.5),
            weights=weights,
            do_averaging=False,
            beta=0,
            eta=0.5,
        )

        geometry = [
            mp.Cylinder(
                center=mp.Vector3(0.35, 0.1),
                radius=0.1,
                height=mp.inf,
                material=mp.Medium(index=1.5),
            ),
            mp.Block(
                center=mp.Vector3(-0.15, -0.2),
                size=mp.Vector3(0.2, 0.24, mp.inf),
                material=SiN,
            ),
            mp.Block(
                center=mp.Vector3(-0.2, 0.2),
                size=mp.Vector3(0.4, 0.4, mp.inf),
                material=matgrid,
            ),
            mp.Prism(
                vertices=[
                    mp.Vector3(0.05, 0.45),
                    mp.Vector3(0.32, 0.22),
                    mp.Vector3(0.15, 0.10),
                ],
                height=0.5,
                material=Co,
            ),
        ]

        self.sim = mp.Simulation(
            resolution=resolution,
            cell_size=self.cell_size,
            geometry=geometry,
            eps_averaging=False,
        )
        self.sim.init_sim()

    @parameterized.parameterized.expand(
        [
            (mp.Vector3(0.2, 0.2), 1.1),
            (mp.Vector3(-0.2, 0.1), 0.7),
            (mp.Vector3(-0.2, -0.25), 0.55),
            (mp.Vector3(0.4, 0.1), 0),
        ]
    )
    def test_get_epsilon_grid(self, pt, freq):
        eps_grid = self.sim.get_epsilon_grid(
            np.array([pt.x]), np.array([pt.y]), np.array([0]), freq
        )
        eps_pt = self.sim.get_epsilon_point(pt, freq)
        print(f"eps:, ({pt.x},{pt.y}), {eps_grid}, {eps_pt}")
        self.assertAlmostEqual(np.real(eps_grid), np.real(eps_pt), places=6)
        self.assertAlmostEqual(np.imag(eps_grid), np.imag(eps_pt), places=6)


if __name__ == "__main__":
    unittest.main()
