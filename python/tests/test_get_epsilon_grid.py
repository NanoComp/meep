import unittest
import numpy as np
import meep as mp
from meep.materials import SiN, Co
from utils import ApproxComparisonTestCase

class TestGetEpsilonGrid(ApproxComparisonTestCase):

    def setUp(self):
        resolution = 100
        self.cell_size = mp.Vector3(1.0,1.0,0)

        matgrid_resolution = 200
        matgrid_size = mp.Vector3(1.0,1.0,mp.inf)
        Nx = int(matgrid_resolution*matgrid_size.x) + 1
        Ny = int(matgrid_resolution*matgrid_size.y) + 1
        x = np.linspace(-0.5*matgrid_size.x,0.5*matgrid_size.x,Nx)
        y = np.linspace(-0.5*matgrid_size.y,0.5*matgrid_size.y,Ny)
        xv, yv = np.meshgrid(x,y)
        rad = 0.201943
        w = 0.104283
        weights = np.logical_and(np.sqrt(np.square(xv) + np.square(yv)) > rad,
                                 np.sqrt(np.square(xv) + np.square(yv)) < rad+w,
                                 dtype=np.double)

        matgrid = mp.MaterialGrid(mp.Vector3(Nx,Ny),
                                  mp.air,
                                  mp.Medium(index=3.5),
                                  weights=weights,
                                  do_averaging=False,
                                  beta=0,
                                  eta=0.5)

        geometry = [mp.Cylinder(center=mp.Vector3(0.35,0.1),
                                radius=0.1,
                                height=mp.inf,
                                material=mp.Medium(index=1.5)),
                    mp.Block(center=mp.Vector3(-0.15,-0.2),
                             size=mp.Vector3(0.2,0.24,mp.inf),
                             material=SiN),
                    mp.Block(center=mp.Vector3(-0.2,0.2),
                             size=mp.Vector3(0.4,0.4,mp.inf),
                             material=matgrid),
                    mp.Prism(vertices=[mp.Vector3(0.05,0.45),
                                       mp.Vector3(0.32,0.22),
                                       mp.Vector3(0.15,0.10)],
                             height=0.5,
                             material=Co)]

        self.sim = mp.Simulation(resolution=resolution,
                                 cell_size=self.cell_size,
                                 geometry=geometry,
                                 eps_averaging=False)

    def test_get_epsilon_grid(self):
        grid_resolution = 3 if mp.is_single_precision() else 10
        Nx = int(grid_resolution*self.cell_size.x) + 1
        Ny = int(grid_resolution*self.cell_size.y) + 1
        xtics = np.linspace(-0.5*self.cell_size.x,0.5*self.cell_size.x,Nx)
        ytics = np.linspace(-0.5*self.cell_size.y,0.5*self.cell_size.y,Ny)

        frequency = [0, 0.3, 1.1, 1.5]
        self.sim.init_sim()

        for frq in frequency:
            eps_grid = self.sim.get_epsilon_grid(xtics,ytics,np.array([0]),frq)
            eps_pt = []
            for xi, xv in enumerate(xtics):
                for yi, yv in enumerate(ytics):
                    eps_pt.append(self.sim.get_epsilon_point(mp.Vector3(xv,yv),frq))
                    print("eps:, {}, ({:.4f}, {:.4f}), {}, {}".format(frq,xv,yv,eps_pt[-1],eps_grid[xi,yi]))

            self.assertClose(eps_grid, eps_pt, epsilon=1e-6)

if __name__ == '__main__':
    unittest.main()
