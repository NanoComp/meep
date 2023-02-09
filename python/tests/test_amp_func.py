import unittest
import meep as mp
import meep.adjoint as mpa
import numpy as np
from autograd import numpy as npa
import multiprocessing

def amp_func():
        Si = mp.Medium(index=3.4)
        SiO2 = mp.Medium(index=1.44)
        resolution = 20
        Sx = 3
        Sy = 3
        Sz = 3
        d3 = False
        si_thick = 1
        cell_size = mp.Vector3(Sx, Sy, Sz if d3 else 0)

        pml_layers = [mp.PML(0.5)]

        fcen = 1 / 1.55
        width = 0.2
        fwidth = width * fcen
        source_center = [-.9, 0, 0]
        source_size = mp.Vector3(0, 2, 1 if d3 else 0)
        src = mp.GaussianSource(frequency=fcen, fwidth=fwidth)

        def amplitude(xyz: mp.Vector3()):
            return np.exp(-15 * xyz.y**2)

        source = [
            mp.Source(
                src,
                component=mp.Ez,
                center = source_center,
                size = source_size,
                amp_func=amplitude
            )
        ]

        design_region_resolution = 10
        Nx = design_region_resolution
        Ny = design_region_resolution

        design_variables = mp.MaterialGrid(mp.Vector3(Nx, Ny), SiO2, Si, grid_type="U_MEAN")
        design_region = mpa.DesignRegion(
            design_variables, volume=mp.Volume(center=mp.Vector3(), size=mp.Vector3(1, 1, si_thick if d3 else 0))
        )

        geometry = [
            mp.Block(
                center=mp.Vector3(x=-Sx / 4), material=Si, size=mp.Vector3(Sx / 2, 0.5, si_thick if d3 else 0)
            ),  # horizontal waveguide
            mp.Block(
                center=mp.Vector3(y=Sy / 4), material=Si, size=mp.Vector3(0.5, Sy / 2, si_thick if d3 else 0)
            ),  # vertical waveguide
            mp.Block(
                center=design_region.center, size=design_region.size, material=design_variables
            ),  # design region
        ]

        sim = mp.Simulation(
            cell_size=cell_size,
            boundary_layers=pml_layers,
            geometry=geometry,
            sources=source,
            eps_averaging=False,
            resolution=resolution,
            force_all_components=True
        )

        TE0 = mpa.EigenmodeCoefficient(
            sim, mp.Volume(center=mp.Vector3(-.75, 0, 0), size=mp.Vector3(y=2, z=2 if d3 else 0)), mode=1
        )

        TE_top = mpa.EigenmodeCoefficient(
            sim, mp.Volume(center=mp.Vector3(0, .75, 0), size=mp.Vector3(x=2, y=0, z = 2 if d3 else 0)), mode=1
        )

        ob_list = [TE0, TE_top]

        def J(source, top):
            return npa.abs(top / source) ** 2

        opt = mpa.OptimizationProblem(
            simulation=sim,
            objective_functions=J,
            objective_arguments=ob_list,
            design_regions=[design_region],
            fcen=fcen,
            df=0,
            nf=1,
        )

        x0 = 0.5 * np.ones((Nx * Ny,))

        opt([x0])
        opt.forward_run() 

class TestAmpFunc(unittest.TestCase):   
    def test_amp_func(self):
        process = multiprocessing.Process(target=amp_func)
        process.start()
        process.join()

        self.assertEqual(process.exitcode, 0)

if __name__ == "__main__":
    unittest.main()
