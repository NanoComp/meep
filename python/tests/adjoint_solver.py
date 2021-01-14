import meep as mp
import meep.adjoint as mpa
import numpy as np
from autograd import numpy as npa
from copy import deepcopy
import unittest

class TestAdjointSolver(unittest.TestCase):

    def test_adjoint_solver(self):
        np.random.seed(9861548)
        mp.verbosity(1)

        resolution = 25

        silicon = mp.Medium(epsilon=12)

        sxy = 5.0
        cell_size = mp.Vector3(sxy,sxy,0)

        dpml = 1.0
        boundary_layers = [mp.PML(thickness=dpml)]

        eig_parity = mp.EVEN_Y+mp.ODD_Z

        design_shape = mp.Vector3(1.5,1.5)
        design_region_resolution = int(2*resolution)
        Nx = int(design_region_resolution*design_shape.x)
        Ny = int(design_region_resolution*design_shape.y)

        w = 1.0 # waveguide width
        waveguide_geometry = [mp.Block(material=silicon,
                                       center=mp.Vector3(),
                                       size=mp.Vector3(mp.inf,w,mp.inf))]

        fcen = 1/1.55
        df = 0.2*fcen
        sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen,fwidth=df),
                                      center=mp.Vector3(-0.5*sxy+dpml,0),
                                      size=mp.Vector3(0,sxy,0),
                                      eig_band=1,
                                      eig_parity=eig_parity,
                                      eig_match_freq=True)]

        def forward_simulation(design_params):
            matgrid = mp.MaterialGrid(mp.Vector3(Nx,Ny),
                                      mp.air,
                                      silicon,
                                      design_parameters=design_params.reshape(Nx,Ny),
                                      grid_type='U_SUM')
            
            matgrid_geometry = [mp.Block(center=mp.Vector3(),
                                         size=mp.Vector3(design_shape.x,design_shape.y,0),
                                         material=matgrid)]

            geometry = waveguide_geometry + matgrid_geometry

            sim = mp.Simulation(resolution=resolution,
                                cell_size=cell_size,
                                boundary_layers=boundary_layers,
                                sources=sources,
                                geometry=geometry)

            mode = sim.add_mode_monitor(fcen, 0, 1,
                                        mp.ModeRegion(center=mp.Vector3(0.5*sxy-dpml),size=mp.Vector3(0,sxy,0)),
                                        yee_grid=True)

            sim.run(until_after_sources=20)

            # mode coefficients
            coeff = sim.get_eigenmode_coefficients(mode,[1],eig_parity).alpha[0,0,0]

            # S parameters
            S12 = abs(coeff)**2

            sim.reset_meep()

            return S12

        def adjoint_solver(design_params):
            matgrid = mp.MaterialGrid(mp.Vector3(Nx,Ny),
                                      mp.air,
                                      silicon,
                                      design_parameters=np.ones((Nx,Ny)))

            matgrid_region = mpa.DesignRegion(matgrid,
                                              volume=mp.Volume(center=mp.Vector3(),
                                                               size=mp.Vector3(design_shape.x,design_shape.y,0)))

            matgrid_geometry = [mp.Block(center=matgrid_region.center,
                                         size=matgrid_region.size,
                                         material=matgrid)]

            geometry = waveguide_geometry + matgrid_geometry

            sim = mp.Simulation(resolution=resolution,
                                cell_size=cell_size,
                                boundary_layers=boundary_layers,
                                sources=sources,
                                geometry=geometry)

            obj_list = [mpa.EigenmodeCoefficient(sim,mp.Volume(center=mp.Vector3(0.5*sxy-dpml),size=mp.Vector3(0,sxy,0)),1)]

            def J(mode_mon):
                return npa.abs(mode_mon)**2

            opt = mpa.OptimizationProblem(
                simulation = sim,
                objective_functions = J,
                objective_arguments = obj_list,
                design_regions = [matgrid_region],
                frequencies=[fcen],
                decay_by = 1e-4,
                decay_fields=[mp.Ez])

            f, dJ_du = opt([design_params])

            sim.reset_meep()

            return dJ_du

        design_params_orig = np.random.rand(Nx*Ny)

        ## compute gradient using adjoint solver
        adjsol_grad = adjoint_solver(design_params_orig)

        ## compute gradient using finite differences
        S12_unperturbed = forward_simulation(design_params_orig)

        ## epsilon perturbation
        deps = 1e-7

        design_params = deepcopy(design_params_orig)
        design_params[2984] = design_params[2984] + deps
        S12_perturbed = forward_simulation(design_params)
        fd_grad1 = (S12_perturbed - S12_unperturbed) / deps

        print("gradient #1:, {:.10f}, {:.10f}".format(fd_grad1,adjsol_grad[2984]))
        self.assertAlmostEqual(fd_grad1,adjsol_grad[2984],places=3)

        design_params = deepcopy(design_params_orig)
        design_params[5044] = design_params[5044] + deps
        S12_perturbed = forward_simulation(design_params)
        fd_grad2 = (S12_perturbed - S12_unperturbed) / deps

        print("gradient #2:, {:.10f}, {:.10f}".format(fd_grad2,adjsol_grad[5044]))
        self.assertAlmostEqual(fd_grad2,adjsol_grad[5044],places=3)

if __name__ == '__main__':
    unittest.main()
