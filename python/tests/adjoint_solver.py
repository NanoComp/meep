import meep as mp
try:
    import meep.adjoint as mpa
except:
    import adjoint as mpa
import numpy as np
from autograd import numpy as npa
from autograd import tensor_jacobian_product
import unittest
from enum import Enum

MonitorObject = Enum('MonitorObject', 'EIGENMODE DFT')

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

## ensure reproducible results
np.random.seed(9861548)

## random design region
p = np.random.rand(Nx*Ny)

## random epsilon perturbation for design region
deps = 1e-5
dp = deps*np.random.rand(Nx*Ny)

w = 1.0
waveguide_geometry = [mp.Block(material=silicon,
                               center=mp.Vector3(),
                               size=mp.Vector3(mp.inf,w,mp.inf))]

fcen = 1/1.55
df = 0.2*fcen
sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen,fwidth=df),
                              center=mp.Vector3(-0.5*sxy+dpml,0),
                              size=mp.Vector3(0,sxy),
                              eig_band=1,
                              eig_parity=eig_parity)]


def forward_simulation(design_params,mon_type):
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

    if mon_type.name == 'eigenmode':
        mode = sim.add_mode_monitor(fcen,
                                    0,
                                    1,
                                    mp.ModeRegion(center=mp.Vector3(0.5*sxy-dpml),size=mp.Vector3(0,sxy,0)),
                                    yee_grid=True)

    elif mon_type.name == 'DFT_fields':
        mode = sim.add_dft_fields([mp.Ez],
                                  fcen,
                                  0,
                                  1,
                                  center=mp.Vector3(0.5*sxy-dpml),
                                  size=mp.Vector3(0,sxy),
                                  yee_grid=False)

    sim.run(until_after_sources=20)

    if mon_type.name == 'eigenmode':
        coeff = sim.get_eigenmode_coefficients(mode,[1],eig_parity).alpha[0,0,0]
        S12 = abs(coeff)**2

    elif mon_type.name == 'DFT_fields':
        Ez_dft = sim.get_dft_array(mode, mp.Ez, 0)
        Ez2 = abs(Ez_dft[47])**2

    sim.reset_meep()

    if mon_type.name == 'eigenmode':
        return S12
    elif mon_type.name == 'DFT_fields':
        return Ez2


def adjoint_solver(design_params, mon_type):
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

    if mon_type.name == 'eigenmode':
        obj_list = [mpa.EigenmodeCoefficient(sim,
                                             mp.Volume(center=mp.Vector3(0.5*sxy-dpml),
                                                       size=mp.Vector3(0,sxy,0)),1)]

        def J(mode_mon):
            return npa.abs(mode_mon)**2

    elif mon_type.name == 'DFT_fields':
        obj_list = [mpa.FourierFields(sim,
                                      mp.Volume(center=mp.Vector3(0.5*sxy-dpml),
                                                size=mp.Vector3(0,sxy,0)),
                                      mp.Ez)]

        def J(mode_mon):
            return npa.abs(mode_mon[0,47])**2

    opt = mpa.OptimizationProblem(
        simulation = sim,
        objective_functions = J,
        objective_arguments = obj_list,
        design_regions = [matgrid_region],
        frequencies=[fcen],
        decay_fields=[mp.Ez])

    f, dJ_du = opt([design_params])

    sim.reset_meep()

    return f, dJ_du


def mapping(x,filter_radius,eta,beta):
    filtered_field = mpa.conic_filter(x,filter_radius,design_shape.x,design_shape.y,design_region_resolution)

    projected_field = mpa.tanh_projection(filtered_field,beta,eta)

    return projected_field.flatten()


class TestAdjointSolver(unittest.TestCase):

    def test_adjoint_solver_DFT_fields(self):
        ## compute gradient using adjoint solver
        adjsol_obj, adjsol_grad = adjoint_solver(p, MonitorObject.DFT)

        ## compute unperturbed Ez2
        Ez2_unperturbed = forward_simulation(p, MonitorObject.DFT)

        print("Ez2:, {:.6f}, {:.6f}".format(adjsol_obj,Ez2_unperturbed))
        self.assertAlmostEqual(adjsol_obj,Ez2_unperturbed,places=3)

        ## compute perturbed Ez2
        Ez2_perturbed = forward_simulation(p+dp, MonitorObject.DFT)

        print("directional_derivative:, {:.6f}, {:.6f}".format(np.dot(dp,adjsol_grad),Ez2_perturbed-Ez2_unperturbed))
        self.assertAlmostEqual(np.dot(dp,adjsol_grad),Ez2_perturbed-Ez2_unperturbed,places=5)


    def test_adjoint_solver_eigenmode(self):
        ## compute gradient using adjoint solver
        adjsol_obj, adjsol_grad = adjoint_solver(p, MonitorObject.EIGENMODE)

        ## compute unperturbed S12
        S12_unperturbed = forward_simulation(p, MonitorObject.EIGENMODE)

        print("S12:, {:.6f}, {:.6f}".format(adjsol_obj,S12_unperturbed))
        self.assertAlmostEqual(adjsol_obj,S12_unperturbed,places=3)

        ## compute perturbed S12
        S12_perturbed = forward_simulation(p+dp, MonitorObject.EIGENMODE)

        print("directional_derivative:, {:.6f}, {:.6f}".format(np.dot(dp,adjsol_grad),S12_perturbed-S12_unperturbed))
        self.assertAlmostEqual(np.dot(dp,adjsol_grad),S12_perturbed-S12_unperturbed,places=5)


    def test_gradient_backpropagation(self):
        ## filter/thresholding parameters
        filter_radius = 0.21985
        eta = 0.49093
        beta = 4.0698

        mapped_p = mapping(p,filter_radius,eta,beta)

        ## compute gradient using adjoint solver
        adjsol_obj, adjsol_grad = adjoint_solver(mapped_p, MonitorObject.EIGENMODE)

        ## backpropagate the gradient
        bp_adjsol_grad = tensor_jacobian_product(mapping,0)(p,filter_radius,eta,beta,adjsol_grad)

        ## compute unperturbed S12
        S12_unperturbed = forward_simulation(mapped_p, MonitorObject.EIGENMODE)

        print("S12:, {:.6f}, {:.6f}".format(adjsol_obj,S12_unperturbed))
        self.assertAlmostEqual(adjsol_obj,S12_unperturbed,places=3)

        ## compute perturbed S12
        S12_perturbed = forward_simulation(mapping(p+dp,filter_radius,eta,beta), MonitorObject.EIGENMODE)

        print("directional_derivative:, {:.6f}, {:.6f}".format(np.dot(dp,bp_adjsol_grad),S12_perturbed-S12_unperturbed))
        self.assertAlmostEqual(np.dot(dp,bp_adjsol_grad),S12_perturbed-S12_unperturbed,places=5)


if __name__ == '__main__':
    unittest.main()
