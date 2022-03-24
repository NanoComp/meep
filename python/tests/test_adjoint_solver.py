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
from utils import ApproxComparisonTestCase

MonitorObject = Enum('MonitorObject', 'EIGENMODE DFT')

resolution = 30

silicon = mp.Medium(epsilon=12)
sapphire = mp.Medium(epsilon_diag=(10.225,10.225,9.95),
                     epsilon_offdiag=(-0.825,-0.55*np.sqrt(3/2),0.55*np.sqrt(3/2)))

sxy = 5.0
cell_size = mp.Vector3(sxy,sxy,0)

dpml = 1.0
pml_xy = [mp.PML(thickness=dpml)]
pml_x = [mp.PML(thickness=dpml,direction=mp.X)]

eig_parity = mp.EVEN_Y + mp.ODD_Z

design_region_size = mp.Vector3(1.5,1.5)
design_region_resolution = int(2*resolution)
Nx, Ny = int(design_region_size.x*design_region_resolution), int(design_region_size.y*design_region_resolution)

## ensure reproducible results
rng = np.random.RandomState(9861548)

## random design region
p = 0.5*rng.rand(Nx*Ny)

## random epsilon perturbation for design region
deps = 1e-5
dp = deps*rng.rand(Nx*Ny)

w = 1.0
waveguide_geometry = [mp.Block(material=silicon,
                               center=mp.Vector3(),
                               size=mp.Vector3(mp.inf,w,mp.inf))]

fcen = 1/1.55
df = 0.23*fcen
wvg_source = [mp.EigenModeSource(src=mp.GaussianSource(fcen,fwidth=df),
                                 center=mp.Vector3(-0.5*sxy+dpml+0.1,0),
                                 size=mp.Vector3(0,sxy-2*dpml),
                                 eig_parity=eig_parity)]

pt_source = [mp.Source(src=mp.GaussianSource(fcen,fwidth=df),
                       center=mp.Vector3(-0.5*sxy+dpml,0),
                       size=mp.Vector3(),
                       component=mp.Ez)]

k_point = mp.Vector3(0.23,-0.38)


def forward_simulation(design_params, mon_type, frequencies=None, mat2=silicon):
    matgrid = mp.MaterialGrid(mp.Vector3(Nx,Ny),
                              mp.air,
                              mat2,
                              weights=design_params.reshape(Nx,Ny))

    matgrid_geometry = [mp.Block(center=mp.Vector3(),
                                 size=mp.Vector3(design_region_size.x,design_region_size.y,0),
                                 material=matgrid)]

    geometry = waveguide_geometry + matgrid_geometry

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=pml_xy,
                        sources=wvg_source,
                        geometry=geometry)

    if not frequencies:
        frequencies = [fcen]

    if mon_type.name == 'EIGENMODE':
        mode = sim.add_mode_monitor(frequencies,
                                    mp.ModeRegion(center=mp.Vector3(0.5*sxy-dpml-0.1),
                                                  size=mp.Vector3(0,sxy-2*dpml,0)),
                                    yee_grid=True,
                                    eig_parity=eig_parity)
    elif mon_type.name == 'DFT':
        mode = sim.add_dft_fields([mp.Ez],
                                  frequencies,
                                  center=mp.Vector3(1.25),
                                  size=mp.Vector3(0.25,1,0),
                                  yee_grid=False)

    sim.run(until_after_sources=mp.stop_when_dft_decayed())

    if mon_type.name == 'EIGENMODE':
        coeff = sim.get_eigenmode_coefficients(mode,[1],eig_parity).alpha[0,:,0]
        S12 = np.power(np.abs(coeff),2)
    elif mon_type.name == 'DFT':
        Ez2 = []
        for f in range(len(frequencies)):
            Ez_dft = sim.get_dft_array(mode, mp.Ez, f)
            Ez2.append(np.power(np.abs(Ez_dft[4,10]),2))
        Ez2 = np.array(Ez2)

    sim.reset_meep()

    if mon_type.name == 'EIGENMODE':
        return S12
    elif mon_type.name == 'DFT':
        return Ez2


def adjoint_solver(design_params, mon_type, frequencies=None, mat2=silicon):
    matgrid = mp.MaterialGrid(mp.Vector3(Nx,Ny),
                              mp.air,
                              mat2,
                              weights=np.ones((Nx,Ny)))

    matgrid_region = mpa.DesignRegion(matgrid,
                                      volume=mp.Volume(center=mp.Vector3(),
                                                       size=mp.Vector3(design_region_size.x,
                                                                       design_region_size.y,
                                                                       0)))

    matgrid_geometry = [mp.Block(center=matgrid_region.center,
                                 size=matgrid_region.size,
                                 material=matgrid)]

    geometry = waveguide_geometry + matgrid_geometry

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=pml_xy,
                        sources=wvg_source,
                        geometry=geometry)

    if not frequencies:
        frequencies = [fcen]

    if mon_type.name == 'EIGENMODE':
        obj_list = [mpa.EigenmodeCoefficient(sim,
                                             mp.Volume(center=mp.Vector3(0.5*sxy-dpml-0.1),
                                                       size=mp.Vector3(0,sxy-2*dpml,0)),
                                             1,
                                             eig_parity=eig_parity)]

        def J(mode_mon):
            return npa.power(npa.abs(mode_mon),2)
    elif mon_type.name == 'DFT':
        obj_list = [mpa.FourierFields(sim,
                                      mp.Volume(center=mp.Vector3(1.25),
                                                size=mp.Vector3(0.25,1,0)),
                                      mp.Ez)]

        def J(mode_mon):
            return npa.power(npa.abs(mode_mon[:,4,10]),2)

    opt = mpa.OptimizationProblem(
        simulation=sim,
        objective_functions=J,
        objective_arguments=obj_list,
        design_regions=[matgrid_region],
        frequencies=frequencies)

    f, dJ_du = opt([design_params])

    sim.reset_meep()

    return f, dJ_du


def forward_simulation_complex_fields(design_params, frequencies=None):
    matgrid = mp.MaterialGrid(mp.Vector3(Nx,Ny),
                              mp.air,
                              silicon,
                              weights=design_params.reshape(Nx,Ny))

    geometry = [mp.Block(center=mp.Vector3(),
                         size=mp.Vector3(design_region_size.x,
                                         design_region_size.y,
                                         0),
                         material=matgrid)]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        k_point=k_point,
                        boundary_layers=pml_x,
                        sources=pt_source,
                        geometry=geometry)

    if not frequencies:
        frequencies = [fcen]

    mode = sim.add_dft_fields([mp.Ez],
                              frequencies,
                              center=mp.Vector3(0.9),
                              size=mp.Vector3(0.2,0.5),
                              yee_grid=False)

    sim.run(until_after_sources=mp.stop_when_dft_decayed())

    Ez2 = []
    for f in range(len(frequencies)):
        Ez_dft = sim.get_dft_array(mode, mp.Ez, f)
        Ez2.append(np.power(np.abs(Ez_dft[3,9]),2))
    Ez2 = np.array(Ez2)

    sim.reset_meep()

    return Ez2


def adjoint_solver_complex_fields(design_params, frequencies=None):
    matgrid = mp.MaterialGrid(mp.Vector3(Nx,Ny),
                              mp.air,
                              silicon,
                              weights=np.ones((Nx,Ny)))

    matgrid_region = mpa.DesignRegion(matgrid,
                                      volume=mp.Volume(center=mp.Vector3(),
                                                       size=mp.Vector3(design_region_size.x,
                                                                       design_region_size.y,
                                                                       0)))

    geometry = [mp.Block(center=matgrid_region.center,
                         size=matgrid_region.size,
                         material=matgrid)]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        k_point=k_point,
                        boundary_layers=pml_x,
                        sources=pt_source,
                        geometry=geometry)

    if not frequencies:
        frequencies = [fcen]

    obj_list = [mpa.FourierFields(sim,
                                  mp.Volume(center=mp.Vector3(0.9),
                                            size=mp.Vector3(0.2,0.5)),
                                  mp.Ez)]

    def J(dft_mon):
        return npa.power(npa.abs(dft_mon[:,3,9]),2)

    opt = mpa.OptimizationProblem(
        simulation=sim,
        objective_functions=J,
        objective_arguments=obj_list,
        design_regions=[matgrid_region],
        frequencies=frequencies)

    f, dJ_du = opt([design_params])

    sim.reset_meep()

    return f, dJ_du
    
def forward_simulation_damping(design_params, frequencies=None, mat2=silicon):
    matgrid = mp.MaterialGrid(mp.Vector3(Nx,Ny),
                              mp.air,
                              mat2,
                              weights=design_params.reshape(Nx,Ny),
                              damping = 3.14*fcen)

    matgrid_geometry = [mp.Block(center=mp.Vector3(),
                                 size=mp.Vector3(design_region_size.x,design_region_size.y,0),
                                 material=matgrid)]

    geometry = waveguide_geometry + matgrid_geometry

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=pml_xy,
                        sources=wvg_source,
                        geometry=geometry)

    if not frequencies:
        frequencies = [fcen]

    mode = sim.add_mode_monitor(frequencies,
                                    mp.ModeRegion(center=mp.Vector3(0.5*sxy-dpml-0.1),
                                                  size=mp.Vector3(0,sxy-2*dpml,0)),
                                    yee_grid=True,
                                    eig_parity=eig_parity)

    sim.run(until_after_sources=mp.stop_when_dft_decayed())


    coeff = sim.get_eigenmode_coefficients(mode,[1],eig_parity).alpha[0,:,0]
    S12 = np.power(np.abs(coeff),2)
    sim.reset_meep()
    return S12

def adjoint_solver_damping(design_params, frequencies=None, mat2=silicon):
    matgrid = mp.MaterialGrid(mp.Vector3(Nx,Ny),
                              mp.air,
                              mat2,
                              weights=np.ones((Nx,Ny)),
                              damping = 3.14*fcen)
    matgrid_region = mpa.DesignRegion(matgrid,
                                      volume=mp.Volume(center=mp.Vector3(), size=mp.Vector3(design_region_size.x, design_region_size.y, 0)))

    matgrid_geometry = [mp.Block(center=matgrid_region.center,
                                 size=matgrid_region.size,
                                 material=matgrid)]

    geometry = waveguide_geometry + matgrid_geometry

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=pml_xy,
                        sources=wvg_source,
                        geometry=geometry)

    if not frequencies:
        frequencies = [fcen]

    obj_list = [mpa.EigenmodeCoefficient(sim, mp.Volume(center=mp.Vector3(0.5*sxy-dpml-0.1),
                                        size=mp.Vector3(0,sxy-2*dpml,0)), 1, eig_parity=eig_parity)]

    def J(mode_mon):
        return npa.power(npa.abs(mode_mon),2)


    opt = mpa.OptimizationProblem(
        simulation=sim,
        objective_functions=J,
        objective_arguments=obj_list,
        design_regions=[matgrid_region],
        frequencies=frequencies,
        minimum_run_time=150)

    f, dJ_du = opt([design_params])

    sim.reset_meep()
    return f, dJ_du

def mapping(x,filter_radius,eta,beta):
    filtered_field = mpa.conic_filter(x,
                                      filter_radius,
                                      design_region_size.x,
                                      design_region_size.y,
                                      design_region_resolution)

    projected_field = mpa.tanh_projection(filtered_field,beta,eta)

    return projected_field.flatten()


class TestAdjointSolver(ApproxComparisonTestCase):

    def test_adjoint_solver_DFT_fields(self):
        print("*** TESTING DFT ADJOINT ***")

        ## test the single frequency and multi frequency cases
        for frequencies in [[fcen], [1/1.58, fcen, 1/1.53]]:
            ## compute gradient using adjoint solver
            adjsol_obj, adjsol_grad = adjoint_solver(p, MonitorObject.DFT, frequencies)

            ## compute unperturbed |Ez|^2
            Ez2_unperturbed = forward_simulation(p, MonitorObject.DFT, frequencies)

            ## compare objective results
            print("|Ez|^2 -- adjoint solver: {}, traditional simulation: {}".format(adjsol_obj,Ez2_unperturbed))
            self.assertClose(adjsol_obj,Ez2_unperturbed,epsilon=1e-6)

            ## compute perturbed Ez2
            Ez2_perturbed = forward_simulation(p+dp, MonitorObject.DFT, frequencies)

            ## compare gradients
            if adjsol_grad.ndim < 2:
                adjsol_grad = np.expand_dims(adjsol_grad,axis=1)
            adj_scale = (dp[None,:]@adjsol_grad).flatten()
            fd_grad = Ez2_perturbed-Ez2_unperturbed
            print("Directional derivative -- adjoint solver: {}, FD: {}".format(adj_scale,fd_grad))
            tol = 0.07 if mp.is_single_precision() else 0.006
            self.assertClose(adj_scale,fd_grad,epsilon=tol)


    def test_adjoint_solver_eigenmode(self):
        print("*** TESTING EIGENMODE ADJOINT ***")

        ## test the single frequency and multi frequency case
        for frequencies in [[fcen], [1/1.58, fcen, 1/1.53]]:
            ## compute gradient using adjoint solver
            adjsol_obj, adjsol_grad = adjoint_solver(p, MonitorObject.EIGENMODE, frequencies)

            ## compute unperturbed S12
            S12_unperturbed = forward_simulation(p, MonitorObject.EIGENMODE, frequencies)

            ## compare objective results
            print("S12 -- adjoint solver: {}, traditional simulation: {}".format(adjsol_obj,S12_unperturbed))
            self.assertClose(adjsol_obj,S12_unperturbed,epsilon=1e-6)

            ## compute perturbed S12
            S12_perturbed = forward_simulation(p+dp, MonitorObject.EIGENMODE, frequencies)

            ## compare gradients
            if adjsol_grad.ndim < 2:
                adjsol_grad = np.expand_dims(adjsol_grad,axis=1)
            adj_scale = (dp[None,:]@adjsol_grad).flatten()
            fd_grad = S12_perturbed-S12_unperturbed
            print("Directional derivative -- adjoint solver: {}, FD: {}".format(adj_scale,fd_grad))
            tol = 0.04 if mp.is_single_precision() else 0.01
            self.assertClose(adj_scale,fd_grad,epsilon=tol)


    def test_gradient_backpropagation(self):
        print("*** TESTING BACKPROP ***")

        for frequencies in [[fcen], [1/1.58, fcen, 1/1.53]]:
            ## filter/thresholding parameters
            filter_radius = 0.21985
            eta = 0.49093
            beta = 4.0698

            mapped_p = mapping(p,filter_radius,eta,beta)

            ## compute gradient using adjoint solver
            adjsol_obj, adjsol_grad = adjoint_solver(mapped_p, MonitorObject.EIGENMODE, frequencies)

            ## backpropagate the gradient
            if len(frequencies) > 1:
                bp_adjsol_grad = np.zeros(adjsol_grad.shape)
                for i in range(len(frequencies)):
                    bp_adjsol_grad[:,i] = tensor_jacobian_product(mapping,0)(p,filter_radius,eta,beta,adjsol_grad[:,i])
            else:
                bp_adjsol_grad = tensor_jacobian_product(mapping,0)(p,filter_radius,eta,beta,adjsol_grad)

            ## compute unperturbed S12
            S12_unperturbed = forward_simulation(mapped_p,MonitorObject.EIGENMODE,frequencies)

            ## compare objective results
            print("S12 -- adjoint solver: {}, traditional simulation: {}".format(adjsol_obj,S12_unperturbed))
            self.assertClose(adjsol_obj,S12_unperturbed,epsilon=1e-6)

            ## compute perturbed S12
            S12_perturbed = forward_simulation(mapping(p+dp,filter_radius,eta,beta),MonitorObject.EIGENMODE,frequencies)

            if bp_adjsol_grad.ndim < 2:
                bp_adjsol_grad = np.expand_dims(bp_adjsol_grad,axis=1)
            adj_scale = (dp[None,:]@bp_adjsol_grad).flatten()
            fd_grad = S12_perturbed-S12_unperturbed
            print("Directional derivative -- adjoint solver: {}, FD: {}".format(adj_scale,fd_grad))
            tol = 0.02 if mp.is_single_precision() else 0.01
            self.assertClose(adj_scale,fd_grad,epsilon=tol)


    def test_complex_fields(self):
        print("*** TESTING COMPLEX FIELDS ***")

        for frequencies in [[fcen], [1/1.58, fcen, 1/1.53]]:
            ## compute gradient using adjoint solver
            adjsol_obj, adjsol_grad = adjoint_solver_complex_fields(p, frequencies)

            ## compute unperturbed |Ez|^2
            Ez2_unperturbed = forward_simulation_complex_fields(p, frequencies)

            ## compare objective results
            print("Ez2 -- adjoint solver: {}, traditional simulation: {}".format(adjsol_obj,Ez2_unperturbed))
            self.assertClose(adjsol_obj,Ez2_unperturbed,epsilon=1e-6)

            ## compute perturbed |Ez|^2
            Ez2_perturbed = forward_simulation_complex_fields(p+dp, frequencies)

            ## compare gradients
            if adjsol_grad.ndim < 2:
                adjsol_grad = np.expand_dims(adjsol_grad,axis=1)
            adj_scale = (dp[None,:]@adjsol_grad).flatten()
            fd_grad = Ez2_perturbed-Ez2_unperturbed
            print("Directional derivative -- adjoint solver: {}, FD: {}".format(adj_scale,fd_grad))
            tol = 0.018 if mp.is_single_precision() else 0.002
            self.assertClose(adj_scale,fd_grad,epsilon=tol)

    def test_damping(self):
        print("*** TESTING CONDUCTIVITIES ***")

        for frequencies in [[1/1.58, fcen, 1/1.53]]:
            ## compute gradient using adjoint solver
            adjsol_obj, adjsol_grad = adjoint_solver_damping(p, frequencies)

            ## compute unperturbed S12
            S12_unperturbed = forward_simulation_damping(p, frequencies)

            ## compare objective results
            print("S12 -- adjoint solver: {}, traditional simulation: {}".format(adjsol_obj,S12_unperturbed))
            self.assertClose(adjsol_obj,S12_unperturbed,epsilon=1e-6)

            ## compute perturbed S12
            S12_perturbed = forward_simulation_damping(p+dp, frequencies)

            ## compare gradients
            if adjsol_grad.ndim < 2:
                adjsol_grad = np.expand_dims(adjsol_grad,axis=1)
            adj_scale = (dp[None,:]@adjsol_grad).flatten()
            fd_grad = S12_perturbed-S12_unperturbed
            print("Directional derivative -- adjoint solver: {}, FD: {}".format(adj_scale,fd_grad))
            tol = 0.06 if mp.is_single_precision() else 0.03
            self.assertClose(adj_scale,fd_grad,epsilon=tol)

    def test_offdiagonal(self):
        print("*** TESTING OFFDIAGONAL COMPONENTS ***")

        ## test the single frequency and multi frequency case
        for frequencies in [[fcen], [1/1.58, fcen, 1/1.53]]:
            ## compute gradient using adjoint solver
            adjsol_obj, adjsol_grad = adjoint_solver(p, MonitorObject.EIGENMODE, frequencies, sapphire)

            ## compute unperturbed S12
            S12_unperturbed = forward_simulation(p, MonitorObject.EIGENMODE, frequencies, sapphire)

            ## compare objective results
            print("S12 -- adjoint solver: {}, traditional simulation: {}".format(adjsol_obj,S12_unperturbed))
            self.assertClose(adjsol_obj,S12_unperturbed,epsilon=1e-6)

            ## compute perturbed S12
            S12_perturbed = forward_simulation(p+dp, MonitorObject.EIGENMODE, frequencies, sapphire)

            ## compare gradients
            if adjsol_grad.ndim < 2:
                adjsol_grad = np.expand_dims(adjsol_grad,axis=1)
            adj_scale = (dp[None,:]@adjsol_grad).flatten()
            fd_grad = S12_perturbed-S12_unperturbed
            print("Directional derivative -- adjoint solver: {}, FD: {}".format(adj_scale,fd_grad))
            tol = 0.05 if mp.is_single_precision() else 0.04
            self.assertClose(adj_scale,fd_grad,epsilon=tol)
if __name__ == '__main__':
    unittest.main()
