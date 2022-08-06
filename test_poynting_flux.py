import meep as mp
import numpy as np
import meep.adjoint as mpa
from autograd import numpy as npa
from autograd import tensor_jacobian_product
import nlopt
from matplotlib import pyplot as plt

resolution = 30  # pixels/μm

silicon = mp.Medium(epsilon=12)
sapphire = mp.Medium(epsilon_diag=(10.225,10.225,9.95),
                            epsilon_offdiag=(-0.825,-0.55*np.sqrt(3/2),0.55*np.sqrt(3/2)))

sxy = 5.0
cell_size = mp.Vector3(sxy,sxy,0)
#mp.integrate_field_function()

dpml = 1.0
pml_xy = [mp.PML(thickness=dpml)]
pml_x = [mp.PML(thickness=dpml,direction=mp.X)]

eig_parity = mp.EVEN_Y+mp.ODD_Z

design_region_size = mp.Vector3(1.5,1.5)
design_region_resolution = int(2*resolution)
Nx = int(design_region_size.x*design_region_resolution)
Ny = int(design_region_size.y*design_region_resolution)

# ensure reproducible results
rng = np.random.RandomState(9861548)

# random design region
#p = 0.5*rng.rand(Nx*Ny)
p = np.ones(Nx*Ny)

# random perturbation for design region
deps = 1e-5
dp = deps*rng.rand(Nx*Ny)

w = 1.0
waveguide_geometry = [mp.Block(material=silicon,
                                    center=mp.Vector3(),
                                    size=mp.Vector3(mp.inf,w,mp.inf))]

fcen = 1/1.55
df = 0.2*fcen
mode_source = [mp.EigenModeSource(src=mp.GaussianSource(fcen,fwidth=df),
                                        center=mp.Vector3(-0.5*sxy+dpml,0),
                                        size=mp.Vector3(0,sxy-2*dpml),
                                        eig_parity=eig_parity)]

pt_source = [mp.Source(src=mp.GaussianSource(fcen,fwidth=df),
                            center=mp.Vector3(-0.5*sxy+dpml,0),
                            size=mp.Vector3(),
                            component=mp.Ez)]

line_source = [mp.Source(src=mp.GaussianSource(fcen,fwidth=df),
                                center=mp.Vector3(-0.85,0),
                                size=mp.Vector3(0,sxy-2*dpml),
                                component=mp.Ez)]

k_point = mp.Vector3(0.23,-0.38)


matgrid = mp.MaterialGrid(mp.Vector3(Nx,Ny),
                                mp.air,
                                silicon,
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
                    sources=mode_source,
                    geometry=geometry)


frequencies = [fcen]
print(fcen)
volume = mp.Volume(center=mp.Vector3(1.25), size=mp.Vector3(0,2,0))
# interesting_components = [mp.Hz,mp.Hy,mp.Ez,mp.Ey]
# dft_fields = sim.add_dft_fields(interesting_components, frequencies, where = volume,yee_grid = False)
# sim.run(until_after_sources=mp.stop_when_fields_decayed(100,mp.Ez,mp.Vector3(1.25,0),1e-12))
# dft_field_data = []
# for comp in interesting_components:
#     dft_field_data.append(sim.get_dft_array(dft_fields,comp,0))
# print(dft_field_data)
# discretization_factor = 1/2*dft_field_data[1].size
# print(discretization_factor)
# manual_flux_from_dft_fields = -np.sum(np.conj(dft_field_data[2])*dft_field_data[1])/(discretization_factor)
# print("manual_flux_from_dft_fields:")
# print(manual_flux_from_dft_fields)
# manual_flux = np.conj(dft_field_data)

# refl_fr = mp.FluxRegion(center=mp.Vector3(1.25,0), size=mp.Vector3(0,2,0))
# refl = sim.add_flux(fcen, df, 1,refl_fr)
# sim.run(until_after_sources=mp.stop_when_fields_decayed(100,mp.Ez,mp.Vector3(1.25,0),1e-2))

# straight_tran_data = sim.get_flux_data(refl)
# manual_flux = np.conj(straight_tran_data.H)*(straight_tran_data.E)
# print("manual flux is:")
# print(np.sum(manual_flux))
# straight_tran_flux = mp.get_fluxes(refl)
# print(mp.get_flux_freqs(refl))
# print(straight_tran_data)
# print(straight_tran_flux)

obj_list = [mpa.PoyntingFlux(sim, volume)]
eval_hist = []
def J(mode_mon):
    eval = npa.sum(mode_mon)
    return eval
opt = mpa.OptimizationProblem(simulation=sim,
                                      objective_functions=J,
                                      objective_arguments=obj_list,
                                      design_regions=[matgrid_region],
                                      frequencies=frequencies)

f, dJ_du = opt([p],need_gradient=False)
print(f)

# algorithm = nlopt.LD_MMA
# n = Nx*Ny
# maxeval = 10

# #all_fouriersrcdata = monitor.swigobj.fourier_sourcedata(self.volume.swigobj, self.component, self.sim.fields, dJ)

# evaluation_history = []
# sensitivity = [0]
# def f_fun(x, grad):
#     f0, dJ_du = opt([x])
#     f0 = f0[0] # f0 is an array of length 1 
#     if grad.size > 0:
#         grad[:] = np.squeeze(dJ_du)
#     evaluation_history.append(np.real(f0))
#     sensitivity[0] = dJ_du
#     return np.real(f0)


# solver = nlopt.opt(algorithm, n)
# solver.set_lower_bounds(0)
# solver.set_upper_bounds(1)
# solver.set_max_objective(f_fun)
# solver.set_maxeval(maxeval)
# x = solver.optimize(p)
# plt.figure()
# plt.plot(evaluation_history,'o-')
# plt.grid(True)
# plt.xlabel('Iteration')
# plt.ylabel('FOM')
# plt.show()