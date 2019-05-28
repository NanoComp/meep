import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt

## Material parameters
sn = 0.0002
fn = 0.5
gn = 1e-5
epsn = 15
psat = 100

def meep_sim(gn, sn, epsn, tmax):

    ## For mu ~ 1 and epsn ~ 15, wavelength will be ~0.26.  Set dx and dt to ~0.01
    resolution = 100

    cellsize = [0.01, 0.01, 5.0] # Don't set to zero size: need 3D fields
    fsrc = 0.8
    src_z = -1.8
    pml_depth = 0.25

    import meep as mp

    cell = mp.Vector3(*cellsize)

    ## Define gyromagnetic material biased along z direction
    susc = [mp.GyrotropicSusceptibility(frequency=fn, gamma=gn, sigma=sn,
                                        bias=mp.Vector3(0, 0, psat))]
    mat = mp.Medium(epsilon=epsn, mu=1, H_susceptibilities=susc)

    pml_layers = [mp.PML(thickness=pml_depth, direction=mp.Z)]

    sources = [mp.Source(mp.ContinuousSource(frequency=fsrc),
                         component=mp.Ex, center=mp.Vector3(0, 0, src_z))]

    sim = mp.Simulation(cell_size=cell, geometry=[], sources=sources,
                        k_point=mp.Vector3(), # Periodic BCs in x and y
                        boundary_layers=pml_layers,
                        default_material=mat, resolution=resolution)

    sim.run(until=tmax)

    ex_data = sim.get_array(center=mp.Vector3(), size=mp.Vector3(0, 0, cellsize[2]), component=mp.Ex)
    ey_data = sim.get_array(center=mp.Vector3(), size=mp.Vector3(0, 0, cellsize[2]), component=mp.Ey)

    L = cellsize[2]/2
    z = np.linspace(-L, L, len(ex_data))
    plt.title('t = {}'.format(tmax))
    plt.plot(z, ex_data, label='Ex')
    plt.plot(z, ey_data, label='Ey')
    plt.xlim(-L, L)
    plt.xlabel('z')
    plt.legend()

plt.figure(2, figsize=(8,8))
plt.subplot(2,1,1)
meep_sim(gn, sn, epsn, 200)
plt.subplot(2,1,2)
meep_sim(gn, sn, epsn, 500)
plt.tight_layout()
