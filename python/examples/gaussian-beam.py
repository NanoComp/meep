## generate a titled Gaussian beam profile by defining the amplitude function of the source

import meep as mp
import math
import cmath
import numpy as np
import matplotlib.pyplot as plt

resolution = 40 # pixels/μm

cell_size = mp.Vector3(20,10,0)

pml_layers = [mp.PML(thickness=1.0,direction=mp.Y)]

fcen = 1.0 # center frequency of CW source (wavelength is 1 μm)

tilt_angle = math.radians(-10) # angle of tilted beam
k = mp.Vector3(y=1).rotate(mp.Vector3(z=1),tilt_angle).scale(fcen)

sigma = 1.5 # beam width

def gaussian_beam(sigma, k, x0):
    def _gaussian_beam(x):
        return cmath.exp(1j*2*math.pi*k.dot(x-x0)-(x-x0).dot(x-x0)/(2*sigma**2))
    return _gaussian_beam

src_pt = mp.Vector3(y=4)
sources = [mp.Source(src=mp.ContinuousSource(fcen, fwidth=0.2*fcen),
                     component=mp.Ez,
                     center=src_pt,
                     size=mp.Vector3(20),
                     amp_func=gaussian_beam(sigma,k,src_pt))]

sim = mp.Simulation(cell_size=cell_size,
                    sources=sources,
                    k_point=k,
                    boundary_layers=pml_layers,
                    resolution=resolution)

non_pml_vol = mp.Volume(center=mp.Vector3(), size=mp.Vector3(20,8,0))
sim.run(until=50)

ez_data = sim.get_array(vol=non_pml_vol, component=mp.Ez)

plt.figure()
plt.imshow(np.flipud(np.transpose(np.real(ez_data))), interpolation='spline36', cmap='RdBu')
plt.axis('off')
plt.show()
