import meep as mp

## Define a gyroelectric medium
f0 = 1.0
gamma = 1e-6
epsn = 1.5
b0 = 0.15
sn = 0.1

susc = [mp.GyrotropicLorentzianSusceptibility(frequency=f0, gamma=gamma, sigma=sn,
                                              bias=mp.Vector3(0, 0, b0))]
mat = mp.Medium(epsilon=epsn, mu=1, E_susceptibilities=susc)

## Set up and run the Meep simulation:
tmax = 100
L = 20.0
cell = mp.Vector3(0, 0, L)
fsrc, src_z = 0.8, -8.5
pml_layers = [mp.PML(thickness=1.0, direction=mp.Z)]

sources = [mp.Source(mp.ContinuousSource(frequency=fsrc),
                     component=mp.Ex, center=mp.Vector3(0, 0, src_z))]

sim = mp.Simulation(cell_size=cell, geometry=[], sources=sources,
                    k_point=mp.Vector3(),   # Periodic boundary conditions
                    boundary_layers=pml_layers,
                    default_material=mat, resolution=100)
sim.run(until=tmax)

## Plot results:
import numpy as np
import matplotlib.pyplot as plt

ex_data = sim.get_array(center=mp.Vector3(), size=mp.Vector3(0, 0, L), component=mp.Ex)
ey_data = sim.get_array(center=mp.Vector3(), size=mp.Vector3(0, 0, L), component=mp.Ey)

z = np.linspace(-L/2, L/2, len(ex_data))
plt.figure(1)
plt.plot(z, ex_data, label='Ex')
plt.plot(z, ey_data, label='Ey')
plt.xlim(-L/2, L/2); plt.xlabel('z')
plt.legend()

## Comparison with analytic result:
dfsq = (f0**2 - 1j*fsrc*gamma - fsrc**2)
eperp = epsn + sn * f0**2 * dfsq / (dfsq**2 + (fsrc*b0)**2)
eta = sn * f0**2 * fsrc * b0 / (dfsq**2 + (fsrc*b0)**2)

k_gyro = 2*np.pi*fsrc * np.sqrt(0.5*(eperp - np.sqrt(eperp**2 - eta**2)))
Ex_theory = 0.37 * np.cos(k_gyro * (z - src_z)).real
Ey_theory = 0.37 * np.sin(k_gyro * (z - src_z)).real

plt.figure(2)
plt.subplot(2,1,1)
plt.plot(z, ex_data, label='Ex (MEEP)')
plt.plot(z, Ex_theory, 'k--')
plt.plot(z, -Ex_theory, 'k--', label='Ex (theory)')
plt.xlim(-L/2, L/2); plt.xlabel('z')
plt.legend()

plt.subplot(2,1,2)
plt.plot(z, ey_data, label='Ey (MEEP)')
plt.plot(z, Ey_theory, 'k--')
plt.plot(z, -Ey_theory, 'k--', label='Ey (theory)')
plt.xlim(-L/2, L/2); plt.xlabel('z')
plt.legend()
plt.tight_layout()
plt.show()
