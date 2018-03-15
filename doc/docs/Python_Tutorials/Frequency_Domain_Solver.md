---
# Frequency Domain Solver
---

This tutorial demonstrates Meep's [frequency-domain solver](../Python_User_Interface/#frequency-domain-solver) which is used to compute the fields produced in a geometry in response to a [constant-frequency source](https://en.wikipedia.org/wiki/Continuous_wave). For more details regarding what this feature is and how it works, refer to this [FAQ](../FAQ/#what-is-meeps-frequency-domain-solver-and-how-does-it-work). This example involves using the frequency-domain solver to compute the fields of a ring resonator which has been described in a [separate tutorial](Basics/#modes-of-a-ring-resonator). We will verify that the error in the computed fields decreases monotonically with decreasing tolerance of the iterative solver.

Only two changes are necessary to the [original simulation](https://github.com/stevengj/meep/blob/master/python/examples/ring.py) which involved time-stepping the fields in response to a broadband source: (1) replace the Gaussian-pulse source with a [continuous-wave (CW) source](../Python_User_Interface/#continuoussource), and (2) turn on complex fields since, by default, real fields are used. Everything else remains unchanged.

Since the frequency-domain solver uses an [iterative method](https://en.wikipedia.org/wiki/Iterative_method), there are two things we can do to improve its convergence: (1) use a non-zero smoothing width (default is 0) for the CW source to reduce the high-frequency oscillations produced by its abrupt turn on (which have slow group velocities and are absorbed poorly by [PML](../Perfectly_Matched_Layer/)), and (2) increase the value of $L$ (default is 2) which is one of the parameters used in the [BiCGSTAB-L](https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method) iterative solver.

We will use the frequency-domain solver to compute the fundamental mode at five different tolerance values on a logarithmic scale. We will plot the L2 norm of the difference of the fields at the first four tolerances with the fields at the lowest tolerance as a function of the tolerance. The simulation script and the resulting plot are shown below.

```py
import meep as mp
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

n = 3.4
w = 1
r = 1
pad = 4
dpml = 2

sxy = 2*(r+w+pad+dpml)

c1 = mp.Cylinder(radius=r+w, material=mp.Medium(index=n))
c2 = mp.Cylinder(radius=r)

fcen = 0.118
df = 0.08
src = mp.Source(mp.ContinuousSource(fcen,fwidth=df), mp.Ez, mp.Vector3(r+0.1))

sim = mp.Simulation(cell_size=mp.Vector3(sxy,sxy),
                    geometry=[c1,c2],
                    sources=[src],
                    resolution=10,
                    force_complex_fields=True,
                    symmetries=[mp.Mirror(mp.Y)],
                    boundary_layers=[mp.PML(dpml)])

num_tols = 5
tols = np.power(10, np.arange(-8.0,-8.0-num_tols,-1.0))
ez_dat = np.zeros((122,122,num_tols), dtype=np.complex_)

for i in range(num_tols):
    sim.init_fields()
    sim.solve_cw(tols[i], 10000, 10)
    ez_dat[:,:,i] = sim.get_array(center=mp.Vector3(), size=mp.Vector3(sxy-2*dpml,sxy-2*dpml), component=mp.Ez)

err_dat = np.zeros(num_tols-1)
for i in range(num_tols-1):
    err_dat[i] = LA.norm(ez_dat[:,:,i]-ez_dat[:,:,num_tols-1])

plt.figure()    
plt.loglog(tols[:num_tols-1], err_dat, 'bo-');
plt.xlabel("frequency-domain solver tolerance");
plt.ylabel("L2 norm of error in fields");
plt.show()

plt.figure()
eps_data = sim.get_array(mp.Vector3(), mp.Vector3(sxy-2*dpml,sxy-2*dpml), mp.Dielectric)
ez_data = ez_dat[:,:,num_tols-1].imag
plt.figure(dpi=100)
plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
plt.imshow(ez_data.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.9)
plt.axis('off')
plt.show()

if np.all(np.diff(err_dat) < 0):
    print("PASSED solve_cw test: error in the fields is decreasing with increasing resolution")
else:
    print("FAILED solve_cw test: error in the fields is NOT decreasing with increasing resolution")
```

The inset of the figure below shows the scalar E$_z$ field (imaginary part), computed using a tolerance of 10$^{-12}$, overlapping the ring-resonator geometry. Note the three-fold mirror symmetry of the field pattern and presence of the point source. The error in the fields decreases monotonically with decreasing tolerance of the frequency-domain solver.

<center>
![](../images/CWsolver-python.png)
</center>