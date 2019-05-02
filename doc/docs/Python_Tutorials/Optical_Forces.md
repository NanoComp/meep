---
# Optical Forces
---

This tutorial demonstrates Meep's ability to compute classical forces via the [Maxwell stress tensor](https://en.wikipedia.org/wiki/Maxwell_stress_tensor) (MST). The geometry consists of two identical, parallel, silicon waveguides with square cross section in vacuum. A schematic of the geometry is shown below. Due to the parallel orientation of the waveguides, the two modes can be chosen to be either symmetric or anti-symmetric with respect to an *x* mirror-symmetry plane between them. As the two waveguides are brought closer and closer together, their modes increasingly couple and give rise to a gradient force that is *transverse* to the waveguide axis (i.e., in the *x* direction). This is different from [radiation pressure](https://en.wikipedia.org/wiki/Radiation_pressure) which involves momentum exchange between photons and is *longitudinal* in nature. An interesting phenomena that occurs for this coupled system is that the force can be tuned to be either attractive or repulsive depending on the relative phase of the modes. This tutorial will demonstrate this effect.

<center>
![](../images/Waveguide_forces.png)
</center>

The gradient force on each waveguide arising from the evanescent coupling of the two waveguide modes can be computed analytically:

$$F=-\frac{1}{\omega}\frac{d\omega}{ds}\Bigg\vert_\vec{k}U,$$

where ω is the mode frequency of the coupled-waveguide system, $s$ is the separation distance between the parallel waveguides, $k$ is the conserved wave vector and $U$ is the total energy of the electromagnetic fields. By convention, negative and positive values correspond to attractive and repulsive forces, respectively. For more details, see [Optics Letters, Vol. 30, pp. 3042-4, 2005](https://www.osapublishing.org/ol/abstract.cfm?uri=ol-30-22-3042). This expression has been shown to be mathematically equivalent to the MST in [Optics Express, Vol. 17, pp. 18116-35, 2009](http://www.opticsinfobase.org/oe/abstract.cfm?URI=oe-17-20-18116). We will verify this result in this tutorial. Note: in this particular example, only the fundamental `ODD_Y` mode shows the bidirectional force.

It is convenient to normalize the force in order to work with dimensionless quantities. Since the total power transmitted through the waveguide is $P=v_gU/L$ where $v_g$ is the group velocity, $L$ is the waveguide length, and $U$ is defined as before, we focus instead on the force per unit length per unit power $(F/L)(ac/P)$ where $a$ is the waveguide width and $c$ is the speed of light. This dimensionless quantity enables us to compute both the flux and the force in a single simulation.

The gradient force can be computed using two different methods: (1) using MPB, compute the frequency and group velocity for a given mode over a range of separation distances and then use a centered [finite-difference](https://en.wikipedia.org/wiki/Finite_difference) scheme to numerically evaluate the formula from above, and (2) using Meep, directly compute both the gradient force and the power transmitted through the waveguide for the guided mode over the same range of separation distances. This tutorial verifies that (1) and (2) produce equivalent results.

The simulation script is in [examples/parallel-wvgs-force.py](https://github.com/NanoComp/meep/blob/master/python/examples/parallel-wvgs-force.py).

The main component of the script is the function `parallel_waveguide(s,xodd)` which computes the Poynting flux and the force given the waveguide separation distance `s` and parity of the waveguide mode `xodd`. Since the eigenmode frequency is not known apriori, a preliminary [`Harminv`](../Python_User_Interface.md#harminv) run is required using a broadband pulsed source. The propagating mode never decays away and the runtime is therefore chosen arbitrarily as 200 time units after the pulsed sources have turned off. Once we have determined the eigenmode frequency, we then replace the `Source` with [`EigenModeSource`](../Python_User_Interface.md#eigenmodesource) to compute: (1) the force on each waveguide due to the mode coupling and (2) the power in the mode. The [eigenmode source](Eigenmode_Source.md) enables a more efficient mode excitation than simply using a constant-amplitude point/area source.

```py
import meep as mp
import numpy as np
import matplotlib.pyplot as plt

resolution = 30   # pixels/μm
    
Si = mp.Medium(index=3.45)

dpml = 1.0
pml_layers = [mp.PML(dpml)]
    
sx = 5
sy = 3
cell = mp.Vector3(sx+2*dpml,sy+2*dpml,0)

a = 1.0     # waveguide width

k_point = mp.Vector3(z=0.5)

fcen = 0.22
df = 0.06

def parallel_waveguide(s,xodd):
    geometry = [mp.Block(center=mp.Vector3(-0.5*(s+a)),
                         size=mp.Vector3(a,a,mp.inf),
                         material=Si),
                mp.Block(center=mp.Vector3(0.5*(s+a)),
                         size=mp.Vector3(a,a,mp.inf),
                         material=Si)]

    symmetries = [mp.Mirror(mp.X, phase=-1.0 if xodd else 1.0),
                  mp.Mirror(mp.Y, phase=-1.0)]

    sources = [mp.Source(src=mp.GaussianSource(fcen, fwidth=df),
                         component=mp.Ey,
                         center=mp.Vector3(-0.5*(s+a)),
                         size=mp.Vector3(a,a)),
               mp.Source(src=mp.GaussianSource(fcen, fwidth=df),
                         component=mp.Ey,
                         center=mp.Vector3(0.5*(s+a)),
                         size=mp.Vector3(a,a),
                         amplitude=-1.0 if xodd else 1.0)]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell,
                        boundary_layers=pml_layers,
                        geometry=geometry,
                        symmetries=symmetries,
                        k_point=k_point,
                        sources=sources)

    h = mp.Harminv(mp.Ey, mp.Vector3(0.5*(s+a)), fcen, df)

    sim.run(mp.after_sources(h), until_after_sources=200)

    f = h.modes[0].freq
    print("freq:, {}, {}".format(s, f))

    sim.reset_meep()

    eig_sources = [mp.EigenModeSource(src=mp.GaussianSource(f, fwidth=df),
                                      size=mp.Vector3(a,a),
                                      center=mp.Vector3(-0.5*(s+a)),
                                      eig_kpoint=k_point,
                                      eig_match_freq=True,
                                      eig_parity=mp.ODD_Y),
                   mp.EigenModeSource(src=mp.GaussianSource(f, fwidth=df),
                                      size=mp.Vector3(a,a),
                                      center=mp.Vector3(0.5*(s+a)),
                                      eig_kpoint=k_point,
                                      eig_match_freq=True,
                                      eig_parity=mp.ODD_Y,
                                      amplitude=-1.0 if xodd else 1.0)]

    sim.change_sources(eig_sources)

    flux_reg = mp.FluxRegion(direction=mp.Z, center=mp.Vector3(), size=mp.Vector3(1.2*(2*a+s),1.2*a))
    wvg_flux = sim.add_flux(f, 0, 1, flux_reg)

    force_reg1 = mp.ForceRegion(mp.Vector3(0.5*s), direction=mp.X, weight=1.0, size=mp.Vector3(y=a))
    force_reg2 = mp.ForceRegion(mp.Vector3(0.5*s+a), direction=mp.X, weight=-1.0, size=mp.Vector3(y=a))
    wvg_force = sim.add_force(f, 0, 1, force_reg1, force_reg2)

    sim.run(until_after_sources=5000)

    flux = mp.get_fluxes(wvg_flux)[0]
    force = mp.get_forces(wvg_force)[0]
    
    sim.reset_meep()
    return flux, force
```

There are two important items to note in `parallel_waveguide`: (1) a single flux surface is used to compute the Poynting flux in $z$ which spans an area slightly larger than both waveguides rather than two separate flux surfaces (one for each waveguide). This is because in the limit of small separation, two flux surfaces overlap whereas the total power through a single flux surface need, by symmetry, only be halved in order to determine the value for just one of the two waveguides. (2) Instead of defining a closed, four-sided "box" surrounding the waveguides, the MST is computed along just two $y$-oriented lines (to obtain the force in the $x$ direction) with different `weight` values to correctly sum the total force. By symmetry, the force in the $y$ direction is zero and need not be computed. Choosing a suitable runtime requires some care. A large runtime is necessary to obtain the steady-state response but this will also lead to large values for the discrete Fourier-transformed fields used to compute both the flux and the MST. Large floating-point numbers may contain [roundoff errors](https://en.wikipedia.org/wiki/Round-off_error).

The simulation is run over the range of separation distances from 0.05 to 1.00 μm in increments of 0.05 μm. The results are compared with those from MPB. This is shown in the figure above. The two methods show good agreement.

```py
s = np.arange(0.05,1.05,0.05)
fluxes_odd = np.zeros(s.size)
forces_odd = np.zeros(s.size)

fluxes_even = np.zeros(s.size)
forces_even = np.zeros(s.size)

for k in range(len(s)):
    fluxes_odd[k], forces_odd[k] = parallel_waveguide(s[k],True)
    fluxes_even[k], forces_even[k] = parallel_waveguide(s[k],False)

plt.figure(dpi=150)
plt.plot(s,-forces_odd/fluxes_odd,'rs',label='anti symmetric')
plt.plot(s,-forces_even/fluxes_even,'bo',label='symmetric')
plt.grid(True)
plt.xlabel('waveguide separation s/a')    
plt.ylabel('optical force (F/L)(ac/P)')
plt.legend(loc='upper right')
plt.show()
```

The MPB simulation is in [examples/parallel-wvgs-mpb.py](https://github.com/NanoComp/meep/blob/master/python/examples/parallel-wvgs-mpb.py). Note: since MPB permits symmetries only in the $y$ and $z$ directions, the coordinate axes used in the MPB script to define the waveguide geometry are different than those in the Meep script. In MPB, the propagating axis is $x$ whereas in Meep it is $z$.

```py
import meep as mp
from meep import mpb
import numpy as np
import matplotlib.pyplot as plt

resolution = 128  # pixels/μm

Si = mp.Medium(index=3.45)

syz = 10
geometry_lattice = mp.Lattice(size=mp.Vector3(0,syz,syz))

k_points = [mp.Vector3(0.5)]

num_bands = 1
tolerance = 1e-9

a = 1.0  # waveguide width

def parallel_waveguide(s,yodd):
    geometry = [mp.Block(center=mp.Vector3(0,-0.5*(s+a),0),
                         size=mp.Vector3(mp.inf,a,a),
                         material=Si),
                mp.Block(center=mp.Vector3(0,0.5*(s+a),0),
                         size=mp.Vector3(mp.inf,a,a),
                         material=Si)]

    ms = mpb.ModeSolver(resolution=resolution,
                        k_points=k_points,
                        geometry_lattice=geometry_lattice,
                        geometry=geometry,
                        num_bands=num_bands,
                        tolerance=tolerance)

    if yodd:
        ms.run_yodd_zodd()
    else:
        ms.run_yeven_zodd()

    f = ms.get_freqs()[0]
    vg = ms.compute_group_velocity_component(mp.Vector3(1,0,0))[0]

    return f,vg

ss = np.arange(0.05,1.05,0.05)

f_odd = np.zeros(len(ss))
vg_odd = np.zeros(len(ss))
f_even = np.zeros(len(ss))
vg_even = np.zeros(len(ss))

for j in range(len(ss)):
    f_odd[j], vg_odd[j] = parallel_waveguide(ss[j],True)
    f_even[j], vg_even[j] = parallel_waveguide(ss[j],False)

ds = ss[1]-ss[0]

def compute_force(f,vg):
    f_avg = 0.5*(f[:-1]+f[1:])
    df = f[1:]-f[:-1]
    vg_avg = 0.5*(vg[:-1]+vg[1:])
    return np.multiply(np.multiply(-1/f_avg,df/ds), 1/vg_avg)

force_odd = compute_force(f_odd,vg_odd)
force_even = compute_force(f_even,vg_even)

plt.plot(ss[:-1],force_odd,'b-',label='antisymmetric')
plt.plot(ss[:-1],force_even,'r-',label='symmetric')
plt.xlabel("waveguide separation s/a")
plt.ylabel("optical force (F/L)(ac/P)")
plt.legend(loc='upper right')
plt.xticks(np.arange(0,1.2,0.2))
plt.yticks(np.arange(-1.5,1.0,0.5))
plt.show()
```
