---
# Mode Decomposition
---

This tutorial demonstrates the mode-decomposition feature which is used to decompose a given mode profile into a superposition of harmonic basis modes. There are examples for two kinds of modes in lossless, dielectric media: (1) localized (i.e., guided) and (2) non-localized (i.e., planewave).

[TOC]

Reflectance of a Waveguide Taper
--------------------------------

This example involves computing the reflectance &mdash; the fraction of the incident power which is reflected &mdash; of the fundamental mode of a linear waveguide taper. The structure and the simulation parameters are shown in the schematic below. We will verify that computing the reflectance using two different methods produces nearly identical results: (1) mode decomposition and (2) directly from the fields via its [Poynting flux](../Introduction.md#transmittancereflectance-spectra). Also, we will demonstrate that the scaling of the reflectance with the taper length is quadratic, consistent with analytical results from [Optics Express, Vol. 16, pp. 11376-92, 2008](http://www.opticsinfobase.org/abstract.cfm?URI=oe-16-15-11376).

<center>
![](../images/waveguide-taper.png)
</center>

The structure, which can be viewed as a [two-port network](https://en.wikipedia.org/wiki/Two-port_network), consists of a single-mode waveguide of width `w1` (1 μm)  coupled to a second waveguide of width `w2` (2 μm) via a linearly-sloped taper of length `Lt`. The structure is homogeneous ε=12 in vacuum. PML absorbing boundaries surround the computational cell. An eigenmode source with E<sub>z</sub> polarization is used to launch the fundamental mode with a wavelength of 6.67 μm. There is an eigenmode-expansion monitor placed at the midpoint of the first waveguide. This is a line monitor which extends beyond the waveguide in order to span the entire mode profile including its evanescent tails. The Fourier-transformed fields along this line monitor are used to compute the basis coefficients of the harmonic modes. These are computed separately via the eigenmode solver [MPB](https://mpb.readthedocs.io). Technical details are described in [Mode Decomposition](../Mode_Decomposition.md) where it is shown that the squared magnitude of the mode coefficient is equivalent to the power (i.e., Poynting flux) in the given eigenmode. The ratio of the complex mode coefficients can be used to compute the [S parameters](https://en.wikipedia.org/wiki/Scattering_parameters). In this example, we are computing |S<sub>11</sub>|<sup>2</sup> which is the reflectance. Another line monitor could have been placed in the second waveguide to compute the transmittance or |S<sub>21</sub>|<sup>2</sup>. Also not considered is the scattering into radiative modes.

The structure has mirror symmetry in the $y$ direction which can be exploited to reduce the computation size by a factor of two. This requires that we use `add_flux` and specify `eig_parity=mp.ODD_Z+mp.EVEN_Y` in the call to `get_eigenmode_coefficients`. This is because `add_mode_monitor`, which is an alias for `add_flux`, is not optimized for symmetry.

At the end of the simulation, the reflectance of the fundamental mode computed using the two methods are displayed. The simulation script is in [examples/mode-decomposition.py](https://github.com/stevengj/meep/blob/master/python/examples/mode-decomposition.py).

```py
import meep as mp

resolution = 61   # pixels/μm

w1 = 1            # width of waveguide 1
w2 = 2            # width of waveguide 2
Lw = 10           # length of waveguide 1/2

dair = 3.0        # length of air region
dpml_x = 6.0      # length of PML in x direction
dpml_y = 2.0      # length of PML in y direction

sy = dpml_y+dair+w2+dair+dpml_y

Si = mp.Medium(epsilon=12.0)

boundary_layers = [mp.PML(dpml_x,direction=mp.X),
                   mp.PML(dpml_y,direction=mp.Y)]

# mode wavelength
lcen = 6.67
# mode frequency
fcen = 1/lcen

symmetries = [mp.Mirror(mp.Y)]

for m in range(5):
    Lt = 2**m
    sx = dpml_x+Lw+Lt+Lw+dpml_x
    cell_size = mp.Vector3(sx,sy,0)

    src_pt = mp.Vector3(-0.5*sx+dpml_x+0.2*Lw,0,0)
    sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen,fwidth=0.2*fcen),
                                  component=mp.Ez,
                                  center=src_pt,
                                  size=mp.Vector3(0,sy-2*dpml_y,0),
                                  eig_match_freq=True,
                                  eig_parity=mp.ODD_Z+mp.EVEN_Y)]

    # straight waveguide
    vertices = [mp.Vector3(-0.5*sx-1,0.5*w1),
                mp.Vector3(0.5*sx+1,0.5*w1),
                mp.Vector3(0.5*sx+1,-0.5*w1),
                mp.Vector3(-0.5*sx-1,-0.5*w1)]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=boundary_layers,
                        geometry=[mp.Prism(vertices,height=mp.inf,material=Si)],
                        sources=sources,
                        symmetries=symmetries)

    mon_pt = mp.Vector3(-0.5*sx+dpml_x+0.7*Lw,0,0)
    flux = sim.add_flux(fcen,0,1,mp.FluxRegion(center=mon_pt,size=mp.Vector3(0,sy-2*dpml_y,0)))

    sim.run(until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,mon_pt,1e-9))

    res = sim.get_eigenmode_coefficients(flux,[1],eig_parity=mp.ODD_Z+mp.EVEN_Y)
    incident_coeffs = res.alpha
    incident_flux = mp.get_fluxes(flux)
    incident_flux_data = sim.get_flux_data(flux)

    sim.reset_meep()    
    
    # linear taper
    vertices = [mp.Vector3(-0.5*sx-1,0.5*w1),
                mp.Vector3(-0.5*Lt,0.5*w1),
                mp.Vector3(0.5*Lt,0.5*w2),
                mp.Vector3(0.5*sx+1,0.5*w2),
                mp.Vector3(0.5*sx+1,-0.5*w2),
                mp.Vector3(0.5*Lt,-0.5*w2),
                mp.Vector3(-0.5*Lt,-0.5*w1),
                mp.Vector3(-0.5*sx-1,-0.5*w1)]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=boundary_layers,
                        geometry=[mp.Prism(vertices,height=mp.inf,material=Si)],
                        sources=sources,
                        symmetries=symmetries)

    refl_flux = sim.add_flux(fcen,0,1,mp.FluxRegion(center=mon_pt,size=mp.Vector3(0,sy-2*dpml_y,0)))
    sim.load_minus_flux_data(refl_flux,incident_flux_data)
    
    sim.run(until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,mon_pt,1e-9))
        
    res = sim.get_eigenmode_coefficients(refl_flux,[1],eig_parity=mp.ODD_Z+mp.EVEN_Y)
    coeffs = res.alpha
    taper_flux = mp.get_fluxes(refl_flux)
    print("refl:, {}, {:.8f}, {:.8f}".format(Lt,abs(coeffs[0,0,1])**2/abs(incident_coeffs[0,0,0])**2,-taper_flux[0]/incident_flux[0]))
```

We compute the reflectance for five different taper lengths: 1, 2, 4, 8, and 16 μm. A quadratic scaling of the reflectance with the taper length appears as a straight line on a log-log plot. The Bash commands to run the simulation and extract the plotting data from the output are:

```sh
mpirun -np 3 python -u mode-decomposition.py |tee taper_data.out;
grep refl: taper_data.out |cut -d , -f2- > taper_data.dat
```

The results are plotted using the Python script below. The plot is shown in the accompanying figure. The reflectance values computed using the two methods of mode decomposition and flux are nearly identical. For reference, a line with quadratic scaling is shown in black. The reflectance of the linear waveguide taper decreases quadratically with the taper length consistent with analytic theory.

```py
import numpy as np
import matplotlib.pyplot as plt

f = np.genfromtxt("taper_data.dat", delimiter=",")
Lt = f[:,0]
Rmode = f[:,1]
Rflux = f[:,2]

plt.figure(dpi=150)
plt.loglog(Lt,Rmode,'bo-',label='mode decomposition')
plt.loglog(Lt,Rflux,'ro-',label='flux')
plt.loglog(Lt,0.005/Lt**2,'k-',label=r'quadratic reference (1/Lt$^2$)')
plt.legend(loc='upper right')
plt.xlabel('taper length Lt (μm)')
plt.ylabel('reflectance')
plt.show()
```

<center>
![](../images/refl_coeff_vs_taper_length.png)
</center>

Diffraction Spectrum of a Binary Grating
----------------------------------------

The mode-decomposition feature can also be applied to planewaves in homogeneous media with scalar permittivity/permeability (i.e., no anisotropy). This will be demonstrated in this example to compute the diffraction spectrum of a binary phase [grating](https://en.wikipedia.org/wiki/Diffraction_grating). The unit cell geometry of the grating is shown in the schematic below. The grating is periodic in the $y$ direction with periodicity `gp` and has a rectangular profile of height `gh` and duty cycle `gdc`. The grating parameters are `gh`=0.5 μm, `gdc`=0.5, and `gp`=10 μm. There is a semi-infinite substrate of thickness `dsub` adjacent to the grating. The substrate and grating are glass with a refractive index of 1.5. The surrounding is air/vacuum. Perfectly matched layers (PML) of thickness `dpml` are used in the $\pm x$ boundaries.

### Transmittance Spectra for Planewave at Normal Incidence

A pulsed planewave with E<sub>z</sub> polarization spanning wavelengths of 0.4 to 0.6 μm is normally incident on the grating from the glass substrate. The eigenmode monitor is placed in the air region. We will use mode decomposition to compute the transmittance &mdash; the ratio of the power in the $+x$ direction of the diffracted mode relative to that of the incident planewave &mdash; for the first ten diffraction orders. Two simulations are required: (1) an *empty* cell of homogeneous glass to obtain the incident power of the source, and (2) the grating structure to obtain the diffraction orders. At the end of the simulation, the wavelength, angle, and transmittance for each diffraction order are computed.

The simulation script is in [examples/binary_grating.py](https://github.com/stevengj/meep/blob/master/python/examples/binary_grating.py).

<center>
![](../images/grating.png)
</center>

```py
import meep as mp
import math

resolution = 60        # pixels/μm

dpml = 1.0             # PML thickness
dsub = 3.0             # substrate thickness
dpad = 3.0             # padding between grating and PML
gp = 10.0              # grating period
gh = 0.5               # grating height
gdc = 0.5              # grating duty cycle

sx = dpml+dsub+gh+dpad+dpml
sy = gp

cell_size = mp.Vector3(sx,sy,0)
pml_layers = [mp.PML(thickness=dpml,direction=mp.X)]

wvl_min = 0.4           # min wavelength
wvl_max = 0.6           # max wavelength
fmin = 1/wvl_max        # min frequency
fmax = 1/wvl_min        # max frequency
fcen = 0.5*(fmin+fmax)  # center frequency
df = fmax-fmin          # frequency width

src_pt = mp.Vector3(-0.5*sx+dpml+0.5*dsub,0,0)
sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ez, center=src_pt, size=mp.Vector3(0,sy,0))]

k_point = mp.Vector3(0,0,0)

glass = mp.Medium(index=1.5)

symmetries=[mp.Mirror(mp.Y)]

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    k_point=k_point,
                    default_material=glass,
                    sources=sources,
                    symmetries=symmetries)

nfreq = 21
mon_pt = mp.Vector3(0.5*sx-dpml-0.5*dpad,0,0)
flux_mon = sim.add_flux(fcen, df, nfreq, mp.FluxRegion(center=mon_pt, size=mp.Vector3(0,sy,0)))

sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mon_pt, 1e-9))

input_flux = mp.get_fluxes(flux_mon)

sim.reset_meep()

geometry = [mp.Block(material=glass, size=mp.Vector3(dpml+dsub,mp.inf,mp.inf), center=mp.Vector3(-0.5*sx+0.5*(dpml+dsub),0,0)),
            mp.Block(material=glass, size=mp.Vector3(gh,gdc*gp,mp.inf), center=mp.Vector3(-0.5*sx+dpml+dsub+0.5*gh,0,0))]

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    k_point=k_point,
                    sources=sources,
                    symmetries=symmetries)

mode_mon = sim.add_flux(fcen, df, nfreq, mp.FluxRegion(center=mon_pt, size=mp.Vector3(0,sy,0)))

sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mon_pt, 1e-9))

freqs = mp.get_eigenmode_freqs(mode_mon)

nmode = 10
res = sim.get_eigenmode_coefficients(mode_mon, range(1,nmode+1), eig_parity=mp.ODD_Z+mp.EVEN_Y)
coeffs = res.alpha
kdom = res.kdom

for nm in range(nmode):
  for nf in range(nfreq):
    mode_wvl = 1/freqs[nf]
    mode_angle = math.degrees(math.acos(kdom[nm*nfreq+nf].x/freqs[nf]))
    mode_tran = abs(coeffs[nm,nf,0])**2/input_flux[nf]
    if nm != 0:
      mode_tran = 0.5*mode_tran
    print("grating{}:, {:.5f}, {:.2f}, {:.8f}".format(nm,mode_wvl,mode_angle,mode_tran))
```

Note the use of the keyword parameter argument `eig_parity=mp.ODD_Z+mp.EVEN_Y` in the call to `get_eigenmode_coefficients`. This is important for specifying **non-degenerate** modes in MPB since the `k_point` is (0,0,0). `ODD_Z` is for modes with E<sub>z</sub> polarization. `EVEN_Y` is necessary since each diffraction order which is based on a given k<sub>x</sub> consists of *two* modes: one going in the +y direction and the other in the -y direction. `EVEN_Y` forces MPB to compute only the +k<sub>y</sub> + -k<sub>y</sub> (cosine) mode. As a result, the total transmittance must be halved in this case to obtain the transmittance for the individual +k<sub>y</sub> or -k<sub>y</sub> mode. For `ODD_Y`, MPB will compute the +k<sub>y</sub> - -k<sub>y</sub> (sine) mode but this will have zero power because the source is even. If the $y$ parity is left out, MPB will return a random superposition of the cosine and sine modes. Alternatively, in this example an input planewave with H<sub>z</sub> instead of E<sub>z</sub> polarization can be used which requires `eig_parity=mp.EVEN_Z+mp.ODD_Y` as well as an odd mirror symmetry plane in *y*. Finally, note the use of `add_flux` instead of `add_mode_monitor` when using symmetries.

The simulation is run and the results piped to a file (the grating data is extracted to a separate file for plotting) using the following shell script:

```sh
#!/bin/bash

python -u binary_grating.py |tee grating.out
grep grating grating.out |cut -d , -f2- > grating.dat
```

The diffraction spectrum is plotted using the following script and shown in the figure below.

```py

import matplotlib.pyplot as plt
import numpy as np

d = np.genfromtxt("grating.dat",delimiter=",")

nmode = 10
nfreq = 21

thetas = np.empty((nmode,nfreq))
wvls = np.empty((nmode,nfreq))
tran = np.empty((nmode,nfreq))

for j in range(nfreq):
    tran[:,j] = d[j::nfreq,2]
    thetas[:,j] = d[j::nfreq,1]
    wvls[:,j] = d[j,0]

tran_max = round(tran.max(),1)

plt.figure(dpi=150)
plt.pcolormesh(wvls, thetas, tran, cmap='Blues', shading='flat', vmin=0, vmax=tran_max)
plt.axis([wvls.min(), wvls.max(), thetas.min(), thetas.max()])
plt.xlabel("wavelength (μm)")
plt.ylabel("diffraction angle (degrees)")
plt.xticks([t for t in np.arange(0.4,0.7,0.1)])
plt.yticks([t for t in range(0,35,5)])
plt.title("transmittance of diffraction orders")
cbar = plt.colorbar()
cbar.set_ticks([t for t in np.arange(0,tran_max+0.1,0.1)])
cbar.set_ticklabels(["{:.1f}".format(t) for t in np.arange(0,tran_max+0.1,0.1)])
plt.show()
```

Each diffraction order corresponds to a single angle. In the figure below, this angle is represented by the *lower* boundary of each labeled region. For example, the m=0 order has a diffraction angle of 0° at all wavelengths. The representation of the diffraction orders as finite angular regions is an artifact of matplotlib's [pcolormesh](https://matplotlib.org/api/_as_gen/matplotlib.pyplot.pcolormesh.html) routine. Note that only the positive diffraction orders are shown as these are equivalent to the negative orders due to the symmetry of the source and the structure. The transmittance of each diffraction order should ideally be a constant for all wavelengths. The slight wavelength dependence shown in the figure is due to numerical discretization which can be mitigated by increasing the resolution.

The diffraction orders/modes are a finite set of propagating planewaves. The wavevector k<sub>x</sub> of these modes can be computed analytically: for a frequency of ω (in c=1 units), these propagating modes are the **real** solutions of sqrt(ω²n²-(k<sub>y</sub>+2πm/Λ)²) where m is the diffraction order (an integer), Λ is the periodicity of the grating, and n is the refractive index of the propagating medium. In this example, n=1, k<sub>y</sub>=0, and Λ=10 μm. Thus, at a wavelength of 0.5 μm there are a total of 20 diffraction orders of which we only computed the first 10. The wavevector k<sub>x</sub> is used to compute the angle of the diffraction order as cos<sup>-1</sup>(k<sub>x</sub>/ω). Evanescent modes, those with an imaginary k<sub>x</sub>, exist for |m|>20 but these modes carry no power. Note that currently Meep does not compute the number of propagating modes for you. If the mode number passed to `get_eigenmode_coefficients` is larger than the number of propagating modes at a given frequency/wavelength, MPB's Newton solver will fail to converge and will return zero for the mode coefficient. It is therefore a good idea to know beforehand the number of propagating modes.

<center>
![](../images/grating_diffraction_spectra.png)
</center>

In the limit where the grating periodicity is much larger than the wavelength and the size of the diffracting element (i.e., more than 10 times), as it is in this example, the [diffraction efficiency](https://en.wikipedia.org/wiki/Diffraction_efficiency) can be computed analytically using scalar theory. This is described in the OpenCourseWare [Optics course](https://ocw.mit.edu/courses/mechanical-engineering/2-71-optics-spring-2009/) in the Lecture 16 (Gratings: Amplitude and Phase, Sinusoidal and Binary) [notes](https://ocw.mit.edu/courses/mechanical-engineering/2-71-optics-spring-2009/video-lectures/lecture-16-gratings-amplitude-and-phase-sinusoidal-and-binary/MIT2_71S09_lec16.pdf) and [video](https://www.youtube.com/watch?v=JmWguqCZRxk). For a review of scalar diffraction theory, see Chapter 3 ("Analysis of Two-Dimensional Signals and Systems") of [Introduction to Fourier Optics (fourth edition)](https://www.amazon.com/Introduction-Fourier-Optics-Joseph-Goodman-ebook/dp/B076TBP48F) by J.W. Goodman. From the scalar theory, the diffraction efficiency of the binary grating is 4/(mπ)<sup>2</sup> when the phase difference between the propagating distance in the glass relative to the same distance in air is π. The phase difference/contrast is (2π/λ)(n-1)s where λ is the wavelength, n is the refractive index of the grating, and s is the propagation distance in the grating (`gh` in the script). A special feature of the binary grating is that the diffraction efficiency is 0 for all *even* orders. This is verified by the diffraction spectrum shown above.

To convert the diffraction efficiency into transmittance in the *x* direction (in order to be able to compare the scalar-theory results with those from Meep), the diffraction efficiency must be multiplied by the Fresnel transmittance from air to glass and by the cosine of the diffraction angle. We compare the analytic and simulated results at a wavelength of 0.5 μm for diffraction orders 1, 3, 5, and 7. The analytic results are 0.3886, 0.0427, 0.0151, and 0.0074. The Meep results are 0.3891, 0.04287, 0.0152, and 0.0076. This corresponds to relative errors of approximately 1.3%, 0.4%, 0.8%, and 2.1% which indicates good agreement.

### Reflectance and Transmittance Spectra for Planewave at Oblique Incidence

As an additional demonstration of the mode-decomposition feature, the reflectance and transmittance of all diffracted orders for any grating with no material absorption and a planewave source incident at any arbitrary angle and wavelength must necessarily sum to unity. Also, the total reflectance and transmittance must be equivalent to values computed using the Poynting flux. This demonstration is somewhat similar to the [single-mode waveguide example](#reflectance-of-a-waveguide-taper).

The following script is adapted from the previous binary-grating example involving a [normally-incident planewave](#transmittance-spectra-for-planewave-at-normal-incidence). The total reflectance, transmittance, and their sum are displayed at the end of the simulation on two separate lines prefixed by `mode-coeff:` and `poynting-flux:`.

Results are computed for a single wavelength of 0.5 μm. The pulsed planewave is incident at an angle of 10.7°. Its spatial profile is defined using the source amplitude function `pw_amp`. This [anonymous function](https://en.wikipedia.org/wiki/Anonymous_function) takes two arguments, the wavevector and a point in space (both `mp.Vector3`s), and returns a function of one argument which defines the planewave amplitude at that point. A narrow bandwidth pulse is used in order to mitigate the intrinsic discretization effects of the [Yee grid](../Yee_Lattice.md) for oblique planewaves. Also, the `stop_when_fields_decayed` termination criteria is replaced with `until_after_sources`. As a general rule of thumb, the more oblique the planewave source, the longer the run time required to ensure accurate results. There is an additional line monitor between the source and the grating for computing the reflectance. The angle of each reflected/transmitted mode, which can be positive or negative, is computed using its dominant planewave vector. Since the oblique source breaks the symmetry in the $y$ direction, each diffracted order must be computed separately. In total, there are 59 reflected and 39 transmitted orders.

The simulation script is in [examples/binary_grating_oblique.py](https://github.com/stevengj/meep/blob/master/python/examples/binary_grating_oblique.py).

```py
import meep as mp
import math
import cmath
import numpy as np

resolution = 50        # pixels/μm

dpml = 1.0             # PML thickness
dsub = 3.0             # substrate thickness
dpad = 3.0             # length of padding between grating and PML
gp = 10.0              # grating period
gh = 0.5               # grating height
gdc = 0.5              # grating duty cycle

sx = dpml+dsub+gh+dpad+dpml
sy = gp

cell_size = mp.Vector3(sx,sy,0)
pml_layers = [mp.PML(thickness=dpml,direction=mp.X)] 

wvl = 0.5              # center wavelength
fcen = 1/wvl           # center frequency
df = 0.05*fcen         # frequency width

ng = 1.5
glass = mp.Medium(index=ng)

use_cw_solver = False  # CW solver or time stepping?
tol = 1e-6             # CW solver tolerance
max_iters = 2000       # CW solver max iterations
L = 10                 # CW solver L

# rotation angle of incident planewave; CCW about Z axis, 0 degrees along +X axis
theta_in = math.radians(10.7)

# k (in source medium) with correct length (plane of incidence: XY)
k = mp.Vector3(math.cos(theta_in),math.sin(theta_in),0).scale(fcen*ng)

symmetries = []
eig_parity = mp.ODD_Z
if theta_in == 0:
  k = mp.Vector3(0,0,0)
  symmetries = [mp.Mirror(mp.Y)]
  eig_parity += mp.EVEN_Y

def pw_amp(k,x0):
  def _pw_amp(x):
    return cmath.exp(1j*2*math.pi*k.dot(x+x0))
  return _pw_amp

src_pt = mp.Vector3(-0.5*sx+dpml+0.3*dsub,0,0)
sources = [mp.Source(mp.ContinuousSource(fcen,fwidth=df) if use_cw_solver else mp.GaussianSource(fcen,fwidth=df),
                     component=mp.Ez,
                     center=src_pt,
                     size=mp.Vector3(0,sy,0),
                     amp_func=pw_amp(k,src_pt))]

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    k_point=k,
                    default_material=glass,
                    sources=sources,
                    symmetries=symmetries)

refl_pt = mp.Vector3(-0.5*sx+dpml+0.5*dsub,0,0)
refl_flux = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=refl_pt, size=mp.Vector3(0,sy,0)))

if use_cw_solver:
  sim.init_sim()
  sim.solve_cw(tol, max_iters, L)
else:
  sim.run(until_after_sources=100)

input_flux = mp.get_fluxes(refl_flux)
input_flux_data = sim.get_flux_data(refl_flux)

sim.reset_meep()

geometry = [mp.Block(material=glass, size=mp.Vector3(dpml+dsub,mp.inf,mp.inf), center=mp.Vector3(-0.5*sx+0.5*(dpml+dsub),0,0)),
            mp.Block(material=glass, size=mp.Vector3(gh,gdc*gp,mp.inf), center=mp.Vector3(-0.5*sx+dpml+dsub+0.5*gh,0,0))]

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    k_point=k,
                    sources=sources,
                    symmetries=symmetries)

refl_flux = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=refl_pt, size=mp.Vector3(0,sy,0)))
sim.load_minus_flux_data(refl_flux,input_flux_data)

tran_pt = mp.Vector3(0.5*sx-dpml-0.5*dpad,0,0)
tran_flux = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=tran_pt, size=mp.Vector3(0,sy,0)))

if use_cw_solver:
  sim.init_sim()
  sim.solve_cw(tol, max_iters, L)
else:
  sim.run(until_after_sources=200)

nm_r = np.floor((fcen*ng-k.y)*gp)-np.ceil((-fcen*ng-k.y)*gp) # number of reflected orders
if theta_in == 0:
  nm_r = nm_r/2 # since eig_parity removes degeneracy in y-direction
nm_r = int(nm_r)

res = sim.get_eigenmode_coefficients(refl_flux, range(1,nm_r+1), eig_parity=eig_parity)
r_coeffs = res.alpha

Rsum = 0
for nm in range(nm_r):
  r_kdom = res.kdom[nm]
  Rmode = abs(r_coeffs[nm,0,1])**2/input_flux[0]
  r_angle = np.sign(r_kdom.y)*math.acos(r_kdom.x/(ng*fcen))
  print("refl:, {}, {:.2f}, {:.8f}".format(nm,math.degrees(r_angle),Rmode))
  Rsum += Rmode

nm_t = np.floor((fcen-k.y)*gp)-np.ceil((-fcen-k.y)*gp)       # number of transmitted orders
if theta_in == 0:
  nm_t = nm_t/2 # since eig_parity removes degeneracy in y-direction
nm_t = int(nm_t)

res = sim.get_eigenmode_coefficients(tran_flux, range(1,nm_t+1), eig_parity=eig_parity)
t_coeffs = res.alpha

Tsum = 0
for nm in range(nm_t):
  t_kdom = res.kdom[nm]
  Tmode = abs(t_coeffs[nm,0,0])**2/input_flux[0]
  t_angle = np.sign(t_kdom.y)*math.acos(t_kdom.x/fcen)
  print("tran:, {}, {:.2f}, {:.8f}".format(nm,math.degrees(t_angle),Tmode))
  Tsum += Tmode

print("mode-coeff:, {:.6f}, {:.6f}, {:.6f}".format(Rsum,Tsum,Rsum+Tsum))

r_flux = mp.get_fluxes(refl_flux)
t_flux = mp.get_fluxes(tran_flux)
Rflux = -r_flux[0]/input_flux[0]
Tflux =  t_flux[0]/input_flux[0]
print("poynting-flux:, {:.6f}, {:.6f}, {:.6f}".format(Rflux,Tflux,Rflux+Tflux))
```

Since this is a single-wavelength calculation, the [frequency-domain solver](../Python_User_Interface.md#frequency-domain-solver) can be used instead of time stepping for a possible performance enhancement. The only changes necessary to the original script are to replace two objects: (1) `GaussianSource` with `ContinuousSource` and (2) `run` with `solve_cw`. Choosing which approach to use is determined by the `use_cw_solver` boolean variable. In this example, mainly because of the oblique source, the frequency-domain solver converges slowly and is less efficient than the time-stepping simulation. The results from both approaches are nearly identical. Time stepping is therefore the default.

The following are several lines of output for eight of the reflected and transmitted orders. The first numerical column is the mode number, the second is the mode angle (in degrees), and the third is the fraction of the input power that is concentrated in the mode. Note that the thirteenth transmitted order at 19.18° contains nearly 38% of the input power.

```
...
refl:, 7, 6.83, 0.00006655
refl:, 8, -8.49, 0.00005703
refl:, 9, 8.76, 0.00015782
refl:, 10, -10.43, 0.00001277
refl:, 11, 10.70, 0.04414104
refl:, 12, -12.38, 0.00005981
refl:, 13, 12.65, 0.00041466
refl:, 14, -14.34, 0.00001991
...
```

```
...
tran:, 12, -18.75, 0.00095295
tran:, 13, 19.18, 0.38261656
tran:, 14, -21.81, 0.00198510
tran:, 15, 22.24, 0.00107184
tran:, 16, -24.93, 0.00098452
tran:, 17, 25.37, 0.04148787
tran:, 18, -28.13, 0.00137329
tran:, 19, 28.59, 0.00113850
...
```

The mode number is equivalent to the band index from the MPB calculation. The ordering of the modes is according to *decreasing* values of k<sub>x</sub>. The first mode has the largest k<sub>x</sub> and thus angle closest to 0°. As a corollary, the first mode has the smallest |k<sub>y</sub>+2πm/Λ|. For a non-zero k<sub>y</sub> (as in the case of an obliquely incident source), this expression will not necessarily be zero. The first seven reflected modes have m values of -3, -4, -2, -5, -1, -6, and 0. These m values are not monotonic. This is because k<sub>x</sub> is a nonlinear function of m as shown earlier. The ordering of the transmitted modes is different since these modes are in vacuum and not glass (recall that the medium's refractive index is also a part of this nonlinear function). In the first example involving a normally incident source with k<sub>y</sub>=0, the ordering of the modes is monotonic: m = 0, &#177;1, &#177;2, ...

The two main lines of the output are:

```
mode-coeff:, 0.061007, 0.937897, 0.998904
poynting-flux:, 0.061063, 0.938384, 0.999447
```

The first numerical column is the total reflectance, the second is the total transmittance, and the third is their sum. Results from the mode coefficients agree with the Poynting flux values to three decimal places. Also, the total reflectance and transmittance sum to unity. These results indicate that approximately 6% of the input power is reflected and the remaining 94% is transmitted.

Diffraction Spectra of Liquid-Crystal Polarization Gratings
-----------------------------------------------------------

As a final demonstration of mode decomposition, we compute the diffraction spectrum of a [liquid-crystal](https://en.wikipedia.org/wiki/Liquid_crystal) polarization grating. These types of beam splitters use [birefringence](https://en.wikipedia.org/wiki/Birefringence) to produce diffraction orders which are [circularly polarized](https://en.wikipedia.org/wiki/Circular_polarization). We will investigate two kinds of polarization gratings: (1) a homogeneous [uniaxial](https://en.wikipedia.org/wiki/Birefringence#Uniaxial_materials) grating (also known as a circular-polarization grating), and (2) a [twisted-nematic](https://en.wikipedia.org/wiki/Liquid_crystal#Chiral_phases) bilayer grating as described in [Optics Letters, Vol. 33, No. 20, pp. 2287-9, 2008](https://www.osapublishing.org/ol/abstract.cfm?uri=ol-33-20-2287) ([pdf](https://www.imagineoptix.com/cms/wp-content/uploads/2017/01/OL_08_Oh-broadband_PG.pdf)). The homogeneous uniaxial grating is just a special case of the twisted-nematic grating with a nematic [director](https://en.wikipedia.org/wiki/Liquid_crystal#Director) rotation angle of φ=0°.

A schematic of the grating geometry is shown below. The grating is a 2d slab in the *xy*-plane with two parameters: birefringence (Δn) and thickness (d). The twisted-nematic grating consists of two layers of thickness d each with equal and opposite rotation angles of φ=70° for the nematic director. Both gratings contain only three diffraction orders: m=0, ±1. The m=0 order is linearly polarized and the m=±1 orders are circularly polarized with opposite chirality. For the uniaxial grating, the diffraction efficiencies for a mode with wavelength λ can be computed analytically: η<sub>0</sub>=cos<sup>2</sup>(πΔnd/λ), η<sub>±1</sub>=0.5sin<sup>2</sup>(πΔnd/λ). The derivation of these formulas is presented in [Optics Letters, Vol, 24, No. 9, pp. 584-6, 1999](https://www.osapublishing.org/ol/abstract.cfm?uri=ol-24-9-584). We will verify these analytic results and also demonstrate that the twisted-nematic grating produces a broader bandwidth response for the ±1 orders than the homogeneous uniaxial grating. An important property of these polarization gratings for e.g. display applications is that for a circular-polarized input planewave and phase delay (Δnd/λ) of nearly 0.5, there is only a single diffraction order with *opposite* chiraity.

<center>
![](../images/polarization_grating_schematic.png)
</center>

In this example, the input is a linear-polarized planewave with center wavelength of λ=0.54 μm at normal incidence. The linear polarization is in the *yz*-plane with a rotation angle of 45° counter clockwise about the *x* axis. Two sets of mode coefficients are computed in the air region adjacent to the grating for each orthogonal polarization: `ODD_Z+EVEN_Y` and `EVEN_Z+ODD_Y`, which correspond to +k<sub>y</sub> + -k<sub>y</sub> (cosine) and +k<sub>y</sub> - -k<sub>y</sub> (sine) modes. From these linear-polarized mode coefficients, the circular-polarized mode coefficients with opposite chirality can be computed: (cosine amplitudes)+i(sine amplitudes) and (cosine amplitudes)-i(sine amplitudes), which correspond to modes with +k<sub>y</sub> and -k<sub>y</sub> separately. The transmittance for the diffraction orders are computed from the mode coefficients. As usual, this requires a separate normalization run to compute the power of the input planewave.

The anisotropic permittivity of the grating is specified using a [material function](../Python_User_Interface.md#medium). The nematic director is oriented along the *z* axis for φ=0°. Thus, the E<sub>z</sub> electric field has a larger permittivity than the E<sub>y</sub> field where the birefringence (Δn) is 0.159. The grating has a periodicity of Λ=6.5 μm in the *y* direction.

The simulation script is in [examples/polarization_grating.py](https://github.com/stevengj/meep/blob/master/python/examples/polarization_grating.py).

```py
import meep as mp
import math
import argparse

resolution = 50        # pixels/μm

dpml = 1.0             # PML thickness
dsub = 1.0             # substrate thickness
dpad = 1.0             # padding thickness

k_point = mp.Vector3(0,0,0)

n_0 = 1.55
delta_n = 0.159
epsilon_diag = mp.Matrix(mp.Vector3(n_0**2,0,0),mp.Vector3(0,n_0**2,0),mp.Vector3(0,0,(n_0+delta_n)**2))

wvl = 0.54             # center wavelength
fcen = 1/wvl           # center frequency
df = 0.05*fcen         # frequency width

def pol_grating(d,ph,gp,nmode):
    sx = dpml+dsub+d+d+dpad+dpml
    sy = gp

    cell_size = mp.Vector3(sx,sy,0)
    pml_layers = [mp.PML(thickness=dpml,direction=mp.X)]

    # twist angle of nematic director; from equation 1b
    def phi(p):
        xx  = p.x-(-0.5*sx+dpml+dsub)
        if (xx >= 0) and (xx <= d):
            return math.pi*p.y/gp + ph*xx/d
        else:
            return math.pi*p.y/gp - ph*xx/d + 2*ph

    # return the anisotropic permittivity tensor for a uniaxial, twisted nematic liquid crystal
    def lc_mat(p):
        # rotation matrix for rotation around x axis
        Rx = mp.Matrix(mp.Vector3(1,0,0),mp.Vector3(0,math.cos(phi(p)),math.sin(phi(p))),mp.Vector3(0,-math.sin(phi(p)),math.cos(phi(p))))
        lc_epsilon = Rx * epsilon_diag * Rx.transpose()
        lc_epsilon_diag = mp.Vector3(lc_epsilon[0].x,lc_epsilon[1].y,lc_epsilon[2].z)
        lc_epsilon_offdiag = mp.Vector3(lc_epsilon[1].x,lc_epsilon[2].x,lc_epsilon[2].y)
        return mp.Medium(epsilon_diag=lc_epsilon_diag,epsilon_offdiag=lc_epsilon_offdiag)

    geometry = [mp.Block(center=mp.Vector3(-0.5*sx+0.5*(dpml+dsub)),size=mp.Vector3(dpml+dsub,mp.inf,mp.inf),material=mp.Medium(index=n_0)),
                mp.Block(center=mp.Vector3(-0.5*sx+dpml+dsub+d),size=mp.Vector3(2*d,mp.inf,mp.inf),material=lc_mat)]

    # linear-polarized planewave pulse source
    src_pt = mp.Vector3(-0.5*sx+dpml+0.3*dsub,0,0)
    sources = [mp.Source(mp.GaussianSource(fcen,fwidth=df), component=mp.Ez, center=src_pt, size=mp.Vector3(0,sy,0)),
               mp.Source(mp.GaussianSource(fcen,fwidth=df), component=mp.Ey, center=src_pt, size=mp.Vector3(0,sy,0))]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=pml_layers,
                        k_point=k_point,
                        sources=sources,
                        default_material=mp.Medium(index=n_0))

    refl_pt = mp.Vector3(-0.5*sx+dpml+0.5*dsub,0,0)
    refl_flux = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=refl_pt, size=mp.Vector3(0,sy,0)))

    sim.run(until_after_sources=100)

    input_flux = mp.get_fluxes(refl_flux)
    input_flux_data = sim.get_flux_data(refl_flux)

    sim.reset_meep()

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=pml_layers,
                        k_point=k_point,
                        sources=sources,
                        geometry=geometry)

    refl_flux = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=refl_pt, size=mp.Vector3(0,sy,0)))
    sim.load_minus_flux_data(refl_flux,input_flux_data)

    tran_pt = mp.Vector3(0.5*sx-dpml-0.5*dpad,0,0)
    tran_flux = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=tran_pt, size=mp.Vector3(0,sy,0)))

    sim.run(until_after_sources=300)

    res1 = sim.get_eigenmode_coefficients(tran_flux, range(1,nmode+1), eig_parity=mp.ODD_Z+mp.EVEN_Y)
    res2 = sim.get_eigenmode_coefficients(tran_flux, range(1,nmode+1), eig_parity=mp.EVEN_Z+mp.ODD_Y)
    angles = [math.degrees(math.acos(kdom.x/fcen)) for kdom in res1.kdom]

    return input_flux[0], angles, res1.alpha[:,0,0], res2.alpha[:,0,0];


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-dd', type=float, default=1.7, help='chiral layer thickness (default: 1.7 μm)')
    parser.add_argument('-ph', type=float, default=70, help='chiral layer twist angle (default: 70°)')
    parser.add_argument('-gp', type=float, default=6.5, help='grating period (default: 6.5 μm)')
    parser.add_argument('-nmode', type=int, default=5, help='number of mode coefficients to compute (default: 5)')
    args = parser.parse_args()
    
    input_flux, angles, coeffs1, coeffs2 = pol_grating(args.dd,math.radians(args.ph),args.gp,args.nmode)

    tran1 = abs(coeffs1+1j*coeffs2)**2/input_flux
    tran2 = abs(coeffs1-1j*coeffs2)**2/input_flux
    print("tran:, 0, 0, {:.5f}".format(0.5*(tran1[0]+tran2[0])))
    for m in range(1,args.nmode):
        print("tran:, +{}, +{:.2f}, {:.5f}".format(m,angles[m],tran1[m]))
        print("tran:, -{}, -{:.2f}, {:.5f}".format(m,angles[m],tran2[m]))
```

The Bash script below runs the grating simulations over a range of grating thicknesses corresponding to a range of phase delays (Δnd/λ) of 0 to 1. The entire output is saved to a file and the transmittance data is extracted from the output and placed in a separate file.

```sh
for d in `seq 0.1 0.1 3.4`; do
    echo "circular polarization grating with d=${d}";
    dd=$(printf "%0.2f" $(echo "scale=2;0.5*${d}" |bc));
    python polarization_grating.py -dd ${dd} -ph 0 |tee -a circ_pol_grating.out;

    echo "bilayer twisted nematic polarization grating with d=${dd}";
    python polarization_grating.py -dd ${d} -ph 70 |tee -a bilayer_pol_grating.out;
done

grep tran: circ_pol_grating.out |cut -d , -f2- > circ_pol_grating.dat
grep tran: bilayer_pol_grating.out |cut -d , -f2- > bilayer_pol_grating.dat
```

The output from the simulation for the homogeneous uniaxial grating is plotted using the script below. The diffraction spectra for the two gratings are shown in the accompanying figures.

```py
import numpy as np
import math
import matplotlib.pyplot as plt

d = np.genfromtxt("circ_pol_grating.dat",delimiter=",")
m0 = d[0::9,2]
m1p = d[1::9,2]
m1m = d[2::9,2]
angles = d[1::9,1]

cos_angles = [math.cos(math.radians(t)) for t in angles]
tran = m0+m1p+m1m
eff_m0 = m0/tran
eff_m1p_m1m = ((m1p+m1m)/tran)/cos_angles

dd = np.arange(0.1,3.5,0.1)
delta_n = 0.159
wvl = 0.54
phase = delta_n*dd/wvl

eff_m0_analytic = [math.cos(math.pi*p)**2 for p in phase]
eff_m1p_m1m_analytic = [math.sin(math.pi*p)**2 for p in phase]

plt.figure(dpi=100)
plt.plot(phase,eff_m0,'bo-',clip_on=False,label='0th order (meep)')
plt.plot(phase,eff_m0_analytic,'b--',clip_on=False,label='0th order (analytic)')
plt.plot(phase,eff_m1p_m1m,'ro-',clip_on=False,label='±1 orders (meep)')
plt.plot(phase,eff_m1p_m1m_analytic,'r--',clip_on=False,label='±1 orders (analytic)')
plt.axis([0, 1.0, 0, 1])
plt.xticks([t for t in np.arange(0,1.2,0.2)])
plt.xlabel("phase delay Δnd/λ")
plt.ylabel("diffraction efficiency @ λ = 0.54 μm")
plt.legend(loc='center')
plt.title("homogeneous uniaxial grating")
plt.show()
```

<center>
![](../images/polarization_grating_diffraction_spectra.png)
</center>

The left figure shows good agreement between the simulation results and analytic theory for the homogeneous uniaxial grating. Approximately 6% of the power in the input planewave is lost due to reflection from the grating. This value is an average over the range of phase delays. The total transmittance is therefore around 94%. The twisted-nematic grating, with results shown in the right figure, produces ±1 diffraction orders with nearly-constant peak transmittance over a broader bandwidth around Δnd/λ=0.5 than the homogeneous uniaxial polarization grating. This is consistent with results from the reference. The average reflectance and transmittance for the twisted-nematic grating are similar to those for the homogeneous uniaxial grating.