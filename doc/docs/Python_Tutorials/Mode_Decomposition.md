---
# Mode Decomposition
---

This tutorial demonstrates the mode-decomposition feature which is used to decompose a given mode profile into a superposition of harmonic basis modes. There are examples for two kinds of modes in lossless, dielectric media: (1) localized (i.e., guided) and (2) non-localized (i.e., planewave).

[TOC]

Reflectance of a Waveguide Taper
--------------------------------

This example involves computing the reflectance &mdash; the fraction of the reflected power over the incident power &mdash; of the fundamental mode of a linear waveguide taper. The structure and the simulation parameters are shown in the schematic below. We will verify that computing the reflectance using two differerent methods produces nearly identical results: (1) mode decomposition and (2) directly from the fields via its [Poynting flux](../Introduction/#transmittancereflectance-spectra). Also, we will demonstrate that the scaling of the reflectance with the taper length is quadratic, consistent with analytical results from [Optics Express, Vol. 16, pp. 11376-92, 2008](http://www.opticsinfobase.org/abstract.cfm?URI=oe-16-15-11376).

<center>
![](../images/waveguide-taper.png)
</center>

The structure, which can be viewed as a [two-port network](https://en.wikipedia.org/wiki/Two-port_network), consists of a single-mode waveguide of width `w1` (1 μm)  coupled to a second waveguide of width `w2` (2 μm) via a linearly-sloped taper of length `Lt`. The structure is homogeneous ε=12 in vacuum. PML absorbing boundaries surround the computational cell. An eigenmode source with E<sub>z</sub> polarization is used to launch the fundamental mode with a wavelength of 6.67 μm. There is an eigenmode-expansion monitor placed at the midpoint of the first waveguide. This is a line monitor which extends beyond the waveguide in order to span the entire mode profile including its evanescent tails. The Fourier-transformed fields along this line monitor are used to compute the basis coefficients of the harmonic modes. These are computed separately via the eigenmode solver [MPB](https://mpb.readthedocs.io). Technical details are described in [Mode Decomposition](../Mode_Decomposition) where it is shown that the squared magnitude of the mode coefficient is equivalent to the power (i.e., Poynting flux) in the given eigenmode. The ratio of the complex mode coefficients can be used to compute the [S parameters](https://en.wikipedia.org/wiki/Scattering_parameters). In this example, we are computing |S<sub>11</sub>|<sup>2</sup> which is the reflectance. Another line monitor could have been placed in the second waveguide to compute the transmittance or |S<sub>21</sub>|<sup>2</sup>. Also not considered in this example is the scattering into radiative modes.

The structure has mirror symmetry in the $y$ direction which can be exploited to reduce the computation size by a factor of two. This requires that we use `add_flux` and specify `eig_parity=mp.ODD_Z+mp.EVEN_Y` in the call to `get_eigenmode_coefficients`. This is because `add_mode_monitor`, which is an alias for `add_flux`, is not optimized for symmetry.

At the end of the simulation, the reflectance of the fundamental mode computed using the two methods are displayed. The simulation script is shown below and in [mode-decomposition.py](https://github.com/stevengj/meep/blob/master/python/examples/mode-decomposition.py).

```py
import meep as mp

resolution = 60   # pixels/μm

w1 = 1            # width of waveguide 1
w2 = 2            # width of waveguide 2
Lw = 10           # length of waveguide 1/2

dair = 3.0        # length of air region
dpml = 5.0        # length of PML

sy = dpml+dair+w2+dair+dpml

half_w1 = 0.5*w1
half_w2 = 0.5*w2

Si = mp.Medium(epsilon=12.0)

boundary_layers = [mp.PML(dpml)]

# mode wavelength
lcen = 6.67
# mode frequency
fcen = 1/lcen

symmetries = [mp.Mirror(mp.Y)]

for m in range(5):    
    Lt = 2**m
    sx = dpml+Lw+Lt+Lw+dpml
    cell_size = mp.Vector3(sx,sy,0)

    prism_x = sx+1
    half_Lt = 0.5*Lt

    src_pt = mp.Vector3(-0.5*sx+dpml+0.2*Lw,0,0)
    sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen,fwidth=0.2*fcen),
                                  component=mp.Ez,
                                  center=src_pt,
                                  size=mp.Vector3(0,sy-2*dpml,0),
                                  eig_match_freq=True,
                                  eig_parity=mp.ODD_Z+mp.EVEN_Y)]
    
    # straight waveguide
    vertices = [mp.Vector3(-prism_x,half_w1),
                mp.Vector3(prism_x,half_w1),
                mp.Vector3(prism_x,-half_w1),
                mp.Vector3(-prism_x,-half_w1)]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=boundary_layers,
                        geometry=[mp.Prism(vertices,height=mp.inf,material=Si)],
                        sources=sources,
                        symmetries=symmetries)

    mon_pt = mp.Vector3(-0.5*sx+dpml+0.5*Lw,0,0)
    flux = sim.add_flux(fcen,0,1,mp.FluxRegion(center=mon_pt,size=mp.Vector3(0,sy-2*dpml,0)))

    sim.run(until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,src_pt,1e-9))

    incident_coeffs, vgrp, kpoints = sim.get_eigenmode_coefficients(flux,[1],eig_parity=mp.ODD_Z+mp.EVEN_Y)
    incident_flux = mp.get_fluxes(flux)
    incident_flux_data = sim.get_flux_data(flux)

    sim.reset_meep()    
    
    # linear taper
    vertices = [mp.Vector3(-prism_x,half_w1),
                mp.Vector3(-half_Lt,half_w1),
                mp.Vector3(half_Lt,half_w2),
                mp.Vector3(prism_x,half_w2),
                mp.Vector3(prism_x,-half_w2),
                mp.Vector3(half_Lt,-half_w2),
                mp.Vector3(-half_Lt,-half_w1),
                mp.Vector3(-prism_x,-half_w1)]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=boundary_layers,
                        geometry=[mp.Prism(vertices,height=mp.inf,material=Si)],
                        sources=sources,
                        symmetries=symmetries)

    refl_flux = sim.add_flux(fcen,0,1,mp.FluxRegion(center=mon_pt,size=mp.Vector3(0,sy-2*dpml,0)))
    sim.load_minus_flux_data(refl_flux,incident_flux_data)
    
    sim.run(until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,src_pt,1e-9))
        
    coeffs, vgrp, kpoints = sim.get_eigenmode_coefficients(refl_flux,[1],eig_parity=mp.ODD_Z+mp.EVEN_Y)
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

The mode-decomposition feature can also be applied to planewaves in homogeneous media with scalar permittivity/permeability (i.e., no anisotropy). This will be demonstrated in this example to compute the diffraction spectrum of a binary phase [grating](https://en.wikipedia.org/wiki/Diffraction_grating). The unit cell geometry of the grating is shown in the schematic below. The grating is periodic in the $y$ direction with periodicity `gp` and has a rectangular profile of height `gh` and duty cycle `gdc`. The grating parameters are `gh`=0.5 μm, `gdc`=0.5, and `gp`=10 μm. There is a semi-infinite substrate of thickness `dsub` adjacent to the grating. The substrate and grating are glass with a constant refractive index of 1.5. The surrounding is air/vacuum. Perfectly matched layers (PML) of thickness `dpml` are used in the $\pm x$ boundaries. A pulsed planewave with E<sub>z</sub> polarization spanning wavelengths of 0.4 to 0.6 μm is normally incident on the grating from the glass substrate. The eigenmode monitor is placed in the air region. We will use mode decomposition to compute the transmittance &mdash; the ratio of the power in the $+x$ direction of the diffracted mode relative to that of the incident planewave &mdash; for the first ten diffraction orders. Two simulations are required: (1) an empty cell of homogeneous glass to obtain the incident power of the source, and (2) the grating structure to obtain the diffraction orders. At the end of the simulation, the wavelength, angle, and transmittance for each diffraction order are computed.

The simulation script is shown below and in [binary_grating.py](https://github.com/stevengj/meep/blob/master/python/examples/binary_grating.py).

<center>
![](../images/grating.png)
</center>

```py
import meep as mp
import math

resolution = 60        # pixels/μm

dsub = 3.0             # substrate thickness
dpad = 3.0             # padding between grating and pml
gp = 10.0              # grating period
gh = 0.5               # grating height
gdc = 0.5              # grating duty cycle

dpml = 1.0             # PML thickness
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
sources = [mp.Source(mp.GaussianSource(fcen,fwidth=df), component=mp.Ez, center=src_pt, size=mp.Vector3(0,sy,0))]

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

sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(xm,0,0), 1e-9))

freqs = mp.get_eigenmode_freqs(mode_mon)

nmode = 10
for nm in range(nmode):
  coeffs, vgrps, kpoints = sim.get_eigenmode_coefficients(mode_mon, [nm+1], eig_parity=mp.ODD_Z+mp.EVEN_Y)
  for nf in range(nfreq):
      mode_wvl = 1/freqs[nf]
      mode_angle = math.degrees(math.acos(kpoints[nf].x/freqs[nf]))
      mode_tran = abs(coeffs[0,nf,0])**2/input_flux[nf]
      if nm != 0:
      	 mode_tran = 0.5*mode_tran
      print("grating{}:, {:.5f}, {:.2f}, {:.8f}".format(nm,mode_wvl,mode_angle,mode_tran))
```

Note the use of the keyword parameter argument `eig_parity=mp.ODD_Z+mp.EVEN_Y` in the call to `get_eigenmode_coefficients`. This is important for specifying **non-degenerate** modes in MPB since the `k_point` is (0,0,0). `ODD_Z` is for modes with E<sub>z</sub> polarzation. `EVEN_Y` is necessary since each diffraction order m (an integer) with |m|>0 and a given k<sub>x</sub> consists of two modes: one going in the +y direction and the other in the -y direction. `EVEN_Y` forces MPB to compute only the +k<sub>y</sub> + -k<sub>y</sub> (cosine) mode. Additionally, the total transmittance must be halved in this case to obtain the transmittance for either the +k<sub>y</sub> or -k<sub>y</sub> mode. For `ODD_Y`, MPB will compute the sine mode but this will have zero power because the source is even. If the $y$ parity is left out, MPB will return a random superposition of the cosine and sine modes. Specifying the `eig_parity` parameter this way ensures that the ordering of the modes corresponds to only the non-degenerate diffraction orders. Finally, note the use of `add_flux` instead of `add_mode_monitor` when using symmetries.

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

Each diffraction order corresponds to a single angle. In the figure below, this angle is represented by the *lower* boundary of each labeled region. For example, the m=0 order has a diffraction angle of 0 degrees at all wavelengths. The representation of the diffraction orders as finite angular regions is an artifact of matplotlib's [pcolormesh](https://matplotlib.org/api/_as_gen/matplotlib.pyplot.pcolormesh.html) routine. Note that only the positive diffraction orders are shown as these are equivalent to the negative orders due to the symmetry of the source and the structure.

The transmittance of each diffraction order should ideally be a constant. The slight wavelength dependence shown in the figure is a numerical artifact which can be mitigated by (1) increasing the resolution or (2) time-stepping for a longer duration to ensure that the fields have sufficiently decayed away.

The diffraction orders/modes are a finite set of propagating planewaves. The wavevector k<sub>x</sub> of these modes can be computed analytically: for a frequency of ω (in c=1 units), these propagating modes are the **real** solutions of sqrt(ω²n²-(k<sub>y</sub>+2πm/Λ)²) where m is the diffraction order (an integer) and Λ is the periodicity of the grating. In this example, n=1, k<sub>y</sub>=0, and Λ=10 μm. Thus, as an example, at a wavelength of 0.5 μm there are 20 diffraction orders (though we only computed the first ten). The wavevector k<sub>x</sub> is used to compute the angle of the diffraction order as cos<sup>-1</sup>(k<sub>x</sub>/ω). Evanescent modes, those with an imaginary k<sub>x</sub>, exist for |m|>20 but these modes carry no power. Note that currently Meep does not compute the number of propagating modes for you. If the mode number passed to `get_eigenmode_coefficients` is larger than the number of propagating modes at a given frequency/wavelength, MPB's Newton solver will fail to converge and will return zero for the mode coefficient. It is therefore a good idea to know beforehand the number of propagating modes.

<center>
![](../images/grating_diffraction_spectra.png)
</center>

In the limit where the grating periodicity is much larger than the wavelength and the size of the diffracting element (i.e., more than 10 times), as it is in this example, the [diffraction efficiency](https://en.wikipedia.org/wiki/Diffraction_efficiency) can be computed analytically using scalar theory. This is described in the OpenCourseWare [Optics course](https://ocw.mit.edu/courses/mechanical-engineering/2-71-optics-spring-2009/) in the Lecture 16 (Gratings: Amplitude and Phase, Sinusoidal and Binary) [notes](https://ocw.mit.edu/courses/mechanical-engineering/2-71-optics-spring-2009/video-lectures/lecture-16-gratings-amplitude-and-phase-sinusoidal-and-binary/MIT2_71S09_lec16.pdf) and [video](https://www.youtube.com/watch?v=JmWguqCZRxk). For a review of scalar diffraction theory, see Chapter 3 ("Analysis of Two-Dimensional Signals and Systems") of [Introduction to Fourier Optics (fourth edition)](https://www.amazon.com/Introduction-Fourier-Optics-Joseph-Goodman-ebook/dp/B076TBP48F) by J.W. Goodman. From the scalar theory, the diffraction efficiency of the binary grating is (2/(mπ))<sup>2</sup> when the phase difference between the propagating distance in the glass relative to the same distance in air is π. The phase differerence/contrast is (2π/λ)(n-1)s where λ is the wavelength, n is the refractive index of the grating, and s is the propagation distance in the grating (`gh` in the simulation script). A special feature of the binary grating is that the diffraction efficiency is 0 for all *even* orders. This is verified by the diffraction spectrum shown above.

To convert the diffraction efficiency into transmittance in the *x* direction (in order to be able to compare the analytic results with those from Meep), the diffraction efficiency must be multiplied by the Fresnel transmittance from air to glass and by the cosine of the diffraction angle. We compare the analytic and simulated results at a wavelength of 0.5 μm (for which the scalar theory is valid) for diffraction orders 1, 3, 5, and 7. The analytic results are 0.3886, 0.0427, 0.0151, and 0.0074. The Meep results are 0.3942, 0.04371, 0.0154, and 0.0077. This corresponds to relative errors of approximately 1.4%, 2.3%, 2.0%, and 3.1% which indicates good agreement.
