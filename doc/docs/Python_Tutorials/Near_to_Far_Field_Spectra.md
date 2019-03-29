---
# Near to Far Field Spectra
---

We demonstrate Meep's [near-to-far field transformation](../Python_User_Interface.md#near-to-far-field-spectra) feature using two examples. There are three steps involved in this type of calculation. First, we need to define the "near" surface(s) as a set of surfaces capturing *all* outgoing radiation in the desired direction(s). Second, we run the simulation using a pulsed source (or alternatively, a CW source via the [frequency-domain solver](../Python_User_Interface.md#frequency-domain-solver)) to allow Meep to accumulate the Fourier transforms on the near surface(s). Third, we have Meep compute the far fields at any desired points with the option to save the far fields to an HDF5 file.

[TOC]

### Radiation Pattern of an Antenna

In this example, we compute the [radiation pattern](https://en.wikipedia.org/wiki/Radiation_pattern) of an antenna. This involves an electric-current point dipole source as the emitter in vacuum. The source is placed at the center of a 2d square cell surrounded by PML. The near fields are obtained on a bounding box defined along the edges of the non-PML region. The far fields are computed in two ways: along the (1) sides of a square box and (2) circumference of a circle, having a length/radius many times larger than the source wavelength and lying beyond the cell. From both the near and far fields, we will also compute the total outgoing Poynting flux and demonstrate that they are equivalent. Results will be shown for three orthogonal polarizations of the input source.

The simulation geometry is shown in the following schematic.

<center>
![](../images/Near2far_simulation_geometry.png)
</center>

In the first part of the simulation, we define the cell and sources as well as the near field and flux regions. Since we are using a pulsed source (with center wavelength of 1 μm), the fields are timestepped until they have sufficiently decayed away.

The simulation script is in [examples/antenna-radiation.py](https://github.com/NanoComp/meep/blob/master/python/examples/antenna-radiation.py). The notebook is [examples/antenna-radiation.ipynb](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/antenna-radiation.ipynb)

```py
import meep as mp
import math
import numpy as np
import matplotlib.pyplot as plt

resolution = 50  # pixels/um

sxy = 4
dpml = 1
cell = mp.Vector3(sxy+2*dpml,sxy+2*dpml,0)

pml_layers = [mp.PML(dpml)]

fcen = 1.0
df = 0.4
src_cmpt = mp.Ez
sources = [mp.Source(src=mp.GaussianSource(fcen,fwidth=df),
                     center=mp.Vector3(),
                     component=src_cmpt)]

if src_cmpt == mp.Ex:
    symmetries = [mp.Mirror(mp.X,phase=-1),
                  mp.Mirror(mp.Y,phase=+1)]
elif src_cmpt == mp.Ey:
    symmetries = [mp.Mirror(mp.X,phase=+1),
                  mp.Mirror(mp.Y,phase=-1)]
elif src_cmpt == mp.Ez:
    symmetries = [mp.Mirror(mp.X,phase=+1),
                  mp.Mirror(mp.Y,phase=+1)]

sim = mp.Simulation(cell_size=cell,
                    resolution=resolution,
                    sources=sources,
                    symmetries=symmetries,
                    boundary_layers=pml_layers)

nearfield_box = sim.add_near2far(fcen, 0, 1,
                                 mp.Near2FarRegion(mp.Vector3(y=0.5*sxy), size=mp.Vector3(sxy)),
                                 mp.Near2FarRegion(mp.Vector3(y=-0.5*sxy), size=mp.Vector3(sxy), weight=-1),
                                 mp.Near2FarRegion(mp.Vector3(0.5*sxy), size=mp.Vector3(y=sxy)),
                                 mp.Near2FarRegion(mp.Vector3(-0.5*sxy), size=mp.Vector3(y=sxy), weight=-1))

flux_box = sim.add_flux(fcen, 0, 1,
                        mp.FluxRegion(mp.Vector3(y=0.5*sxy), size=mp.Vector3(sxy)),
                        mp.FluxRegion(mp.Vector3(y=-0.5*sxy), size=mp.Vector3(sxy), weight=-1),
                        mp.FluxRegion(mp.Vector3(0.5*sxy), size=mp.Vector3(y=sxy)),
                        mp.FluxRegion(mp.Vector3(-0.5*sxy), size=mp.Vector3(y=sxy), weight=-1))

sim.run(until_after_sources=mp.stop_when_fields_decayed(50, src_cmpt, mp.Vector3(), 1e-8))
```

After the time stepping, the flux of the near fields is computed using `get_fluxes`:

```py
near_flux = mp.get_fluxes(flux_box)[0]
```

In the first of two cases, the flux of the far fields is computed using the `flux` routine for a square box of side length 2 mm which is 2000 times larger than the source wavelength. This requires computing the outgoing flux on each of the four sides of the box separately and summing the values. The resolution of the far fields is chosen arbitrarily as 1 point/μm. This means there are 2x10<sup>6</sup> points per side length.

```py
r = 1000/fcen      # half side length of far-field square box OR radius of far-field circle
res_ff = 1         # resolution of far fields (points/μm)
far_flux_box = (nearfield_box.flux(mp.Y, mp.Volume(center=mp.Vector3(y=r), size=mp.Vector3(2*r)), res_ff)[0]
               - nearfield_box.flux(mp.Y, mp.Volume(center=mp.Vector3(y=-r), size=mp.Vector3(2*r)), res_ff)[0]
               + nearfield_box.flux(mp.X, mp.Volume(center=mp.Vector3(r), size=mp.Vector3(y=2*r)), res_ff)[0]
               - nearfield_box.flux(mp.X, mp.Volume(center=mp.Vector3(-r), size=mp.Vector3(y=2*r)), res_ff)[0])
```

For the second of two cases, we use the `get_farfield` routine to compute the far fields by looping over a set of 100 equally-spaced points along the circumference of a circle with radius of 1 mm. The six far field components (E$_x$, E$_y$, E$_z$, H$_x$, H$_y$, H$_z$) are stored as separate arrays of complex numbers. From the far fields at each point $\mathbf{r}$, we compute the outgoing or radial flux: $\sqrt{P_x^2+P_y^2}$, where P$_x$ and P$_y$ are the components of the Poynting vector $\mathbf{P}(\mathbf{r})=(P_x,P_y,P_z)=\mathrm{Re}\, \mathbf{E}(\mathbf{r})^*\times\mathbf{H}(\mathbf{r})$. Note that $P_z$ is always 0 since this is a 2d simulation. The total flux is computed and the three flux values are displayed.

```py
npts = 100         # number of points in [0,2*pi) range of angles
angles = 2*math.pi/npts*np.arange(npts)

E = np.zeros((npts,3),dtype=np.complex128)
H = np.zeros((npts,3),dtype=np.complex128)
for n in range(npts):
    ff = sim.get_farfield(nearfield_box,
                          mp.Vector3(r*math.cos(angles[n]),
                                     r*math.sin(angles[n])))
    E[n,:] = [np.conj(ff[j]) for j in range(3)]
    H[n,:] = [ff[j+3] for j in range(3)]

Px = np.real(np.multiply(E[:,1],H[:,2])-np.multiply(E[:,2],H[:,1]))
Py = np.real(np.multiply(E[:,2],H[:,0])-np.multiply(E[:,0],H[:,2]))
Pr = np.sqrt(np.square(Px)+np.square(Py))

far_flux_circle = np.sum(Pr)*2*np.pi*r/len(Pr)

print("flux:, {:.6f}, {:.6f}, {:.6f}".format(near_flux,far_flux_box,far_flux_circle))
```

By [Poynting's theorem](https://en.wikipedia.org/wiki/Poynting%27s_theorem), the total outgoing flux obtained by integrating around a *closed* surface should be the same whether it is calculated from the near or far fields (unless there are sources or absorbers in between). The flux of the near fields for the J$_z$ source is `2.456196` and that for the far fields is `2.458030` (box) and `2.457249` (circle). The ratio of near- to far-field (circle) flux is `0.999571`. Similarly, for the J$_x$ source, the values are `1.227786` (near-field), `1.227651` (far-field box), and `1.227260` (far-field circle). The ratio of near- to far-field (circle) flux is `1.000429`. The slight differences in the flux values are due to discretization effects and will decrease as the resolution is increased.

Finally, we plot the radial flux normalized by its maximum value over the entire interval to obtain a range of values between 0 and 1. These are shown below in the linearly-scaled, polar-coordinate plots. The three figures are obtained using separate runs involving a `src_cmpt` of E$_x$, E$_y$, and E$_z$. As expected, the J$_x$ and J$_y$ sources produce [dipole](https://en.wikipedia.org/wiki/Electric_dipole_moment) radiation patterns while J$_z$ has a monopole pattern.

```py
ax = plt.subplot(111, projection='polar')
ax.plot(angles,Pr/max(Pr),'b-')
ax.set_rmax(1)
ax.set_rticks([0,0.5,1])
ax.grid(True)
ax.set_rlabel_position(22)
plt.show()
```

<center>
![](../images/Source_radiation_pattern.png)
</center>

### Far-Field Intensity of a Cavity

For this demonstration, we will compute the far-field spectra of a resonant cavity mode in a holey waveguide; a structure we had explored in [Tutorial/Resonant Modes and Transmission in a Waveguide Cavity](Resonant_Modes_and_Transmission_in_a_Waveguide_Cavity.md). The script is in [examples/cavity-farfield.py](https://github.com/NanoComp/meep/blob/master/python/examples/cavity-farfield.py). The notebook is [examples/cavity-farfield.ipynb](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/cavity-farfield.ipynb). The structure is shown at the bottom of the left image below.

![center|Schematic of the computational cell for a holey waveguide with cavity showing the location of the "near" boundary surface and the far-field region.](../images/N2ff_comp_cell.png)

To set this up, we simply remove the last portion of [examples/holey-wvg-cavity.py](https://github.com/NanoComp/meep/blob/master/python/examples/holey-wvg-cavity.py), beginning right after the line:

```py
sim.symmetries.append(mp.Mirror(mp.Y, phase=-1))
sim.symmetries.append(mp.Mirror(mp.X, phase=-1))
```

and insert the following lines:

```py
d1 = 0.2

sim = mp.Simulation(cell_size=cell,
                    geometry=geometry,
                    sources=[sources],
                    symmetries=symmetries,
                    boundary_layers=[pml_layers],
                    resolution=resolution)

nearfield = sim.add_near2far(
    fcen, 0, 1,
    mp.Near2FarRegion(mp.Vector3(0, 0.5 * w + d1), size=mp.Vector3(2 * dpml - sx)),
    mp.Near2FarRegion(mp.Vector3(-0.5 * sx + dpml, 0.5 * w + 0.5 * d1), size=mp.Vector3(0, d1), weight=-1.0),
    mp.Near2FarRegion(mp.Vector3(0.5 * sx - dpml, 0.5 * w + 0.5 * d1), size=mp.Vector3(0, d1))
)
```

We are creating a "near" bounding surface, consisting of three separate regions surrounding the cavity, that captures <i>all</i> outgoing waves in the top-half of the cell. Note that the *x*-normal surface on the left has a `weight` of -1 corresponding to the direction of the *outward normal* vector relative to the *x* direction so that the far-field spectra is correctly computed from the outgoing fields, similar to the flux and force features. The parameter `d1` is the distance between the edge of the waveguide and the bounding surface, as shown in the schematic above, and we will demonstrate that changing this parameter does not change the far-field spectra which we compute at a single frequency corresponding to the cavity mode.

We then time step the fields until, at a random point, they have sufficiently decayed away as the cell is surrounded by PMLs, and output the far-field spectra over a rectangular area that lies <i>outside</i> of the cell:

```py
sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Hz, mp.Vector3(0.12, -0.37), 1e-8))

d2 = 20
h = 4

sim.output_farfields(nearfield, "spectra-{}-{}-{}".format(d1, d2, h),
                     mp.Volume(mp.Vector3(0, (0.5 * w) + d2 + (0.5 * h)), size=mp.Vector3(sx - 2 * dpml, h)),
                     resolution)
```

The first item to note is that the far-field region is located <i>outside</i> of the cell, although in principle it can be located anywhere. The second is that the far-field spectra can be interpolated onto a spatial grid that has any given resolution but in this example we used the same resolution as the simulation. Note that the simulation itself used purely real fields but the output, given its analytical nature, contains complex fields. Finally, given that the far-field spectra is derived from the Fourier-transformed fields which includes an arbitrary constant factor, we should expect an overall scale and phase difference in the results obtained using the near-to-far-field feature with those from a corresponding simulation involving the full computational volume. The key point is that the results will be qualitatively but not quantitatively identical. The data will be written out to an HDF5 file having a filename prefix with the values of the three main parameters. This file will includes the far-field spectra for all six field components, including real and imaginary parts.

We run the above modified control file and in post-processing create an image of the real and imaginary parts of H$_z$ over the far-field region which is shown in insets (a) above. For comparison, we compute the steady-state fields using a larger cell that contains within it the far-field region. This involves a continuous source and complex fields. Results are shown in figure (b) above. The difference in the relative phases among any two points within each of the two field spectra is zero, which can be confirmed numerically. Also, as would be expected, it can be shown that increasing `d1` does not change the far-field spectra as long as the results are sufficiently converged. This indicates that discretization effects are irrelevant.

In general, it is tricky to interpret the overall scale and phase of the far fields, because it is related to the scaling of the Fourier transforms of the near fields. It is simplest to use the `near2far` feature in situations where the overall scaling is irrelevant, e.g. when you are computing a ratio of fields in two simulations, or a fraction of the far field in some region, etcetera.