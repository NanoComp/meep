import meep as mp
import numpy as np
import math
import matplotlib.pyplot as plt

resolution = 50         # pixels/μm

dpml = 1.0              # PML thickness
dsub = 3.0              # substrate thickness
dpad = 3.0              # padding between grating and PML
gp = 1.0                # grating periodicity
gh = 0.5                # grating height
gdc = 0.5               # grating duty cycle

num_cells = 20
gdc_list = [gdc for _ in range(num_cells)]

wvl = 0.5               # center wavelength
fcen = 1/wvl            # center frequency

k_point = mp.Vector3()

glass = mp.Medium(index=1.5)

pml_layers = [mp.PML(thickness=dpml)]

symmetries=[mp.Mirror(mp.Y)]

sx = dpml+dsub+gh+dpad+dpml
sy = dpml+dpad+num_cells*gp+dpad+dpml
cell_size = mp.Vector3(sx,sy)

src_pt = mp.Vector3(-0.5*sx+dpml+0.5*dsub)
sources = [mp.Source(mp.GaussianSource(fcen,fwidth=0.2*fcen),
                     component=mp.Ez,
                     center=src_pt,
                     size=mp.Vector3(y=sy-2*dpml))]

geometry = [mp.Block(material=glass,
                     size=mp.Vector3(dpml+dsub,mp.inf,mp.inf),
                     center=mp.Vector3(-0.5*sx+0.5*(dpml+dsub)))]

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    k_point=k_point,
                    sources=sources,
                    symmetries=symmetries)

mon_pt = mp.Vector3(0.5*sx-dpml-0.5*dpad)
near_fields = sim.add_dft_fields([mp.Ez], fcen, fcen, 1, center=mon_pt, size=mp.Vector3(y=sy-2*dpml))

sim.run(until_after_sources=100)

flat_dft = sim.get_dft_array(near_fields, mp.Ez, 0)

sim.reset_meep()

for j in range(num_cells):
  geometry.append(mp.Block(material=glass,
                           size=mp.Vector3(gh,gdc_list[j]*gp,mp.inf),
                           center=mp.Vector3(-0.5*sx+dpml+dsub+0.5*gh,-0.5*sy+dpml+dpad+(j+0.5)*gp)))

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    k_point=k_point,
                    sources=sources,
                    symmetries=symmetries)

near_fields = sim.add_dft_fields([mp.Ez], fcen, fcen, 1, center=mon_pt, size=mp.Vector3(y=sy-2*dpml))

sim.run(until_after_sources=100)

grating_dft = sim.get_dft_array(near_fields, mp.Ez, 0)

scattered_field = grating_dft-flat_dft
FT_scattered_field = np.fft.fftshift(np.fft.fft(scattered_field))

[x,y,z,w] = sim.get_array_metadata(dft_cell=near_fields)
ky = np.fft.fftshift(np.fft.fftfreq(len(scattered_field), 1/resolution))

if mp.am_master():
  plt.subplot(2,1,1)
  plt.title("finite grating with {} periods".format(num_cells))
  plt.plot(y,np.abs(scattered_field)**2,'bo-')
  plt.gca().get_yaxis().set_ticks([])
  plt.xlabel("y (μm)")
  plt.ylabel("scattered field\namplitude (a.u.)")

  plt.subplot(2,1,2)
  plt.plot(ky,np.abs(FT_scattered_field)**2,'ro-')
  plt.gca().get_yaxis().set_ticks([])
  plt.xlabel(r'wavevector k$_y$, 2π (μm)$^{-1}$')
  plt.ylabel("diffraction spectra of\nscattered field (a.u.)")
  plt.gca().set_xlim([-3, 3])

  plt.subplots_adjust(hspace=0.3)
  plt.tight_layout(pad=1.0)
  plt.show()
