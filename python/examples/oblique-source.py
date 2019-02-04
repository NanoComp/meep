import meep as mp
import numpy as np

resolution = 50 # pixels/μm

cell_size = mp.Vector3(14,14,0)
    
pml_layers = [mp.PML(thickness=2)]

rot_angle = np.radians(0)  # rotation angle (in degrees) of waveguide, CCW around z-axis

geometry = [mp.Block(center=mp.Vector3(),
                     size=mp.Vector3(mp.inf,1,mp.inf),
                     e1=mp.Vector3(1,0,0).rotate(mp.Vector3(0,0,1), rot_angle),
                     e2=mp.Vector3(0,1,0).rotate(mp.Vector3(0,0,1), rot_angle),
                     material=mp.Medium(epsilon=12))]

fsrc = 0.15 # frequency of eigenmode or continuous-wave (CW) source
kx = 0.4    # initial guess for wavevector in x-direction of eigenmode
bnum = 1    # band index of eigenmode

eig_src = True

if eig_src:
    sources = [mp.EigenModeSource(src=mp.ContinuousSource(fsrc),
                                  center=mp.Vector3(),
                                  size=mp.Vector3(y=14),
                                  direction=mp.AUTOMATIC if rot_angle == 0 else mp.NO_DIRECTION,
                                  eig_kpoint=mp.Vector3(kx,0,0).rotate(mp.Vector3(0,0,1), rot_angle),
                                  eig_band=bnum,
                                  eig_parity=mp.EVEN_Y+mp.ODD_Z if rot_angle == 0 else mp.ODD_Z,
                                  eig_match_freq=True)]
else:
    sources = [mp.Source(src=mp.ContinuousSource(fsrc),
                         center=mp.Vector3(),
                         size=mp.Vector3(y=2),
                         component=mp.Ez)]

sim = mp.Simulation(cell_size=cell_size,
                    resolution=resolution,
                    boundary_layers=pml_layers,
                    sources=sources,
                    geometry=geometry,
                    symmetries=[mp.Mirror(mp.Y)] if rot_angle == 0 else [])

sim.run(until=100)

nonpml_vol = mp.Volume(center=mp.Vector3(), size=mp.Vector3(10,10,0))
eps_data = sim.get_array(vol=nonpml_vol, component=mp.Dielectric)
ez_data = sim.get_array(vol=nonpml_vol, component=mp.Ez)

import matplotlib.pyplot as plt

plt.figure()
plt.imshow(np.flipud(np.transpose(eps_data)), interpolation='spline36', cmap='binary')
plt.imshow(np.flipud(np.transpose(ez_data)), interpolation='spline36', cmap='RdBu', alpha=0.9)
plt.axis('off')
plt.show()
