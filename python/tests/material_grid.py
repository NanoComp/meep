import meep as mp
import numpy as np
import matplotlib.pyplot as plt

resolution = 50 # pixels/μm

cell_size = mp.Vector3(14,14)

pml_layers = [mp.PML(thickness=2)]

# rotation angle (in degrees) of waveguide, counter clockwise (CCW) around z-axis
rot_angle = np.radians(20)

w = 1.0 # width of waveguide

m1 = mp.Medium(epsilon=2)
m2 = mp.Medium(epsilon=12)
n = 100
gs = mp.Vector3(n,n)
np.random.seed(1)
dp = np.random.rand(n*n)
mg = mp.MaterialGrid(gs,m1,m2,dp)

geometry = [mp.Block(center=mp.Vector3(),
                     size=mp.Vector3(2,2,mp.inf),
                     e1=mp.Vector3(x=1).rotate(mp.Vector3(z=1), rot_angle),
                     e2=mp.Vector3(y=1).rotate(mp.Vector3(z=1), rot_angle),
                     material=mg)]

fsrc = 0.15 # frequency of eigenmode or constant-amplitude source
bnum = 1    # band number of eigenmode

kpoint = mp.Vector3(x=1).rotate(mp.Vector3(z=1), rot_angle)

compute_flux = True # compute flux (True) or plot the field profile (False)

eig_src = True # eigenmode (True) or constant-amplitude (False) source

if eig_src:
    sources = [mp.EigenModeSource(src=mp.GaussianSource(fsrc,fwidth=0.2*fsrc) if compute_flux else mp.ContinuousSource(fsrc),
                                  center=mp.Vector3(),
                                  size=mp.Vector3(y=3*w),
                                  direction=mp.NO_DIRECTION,
                                  eig_kpoint=kpoint,
                                  eig_band=bnum,
                                  eig_parity=mp.EVEN_Y+mp.ODD_Z if rot_angle == 0 else mp.ODD_Z,
                                  eig_match_freq=True)]
else:
    sources = [mp.Source(src=mp.GaussianSource(fsrc,fwidth=0.2*fsrc) if compute_flux else mp.ContinuousSource(fsrc),
                         center=mp.Vector3(),
                         size=mp.Vector3(y=3*w),
                         component=mp.Ez)]

sim = mp.Simulation(cell_size=cell_size,
                    resolution=resolution,
                    boundary_layers=pml_layers,
                    sources=sources,
                    #eps_averaging=False,
                    geometry=geometry
                    #symmetries=[mp.Mirror(mp.Y)] if rot_angle == 0 else []
                    )

sim.plot2D()
plt.show()
sim.reset_meep()
mg.update_parameters(dp*0)
sim.init_sim()
sim.plot2D()
plt.show()