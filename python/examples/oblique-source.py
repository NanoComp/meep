import meep as mp
import numpy as np
import matplotlib.pyplot as plt

resolution = 50 # pixels/μm

cell_size = mp.Vector3(14,14,0)

pml_layers = [mp.PML(thickness=2)]

# rotation angle (in degrees) of waveguide, counter clockwise (CCW) around z-axis
rot_angle = np.radians(20)

geometry = [mp.Block(center=mp.Vector3(),
                     size=mp.Vector3(mp.inf,1,mp.inf),
                     e1=mp.Vector3(1).rotate(mp.Vector3(z=1), rot_angle),
                     e2=mp.Vector3(y=1).rotate(mp.Vector3(z=1), rot_angle),
                     material=mp.Medium(epsilon=12))]

fsrc = 0.15 # frequency of eigenmode or constant-amplitude source
kx = 0.4    # initial guess for wavevector in x-direction of eigenmode
bnum = 1    # band number of eigenmode

compute_flux = True # compute flux (True) or plot the field profile (False)

eig_src = True # eigenmode (True) or constant-amplitude (False) source

if eig_src:
    sources = [mp.EigenModeSource(src=mp.GaussianSource(fsrc,fwidth=0.2*fsrc) if compute_flux else mp.ContinuousSource(fsrc),
                                  center=mp.Vector3(),
                                  size=mp.Vector3(y=14),
                                  direction=mp.AUTOMATIC if rot_angle == 0 else mp.NO_DIRECTION,
                                  eig_kpoint=mp.Vector3(kx).rotate(mp.Vector3(z=1), rot_angle),
                                  eig_band=bnum,
                                  eig_parity=mp.EVEN_Y+mp.ODD_Z if rot_angle == 0 else mp.ODD_Z,
                                  eig_match_freq=True)]
else:
    sources = [mp.Source(src=mp.GaussianSource(fsrc,fwidth=0.2*fsrc) if compute_flux else mp.ContinuousSource(fsrc),
                         center=mp.Vector3(),
                         size=mp.Vector3(y=2),
                         component=mp.Ez)]

sim = mp.Simulation(cell_size=cell_size,
                    resolution=resolution,
                    boundary_layers=pml_layers,
                    sources=sources,
                    geometry=geometry,
                    symmetries=[mp.Mirror(mp.Y)] if rot_angle == 0 else [])

if compute_flux:
    tran = sim.add_flux(fsrc, 0, 1, mp.FluxRegion(center=mp.Vector3(x=5), size=mp.Vector3(y=14)))
    sim.run(until_after_sources=50)
    print("flux:, {:.6f}".format(mp.get_fluxes(tran)[0]))
else:
    sim.run(until=100)
    nonpml_vol = mp.Volume(center=mp.Vector3(), size=mp.Vector3(10,10,0))
    eps_data = sim.get_array(vol=nonpml_vol, component=mp.Dielectric)
    ez_data = sim.get_array(vol=nonpml_vol, component=mp.Ez)
    plt.figure()
    plt.imshow(np.flipud(np.transpose(eps_data)), interpolation='spline36', cmap='binary')
    plt.imshow(np.flipud(np.transpose(ez_data)), interpolation='spline36', cmap='RdBu', alpha=0.9)
    plt.axis('off')
    plt.show()
