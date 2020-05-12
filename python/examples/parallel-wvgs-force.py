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

a = 1.0     # waveguide width/height

k_point = mp.Vector3(z=0.5)

fcen = 0.2

def parallel_waveguide(s,xodd):
    geometry = [mp.Block(center=mp.Vector3(-0.5*(s+a)),
                         size=mp.Vector3(a,a,mp.inf),
                         material=Si),
                mp.Block(center=mp.Vector3(0.5*(s+a)),
                         size=mp.Vector3(a,a,mp.inf),
                         material=Si)]

    symmetries = [mp.Mirror(mp.X, phase=-1 if xodd else 1),
                  mp.Mirror(mp.Y, phase=-1)]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell,
                        geometry=geometry,
                        symmetries=symmetries,
                        k_point=k_point)

    sim.init_sim()
    EigenmodeData = sim.get_eigenmode(fcen,
                                      mp.Z,
                                      mp.Volume(center=mp.Vector3(), size=mp.Vector3(sx,sy)),
                                      2 if xodd else 1,
                                      k_point,
                                      match_frequency=False,
                                      parity=mp.ODD_Y)

    f = EigenmodeData.freq
    print("freq-{}:, {}, {}".format("xodd" if xodd else "xeven", s, f))

    sim.reset_meep()

    eig_sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen, fwidth=0.1*fcen),
                                      size=mp.Vector3(sx,sy),
                                      center=mp.Vector3(),
                                      eig_band=2 if xodd else 1,
                                      eig_kpoint=k_point,
                                      eig_match_freq=False,
                                      eig_parity=mp.ODD_Y)]

    sim.change_sources(eig_sources)

    flux_reg = mp.FluxRegion(direction=mp.Z, center=mp.Vector3(), size=mp.Vector3(sx,sy))
    wvg_flux = sim.add_flux(f, 0, 1, flux_reg)

    force_reg1 = mp.ForceRegion(mp.Vector3(0.5*s), direction=mp.X, weight=1, size=mp.Vector3(y=a))
    force_reg2 = mp.ForceRegion(mp.Vector3(0.5*s+a), direction=mp.X, weight=-1, size=mp.Vector3(y=a))
    wvg_force = sim.add_force(f, 0, 1, force_reg1, force_reg2)

    sim.run(until_after_sources=5000)

    flux = mp.get_fluxes(wvg_flux)[0]
    force = mp.get_forces(wvg_force)[0]
    print("flux-{}:, {}, {}".format("xodd" if xodd else "xeven", s, flux))
    print("force-{}:, {}, {}".format("xodd" if xodd else "xeven", s, force))

    sim.reset_meep()
    return flux, force


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
