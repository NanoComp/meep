import numpy as np
import meep as mp

resolution = 50 # pixels/um

sx = 1.5
sy = 6.0

cell_size = mp.Vector3(sx,sy)

dpml = 1.0
pml_layers = [mp.PML(direction=mp.Y,
                     thickness=dpml)]

fcen = 1.0

def src_amp_func(m):
    def _src_amp_func(p):
        return np.cos(2*np.pi*m*(p.x+0.5*sx)/sx)
    return _src_amp_func

sources = [mp.Source(mp.GaussianSource(fcen,fwidth=0.2*fcen),
                     component=mp.Ez,
                     center=mp.Vector3(0,-0.5*sy+dpml),
                     size=mp.Vector3(sx,0),
                     amp_func=src_amp_func(6))]

sim = mp.Simulation(cell_size=cell_size,
                    resolution=resolution,
                    k_point=mp.Vector3(),
                    boundary_layers=pml_layers,
                    sources=sources,
                    default_material=mp.Medium(index=3.5))

flux_mon = sim.add_flux(fcen, 0, 1,
                        mp.FluxRegion(center=mp.Vector3(0,0.5*sy-dpml),size=mp.Vector3(sx)))

sim.run(until_after_sources=20)

flux = mp.get_fluxes(flux_mon)
print("flux:, {}".format(flux[0]))
