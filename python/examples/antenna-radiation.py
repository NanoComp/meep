from __future__ import division, print_function

import meep as mp
import math

resolution = 50

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
    symmetries = [mp.Mirror(mp.Y)]
elif src_cmpt == mp.Ey:
    symmetries = [mp.Mirror(mp.X)]
elif src_cmpt == mp.Ez:
    symmetries = [mp.Mirror(mp.X), mp.Mirror(mp.Y)]

sim = mp.Simulation(cell_size=cell,
                    resolution=resolution,
                    sources=sources,
                    symmetries=symmetries,
                    boundary_layers=pml_layers)

nearfield = sim.add_near2far(fcen, 0, 1,
                             mp.Near2FarRegion(mp.Vector3(0,0.5*sxy), size=mp.Vector3(sxy)),
                             mp.Near2FarRegion(mp.Vector3(0,-0.5*sxy), size=mp.Vector3(sxy), weight=-1),
                             mp.Near2FarRegion(mp.Vector3(0.5*sxy), size=mp.Vector3(0,sxy)),
                             mp.Near2FarRegion(mp.Vector3(-0.5*sxy), size=mp.Vector3(0,sxy), weight=-1))

sim.run(until_after_sources=mp.stop_when_fields_decayed(50, src_cmpt, mp.Vector3(), 1e-8))

r = 1000/fcen      # 1000 wavelengths out from the source
npts = 100         # number of points in [0,2*pi) range of angles

for n in range(npts):
    ff = sim.get_farfield(nearfield, mp.Vector3(r*math.cos(2*math.pi*n/npts),
                                                r*math.sin(2*math.pi*n/npts)))
    print("farfield:, {}, {}, ".format(n, 2*math.pi*n/npts), end='')
    print(", ".join([str(f).strip('()').replace('j','i') for f in ff]))
