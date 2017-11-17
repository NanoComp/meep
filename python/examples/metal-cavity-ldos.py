from __future__ import division

import math
import meep as mp


resolution = 200
sxy = 2
dpml = 1
sxy = sxy + 2 * dpml
cell = mp.Vector3(sxy, sxy, 0)

pml_layers = [mp.PML(dpml)]
a = 1
t = 0.1
geometry = [mp.Block(mp.Vector3(a + 2 * t, a + 2 * t, 1e20)),
            mp.Block(mp.Vector3(a, a, 1e20), material=mp.Medium(epsilon=1.0))]

w = 0
if w > 0:
    geometry.append(mp.Block(center=mp.Vector3(a / 2), size=mp.Vector3(2 * t, w, 1e20),
                             material=mp.Medium(epsilon=1.0)))

fcen = math.sqrt(0.5) / a
df = 0.2
sources = [mp.Source(src=mp.GaussianSource(fcen, fwidth=df), component=mp.Ez,
                     center=mp.Vector3())]

symmetries = [mp.Mirror(mp.Y)]

Th = 500

sim = mp.Simulation(cell_size=cell,
                    geometry=geometry,
                    boundary_layers=pml_layers,
                    sources=sources,
                    symmetries=symmetries,
                    resolution=resolution)

h = mp.Harminv(mp.Ez, mp.Vector3(), fcen, df)
sim.run(mp.after_sources(h), until_after_sources=Th)

m = h.modes[0]
f = m.freq
Q = m.Q
Vmode = 0.25 * a * a
print("ldos0:, {}".format(Q / Vmode / (2 * math.pi * f * math.pi * 0.5)))

sim.reset_meep()
T = 2 * Q * (1 / f)
sim.run(mp.dft_ldos(f, 0, 1), until_after_sources=T)
