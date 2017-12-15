import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# set up the geometry
eps = 13
w = 1.2
r = 0.36
d = 1.4
N = 3
sy = 6
pad = 2
dpml = 1
sx = (2 * (pad + dpml + N)) + d - 1
fcen = 0.25
df = 0.2
nfreq = 500

cell = mp.Vector3(sx, sy, 0)

blk = mp.Block(size=mp.Vector3(mp.inf, w, mp.inf),
               material=mp.Medium(epsilon=eps))

geometry = [blk]

for i in range(3):
    geometry.append(mp.Cylinder(r, center=mp.Vector3(d / 2 + i)))

for i in range(3):
    geometry.append(mp.Cylinder(r, center=mp.Vector3(d / -2 - i)))

sim = mp.Simulation(cell_size=cell,
                    geometry=geometry,
                    sources=[],
                    boundary_layers=[mp.PML(dpml)],
                    resolution=20)

# add sources
sim.sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df),
               mp.Hz, mp.Vector3())]

# run until sources are finished (and no later)
sim._run_sources_until(0, [])

# get 1D and 2D array slices
xMin = -0.25 * sx
xMax = +0.25 * sx
yMin = -0.15 * sy
yMax = +0.15 * sy

# 1D slice of Hz data
size_1d = mp.Vector3(xMax - xMin)
center_1d = mp.Vector3((xMin + xMax) / 2)
slice1d = sim.get_array(center=center_1d, size=size_1d, component=mp.Hz)

# 2D slice of Hz data
size_2d = mp.Vector3(xMax - xMin, yMax - yMin)
center_2d = mp.Vector3((xMin + xMax) / 2, (yMin + yMax) / 2)
slice2d = sim.get_array(center=center_2d, size=size_2d, component=mp.Hz)

# plot 1D slice
plt.subplot(1, 2, 1)
x1d = np.linspace(xMin, xMax, len(slice1d))
plt.plot(x1d, slice1d)

# plot 2D slice
plt.subplot(1, 2, 2)
dy = (yMax - yMin) / slice2d.shape[1]
dx = (xMax - xMin) / slice2d.shape[0]
(x2d, y2d) = np.mgrid[slice(xMin, xMax, dx), slice(yMin, yMax, dy)]
plt.contourf(x2d, y2d, slice2d)
plt.colorbar()
plt.show()
