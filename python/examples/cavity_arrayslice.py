import unittest
import meep as mp
import numpy as np
import matplotlib.pyplot as plt

##################################################
# set up the geometry
##################################################
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

blk = mp.Block(size=mp.Vector3(1e20, w, 1e20),
               material=mp.Medium(epsilon=eps))

geometry = [blk]

for i in range(3):
    geometry.append(mp.Cylinder(r, center=mp.Vector3(d / 2 + i)))

for i in range(3):
    geometry.append(mp.Cylinder(r, center=mp.Vector3(d / -2 - i)))

sim = mp.Simulation(cell_size=cell,
                         geometry=geometry,
                         sources=[],
                         boundary_layers=[mp.pml(dpml)],
                         resolution=20)

##################################################
# add sources 
##################################################
sim.sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df),
                    mp.Hz, mp.Vector3())]

#sim.symmetries = [mp.Mirror(mp.Y, phase=-1), mp.Mirror(mp.X, phase=-1)]

##################################################
# run until sources are finished (and no later)
##################################################
sim._run_sources_until(0,[])

##################################################
# get 1D and 2D array slices
##################################################
xMin = -0.25*sx;
xMax = +0.25*sx;
yMin = -0.15*sy;
yMax = +0.15*sy;

# 1D slice of Hz data
dims1d=np.zeros(3,dtype=np.int32)
v1d=mp.volume( mp.vec(xMin, 0.0), mp.vec(xMax, 0.0) )
rank1d  = sim.fields.get_array_slice_dimensions(v1d, dims1d);
NX1 = dims1d[0];
slice1d = np.zeros(NX1, dtype=np.double);
sim.fields.get_array_slice(v1d, mp.Hz, slice1d);

# 2D slice of Hz data
dims2d=np.zeros(3,dtype=np.int32)
v2d=mp.volume( mp.vec(xMin, yMin), mp.vec(xMax, yMax) )
rank2d  = sim.fields.get_array_slice_dimensions(v2d, dims2d);
NX2 = dims2d[0];
NY2 = dims2d[1];
N2 = NX2*NY2;
slice2d = np.zeros(N2, dtype=np.double);
sim.fields.get_array_slice(v2d, mp.Hz, slice2d);

# plot 1D slice
plt.subplot(1,2,1);
x1d=np.linspace(xMin, xMax, dims1d[0]);
plt.plot(x1d, slice1d);

# plot 2D slice
plt.subplot(1,2,2);
dy = (yMax-yMin) / dims2d[1];
dx = (xMax-xMin) / dims2d[0];
(x2d,y2d)=np.mgrid[ slice(xMin, xMax, dx), slice(yMin, yMax, dy)];
plt.contourf( x2d, y2d, np.reshape(slice2d, (dims2d[0], dims2d[1])) )
plt.colorbar();
plt.show()
