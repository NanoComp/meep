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
                         boundary_layers=[mp.Pml(dpml)],
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

dims1d=np.zeros(3,dtype=np.int32)
v1d=mp.volume( mp.vec(xMin, 0.0), mp.vec(xMax, 0.0) )

dims2d=np.zeros(3,dtype=np.int32)
v2d=mp.volume( mp.vec(xMin, yMin), mp.vec(xMax, yMax) )

###################################################
# TEMPORARY HACK: the following calls should like this ...
# 
# rank1d  = sim.fields.get_array_slice_dimensions(v1d, dims1d);
# slice1d = sim.fields.get_array_slice(v1d, mp.Hz);
# rank2d  = sim.fields.get_array_slice_dimensions(v2d, dims2d);
# slice2d = sim.fields.get_array_slice(v2d, mp.Hz);
#
# ...but I can't figure out how to write typemaps for C++ functions 
# that accept both array arguments and other arguments, so use 
# the following hack versions for now
###################################################
sim.fields.py_set_array_slice_volume(v1d);
sim.fields.py_set_array_slice_component(mp.Hz);
rank1d = sim.fields.py_get_array_slice_dimensions(dims1d);
slice1d = np.zeros(dims1d[0]);
sim.fields.py_get_array_slice(slice1d);

sim.fields.py_set_array_slice_volume(v2d);
sim.fields.py_set_array_slice_component(mp.Hz);
rank2d = sim.fields.py_get_array_slice_dimensions(dims2d);
slice2d = np.zeros(dims2d[0]*dims2d[1]);
sim.fields.py_get_array_slice(slice2d);
###################################################
# END TEMPORARY HACK
###################################################

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
