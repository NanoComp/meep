################################################################
#
################################################################
import argparse
import meep as mp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py

##################################################
# parse arguments 
##################################################
parser = argparse.ArgumentParser()
parser.add_argument('--res',           type=int,   default='20',        help='resolution')
parser.add_argument('--use_symmetry',             default=False, action='store_true')

args         = parser.parse_args()
res          = args.res
use_symmetry = args.use_symmetry

##################################################
# define geometry
##################################################
w    = 1.5      # waveguide thickness
L    = 10       # waveguide length
dair = 1.0
dpml = 1.0

Si = mp.Medium(epsilon=12.0)

sx = dpml+L+dpml
sy = dpml+dair+w+dair+dpml
origin    = mp.Vector3(0,0,0)
cell_size = mp.Vector3(sx,sy,0)

geometry=[ mp.Block( material=mp.Medium(epsilon=10),
                     size=mp.Vector3(mp.inf,w,mp.inf),
                     center=mp.Vector3())]

##################################################
# define sources #################################
##################################################
# temporal envelope 
fcen              = 0.2
df                = 0.2*fcen
temporal_envelope = mp.GaussianSource(fcen, fwidth=df)

# source 1: constant-amplitude Ez source over conventional hypercubic region
src1_center = mp.Vector3(-0.4*L,+0.3*w,0.0);
src1_size   = mp.Vector3(1.7,0.8,0.0)
source1     = mp.Source(src=temporal_envelope, component=mp.Ez,
                            center=src1_center, size=src1_size)

# source 2:
src2_center = mp.Vector3(-0.5,-0.7,0.0)
src2_radius = 1.3
src2_object = mp.Cylinder(center=src2_center, radius=src2_radius)
def radial_phase(p):
    return (       np.cos(4.0*np.pi*np.sqrt(p.x*p.x+p.y*p.y)/src2_radius)
             +1.0j*np.sin(8.0*np.pi*np.sqrt(p.x*p.x+p.y*p.y)/src2_radius)
           )

def vertical_phase(p):
    return (       np.sin(4.0*np.pi*p.y/src2_radius)
             +1.0j*np.sin(1.0*np.pi*p.y/src2_radius)
           )

def horizontal_ramp(p):
    return  np.exp(1.0j*np.radians(36.0))*(p.x-1.0)

source2a = mp.Source(src=temporal_envelope, obj=src2_object,
                     component=mp.Ex, amp_func=radial_phase)
source2b = mp.Source(src=temporal_envelope, obj=src2_object,
                     component=mp.Ez, amp_func=vertical_phase)

# source 3: defined by two prisms
offset=mp.Vector3(0.5,1.0)
src3a_vertices=[
offset+mp.Vector3(2.09400, -2.27500),
offset+mp.Vector3(1.89600, -2.17500),
offset+mp.Vector3(2.13400, -1.55600),
offset+mp.Vector3(2.30900, -1.27500),
offset+mp.Vector3(2.80700, -0.53500),
offset+mp.Vector3(3.03200, 0.91500),
offset+mp.Vector3(3.24700, 1.53500),
offset+mp.Vector3(3.38400, 1.47500),
offset+mp.Vector3(3.17400, 0.89500),
offset+mp.Vector3(2.98400, -0.57500),
offset+mp.Vector3(2.33400, -1.59500)]

src3b_vertices=[
offset+mp.Vector3(3.09400, -2.61100),
offset+mp.Vector3(2.87900, -2.56000),
offset+mp.Vector3(3.13200, -1.88500),
offset+mp.Vector3(3.29000, -0.68600),
offset+mp.Vector3(4.08900, 0.56200),
offset+mp.Vector3(4.30100, 1.14100),
offset+mp.Vector3(4.44500, 1.09900),
offset+mp.Vector3(4.21100, 0.48600),
offset+mp.Vector3(3.45100, -0.76900),
offset+mp.Vector3(3.35700, -1.65600),
offset+mp.Vector3(3.31000, -1.98400)]

src3a_object = mp.Prism(src3a_vertices, height=1)
src3b_object = mp.Prism(src3b_vertices, height=1)
source3a  = mp.Source(src=temporal_envelope, obj=src3a_object,
                      component=mp.Ez,amp_func=radial_phase)
source3b  = mp.Source(src=temporal_envelope, obj=src3b_object,
                      component=mp.Ez,amp_func=horizontal_ramp)

##################################################
# instantiate simulation 
##################################################
boundary_layers = [mp.PML(dpml)]
symmetries = [mp.Mirror(mp.Y)] if use_symmetry else []
sim = mp.Simulation(resolution=res, cell_size=cell_size,
                    boundary_layers=boundary_layers, geometry=geometry,
                    sources=[source1, source2a, source2b, source3a, source3b],
                    symmetries=symmetries)

sim.run(until=1.0)
mp.output_epsilon(sim)

plt.clf();

extents=[-0.5*sx,+0.5*sx,-0.5*sy,+0.5*sy]

matplotlib.rcParams.update({'font.size': 18})
plt.rcParams["font.family"] = "monospace"

plt.subplot(2,2,1);
slice=sim.get_source_slice(type='indicator');
plt.imshow(np.flipud(np.transpose(slice))) 
plt.title("source_slice_type='indicator'",fontsize=22)
plt.colorbar();

plt.subplot(2,2,2);
slice=sim.get_source_slice(type='norm');
plt.imshow(np.flipud(np.transpose(slice)))
plt.title("source_slice_type='norm'",fontsize=22)
plt.colorbar();

plt.subplot(2,2,3);
slice=sim.get_source_slice(type='re ex');
plt.imshow(np.flipud(np.transpose(slice)))
plt.title("source_slice_type='re ex'",fontsize=22)
plt.colorbar();

plt.subplot(2,2,4);
slice=sim.get_source_slice(type='im ez');
plt.imshow(np.flipud(np.transpose(slice)))
plt.title("source_slice_type='im ez'",fontsize=22)
plt.colorbar();

plt.show(False);
