
# visualization.py - Tests the visualization module. Checks 2D 
# plotting of a waveguide with several sources, monitors, and
# boundary conditions. Checks for subdomain plots.
#
# Also tests the animation run function, mp4 output, jshtml output, and git output.

from __future__ import division

import unittest
from subprocess import call

import meep as mp
import numpy as np

# Make sure we have matplotlib installed
import matplotlib
matplotlib.use('agg') # Set backend for consistency and to pull pixels quickly
from matplotlib import pyplot as plt
import io

def hash_figure(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='raw')
    buf.seek(0)
    data = np.frombuffer(buf.getvalue(), dtype=np.uint8)
    return np.sum((data > np.mean(data)) + data)

def setup_sim(zDim=0):
    cell = mp.Vector3(16,8,zDim)

    # A simple waveguide
    geometry = [mp.Block(mp.Vector3(mp.inf,1,1),
                     center=mp.Vector3(),
                     material=mp.Medium(epsilon=12))]
    
    # Add point sources
    sources = [mp.Source(mp.ContinuousSource(frequency=0.15),
                        component=mp.Ez,
                        center=mp.Vector3(-5,0)),
                mp.Source(mp.ContinuousSource(frequency=0.15),
                        component=mp.Ez,
                        center=mp.Vector3(0,2))
                        ]

    # Add line sources
    sources += [mp.Source(mp.ContinuousSource(frequency=0.15),
                        component=mp.Ez,
                        size=mp.Vector3(0,2,0),
                        center=mp.Vector3(-6,0)),
                mp.Source(mp.ContinuousSource(frequency=0.15),
                        component=mp.Ez,
                        size=mp.Vector3(2,0,0),
                        center=mp.Vector3(0,1))]
        
    # Add plane sources
    sources += [mp.Source(mp.ContinuousSource(frequency=0.15),
                        component=mp.Ez,
                        size=mp.Vector3(2,2,0),
                        center=mp.Vector3(-3,0)),
                mp.Source(mp.ContinuousSource(frequency=0.15),
                        component=mp.Ez,
                        size=mp.Vector3(2,2,0),
                    center=mp.Vector3(0,-2))]
    
    # Different pml layers
    pml_layers = [mp.PML(2.0,mp.X),mp.PML(1.0,mp.Y,mp.Low),mp.PML(1.5,mp.Y,mp.High)]

    resolution = 10

    sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=sources,
                    resolution=resolution)
    # Line monitor
    sim.add_flux(1,0,1,mp.FluxRegion(center=mp.Vector3(5,0,0),size=mp.Vector3(0,4), direction=mp.X))

    # Plane monitor
    sim.add_flux(1,0,1,mp.FluxRegion(center=mp.Vector3(2,0,0),size=mp.Vector3(4,4), direction=mp.X))
    
    return sim

class TestVisualization(unittest.TestCase):
    
    def test_plot2D(self):
        
        # Check plotting of geometry with several sources, monitors, and PMLs
        f = plt.figure()
        ax = f.gca()
        sim = setup_sim()
        ax = sim.plot2D(ax=ax)
        if mp.am_master():
            hash_figure(f)
            #self.assertAlmostEqual(hash_figure(f),10231488)

        # Check plotting of fields after timestepping
        f = plt.figure()
        ax = f.gca()
        sim.run(until=200)
        ax = sim.plot2D(ax=ax,fields=mp.Ez)
        if mp.am_master():
            hash_figure(f)
            #self.assertAlmostEqual(hash_figure(f),79786722)

        # Check output_plane feature
        f = plt.figure()
        ax = f.gca()
        vol = mp.Volume(center=mp.Vector3(),size=mp.Vector3(2,2))
        ax = sim.plot2D(ax=ax,fields=mp.Ez,output_plane=vol)
        if mp.am_master():
            hash_figure(f)
            #self.assertAlmostEqual(hash_figure(f),68926258)
    
    @unittest.skipIf(call(['which', 'ffmpeg']) != 0, "ffmpeg is not installed")
    def test_animation_output(self):
        # Check without normalization
        sim = setup_sim()
        Animate = mp.Animate2D(sim,fields=mp.Ez, realtime=False, normalize=False)
        sim.run(mp.at_every(1,Animate),until=5)

        # Check with normalization
        animation = mp.Animate2D(sim,mp.Ez,realtime=False,normalize=True)
        sim.run(mp.at_every(1),until=25)

        # Check mp4 output
        Animate.to_mp4(10,'test.mp4')

        # Check gif output
        Animate.to_gif(10,'test.gif')

        # Check jshtml output
        Animate.to_jshtml(10)
    '''
    Travis does not play well with Mayavi
    def test_3D_mayavi(self):
        sim = setup_sim(4)
        sim.plot3D()
    '''

if __name__ == '__main__':
    unittest.main()
    
