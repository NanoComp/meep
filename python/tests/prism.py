from __future__ import division

import os
import unittest
import numpy as np
import meep as mp

class TestPrism(unittest.TestCase):

  def prism_marching_squares(self,npts):
    resolution = 50

    cell = mp.Vector3(3,3)

    data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))
    vertices_file = os.path.join(data_dir, 'prism_vertices.npz')
    vertices_obj = np.load(vertices_file)

    ## prism vertices precomputed for a circle of radius 1 using
    ## marching squares algorithm of skimage.measure.find_contours
    ## ref: https://github.com/NanoComp/meep/issues/1060
    vertices_data = vertices_obj["N{}".format(npts)]
    vertices = [mp.Vector3(v[0],v[1],0) for v in vertices_data]

    geometry = [mp.Prism(vertices,
                         height=mp.inf,
                         material=mp.Medium(epsilon=12))]

    sim = mp.Simulation(cell_size=cell,
                        geometry=geometry,
                        resolution=resolution)

    sim.init_sim()

    prism_eps = sim.integrate_field_function([mp.Dielectric], lambda r,eps: eps)

    sim.reset_meep()

    geometry = [mp.Cylinder(radius=1.0,
                            center=mp.Vector3(),
                            height=mp.inf,
                            material=mp.Medium(epsilon=12))]

    sim = mp.Simulation(cell_size=cell,
                        geometry=geometry,
                        resolution=resolution)

    sim.init_sim()

    cyl_eps = sim.integrate_field_function([mp.Dielectric], lambda r,eps: eps)

    print("epsilon-sum:, {} (prism-msq), {} (cylinder), {} (relative error)".format(abs(prism_eps),abs(cyl_eps),abs(prism_eps-cyl_eps)/abs(cyl_eps)))

    return abs(prism_eps-cyl_eps)/abs(cyl_eps)


  def prism_circle(self,npts,r):
    resolution = 50

    cell = mp.Vector3(3,3)

    ### prism vertices computed as npts equally-spaced points
    ### along the circumference of a circle with radius r
    angles = 2*np.pi/npts * np.arange(npts)
    vertices = [mp.Vector3(r*np.cos(ang),r*np.sin(ang)) for ang in angles]
    geometry = [mp.Prism(vertices,
                         height=mp.inf,
                         material=mp.Medium(epsilon=12))]

    sim = mp.Simulation(cell_size=cell,
                        geometry=geometry,
                        resolution=resolution)

    sim.init_sim()

    prism_eps = sim.integrate_field_function([mp.Dielectric], lambda r,eps: eps)

    sim.reset_meep()

    geometry = [mp.Cylinder(radius=r,
                            center=mp.Vector3(),
                            height=mp.inf,
                            material=mp.Medium(epsilon=12))]

    sim = mp.Simulation(cell_size=cell,
                        geometry=geometry,
                        resolution=resolution)

    sim.init_sim()

    cyl_eps = sim.integrate_field_function([mp.Dielectric], lambda r,eps: eps)

    print("epsilon-sum:, {} (prism-cyl), {} (cylinder), {} (relative error)".format(abs(prism_eps),abs(cyl_eps),abs(prism_eps-cyl_eps)/abs(cyl_eps)))

    return abs(prism_eps-cyl_eps)/abs(cyl_eps)


  def test_prism(self):
    print("Testing Prism object using marching squares algorithm...")
    d = [self.prism_marching_squares(92),
         self.prism_marching_squares(192),
         self.prism_marching_squares(392)]

    self.assertLess(d[1],d[0])
    self.assertLess(d[2],d[1])

    print("Testing Prism object using circle formula...")
    r = 1.0458710786934182
    d = [self.prism_circle(51,r),
         self.prism_circle(101,r),
         self.prism_circle(201,r)]

    self.assertLess(d[1],d[0])
    self.assertLess(d[2],d[1])

    r = 1.2896871096581341
    d = [self.prism_circle(31,r),
         self.prism_circle(61,r),
         self.prism_circle(121,r)]

    self.assertLess(d[1],d[0])
    self.assertLess(d[2],d[1])

if __name__ == '__main__':
  unittest.main()
