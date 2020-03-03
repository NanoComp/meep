from __future__ import division

import os
import unittest
import numpy as np
import meep as mp

class TestPrism(unittest.TestCase):

  def nonconvex_marching_squares(self,idx,npts):
    resolution = 15

    cell = mp.Vector3(10,10)

    data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))
    vertices_file = os.path.join(data_dir, 'nonconvex_prism_vertices{}.npz'.format(idx))
    vertices_obj = np.load(vertices_file)

    ## prism verticies precomputed from analytic blob shape using
    ## marching squares algorithm of skimage.measure.find_contours
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

    print("epsilon-sum:, {} (prism-msq)".format(abs(prism_eps)))

    return prism_eps


  def convex_marching_squares(self,npts):
    resolution = 50

    cell = mp.Vector3(3,3)

    data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))
    vertices_file = os.path.join(data_dir, 'convex_prism_vertices.npz')
    vertices_obj = np.load(vertices_file)

    ## prism vertices precomputed for a circle of radius 1.0 using
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

    print("epsilon-sum:, {} (prism-msq), {} (cylinder), {} (relative error)".format(abs(prism_eps),abs(cyl_eps),abs((prism_eps-cyl_eps)/cyl_eps)))

    return abs((prism_eps-cyl_eps)/cyl_eps)


  def convex_circle(self,npts,r):
    resolution = 50

    cell = mp.Vector3(3,3)

    ### prism vertices computed as equally-spaced points
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

    print("epsilon-sum:, {} (prism-cyl), {} (cylinder), {} (relative error)".format(abs(prism_eps),abs(cyl_eps),abs((prism_eps-cyl_eps)/cyl_eps)))

    return abs((prism_eps-cyl_eps)/cyl_eps)


  def test_prism(self):
    print("Testing Non-Convex Prism #1 using marching squares algorithm...")
    d1_a = self.nonconvex_marching_squares(1,208)
    d1_b = self.nonconvex_marching_squares(1,448)
    d1_ref = self.nonconvex_marching_squares(1,904)

    self.assertLess(abs((d1_b-d1_ref)/d1_ref),abs((d1_a-d1_ref)/d1_ref))

    print("Testing Non-Convex Prism #2 using marching squares algorithm...")
    d2_a = self.nonconvex_marching_squares(2,128)
    d2_b = self.nonconvex_marching_squares(2,256)
    d2_ref = self.nonconvex_marching_squares(2,516)

    self.assertLess(abs((d2_b-d2_ref)/d2_ref),abs((d2_a-d2_ref)/d2_ref))

    print("Testing Non-Convex Prism #3 using marching squares algorithm...")
    d3_a = self.nonconvex_marching_squares(3,164)
    d3_b = self.nonconvex_marching_squares(3,336)
    d3_ref = self.nonconvex_marching_squares(3,672)

    self.assertLess(abs((d3_b-d3_ref)/d3_ref),abs((d3_a-d3_ref)/d3_ref))

    print("Testing Convex Prism using marching squares algorithm...")
    d = [self.convex_marching_squares(92),
         self.convex_marching_squares(192),
         self.convex_marching_squares(392)]

    self.assertLess(d[1],d[0])
    self.assertLess(d[2],d[1])

    print("Testing Convex Prism #1 using circle formula...")
    r = 1.0458710786934182
    d = [self.convex_circle(51,r),
         self.convex_circle(101,r),
         self.convex_circle(201,r)]

    self.assertLess(d[1],d[0])
    self.assertLess(d[2],d[1])

    print("Testing Convex Prism #2 using circle formula...")
    r = 1.2896871096581341
    d = [self.convex_circle(31,r),
         self.convex_circle(61,r),
         self.convex_circle(121,r)]

    self.assertLess(d[1],d[0])
    self.assertLess(d[2],d[1])

if __name__ == '__main__':
  unittest.main()
