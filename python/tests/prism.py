from __future__ import division

import unittest
import numpy as np
import meep as mp

class TestPrism(unittest.TestCase):

  def prism_cyl_eps(self,npts,r):
    resolution = 70

    cell = mp.Vector3(3,3)

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

    if mp.am_master():
      print("epsilon:, {} (prism), {} (cylinder), {} (relative error)".format(abs(prism_eps),abs(cyl_eps),abs(prism_eps-cyl_eps)/abs(cyl_eps)))

    return abs(prism_eps-cyl_eps)/abs(cyl_eps)


  def test_prism(self):
    r = 1.1958710786934182
    d = [self.prism_cyl_eps(50,r),
         self.prism_cyl_eps(100,r),
         self.prism_cyl_eps(200,r)]

    self.assertLess(d[1],d[0])
    self.assertLess(d[2],d[1])

    r = 1.2896871096581341
    d = [self.prism_cyl_eps(22,r),
         self.prism_cyl_eps(44,r),
         self.prism_cyl_eps(88,r)]

    self.assertLess(d[1],d[0])
    self.assertLess(d[2],d[1])    


if __name__ == '__main__':
  unittest.main()
