import unittest
import meep as mp
import math
import cmath
import numpy as np


class TestDiffractedPlanewave(unittest.TestCase):

  @classmethod
  def setUp(cls):
    cls.resolution = 50        # pixels/μm

    cls.dpml = 1.0             # PML thickness
    cls.dsub = 3.0             # substrate thickness
    cls.dpad = 3.0             # length of padding between grating and PML

    cls.wvl = 0.5              # center wavelength
    cls.fcen = 1/cls.wvl       # center frequency

    cls.ng = 1.5
    cls.glass = mp.Medium(index=cls.ng)

    cls.pml_layers = [mp.PML(thickness=cls.dpml,direction=mp.X)]


  def run_binary_grating_diffraction(self,gp,gh,gdc,theta):
    """Computes the mode coefficient of the transmitted orders of
       a binary grating given an incident planewave using the two
       approaches of band number and `DiffractedPlanewave` object in
       `get_eigenmode_coefficients` and verifies that they produce
       the same result.
    """
    sx = self.dpml+self.dsub+gh+self.dpad+self.dpml
    sy = gp
    cell_size = mp.Vector3(sx,sy,0)

    # rotation angle of incident planewave
    # counter clockwise (CCW) about Z axis, 0 degrees along +X axis
    theta_in = math.radians(theta)

    # k (in source medium) with correct length (plane of incidence: XY)
    k = mp.Vector3(self.fcen*self.ng).rotate(mp.Vector3(z=1), theta_in)

    eig_parity = mp.ODD_Z
    if theta == 0:
      k = mp.Vector3()
      eig_parity += mp.EVEN_Y
      symmetries = [mp.Mirror(direction=mp.Y)]
    else:
      symmetries = []

    def pw_amp(k,x0):
      def _pw_amp(x):
        return cmath.exp(1j*2*math.pi*k.dot(x+x0))
      return _pw_amp

    src_pt = mp.Vector3(-0.5*sx+self.dpml,0,0)
    sources = [mp.Source(mp.GaussianSource(self.fcen,fwidth=0.1*self.fcen),
                         component=mp.Ez,
                         center=src_pt,
                         size=mp.Vector3(0,sy,0),
                         amp_func=pw_amp(k,src_pt))]

    geometry = [mp.Block(material=self.glass,
                         size=mp.Vector3(self.dpml+self.dsub,mp.inf,mp.inf),
                         center=mp.Vector3(-0.5*sx+0.5*(self.dpml+self.dsub),0,0)),
                mp.Block(material=self.glass,
                         size=mp.Vector3(gh,gdc*gp,mp.inf),
                         center=mp.Vector3(-0.5*sx+self.dpml+self.dsub+0.5*gh,0,0))]

    sim = mp.Simulation(resolution=self.resolution,
                        cell_size=cell_size,
                        boundary_layers=self.pml_layers,
                        geometry=geometry,
                        k_point=k,
                        sources=sources,
                        symmetries=symmetries)

    tran_pt = mp.Vector3(0.5*sx-self.dpml,0,0)
    tran_flux = sim.add_mode_monitor(self.fcen,
                                     0,
                                     1,
                                     mp.FluxRegion(center=tran_pt,
                                                   size=mp.Vector3(0,sy,0)))

    sim.run(until_after_sources=mp.stop_when_fields_decayed(20,mp.Ez,src_pt,1e-6))

    m_plus = int(np.floor((self.fcen-k.y)*gp))
    m_minus = int(np.ceil((-self.fcen-k.y)*gp))

    if theta == 0:
      orders = range(m_plus+1)
    else:
      # ordering of the modes computed by MPB is according to *decreasing*
      # values of kx (i.e., closest to propagation direction of 0° or +x)
      ms = range(m_minus,m_plus+1)
      kx = lambda m: np.power(self.fcen,2) - np.power(k.y+m/gp,2)
      kxs = [kx(m) for m in ms]
      ids = np.flip(np.argsort(kxs))
      orders = [ms[d] for d in ids]

    for band,order in enumerate(orders):
      res = sim.get_eigenmode_coefficients(tran_flux,
                                           [band+1],
                                           eig_parity=eig_parity)
      tran_ref = abs(res.alpha[0,0,0])**2
      if theta == 0:
        tran_ref = 0.5*tran_ref
      vg_ref = res.vgrp[0]
      kdom_ref = res.kdom[0]

      res = sim.get_eigenmode_coefficients(tran_flux,
                                           mp.DiffractedPlanewave((0,order,0),
                                                                  mp.Vector3(0,1,0),
                                                                  1,
                                                                  0))
      tran_dp = abs(res.alpha[0,0,0])**2
      if (theta == 0) and (order == 0):
        tran_dp = 0.5*tran_dp
      vg_dp = res.vgrp[0]
      kdom_dp = res.kdom[0]

      err = abs(tran_ref-tran_dp)/tran_ref
      print("tran:, {:2d} (band), {:2d} (order), "
            "{:10.8f} (band num.), {:10.8f} (diff. pw.), "
            "{:10.8f} (error)".format(band,order,tran_ref,tran_dp,err))

      self.assertAlmostEqual(vg_ref,vg_dp,places=4)
      self.assertAlmostEqual(tran_ref,tran_dp,places=4)
      if theta == 0:
        self.assertAlmostEqual(abs(kdom_ref.x),kdom_dp.x,places=5)
        self.assertAlmostEqual(abs(kdom_ref.y),kdom_dp.y,places=5)
        self.assertAlmostEqual(abs(kdom_ref.z),kdom_dp.z,places=5)
      else:
        self.assertAlmostEqual(kdom_ref.x,kdom_dp.x,places=5)
        self.assertAlmostEqual(kdom_ref.y,kdom_dp.y,places=5)
        self.assertAlmostEqual(kdom_ref.z,kdom_dp.z,places=5)

    print("PASSED.")


  def test_diffracted_planewave(self):
    self.run_binary_grating_diffraction(2.6,0.4,0.6,0)
    self.run_binary_grating_diffraction(2.6,0.4,0.6,13.4)

    # self.run_binary_grating_diffraction(10.0,0.5,0.5,0)
    # self.run_binary_grating_diffraction(10.0,0.5,0.5,10.7)


  def run_mode_source(self,gp,gh,gdc,m,diffpw):
    """Computes the transmitted flux of a
       binary grating given an incident planewave
       specified by the diffraction order `m` in the
       y direction. The incident planewave is defined
       using a mode source with either a band number
       or `DiffractedPlanewave` object specified by
       the boolean flag `diffpw`.
    """
    sx = self.dpml+self.dsub+gh+self.dpad+self.dpml
    sy = gp
    cell_size = mp.Vector3(sx,sy,0)

    ky = m/gp
    theta = math.asin(ky/(self.fcen*self.ng))

    # k (in source medium) with correct length (plane of incidence: XY)
    k = mp.Vector3(self.fcen*self.ng).rotate(mp.Vector3(z=1), theta)

    eig_parity = mp.ODD_Z
    if theta == 0:
      k = mp.Vector3()
      eig_parity += mp.EVEN_Y
      symmetries = [mp.Mirror(direction=mp.Y)]
    else:
      symmetries = []

    src_pt = mp.Vector3(-0.5*sx+self.dpml,0,0)
    if diffpw:
      # the *zeroth* diffraction order defines a planewave with a
      # wavevector equal to the `k_point` of the `Simulation` object
      sources = [mp.EigenModeSource(mp.GaussianSource(self.fcen,fwidth=0.1*self.fcen),
                                    center=src_pt,
                                    size=mp.Vector3(0,sy,0),
                                    eig_band=mp.DiffractedPlanewave((0,0,0),
                                                                    mp.Vector3(0,1,0),
                                                                    1,
                                                                    0))]
    else:
      sources = [mp.EigenModeSource(mp.GaussianSource(self.fcen,fwidth=0.1*self.fcen),
                                    center=src_pt,
                                    size=mp.Vector3(0,sy,0),
                                    direction=mp.NO_DIRECTION,
                                    eig_kpoint=k,
                                    eig_band=1,
                                    eig_parity=eig_parity)]

    geometry = [mp.Block(material=self.glass,
                         size=mp.Vector3(self.dpml+self.dsub,mp.inf,mp.inf),
                         center=mp.Vector3(-0.5*sx+0.5*(self.dpml+self.dsub),0,0)),
                mp.Block(material=self.glass,
                         size=mp.Vector3(gh,gdc*gp,mp.inf),
                         center=mp.Vector3(-0.5*sx+self.dpml+self.dsub+0.5*gh,0,0))]

    sim = mp.Simulation(resolution=self.resolution,
                        cell_size=cell_size,
                        boundary_layers=self.pml_layers,
                        k_point=k,
                        geometry=geometry,
                        sources=sources,
                        symmetries=symmetries)

    tran_pt = mp.Vector3(0.5*sx-self.dpml,0,0)
    tran_flux = sim.add_flux(self.fcen,
                             0,
                             1,
                             mp.FluxRegion(center=tran_pt,
                                           size=mp.Vector3(0,sy,0)))

    sim.run(until_after_sources=mp.stop_when_fields_decayed(20,mp.Ez,src_pt,1e-6))

    tran = mp.get_fluxes(tran_flux)[0]

    return tran


  def test_mode_source(self):
    tran_bn = self.run_mode_source(1.5,0.5,0.3,2,False)
    tran_dp = self.run_mode_source(1.5,0.5,0.3,2,True)
    print("mode-source:, "
          "{:.5f} (band number), "
          "{:.5f} (diffraction order)".format(tran_bn,
                                              tran_dp))
    self.assertAlmostEqual(tran_bn,
                           tran_dp,
                           places=3 if mp.is_single_precision() else 4)


if __name__ == '__main__':
  unittest.main()
