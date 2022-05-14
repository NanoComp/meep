import unittest
import meep as mp
import math
import cmath
import numpy as np


class TestDiffractedPlanewave(unittest.TestCase):

  @classmethod
  def setUp(cls):
    cls.resolution = 50   # pixels/μm

    cls.dpml = 1.0        # PML thickness
    cls.dsub = 3.0        # substrate thickness
    cls.dpad = 3.0        # length of padding between grating and PML
    cls.gp = 2.6          # grating period
    cls.gh = 0.4          # grating height
    cls.gdc = 0.5         # grating duty cycle

    cls.sx = cls.dpml+cls.dsub+cls.gh+cls.dpad+cls.dpml
    cls.sy = cls.gp

    cls.cell_size = mp.Vector3(cls.sx,cls.sy,0)
    cls.absorber_layers = [mp.PML(thickness=cls.dpml,direction=mp.X)]

    cls.wvl = 0.5             # center wavelength
    cls.fcen = 1/cls.wvl      # center frequency
    cls.df = 0.05*cls.fcen    # frequency width

    cls.ng = 1.5
    cls.glass = mp.Medium(index=cls.ng)

    cls.eig_parity = mp.ODD_Z

    cls.src_pt = mp.Vector3(-0.5*cls.sx+cls.dpml,0,0)
    cls.tran_pt = mp.Vector3(0.5*cls.sx-cls.dpml,0,0)

    cls.geometry = [mp.Block(material=cls.glass,
                             size=mp.Vector3(cls.dpml+cls.dsub,mp.inf,mp.inf),
                             center=mp.Vector3(-0.5*cls.sx+0.5*(cls.dpml+cls.dsub),0,0)),
                    mp.Block(material=cls.glass,
                             size=mp.Vector3(cls.gh,cls.gdc*cls.gp,mp.inf),
                             center=mp.Vector3(-0.5*cls.sx+cls.dpml+cls.dsub+0.5*cls.gh,0,0))]

    cls.tol = 1e-8  # tolerance for mp.stop_dft_decayed termination criteria


  def run_mode_decomposition(self,theta,bands,orders):
    """Computes the transmittance of the diffraction orders
       of a binary grating given an incident planewave at
       angle `theta`. The transmittance is computed using
       mode decomposition in two different ways:
       (1) band number (`bands`) and (2) diffraction order (`orders`)
       via a `DiffractedPlanewave` object. The test verifies
       that these two methods produce equivalent results.
    """
    # rotation angle of incident planewave
    # counter clockwise (CCW) about Z axis, 0 degrees along +X axis
    theta_in = math.radians(theta)

    # k (in source medium) with correct length (plane of incidence: XY)
    k = mp.Vector3(self.fcen*self.ng).rotate(mp.Vector3(z=1),theta_in)

    symmetries = []
    if theta_in == 0:
      k = mp.Vector3()
      self.eig_parity += mp.EVEN_Y
      symmetries = [mp.Mirror(direction=mp.Y)]

    def pw_amp(k,x0):
      def _pw_amp(x):
        return cmath.exp(1j*2*math.pi*k.dot(x+x0))
      return _pw_amp

    sources = [mp.Source(mp.GaussianSource(self.fcen,fwidth=self.df),
                         component=mp.Ez,
                         center=self.src_pt,
                         size=mp.Vector3(0,self.sy,0),
                         amp_func=pw_amp(k,self.src_pt))]

    sim = mp.Simulation(resolution=self.resolution,
                        cell_size=self.cell_size,
                        boundary_layers=self.absorber_layers,
                        k_point=k,
                        default_material=self.glass,
                        sources=sources,
                        symmetries=symmetries)

    tran_flux = sim.add_flux(self.fcen,
                             0,
                             1,
                             mp.FluxRegion(center=self.tran_pt,
                                           size=mp.Vector3(0,self.sy,0)))

    sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=self.tol))

    input_flux = mp.get_fluxes(tran_flux)

    sim.reset_meep()

    sim = mp.Simulation(resolution=self.resolution,
                        cell_size=self.cell_size,
                        boundary_layers=self.absorber_layers,
                        geometry=self.geometry,
                        k_point=k,
                        sources=sources,
                        symmetries=symmetries)

    tran_flux = sim.add_mode_monitor(self.fcen,
                                     0,
                                     1,
                                     mp.FluxRegion(center=self.tran_pt,
                                                   size=mp.Vector3(0,self.sy,0)))

    sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=self.tol))

    for band,order in zip(bands,orders):
      res = sim.get_eigenmode_coefficients(tran_flux,
                                           [band],
                                           eig_parity=self.eig_parity)
      tran_ref = abs(res.alpha[0,0,0])**2/input_flux[0]
      if (theta_in == 0):
        tran_ref = 0.5*tran_ref
      vg_ref = res.vgrp[0]

      res = sim.get_eigenmode_coefficients(tran_flux,
                                           mp.DiffractedPlanewave((0,order,0),
                                                                  mp.Vector3(0,1,0),
                                                                  1,
                                                                  0))
      if res is not None:
        tran_dp = abs(res.alpha[0,0,0])**2/input_flux[0]
        if ((theta_in == 0) and (order == 0)):
          tran_dp = 0.5*tran_dp
      else:
        tran_dp = 0
      vg_dp = res.vgrp[0]

      err = abs(tran_ref-tran_dp)/tran_ref
      print("info:, {:2d} (band number), {:2d} (diffraction order), "
            "{:10.8f} [trans. (band number)], "
            "{:10.8f} [trans. (diffracted order)], "
            "{:10.8f} (error)".format(band,
                                      order,
                                      tran_ref,
                                      tran_dp,
                                      err))

      self.assertAlmostEqual(vg_ref,vg_dp,places=5)
      self.assertAlmostEqual(tran_ref,tran_dp,places=5)


  def run_mode_source(self,m,diffpw):
    """Computes the transmitted Poynting flux of a
       binary grating given an incident planewave
       specified by the diffraction order `m` in the
       y direction. The incident planewave is defined
       using a mode source with either a band number
       or `DiffractedPlanewave` object specified by
       the boolean flag `diffpw`.
    """
    ky = m/self.gp
    theta = math.asin(ky/(self.fcen*self.ng))

    # k (in source medium) with correct length (plane of incidence: XY)
    k = mp.Vector3(self.fcen*self.ng).rotate(mp.Vector3(z=1), theta)

    symmetries = []
    if theta == 0:
      k = mp.Vector3()
      self.eig_parity += mp.EVEN_Y
      symmetries = [mp.Mirror(direction=mp.Y)]

    if diffpw:
      # the *zeroth* diffraction order defines a planewave with a
      # wavevector equal to the `k_point` of the `Simulation` object
      sources = [mp.EigenModeSource(mp.GaussianSource(self.fcen,fwidth=self.df),
                                    center=self.src_pt,
                                    size=mp.Vector3(0,self.sy,0),
                                    eig_band=mp.DiffractedPlanewave((0,0,0),
                                                                    mp.Vector3(0,1,0),
                                                                    1,
                                                                    0))]
    else:
      sources = [mp.EigenModeSource(mp.GaussianSource(self.fcen,fwidth=self.df),
                                    center=self.src_pt,
                                    size=mp.Vector3(0,self.sy,0),
                                    direction=mp.NO_DIRECTION,
                                    eig_kpoint=k,
                                    eig_band=1,
                                    eig_parity=self.eig_parity)]

    sim = mp.Simulation(resolution=self.resolution,
                        cell_size=self.cell_size,
                        boundary_layers=self.absorber_layers,
                        k_point=k,
                        geometry=self.geometry,
                        sources=sources,
                        symmetries=symmetries)

    tran_flux = sim.add_flux(self.fcen,
                             0,
                             1,
                             mp.FluxRegion(center=self.tran_pt,
                                           size=mp.Vector3(0,self.sy,0)))

    sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=self.tol))

    tran = mp.get_fluxes(tran_flux)[0]

    sim.reset_meep()

    return tran


  def test_mode_decomposition(self):
    self.run_mode_decomposition(0,range(1,6),range(0,5))
    self.run_mode_decomposition(13.4,range(1,6),[-2,-1,-3,0,-4])


  def test_mode_source(self):
    m = 3
    tran_bandnum = self.run_mode_source(m,False)
    tran_diffpw = self.run_mode_source(m,True)
    print("info:, {} (diffraction order), "
          "{:.5f} [trans. (band number)],"
          "{:.5f} [trans. (diffraction order)]".format(m,
                                                       tran_bandnum,
                                                       tran_diffpw))
    self.assertAlmostEqual(tran_bandnum,
                           tran_diffpw,
                           places=3 if mp.is_single_precision() else 4)


if __name__ == '__main__':
  unittest.main()
