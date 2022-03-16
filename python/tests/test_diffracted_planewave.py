import unittest
import meep as mp
import math
import cmath
import numpy as np


class TestDiffractedPlanewave(unittest.TestCase):

  def run_binary_grating_diffraction(self, gp, gh, gdc, theta, bands, orders):
    resolution = 50        # pixels/um

    dpml = 1.0             # PML thickness
    dsub = 3.0             # substrate thickness
    dpad = 3.0             # length of padding between grating and PML

    sx = dpml+dsub+gh+dpad+dpml
    sy = gp

    cell_size = mp.Vector3(sx,sy,0)
    absorber_layers = [mp.Absorber(thickness=dpml,direction=mp.X)]

    wvl = 0.5              # center wavelength
    fcen = 1/wvl           # center frequency
    df = 0.05*fcen         # frequency width

    ng = 1.5
    glass = mp.Medium(index=ng)

    # rotation angle of incident planewave
    # counter clockwise (CCW) about Z axis, 0 degrees along +X axis
    theta_in = math.radians(theta)

    eig_parity = mp.ODD_Z

    # k (in source medium) with correct length (plane of incidence: XY)
    k = mp.Vector3(fcen*ng).rotate(mp.Vector3(z=1), theta_in)

    symmetries = []
    if theta_in == 0:
      k = mp.Vector3()
      eig_parity += mp.EVEN_Y
      symmetries = [mp.Mirror(direction=mp.Y)]

    def pw_amp(k,x0):
      def _pw_amp(x):
        return cmath.exp(1j*2*math.pi*k.dot(x+x0))
      return _pw_amp

    src_pt = mp.Vector3(-0.5*sx+dpml,0,0)
    sources = [mp.Source(mp.GaussianSource(fcen,fwidth=df),
                         component=mp.Ez,
                         center=src_pt,
                         size=mp.Vector3(0,sy,0),
                         amp_func=pw_amp(k,src_pt))]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=absorber_layers,
                        k_point=k,
                        default_material=glass,
                        sources=sources,
                        symmetries=symmetries)

    tran_pt = mp.Vector3(0.5*sx-dpml,0,0)
    tran_flux = sim.add_flux(fcen,
                             0,
                             1,
                             mp.FluxRegion(center=tran_pt,
                                           size=mp.Vector3(0,sy,0)))

    sim.run(until_after_sources=10)

    input_flux = mp.get_fluxes(tran_flux)

    sim.reset_meep()

    geometry = [mp.Block(material=glass,
                         size=mp.Vector3(dpml+dsub,mp.inf,mp.inf),
                         center=mp.Vector3(-0.5*sx+0.5*(dpml+dsub),0,0)),
                mp.Block(material=glass,
                         size=mp.Vector3(gh,gdc*gp,mp.inf),
                         center=mp.Vector3(-0.5*sx+dpml+dsub+0.5*gh,0,0))]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=absorber_layers,
                        geometry=geometry,
                        k_point=k,
                        sources=sources,
                        symmetries=symmetries)

    tran_flux = sim.add_mode_monitor(fcen,
                                     0,
                                     1,
                                     mp.FluxRegion(center=tran_pt,
                                                   size=mp.Vector3(0,sy,0)))

    sim.run(until_after_sources=100)

    for band,order in zip(bands,orders):
      res = sim.get_eigenmode_coefficients(tran_flux,
                                           [band],
                                           eig_parity=eig_parity)
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
      print("tran:, {:2d} (band), {:2d} (order), {:10.8f} (eigensolver), {:10.8f} (planewave), {:10.8f} (error)".format(band,order,tran_ref,tran_dp,err))

      self.assertAlmostEqual(vg_ref,vg_dp,places=5)
      self.assertAlmostEqual(tran_ref,tran_dp,places=5)

  def test_diffracted_planewave(self):
    self.run_binary_grating_diffraction(2.6,0.4,0.6,0,range(1,6),range(0,5))
    self.run_binary_grating_diffraction(2.6,0.4,0.6,13.4,range(1,6),[-2,-1,-3,0,-4])

    # self.run_binary_grating_diffraction(10.0,0.5,0.5,0,[2,4,6],[1,3,5])
    # self.run_binary_grating_diffraction(10.0,0.5,0.5,10.7,[1,4,8],[-6,-4,-2])

if __name__ == '__main__':
  unittest.main()
