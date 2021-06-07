from __future__ import division

import unittest
import meep as mp
import math
import cmath
import numpy as np

class TestEigCoeffs(unittest.TestCase):

  def run_binary_grating_oblique(self, theta):
  
    resolution = 30        # pixels/um

    dpml = 1.0             # PML thickness
    dsub = 1.0             # substrate thickness
    dpad = 1.0             # length of padding between grating and pml
    gp = 6.0               # grating period
    gh = 0.5               # grating height
    gdc = 0.5              # grating duty cycle

    sx = dpml+dsub+gh+dpad+dpml
    sy = gp

    cell_size = mp.Vector3(sx,sy,0)

    # replace anisotropic PML with isotropic Absorber to attenuate parallel-directed fields of oblique source
    abs_layers = [mp.Absorber(thickness=dpml,direction=mp.X)] 

    wvl = 0.5              # center wavelength
    fcen = 1/wvl           # center frequency
    df = 0.05*fcen         # frequency width

    ng = 1.5
    glass = mp.Medium(index=ng)

    # rotation angle of incident planewave; CCW about Y axis, 0 degrees along +X axis
    theta_in = math.radians(theta)

    # k (in source medium) with correct length (plane of incidence: XY)
    k = mp.Vector3(math.cos(theta_in),math.sin(theta_in),0).scale(fcen*ng)

    symmetries = []
    eig_parity = mp.ODD_Z
    if theta_in == 0:
      k = mp.Vector3(0,0,0)
      symmetries = [mp.Mirror(mp.Y)]
      eig_parity += mp.EVEN_Y
  
    def pw_amp(k,x0):
      def _pw_amp(x):
        return cmath.exp(1j*2*math.pi*k.dot(x+x0))
      return _pw_amp

    src_pt = mp.Vector3(-0.5*sx+dpml+0.3*dsub,0,0)
    sources = [mp.Source(mp.GaussianSource(fcen,fwidth=df),
                         component=mp.Ez,
                         center=src_pt,
                         size=mp.Vector3(0,sy,0),
                         amp_func=pw_amp(k,src_pt))]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=abs_layers,
                        k_point=k,
                        default_material=glass,
                        sources=sources,
                        symmetries=symmetries)

    refl_pt = mp.Vector3(-0.5*sx+dpml+0.5*dsub,0,0)
    refl_flux = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=refl_pt, size=mp.Vector3(0,sy,0)))

    sim.run(until_after_sources=100)
  
    input_flux = mp.get_fluxes(refl_flux)
    input_flux_data = sim.get_flux_data(refl_flux)

    sim.reset_meep()

    geometry = [mp.Block(material=glass, size=mp.Vector3(dpml+dsub,mp.inf,mp.inf), center=mp.Vector3(-0.5*sx+0.5*(dpml+dsub),0,0)),
                mp.Block(material=glass, size=mp.Vector3(gh,gdc*gp,mp.inf), center=mp.Vector3(-0.5*sx+dpml+dsub+0.5*gh,0,0))]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=abs_layers,
                        geometry=geometry,
                        k_point=k,
                        sources=sources,
                        symmetries=symmetries)

    refl_flux = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=refl_pt, size=mp.Vector3(0,sy,0)))
    sim.load_minus_flux_data(refl_flux,input_flux_data)

    tran_pt = mp.Vector3(0.5*sx-dpml-0.5*dpad,0,0)
    tran_flux = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=tran_pt, size=mp.Vector3(0,sy,0)))

    sim.run(until_after_sources=100)

    nm_r = np.floor((fcen*ng-k.y)*gp)-np.ceil((-fcen*ng-k.y)*gp) # number of reflected orders
    if theta_in == 0:
      nm_r = nm_r/2 # since eig_parity removes degeneracy in y-direction
    nm_r = int(nm_r)

    res = sim.get_eigenmode_coefficients(refl_flux, range(1,nm_r+1), eig_parity=eig_parity)
    r_coeffs = res.alpha

    Rsum = 0
    for nm in range(nm_r):
      Rsum += abs(r_coeffs[nm,0,1])**2/input_flux[0]

    nm_t = np.floor((fcen-k.y)*gp)-np.ceil((-fcen-k.y)*gp)       # number of transmitted orders
    if theta_in == 0:
      nm_t = nm_t/2 # since eig_parity removes degeneracy in y-direction
    nm_t = int(nm_t)

    res = sim.get_eigenmode_coefficients(tran_flux, range(1,nm_t+1), eig_parity=eig_parity)
    t_coeffs = res.alpha

    Tsum = 0
    for nm in range(nm_t):
      Tsum += abs(t_coeffs[nm,0,0])**2/input_flux[0]

    r_flux = mp.get_fluxes(refl_flux)
    t_flux = mp.get_fluxes(tran_flux)
    Rflux = -r_flux[0]/input_flux[0]
    Tflux =  t_flux[0]/input_flux[0]

    self.assertAlmostEqual(Rsum,Rflux,places=2)
    self.assertAlmostEqual(Tsum,Tflux,places=2)

  def test_binary_grating(self):
    self.run_binary_grating_oblique(0)
    self.run_binary_grating_oblique(10.7)

if __name__ == '__main__':
  unittest.main()
