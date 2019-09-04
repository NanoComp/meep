from __future__ import division

import unittest
import meep as mp
import math
from time import time

class TestSpecialKz(unittest.TestCase):

    def refl_planar(self, theta, special_kz):
        resolution = 100  # pixels/um

        dpml = 1.0
        sx = 3+2*dpml
        sy = 1/resolution
        cell_size = mp.Vector3(sx,sy)
        pml_layers = [mp.PML(dpml,direction=mp.X)]

        fcen = 1.0 # source wavelength = 1 um

        k_point = mp.Vector3(z=math.sin(theta)).scale(fcen)

        sources = [mp.Source(mp.GaussianSource(fcen,fwidth=0.2*fcen),
                             component=mp.Ez,
                             center=mp.Vector3(-0.5*sx+dpml),
                             size=mp.Vector3(y=sy))]

        sim = mp.Simulation(cell_size=cell_size,
                            boundary_layers=pml_layers,
                            sources=sources,
                            k_point=k_point,
                            special_kz=special_kz,
                            resolution=resolution)

        refl_fr = mp.FluxRegion(center=mp.Vector3(-0.25*sx),
                                size=mp.Vector3(y=sy))
        refl = sim.add_flux(fcen,0,1,refl_fr)

        sim.run(until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,mp.Vector3(),1e-9))

        empty_flux = mp.get_fluxes(refl)
        empty_data = sim.get_flux_data(refl)
        sim.reset_meep()

        geometry = [mp.Block(material=mp.Medium(index=3.5),
                             size=mp.Vector3(0.5*sx,mp.inf,mp.inf),
                             center=mp.Vector3(0.25*sx))]

        sim = mp.Simulation(cell_size=cell_size,
                            boundary_layers=pml_layers,
                            geometry=geometry,
                            sources=sources,
                            k_point=k_point,
                            special_kz=special_kz,
                            resolution=resolution)

        refl = sim.add_flux(fcen,0,1,refl_fr)
        sim.load_minus_flux_data(refl,empty_data)

        sim.run(until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,mp.Vector3(),1e-9))

        refl_flux = mp.get_fluxes(refl)

        Rmeep = -refl_flux[0]/empty_flux[0]
        return Rmeep

    def test_special_kz(self):
        n1 = 1
        n2 = 3.5

        # compute angle of refracted planewave in medium n2
        # for incident planewave in medium n1 at angle theta_in
        theta_out = lambda theta_in: math.asin(n1*math.sin(theta_in)/n2)

        # compute Fresnel reflectance for P-polarization in medium n2
        # for incident planewave in medium n1 at angle theta_in
        Rfresnel = lambda theta_in: math.fabs((n1*math.cos(theta_out(theta_in))-n2*math.cos(theta_in))/(n1*math.cos(theta_out(theta_in))+n2*math.cos(theta_in)))**2

        theta = math.radians(23)

        start = time()
        Rmeep_no_kz = self.refl_planar(theta, False)
        t_no_kz = time() - start

        start = time()
        Rmeep_kz = self.refl_planar(theta, True)
        t_kz = time() - start

        Rfres = Rfresnel(theta)

        self.assertAlmostEqual(Rmeep_no_kz,Rfres,places=2)
        self.assertAlmostEqual(Rmeep_kz,Rfres,places=2)
        self.assertLess(t_kz,t_no_kz)

if __name__ == '__main__':
    unittest.main()
