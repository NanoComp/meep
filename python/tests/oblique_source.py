from __future__ import division

import meep as mp
import math
import unittest

class TestEigenmodeSource(unittest.TestCase):

    def test_waveguide_flux(self):
        cell_size = mp.Vector3(10,10,0)

        pml_layers = [mp.PML(thickness=2.0)]

        rot_angles = range(0,60,20) # rotation angle of waveguide, CCW around z-axis

        fluxes = []
        for t in rot_angles:
            rot_angle = math.radians(t)
            sources = [mp.EigenModeSource(src=mp.GaussianSource(1.0,fwidth=0.1),
                                          size=mp.Vector3(y=10),
                                          center=mp.Vector3(x=-3),
                                          direction=mp.NO_DIRECTION,
                                          eig_kpoint=mp.Vector3(math.cos(rot_angle),math.sin(rot_angle),0),
                                          eig_band=1,
                                          eig_parity=mp.ODD_Z,
                                          eig_match_freq=True)]

            geometry = [mp.Block(center=mp.Vector3(),
                                 size=mp.Vector3(mp.inf,1,mp.inf),
                                 e1 = mp.Vector3(1,0,0).rotate(mp.Vector3(0,0,1), rot_angle),
                                 e2 = mp.Vector3(0,1,0).rotate(mp.Vector3(0,0,1), rot_angle),
                                 material=mp.Medium(index=1.5))]

            sim = mp.Simulation(cell_size=cell_size,
                                resolution=50,
                                boundary_layers=pml_layers,
                                sources=sources,
                                geometry=geometry)

            tran = sim.add_flux(1.0, 0, 1, mp.FluxRegion(center=mp.Vector3(x=3), size=mp.Vector3(y=10)))

            sim.run(until_after_sources=100)

            fluxes.append(mp.get_fluxes(tran)[0])
            print("flux:, {:.2f}, {:.6f}".format(t,fluxes[-1]))

        self.assertAlmostEqual(fluxes[0], fluxes[1], places=0)
        self.assertAlmostEqual(fluxes[1], fluxes[2], places=0)
        self.assertAlmostEqual(fluxes[0], fluxes[2], places=0)

if __name__ == '__main__':
    unittest.main()
