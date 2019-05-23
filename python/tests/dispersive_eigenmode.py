
# dispersive_eigenmode.py - Tests the meep eigenmode features (eigenmode source,
# eigenmode decomposition, and get_eigenmode) with dispersive materials.
# TODO: 
#  * check materials with off diagonal components
#  * check magnetic profiles
#  * once imaginary component is supported, check that

from __future__ import division

import unittest
import meep as mp
import numpy as np
from meep import mpb


class TestDispersiveEigenmode(unittest.TestCase):

    def call_chi1(self,material,component,direction,omega):

        sim = mp.Simulation(cell_size=mp.Vector3(1,1,1),
                    default_material=material,
                    resolution=10)

        sim.init_sim()
        v3 = mp.py_v3_to_vec(sim.dimensions, mp.Vector3(0,0,0), sim.is_cylindrical)
        n = 1/np.sqrt(sim.structure.get_chi1inv(int(component),int(direction),v3,omega))
        return n

    def simulate_waveguide_meep(self,material,resolution):
        sx = 1
        sy = 4
        sz = 0
        dpml = 1
        ww = 0.5
        
        cell = mp.Vector3(sx,sy,sz)
        


        geometry = [mp.Block(size=mp.Vector3(mp.inf,ww,mp.inf),
                            center=mp.Vector3(),
                            material=material
                            )]

        fcen = 1/1.55
        df = 0.1*fcen
        numFreqs = 25
        
        sim = mp.Simulation(cell_size=cell,
                            geometry=geometry,
                            resolution=resolution,
                            sources = [mp.Source(mp.GaussianSource(fcen,df),center=mp.Vector3(),component=mp.Ez)]
                            )

        flux1 = sim.add_flux(fcen,df,numFreqs,mp.FluxRegion(center=mp.Vector3(),size=mp.Vector3(y=sy)))

        sim.init_sim()

        res1 = sim.get_eigenmode_coefficients(flux1,[1])

        freqs = np.array(mp.get_flux_freqs(flux1))

        neff_meep = np.squeeze([a.x for a in res1.kpoints]) / freqs

        return freqs, neff_meep
    
    def simulate_waveguide_mpb(self,resolution,material,freqs):
        sx = 0
        sy = 4
        sz = 0
        dpml = 1
        ww = 0.5

        neff = []
        for ifreq in range(freqs.size):
            omega = freqs[ifreq]
            num_bands = 1
            geometry_lattice = mp.Lattice(size=mp.Vector3(0, sy, sz))
            geometry = [mp.Block(size=mp.Vector3(mp.inf,ww,mp.inf),
                                center=mp.Vector3(),
                                material=mp.Medium(index=np.real(np.sqrt(material.epsilon(omega)[0,0]))))
                                ]
            ms = mpb.ModeSolver(
                geometry_lattice=geometry_lattice,
                geometry=geometry,
                resolution=resolution,
                num_bands=num_bands
            )

            k = ms.find_k(mp.NO_PARITY, omega, 1, num_bands, mp.Vector3(1), 1e-3, omega * 3.45,
                    omega * 0.1, omega * 4)
            
            neff.append(k/omega)

        neff_mpb = np.squeeze(neff)

        return neff_mpb
    
    def compare_meep_mpb(self,material):
        resolution = 50
        freqs, neff_meep = self.simulate_waveguide_meep(material,resolution)
        neff_mpb = self.simulate_waveguide_mpb(resolution,material,freqs)

        print('================================')
        print(neff_meep)
        print(neff_mpb)

        from matplotlib import pyplot as plt
        plt.figure()
        plt.plot(1/freqs,neff_meep)
        plt.plot(1/freqs,neff_mpb)
        plt.show()

        self.assertAlmostEqual(neff_meep[0],neff_mpb[0], places=2)
        self.assertAlmostEqual(neff_meep[-1],neff_mpb[-1], places=2)
        

    def test_chi1_routine(self):
        # This test chceks the newly implemented chi1inv routines within the 
        # fields and structure classes by comparing their output to the 
        # python epsilon output.

        from meep.materials import Si, Ag, LiNbO3
        
        # Check Silicon
        w0 = Si.valid_freq_range.min
        w1 = Si.valid_freq_range.max
        self.assertAlmostEqual(np.real(np.sqrt(Si.epsilon([w0,w1])[0,0,0])), self.call_chi1(Si,mp.Ex,mp.X,w0), places=6)
        self.assertAlmostEqual(np.real(np.sqrt(Si.epsilon([w0,w1])[1,0,0])), self.call_chi1(Si,mp.Ex,mp.X,w1), places=6)

        # Check Silver
        w0 = Ag.valid_freq_range.min
        w1 = Ag.valid_freq_range.max
        self.assertAlmostEqual(np.real(np.sqrt(Ag.epsilon([w0,w1])[0,0,0])), self.call_chi1(Ag,mp.Ex,mp.X,w0), places=6)
        self.assertAlmostEqual(np.real(np.sqrt(Ag.epsilon([w0,w1])[1,0,0])), self.call_chi1(Ag,mp.Ex,mp.X,w1), places=6)

        # Check Lithium Niobate (X,X)
        w0 = LiNbO3.valid_freq_range.min
        w1 = LiNbO3.valid_freq_range.max
        self.assertAlmostEqual(np.real(np.sqrt(LiNbO3.epsilon([w0,w1])[0,0,0])), self.call_chi1(LiNbO3,mp.Ex,mp.X,w0), places=6)
        self.assertAlmostEqual(np.real(np.sqrt(LiNbO3.epsilon([w0,w1])[1,0,0])), self.call_chi1(LiNbO3,mp.Ex,mp.X,w1), places=6)

        # Check Lithium Niobate (Z,Z)
        self.assertAlmostEqual(np.real(np.sqrt(LiNbO3.epsilon([w0,w1])[0,2,2])), self.call_chi1(LiNbO3,mp.Ez,mp.Z,w0), places=6)
        self.assertAlmostEqual(np.real(np.sqrt(LiNbO3.epsilon([w0,w1])[1,2,2])), self.call_chi1(LiNbO3,mp.Ez,mp.Z,w1), places=6)

    def test_waveguide(self):
        # This test simultaneously tests the eigenmode source, the eigenmode monitors, 
        # and the mode decomposition routines.

        from meep.materials import Si, Ag, LiNbO3

        # Check Silicon
        self.compare_meep_mpb(Si)
    
        

if __name__ == '__main__':
    unittest.main()
