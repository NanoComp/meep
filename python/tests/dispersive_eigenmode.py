
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
import h5py
class TestDispersiveEigenmode(unittest.TestCase):
    # ----------------------------------------- #
    # ----------- Helper Functions ------------ #
    # ----------------------------------------- #
    # Directly cals the C++ chi1 routine
    def call_chi1(self,material,omega):

        sim = mp.Simulation(cell_size=mp.Vector3(1,1,1),
                    default_material=material,
                    resolution=20)

        sim.init_sim()
        v3 = mp.py_v3_to_vec(sim.dimensions, mp.Vector3(0,0,0), sim.is_cylindrical)
        chi1inv = np.zeros((3,3),dtype=np.float64)
        for i, com in enumerate([mp.Ex,mp.Ey,mp.Ez]):
            for k, dir in enumerate([mp.X,mp.Y,mp.Z]):
                chi1inv[i,k] = sim.structure.get_chi1inv(com,dir,v3,omega)
        n = np.real(np.sqrt(np.linalg.inv(chi1inv.astype(np.complex128))))

        n_actual = np.real(np.sqrt(material.epsilon(omega).astype(np.complex128)))
        
        np.testing.assert_allclose(n,n_actual)
    
    def verify_output_and_slice(self,material,omega):
        # Since the slice routines average the diagonals, we need to do that too:
        chi1 = np.square(np.real(np.sqrt(material.epsilon(omega))))
        chi1inv = (np.linalg.inv(chi1))
        chi1inv = np.diag(chi1inv)
        N = chi1inv.size
        n = np.sqrt(N/np.sum(chi1inv))
        
        sim = mp.Simulation(cell_size=mp.Vector3(2,2,2),
                            default_material=material,
                            resolution=20,
                            eps_averaging=False
                            )
        sim.init_sim()
        
        # Check to make sure the get_slice routine is working with omega
        n_slice = np.sqrt(np.max(sim.get_epsilon(omega)))
        self.assertAlmostEqual(n,n_slice, places=4)

        # Check to make sure h5 output is working with omega
        filename = 'dispersive_eigenmode-eps-000000.00.h5'
        mp.output_epsilon(sim,omega=omega)
        n_h5 = 0
        with h5py.File(filename, 'r') as f:
            n_h5 = np.sqrt(np.mean(f['eps'][()]))
        self.assertAlmostEqual(n,n_h5, places=4)
    
    # ----------------------------------------- #
    # ----------- Test Routines --------------- #
    # ----------------------------------------- #
    def test_chi1_routine(self):
        # Checks the newly implemented get_chi1inv routines within the 
        # fields and structure classes by comparing their output to the 
        # python epsilon output.

        from meep.materials import Si, Ag, LiNbO3, Au
        
        # Check Silicon
        w0 = Si.valid_freq_range.min
        w1 = Si.valid_freq_range.max
        self.call_chi1(Si,w0)
        self.call_chi1(Si,w1)

        # Check Silver
        w0 = Ag.valid_freq_range.min
        w1 = Ag.valid_freq_range.max
        self.call_chi1(Ag,w0)
        self.call_chi1(Ag,w1)

        # Check Gold
        w0 = Au.valid_freq_range.min
        w1 = Au.valid_freq_range.max
        self.call_chi1(Au,w0)
        self.call_chi1(Au,w1)   

        # Check Lithium Niobate (X,X)
        w0 = LiNbO3.valid_freq_range.min
        w1 = LiNbO3.valid_freq_range.max
        self.call_chi1(LiNbO3,w0)
        self.call_chi1(LiNbO3,w1)

        # Now let's rotate LN
        import copy
        rotLiNbO3 = copy.deepcopy(LiNbO3)
        rotLiNbO3.rotate([mp.Vector3(x=np.radians(28)),mp.Vector3(np.radians(45))])
        self.call_chi1(rotLiNbO3,w0)
        self.call_chi1(rotLiNbO3,w1)
    
    def test_get_with_dispersion(self):
        # Checks the get_array_slice and output_fields method
        # with dispersive materials.

        from meep.materials import Si, Ag, LiNbO3, Au
        
        # Check Silicon
        w0 = Si.valid_freq_range.min
        w1 = Si.valid_freq_range.max
        self.verify_output_and_slice(Si,w0)
        self.verify_output_and_slice(Si,w1)

        # Check Silver
        w0 = Ag.valid_freq_range.min
        w1 = Ag.valid_freq_range.max
        self.verify_output_and_slice(Ag,w0)
        self.verify_output_and_slice(Ag,w1)

        # Check Gold
        w0 = Au.valid_freq_range.min
        w1 = Au.valid_freq_range.max
        self.verify_output_and_slice(Au,w0)
        self.verify_output_and_slice(Au,w1)       

        # Check Lithium Niobate
        w0 = LiNbO3.valid_freq_range.min
        w1 = LiNbO3.valid_freq_range.max
        self.verify_output_and_slice(LiNbO3,w0)
        self.verify_output_and_slice(LiNbO3,w1)

        # Now let's rotate LN
        import copy
        rotLiNbO3 = copy.deepcopy(LiNbO3)
        rotLiNbO3.rotate([mp.Vector3(x=np.radians(28)),mp.Vector3(np.radians(45))])
        self.verify_output_and_slice(rotLiNbO3,w0)
        self.verify_output_and_slice(rotLiNbO3,w1)
        

if __name__ == '__main__':
    unittest.main()
