
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


class TestDispersiveEigenmode(unittest.TestCase):

    def call_chi1(self,material,component,direction,omega):

        sim = mp.Simulation(cell_size=mp.Vector3(1,1,1),
                    default_material=material,
                    resolution=10)

        sim.init_sim()
        print(component, direction)
        v3 = mp.py_v3_to_vec(sim.dimensions, mp.Vector3(0,0,0), sim.is_cylindrical)
        n = 1/np.sqrt(sim.fields.get_chi1inv(int(component),int(direction),v3,omega))
        return n

    def test_chi1_routine(self):
        from meep.materials import Si, Ag, LiNbO3, fused_quartz
        # Check Silicon
        w0 = Si.valid_freq_range.min
        w1 = Si.valid_freq_range.max
        self.assertAlmostEqual(np.real(np.sqrt(Si.epsilon([w0,w1])[0,0,0])), self.call_chi1(Si,0,0,w0), places=6)
        self.assertAlmostEqual(np.real(np.sqrt(Si.epsilon([w0,w1])[1,0,0])), self.call_chi1(Si,0,0,w1), places=6)

        # Check Silver
        w0 = Ag.valid_freq_range.min
        w1 = Ag.valid_freq_range.max
        self.assertAlmostEqual(np.real(np.sqrt(Ag.epsilon([w0,w1])[0,0,0])), self.call_chi1(Ag,0,0,w0), places=6)
        self.assertAlmostEqual(np.real(np.sqrt(Ag.epsilon([w0,w1])[1,0,0])), self.call_chi1(Ag,0,0,w1), places=6)

        # Check Lithium Niobate
        w0 = LiNbO3.valid_freq_range.min
        w1 = LiNbO3.valid_freq_range.max
        self.assertAlmostEqual(np.real(np.sqrt(LiNbO3.epsilon([w0,w1])[0,0,0])), self.call_chi1(LiNbO3,0,0,w0), places=6)
        self.assertAlmostEqual(np.real(np.sqrt(LiNbO3.epsilon([w0,w1])[1,0,0])), self.call_chi1(LiNbO3,0,0,w1), places=6)
        self.assertAlmostEqual(np.real(np.sqrt(LiNbO3.epsilon([w0,w1])[0,2,2])), self.call_chi1(LiNbO3,mp.Z,mp.Z,w0), places=6)
        self.assertAlmostEqual(np.real(np.sqrt(LiNbO3.epsilon([w0,w1])[1,2,2])), self.call_chi1(LiNbO3,mp.Z,mp.Z,w1), places=6)

    def test_eigenmode_source(self):
        from meep.materials import Si, Ag, LiNbO3, fused_quartz
    
    def test_eigenmode_decomposition(self):
        from meep.materials import Si, Ag, LiNbO3, fused_quartz
    
    def test_get_eigenmode(self):
        from meep.materials import Si, Ag, LiNbO3, fused_quartz
    
    def test_everything(self):
        from meep.materials import Si, Ag, LiNbO3, fused_quartz
        

if __name__ == '__main__':
    unittest.main()
