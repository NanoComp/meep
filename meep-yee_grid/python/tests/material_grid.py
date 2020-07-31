import meep as mp
import numpy as np
from meep.materials import Si, SiO2
import unittest

class TestMaterialGrid(unittest.TestCase):
    def gen_sim(self):
        resolution = 10 # pixels/um
        cell_size = mp.Vector3(14,14)
        pml_layers = [mp.PML(thickness=2)]
        w = 3.0 # width of waveguide
        m1 = SiO2
        m2 = Si
        n = 10
        gs = mp.Vector3(n,n)
        np.random.seed(1)
        dp = np.random.rand(n*n)
        mg = mp.MaterialGrid(gs,m1,m2,dp,"U_SUM")
        geometry = []
        rot_angle = np.radians(0)
        geometry += [mp.Block(center=mp.Vector3(),
                            size=mp.Vector3(w,w,mp.inf),
                            e1=mp.Vector3(x=1).rotate(mp.Vector3(z=1), rot_angle),
                            e2=mp.Vector3(y=1).rotate(mp.Vector3(z=1), rot_angle),
                            material=mg)]
        geometry += [mp.Block(center=mp.Vector3(),
                            size=mp.Vector3(w,w,mp.inf),
                            e1=mp.Vector3(x=-1).rotate(mp.Vector3(z=1), np.pi/2),
                            e2=mp.Vector3(y=1).rotate(mp.Vector3(z=1), np.pi/2),
                            material=mg)]
        fsrc = 1/1.55 # frequency of eigenmode or constant-amplitude source
        bnum = 1    # band number of eigenmode
        kpoint = mp.Vector3(x=1).rotate(mp.Vector3(z=1), rot_angle)
        compute_flux = True # compute flux (True) or plot the field profile (False)
        eig_src = True # eigenmode (True) or constant-amplitude (False) source
        if eig_src:
            sources = [mp.EigenModeSource(src=mp.GaussianSource(fsrc,fwidth=0.2*fsrc) if compute_flux else mp.ContinuousSource(fsrc),
                                        center=mp.Vector3(),
                                        size=mp.Vector3(y=3*w),
                                        direction=mp.NO_DIRECTION,
                                        eig_kpoint=kpoint,
                                        eig_band=bnum,
                                        eig_parity=mp.EVEN_Y+mp.ODD_Z if rot_angle == 0 else mp.ODD_Z,
                                        eig_match_freq=True)]
        else:
            sources = [mp.Source(src=mp.GaussianSource(fsrc,fwidth=0.2*fsrc) if compute_flux else mp.ContinuousSource(fsrc),
                                center=mp.Vector3(),
                                size=mp.Vector3(y=3*w),
                                component=mp.Ez)]

        sim = mp.Simulation(cell_size=cell_size,
                            resolution=resolution,
                            boundary_layers=pml_layers,
                            sources=sources,
                            extra_materials = [m1,m2],
                            geometry=geometry
                            )
        return sim
    
    def test_eval(self):
        sim = self.gen_sim()
if __name__ == '__main__':
    unittest.main()
