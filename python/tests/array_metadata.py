import unittest
import meep as mp
import numpy as np

class TestArrayMetadata(unittest.TestCase):

    def test_array_metadata(self):
        resolution = 25

        n = 3.4
        w = 1
        r = 1
        pad = 4
        dpml = 2

        sxy = 2*(r+w+pad+dpml)
        nonpml_vol = mp.Volume(mp.Vector3(), size=mp.Vector3(sxy-2*dpml,sxy-2*dpml))

        c1 = mp.Cylinder(radius=r+w, material=mp.Medium(index=n))
        c2 = mp.Cylinder(radius=r)

        fcen = 0.118
        df = 0.08

        # CW source
        src = mp.Source(mp.ContinuousSource(fcen,fwidth=df), mp.Ez, mp.Vector3(r+0.1))

        sim = mp.Simulation(cell_size=mp.Vector3(sxy,sxy),
                            geometry=[c1,c2],
                            sources=[src],
                            resolution=resolution,
                            force_complex_fields=True,
                            symmetries=[mp.Mirror(mp.Y)],
                            boundary_layers=[mp.PML(dpml)])

        sim.init_sim()
        sim.solve_cw(1e-6, 1000, 10)
        cw_modal_volume = sim.modal_volume_in_box(box=nonpml_vol)

        sim.reset_meep()

        # pulsed source
        src = mp.Source(mp.GaussianSource(fcen,fwidth=df), mp.Ez, mp.Vector3(r+0.1))

        sim = mp.Simulation(cell_size=mp.Vector3(sxy,sxy),
                            geometry=[c1,c2],
                            k_point=mp.Vector3(),
                            sources=[src],
                            resolution=resolution,
                            symmetries=[mp.Mirror(mp.Y)],
                            boundary_layers=[mp.PML(dpml)])

        dft_obj = sim.add_dft_fields([mp.Ez], fcen, fcen, 1, where=nonpml_vol)
        sim.run(until_after_sources=100)

        Ez = sim.get_dft_array(dft_obj,mp.Ez,0)
        (X,Y,Z,W) = sim.get_array_metadata(dft_cell=dft_obj)
        Eps = sim.get_array(vol=nonpml_vol,component=mp.Dielectric)
        EpsE2 = np.real(Eps*np.conj(Ez)*Ez)
        pulse_modal_volume = np.sum(W*EpsE2)/np.max(EpsE2)

        self.assertAlmostEqual(cw_modal_volume, pulse_modal_volume, places=1)

if __name__ == '__main__':
    unittest.main()
