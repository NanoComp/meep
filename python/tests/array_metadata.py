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
        cell_size = mp.Vector3(sxy,sxy)

        nonpml_vol = mp.Volume(mp.Vector3(), size=mp.Vector3(sxy-2*dpml,sxy-2*dpml))

        geometry = [mp.Cylinder(radius=r+w, material=mp.Medium(index=n)),
                    mp.Cylinder(radius=r)]

        fcen = 0.118
        df = 0.08

        symmetries = [mp.Mirror(mp.X,phase=-1),
                      mp.Mirror(mp.Y,phase=+1)]

        pml_layers = [mp.PML(dpml)]

        # CW source
        src = [mp.Source(mp.ContinuousSource(fcen,fwidth=df), mp.Ez, mp.Vector3(r+0.1)),
               mp.Source(mp.ContinuousSource(fcen,fwidth=df), mp.Ez, mp.Vector3(-(r+0.1)), amplitude=-1)]

        sim = mp.Simulation(cell_size=cell_size,
                            geometry=geometry,
                            sources=src,
                            resolution=resolution,
                            force_complex_fields=True,
                            symmetries=symmetries,
                            boundary_layers=pml_layers)

        sim.init_sim()
        sim.solve_cw(1e-6, 1000, 10)

        def electric_energy(r, ex, ey, ez, eps):
            return np.real(eps * (np.conj(ex)*ex + np.conj(ey)*ey + np.conj(ez)*ez))

        electric_energy_total = sim.integrate_field_function([mp.Ex,mp.Ey,mp.Ez,mp.Dielectric],electric_energy,nonpml_vol)
        electric_energy_max = sim.max_abs_field_function([mp.Ex,mp.Ey,mp.Ez,mp.Dielectric],electric_energy,nonpml_vol)
        cw_modal_volume = electric_energy_total / electric_energy_max

        sim.reset_meep()

        # pulsed source
        src = [mp.Source(mp.GaussianSource(fcen,fwidth=df), mp.Ez, mp.Vector3(r+0.1)),
               mp.Source(mp.GaussianSource(fcen,fwidth=df), mp.Ez, mp.Vector3(-(r+0.1)), amplitude=-1)]

        sim = mp.Simulation(cell_size=cell_size,
                            geometry=geometry,
                            k_point=mp.Vector3(),
                            sources=src,
                            resolution=resolution,
                            symmetries=symmetries,
                            boundary_layers=pml_layers)

        dft_obj = sim.add_dft_fields([mp.Ez], fcen, fcen, 1, where=nonpml_vol)
        sim.run(until_after_sources=100)

        (Ex,Ey,Ez) = [sim.get_dft_array(dft_obj, c, 0) for c in [mp.Ex, mp.Ey, mp.Ez]]
        (X,Y,Z,W) = sim.get_array_metadata(dft_cell=dft_obj)
        Eps = sim.get_array(vol=nonpml_vol,component=mp.Dielectric)
        EpsE2 = np.real(Eps*(np.conj(Ex)*Ex + np.conj(Ey)*Ey + np.conj(Ez)*Ez))
        pulse_modal_volume = np.sum(W*EpsE2)/np.max(EpsE2)

        self.assertAlmostEqual(cw_modal_volume, pulse_modal_volume, places=1)

if __name__ == '__main__':
    unittest.main()
