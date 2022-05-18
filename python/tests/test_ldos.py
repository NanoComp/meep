import unittest
import numpy as np
import meep as mp

# Computes the Purcell enhancement factor of a horizontal dipole in a 3D
# homogeneous dielectric cavity with lossless metallic walls on two sides.
# The simulated result is validated using analytic theory from
# I. Abram et al., IEEE J. Quantum Electronics, Vol. 34, pp. 71-76 (1998).

class TestLDOS(unittest.TestCase):

    @classmethod
    def setUp(cls):
        cls.resolution = 20  # pixels/μm
        cls.dpml = 0.5       # thickness of PML
        cls.L = 6.0          # length of non-PML region
        cls.n = 2.4          # refractive index of surrounding medium
        cls.wvl = 1.0        # wavelength (in vacuum)

        cls.fcen = 1/cls.wvl
        cls.sources = [mp.Source(src=mp.GaussianSource(cls.fcen,fwidth=0.2*cls.fcen),
                                 component=mp.Ex,
                                 center=mp.Vector3())]

        cls.symmetries = [mp.Mirror(direction=mp.X,phase=-1),
                          mp.Mirror(direction=mp.Y),
                          mp.Mirror(direction=mp.Z)]


    def bulk_ldos(self):
        s = self.L+2*self.dpml
        cell_size = mp.Vector3(s,s,s)

        pml_layers = [mp.PML(self.dpml)]

        sim = mp.Simulation(resolution=self.resolution,
                            cell_size=cell_size,
                            boundary_layers=pml_layers,
                            sources=self.sources,
                            symmetries=self.symmetries,
                            default_material=mp.Medium(index=self.n))

        sim.run(mp.dft_ldos(self.fcen,0,1),
                until_after_sources=mp.stop_when_fields_decayed(20,
                                                                mp.Ex,
                                                                mp.Vector3(),
                                                                1e-6))

        return sim.ldos_data[0]


    def cavity_ldos(self,sz):
        sxy = self.L+2*self.dpml
        cell_size = mp.Vector3(sxy,sxy,sz)

        boundary_layers = [mp.PML(self.dpml,direction=mp.X),
                           mp.PML(self.dpml,direction=mp.Y)]

        sim = mp.Simulation(resolution=self.resolution,
                            cell_size=cell_size,
                            boundary_layers=boundary_layers,
                            sources=self.sources,
                            symmetries=self.symmetries,
                            default_material=mp.Medium(index=self.n))

        sim.run(mp.dft_ldos(ldos=mp.Ldos(self.fcen,0,1)),
                until_after_sources=mp.stop_when_fields_decayed(20,
                                                                mp.Ex,
                                                                mp.Vector3(),
                                                                1e-6))

        return sim.ldos_data[0]


    def test_ldos(self):
        ldos_bulk = self.bulk_ldos()
        print("ldos_bulk:, {:.6f}".format(ldos_bulk))

        cavity_thickness = 0.75
        gap = cavity_thickness*self.wvl/self.n

        ldos_cavity = self.cavity_ldos(gap)

        # Purcell enhancement factor (relative to bulk medium)
        pe_meep = ldos_cavity/ldos_bulk

        # equation 7 of reference
        pe_theory = (3*np.fix(cavity_thickness+0.5)/(4*cavity_thickness) +
                     (4*np.power(np.fix(cavity_thickness+0.5),3) -
                      np.fix(cavity_thickness+0.5))/(16*np.power(cavity_thickness,3)))

        rel_err = abs(pe_meep-pe_theory)/pe_theory

        print("ldos:, {:.6f} (Meep), {:.6f} (theory), "
              "{:.6f} (error)".format(pe_meep,pe_theory,rel_err))

        self.assertAlmostEqual(pe_meep, pe_theory, delta=0.1)


if __name__ == '__main__':
    unittest.main()
