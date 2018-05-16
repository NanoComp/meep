import unittest
import meep as mp


class TestFragmentStats(unittest.TestCase):

    def test_1d_stats(self):

        mat = mp.Medium(
            epsilon=12,
            epsilon_offdiag=mp.Vector3(z=1),
            mu_offdiag=mp.Vector3(x=20),
            E_chi2_diag=mp.Vector3(1, 1),
            H_chi3_diag=mp.Vector3(z=1),
            E_susceptibilities=[mp.LorentzianSusceptibility(), mp.NoisyLorentzianSusceptibility()],
            H_susceptibilities=[mp.DrudeSusceptibility()],
            D_conductivity_diag=mp.Vector3(y=1),
            B_conductivity_diag=mp.Vector3(x=1, z=1)
        )

        geom = [mp.Block(size=mp.Vector3(z=10), material=mat)]
        sim = mp.Simulation(cell_size=mp.Vector3(z=30), resolution=10, geometry=geom, dimensions=1)
        gv = sim._create_grid_volume(False)
        fs = mp.compute_fragment_stats(geom, gv, sim.cell_size, sim.default_material,
                                       sim.subpixel_tol, sim.subpixel_maxeval, sim.ensure_periodicity)

        self.assertEqual(len(fs), 3)

        # First and last boxes have no geometry, only default_material
        for i in [0, 2]:
            self.assertEqual(fs[i].num_anisotropic_eps_pixels, 0)
            self.assertEqual(fs[i].num_anisotropic_mu_pixels, 0)
            self.assertEqual(fs[i].num_nonlinear_pixels, 0)
            self.assertEqual(fs[i].num_susceptibility_pixels, 0)
            self.assertEqual(fs[i].num_nonzero_conductivity_pixels, 0)

        # Second box contains entire block
        self.assertEqual(fs[1].num_anisotropic_eps_pixels, 100)
        self.assertEqual(fs[1].num_anisotropic_mu_pixels, 100)
        self.assertEqual(fs[1].num_nonlinear_pixels, 300)
        self.assertEqual(fs[1].num_susceptibility_pixels, 300)
        self.assertEqual(fs[1].num_nonzero_conductivity_pixels, 300)

    def test_2d_stats(self):
        pass

    def test_3d_stats(self):
        pass

    def test_cyl_stats(self):
        pass


if __name__ == '__main__':
    unittest.main()
