import unittest
import meep as mp


class TestFragmentStats(unittest.TestCase):

    def get_fragment_stats(self, block_size, cell_size, dims):
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

        geom = [mp.Block(size=block_size, material=mat)]
        sim = mp.Simulation(cell_size=cell_size, resolution=10, geometry=geom, dimensions=dims)
        gv = sim._create_grid_volume(False)

        stats = mp.compute_fragment_stats(
            geom,
            gv,
            sim.cell_size,
            sim.default_material,
            sim.subpixel_tol,
            sim.subpixel_maxeval,
            sim.ensure_periodicity
        )
        return stats

    def test_1d_stats(self):
        fs = self.get_fragment_stats(mp.Vector3(z=10), mp.Vector3(z=30), 1)

        self.assertEqual(len(fs), 3)

        # First and last boxes have no geometry, only default_material
        for i in [0, 2]:
            # TODO: Test default_material
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
        # TODO: DFT

    def test_2d_stats(self):
        fs = self.get_fragment_stats(mp.Vector3(10, 10), mp.Vector3(30, 30), 2)

        self.assertEqual(len(fs), 9)

        # All boxes besides the middle one have no geometry, only default_material
        for i in [0, 1, 2, 3, 5, 6, 7, 8]:
            # TODO: Test default_material
            self.assertEqual(fs[i].num_anisotropic_eps_pixels, 0)
            self.assertEqual(fs[i].num_anisotropic_mu_pixels, 0)
            self.assertEqual(fs[i].num_nonlinear_pixels, 0)
            self.assertEqual(fs[i].num_susceptibility_pixels, 0)
            self.assertEqual(fs[i].num_nonzero_conductivity_pixels, 0)

        # Middle box contains entire block
        idx = 4
        self.assertEqual(fs[idx].num_anisotropic_eps_pixels, 1000)
        self.assertEqual(fs[idx].num_anisotropic_mu_pixels, 1000)
        self.assertEqual(fs[idx].num_nonlinear_pixels, 3000)
        self.assertEqual(fs[idx].num_susceptibility_pixels, 3000)
        self.assertEqual(fs[idx].num_nonzero_conductivity_pixels, 3000)
        # TODO: DFT

    def test_3d_stats(self):
        fs = self.get_fragment_stats(mp.Vector3(10, 10, 10), mp.Vector3(30, 30, 30), 3)

        self.assertEqual(len(fs), 27)

        # All boxes besides the middle one have no geometry, only default_material
        for i in range(27):
            if i == 13:
                continue
            # TODO: Test default_material
            self.assertEqual(fs[i].num_anisotropic_eps_pixels, 0)
            self.assertEqual(fs[i].num_anisotropic_mu_pixels, 0)
            self.assertEqual(fs[i].num_nonlinear_pixels, 0)
            self.assertEqual(fs[i].num_susceptibility_pixels, 0)
            self.assertEqual(fs[i].num_nonzero_conductivity_pixels, 0)

        # Middle box contains entire block
        idx = 13
        self.assertEqual(fs[idx].num_anisotropic_eps_pixels, 10000)
        self.assertEqual(fs[idx].num_anisotropic_mu_pixels, 10000)
        self.assertEqual(fs[idx].num_nonlinear_pixels, 30000)
        self.assertEqual(fs[idx].num_susceptibility_pixels, 30000)
        self.assertEqual(fs[idx].num_nonzero_conductivity_pixels, 30000)
        # TODO: DFT

    def test_cyl_stats(self):
        pass


if __name__ == '__main__':
    unittest.main()
