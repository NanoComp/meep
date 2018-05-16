import unittest
import meep as mp


class TestFragmentStats(unittest.TestCase):

    def check_stats(self, fragment, a_eps, a_mu, nonlin, susc, cond):
        self.assertEqual(fragment.num_anisotropic_eps_pixels, a_eps)
        self.assertEqual(fragment.num_anisotropic_mu_pixels, a_mu)
        self.assertEqual(fragment.num_nonlinear_pixels, nonlin)
        self.assertEqual(fragment.num_susceptibility_pixels, susc)
        self.assertEqual(fragment.num_nonzero_conductivity_pixels, cond)

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
        # A z=30 cell, split into three fragments of size 10 each, with a block
        # covering the middle fragment.

        fs = self.get_fragment_stats(mp.Vector3(z=10), mp.Vector3(z=30), 1)

        self.assertEqual(len(fs), 3)

        # First and last fragments have no geometry, only default_material
        for i in [0, 2]:
            # TODO: Test default_material
            self.check_stats(fs[i], 0, 0, 0, 0, 0)

        # Second fragment contains entire block
        self.check_stats(fs[1], a_eps=100, a_mu=100, nonlin=300, susc=300, cond=300)
        # TODO: DFT

    def test_1d_stats_with_overlap(self):
        # A z=30 cell split into three fragments of size 10 each, with a block
        # covering the middle fragment, and half of the two outer fragments.

        fs = self.get_fragment_stats(mp.Vector3(z=20), mp.Vector3(z=30), 1)

        self.assertEqual(len(fs), 3)

        # Middle fragment is completely covered by the block
        self.check_stats(fs[1], a_eps=100, a_mu=100, nonlin=300, susc=300, cond=300)

        # Outer two fragments are half covered by the block
        for i in [0, 2]:
            self.check_stats(fs[i], a_eps=50, a_mu=50, nonlin=150, susc=150, cond=150)

    def test_2d_stats(self):
        # A 30 x 30 cell, with a 10 x 10 block in the middle, split into 9 10 x 10 fragments.

        fs = self.get_fragment_stats(mp.Vector3(10, 10), mp.Vector3(30, 30), 2)

        self.assertEqual(len(fs), 9)

        # All fragments besides the middle one have no geometry, only default_material
        for i in [0, 1, 2, 3, 5, 6, 7, 8]:
            # TODO: Test default_material
            self.check_stats(fs[i], 0, 0, 0, 0, 0)

        # Middle fragment contains entire block
        idx = 4
        self.check_stats(fs[idx], a_eps=1000, a_mu=1000, nonlin=3000, susc=3000, cond=3000)
        # TODO: DFT

    def test_2d_stats_with_overlap(self):
        # A 30 x 30 cell, with a 20 x 20 block in the middle, split into 9 10 x 10 fragments.

        fs = self.get_fragment_stats(mp.Vector3(20, 20), mp.Vector3(30, 30), 2)

        self.assertEqual(len(fs), 9)

        # Middle fragment contains entire block
        idx = 4
        self.check_stats(fs[idx], a_eps=1000, a_mu=1000, nonlin=3000, susc=3000, cond=3000)

        # Top-middle, bottom-middle, left-middle, and right-middle fragments are half
        # covered by the block.
        for i in [1, 3, 5, 7]:
            self.check_stats(fs[i], a_eps=500, a_mu=500, nonlin=1500, susc=1500, cond=1500)

        # The four corner fragments are quarter-filled by the block
        for i in [0, 2, 6, 8]:
            self.check_stats(fs[i], a_eps=250, a_mu=250, nonlin=750, susc=750, cond=750)

    def test_3d_stats(self):
        # A 30 x 30 x 30 cell with a 10 x 10 x 10 block placed at the center, split
        # into 27 10 x 10 x 10 fragments

        fs = self.get_fragment_stats(mp.Vector3(10, 10, 10), mp.Vector3(30, 30, 30), 3)

        self.assertEqual(len(fs), 27)

        # All fragments besides the middle one have no geometry, only default_material
        for i in range(27):
            if i == 13:
                continue
            # TODO: Test default_material
            self.check_stats(fs[i], 0, 0, 0, 0, 0)

        # Middle fragments contains entire block
        idx = 13
        self.check_stats(fs[idx], a_eps=10000, a_mu=10000, nonlin=30000, susc=30000, cond=30000)
        # TODO: DFT

    def test_3d_stats_with_overlap(self):
        # A 30 x 30 x 30 cell with a 20 x 20 x 20 block placed at the center, split
        # into 27 10 x 10 x 10 fragments

        fs = self.get_fragment_stats(mp.Vector3(20, 20, 20), mp.Vector3(30, 30, 30), 3)

        self.assertEqual(len(fs), 27)

        # Middle fragment contains entire block
        idx = 13
        self.check_stats(fs[idx], a_eps=10000, a_mu=10000, nonlin=30000, susc=30000, cond=30000)

        # Six fragments adjacent to the middle fragment faces will be half covered by the block
        for i in [4, 10, 12, 14, 16, 22]:
            self.check_stats(fs[i], a_eps=5000, a_mu=5000, nonlin=15000, susc=15000, cond=15000)

        # The corners will be 1/8 covered
        for i in [0, 2, 6, 8, 18, 20, 24, 26]:
            self.check_stats(fs[i], a_eps=1250, a_mu=1250, nonlin=3750, susc=3750, cond=3750)

        # The rest will be 1/4 covered
        for i in [1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25]:
            self.check_stats(fs[i], a_eps=2500, a_mu=2500, nonlin=7500, susc=7500, cond=7500)

    def test_cyl_stats(self):
        pass


if __name__ == '__main__':
    unittest.main()
