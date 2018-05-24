import unittest
import meep as mp


def make_dft_vecs(flx_reg=None, n2f_reg=None, frc_reg=None, fldc=None, flds=None, fldw=None, fld_cmp=None):
    dft_vecs = {
        'flux_regions': flx_reg,
        'n2f_regions': n2f_reg,
        'force_regions': frc_reg,
        'fields_center': fldc,
        'fields_size': flds,
        'fields_where': fldw,
        'fields_components': fld_cmp
    }
    return dft_vecs


class TestFragmentStats(unittest.TestCase):

    def check_stats(self, fragment, a_eps, a_mu, nonlin, susc, cond):
        self.assertEqual(fragment.num_anisotropic_eps_pixels, a_eps)
        self.assertEqual(fragment.num_anisotropic_mu_pixels, a_mu)
        self.assertEqual(fragment.num_nonlinear_pixels, nonlin)
        self.assertEqual(fragment.num_susceptibility_pixels, susc)
        self.assertEqual(fragment.num_nonzero_conductivity_pixels, cond)

    def get_fragment_stats(self, block_size, cell_size, dims, box_center=mp.Vector3(), dft_vecs=None,
                           def_mat=mp.air):
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

        geom = [mp.Block(size=block_size, center=box_center, material=mat)]
        sim = mp.Simulation(cell_size=cell_size, resolution=10, geometry=geom, dimensions=dims,
                            default_material=def_mat)

        if dft_vecs:
            if dft_vecs['flux_regions']:
                sim.add_flux(1, 0.5, 5, *dft_vecs['flux_regions'])
            if dft_vecs['n2f_regions']:
                sim.add_near2far(1, 0.5, 7, *dft_vecs['n2f_regions'])
            if dft_vecs['force_regions']:
                sim.add_force(1, 0.5, 9, *dft_vecs['force_regions'])
            if dft_vecs['fields_components']:
                sim.add_dft_fields(dft_vecs['fields_components'], 0, 1, 5, where=dft_vecs['fields_where'],
                                   center=dft_vecs['fields_center'], size=dft_vecs['fields_size'])

        gv = sim._create_grid_volume(False)
        stats = sim._compute_fragment_stats(gv)

        return stats

    def test_1d(self):
        # A z=30 cell, split into three fragments of size 10 each, with a block
        # covering the middle fragment.

        # flux covering first fragment, near2far covering second, force covering third
        dft_vecs = make_dft_vecs(
            [mp.FluxRegion(mp.Vector3(z=-10), size=mp.Vector3(z=10))],
            [mp.Near2FarRegion(mp.Vector3(), size=mp.Vector3(z=10))],
            [mp.ForceRegion(mp.Vector3(z=10), mp.X, size=mp.Vector3(z=10))]
        )

        fs = self.get_fragment_stats(mp.Vector3(z=10), mp.Vector3(z=30), 1, dft_vecs=dft_vecs)

        self.assertEqual(len(fs), 3)

        # First and last fragments have no geometry, only default_material
        for i in [0, 2]:
            self.check_stats(fs[i], 0, 0, 0, 0, 0)

        # Second fragment contains entire block
        self.check_stats(fs[1], a_eps=100, a_mu=100, nonlin=300, susc=300, cond=300)

        # Check DFT regions
        self.assertEqual(fs[0].num_dft_pixels, 2000)
        self.assertEqual(fs[1].num_dft_pixels, 2800)
        self.assertEqual(fs[2].num_dft_pixels, 5400)

    def test_1d_with_overlap(self):
        # A z=30 cell split into three fragments of size 10 each, with a block
        # covering the middle fragment, and half of the two outer fragments.

        mat = mp.Medium(H_susceptibilities=[mp.DrudeSusceptibility()])
        fs = self.get_fragment_stats(mp.Vector3(z=20), mp.Vector3(z=30), 1, def_mat=mat)

        self.assertEqual(len(fs), 3)

        # Middle fragment is completely covered by the block
        self.check_stats(fs[1], a_eps=100, a_mu=100, nonlin=300, susc=300, cond=300)

        # Outer two fragments are half covered by the block, and half covered by default_material 'mat'
        for i in [0, 2]:
            self.check_stats(fs[i], a_eps=50, a_mu=50, nonlin=150, susc=200, cond=150)

    def test_1d_with_partial_fragment(self):
        # A cell with z=26, with a 16 unit block in the center, split into 3 fragments,
        # with the first and last fragment of length 8, and 3/8 covered by the block,
        # and the middle fragment completely covered.

        # dft_flux with 2 volumes, 1 covering the first fragment and one covering
        # half of the second fragment
        dft_vecs = make_dft_vecs(flx_reg=[
            mp.FluxRegion(mp.Vector3(z=-9), mp.Vector3(z=8)),
            mp.FluxRegion(mp.Vector3(z=-2.5), mp.Vector3(z=5))
        ])
        fs = self.get_fragment_stats(mp.Vector3(z=16), mp.Vector3(z=26), 1, dft_vecs=dft_vecs)

        self.assertEqual(len(fs), 3)
        # Check first and last box sizes
        self.assertEqual(fs[0].box.low.z, -13)
        self.assertEqual(fs[0].box.high.z, -5)
        self.assertEqual(fs[2].box.low.z, 5)
        self.assertEqual(fs[2].box.high.z, 13)

        # Middle fragment is completely covered by block
        self.check_stats(fs[1], a_eps=100, a_mu=100, nonlin=300, susc=300, cond=300)
        # Outer fragments are 3/8 covered
        for i in [0, 2]:
            self.check_stats(fs[i], a_eps=30, a_mu=30, nonlin=90, susc=90, cond=90)

        # Check dft stats
        self.assertEqual(fs[0].num_dft_pixels, 1600)
        self.assertEqual(fs[1].num_dft_pixels, 1000)
        self.assertEqual(fs[2].num_dft_pixels, 0)

    def test_1d_with_shifted_center(self):
        # A cell with z=26, with a 16 unit block shifted so that the right side is flush,
        # with the right side of the cell, split into 3 fragments, with the first and last
        # fragments of length 8, the first uncovered and the last completely covered by the
        # block, and the middle fragment 80% covered by the block.
        fs = self.get_fragment_stats(mp.Vector3(z=16), mp.Vector3(z=26), 1, mp.Vector3(z=5))

        self.assertEqual(len(fs), 3)
        # Check first and last box sizes
        self.assertEqual(fs[0].box.low.z, -13)
        self.assertEqual(fs[0].box.high.z, -5)
        self.assertEqual(fs[2].box.low.z, 5)
        self.assertEqual(fs[2].box.high.z, 13)

        # First fragment is uncovered
        self.check_stats(fs[0], 0, 0, 0, 0, 0)

        # Middel fragment is 80% covered
        self.check_stats(fs[1], a_eps=80, a_mu=80, nonlin=240, susc=240, cond=240)

        # Last fragment is completely covered, but only 8 units long
        self.check_stats(fs[2], a_eps=80, a_mu=80, nonlin=240, susc=240, cond=240)

    def test_1d_dft_fields(self):
        # A z=30 cell, split into three fragments of size 10 each, with a block
        # covering the middle fragment.

        # dft_fields covering first fragment
        dft_vecs = make_dft_vecs(fldc=mp.Vector3(z=-10), flds=mp.Vector3(z=10), fld_cmp=[mp.X, mp.Y])
        fs = self.get_fragment_stats(mp.Vector3(z=10), mp.Vector3(z=30), 1, dft_vecs=dft_vecs)

        self.assertEqual(len(fs), 3)

        for i in [1, 2]:
            self.assertEqual(fs[i].num_dft_pixels, 0)

        self.assertEqual(fs[0].num_dft_pixels, 1000)

        # Same test with volume instead of center and size
        dft_vecs = make_dft_vecs(fldw=mp.Volume(mp.Vector3(z=-10), mp.Vector3(z=10)), fld_cmp=[mp.X, mp.Y])
        fs = self.get_fragment_stats(mp.Vector3(z=10), mp.Vector3(z=30), 1, dft_vecs=dft_vecs)

        for i in [1, 2]:
            self.assertEqual(fs[i].num_dft_pixels, 0)

        self.assertEqual(fs[0].num_dft_pixels, 1000)

    def test_2d(self):
        # A 30 x 30 cell, with a 10 x 10 block in the middle, split into 9 10 x 10 fragments.

        # flux covering top-left fragment, near2far covering top-middle, force covering top-right
        dft_vecs = make_dft_vecs(
            [mp.FluxRegion(mp.Vector3(-10, 10), size=mp.Vector3(10, 10))],
            [mp.Near2FarRegion(mp.Vector3(0, 10), size=mp.Vector3(10, 10))],
            [mp.ForceRegion(mp.Vector3(10, 10), mp.X, size=mp.Vector3(10, 10))]
        )
        fs = self.get_fragment_stats(mp.Vector3(10, 10), mp.Vector3(30, 30), 2, dft_vecs=dft_vecs)

        self.assertEqual(len(fs), 9)

        # Check fragment boxes
        self.assertEqual(fs[0].box.low.x, -15)
        self.assertEqual(fs[0].box.low.y, -15)
        self.assertEqual(fs[0].box.high.x, -5)
        self.assertEqual(fs[0].box.high.y, -5)

        self.assertEqual(fs[1].box.low.x, -15)
        self.assertEqual(fs[1].box.low.y, -5)
        self.assertEqual(fs[1].box.high.x, -5)
        self.assertEqual(fs[1].box.high.y, 5)

        self.assertEqual(fs[2].box.low.x, -15)
        self.assertEqual(fs[2].box.low.y, 5)
        self.assertEqual(fs[2].box.high.x, -5)
        self.assertEqual(fs[2].box.high.y, 15)

        # All fragments besides the middle one have no geometry, only default_material
        for i in [0, 1, 2, 3, 5, 6, 7, 8]:
            self.check_stats(fs[i], 0, 0, 0, 0, 0)

        # Middle fragment contains entire block
        idx = 4
        self.check_stats(fs[idx], a_eps=1000, a_mu=1000, nonlin=3000, susc=3000, cond=3000)

        # Check DFT regions
        for i in [0, 1, 3, 4, 6, 7]:
            self.assertEqual(fs[i].num_dft_pixels, 0)

        self.assertEqual(fs[2].num_dft_pixels, 20000)
        self.assertEqual(fs[5].num_dft_pixels, 28000)
        self.assertEqual(fs[8].num_dft_pixels, 54000)

    def test_2d_with_overlap(self):
        # A 30 x 30 cell, with a 20 x 20 block in the middle, split into 9 10 x 10 fragments.

        mat = mp.Medium(H_susceptibilities=[mp.DrudeSusceptibility()])
        fs = self.get_fragment_stats(mp.Vector3(20, 20), mp.Vector3(30, 30), 2, def_mat=mat)

        self.assertEqual(len(fs), 9)

        # Middle fragment contains entire block
        idx = 4
        self.check_stats(fs[idx], a_eps=1000, a_mu=1000, nonlin=3000, susc=3000, cond=3000)

        # Top-middle, bottom-middle, left-middle, and right-middle fragments are half
        # covered by the block, and half covered by default_material 'mat'.
        for i in [1, 3, 5, 7]:
            self.check_stats(fs[i], a_eps=500, a_mu=500, nonlin=1500, susc=2000, cond=1500)

        # The four corner fragments are quarter-filled by the block, and 3/4 filled by
        # default_material 'mat'
        for i in [0, 2, 6, 8]:
            self.check_stats(fs[i], a_eps=250, a_mu=250, nonlin=750, susc=1500, cond=750)

    def test_2d_with_partial_fragments_and_shifted_center(self):
        # A 26 x 26 cell with a 18 x 18 Block in the lower right corner

        fs = self.get_fragment_stats(mp.Vector3(18, 18), mp.Vector3(26, 26), 2, mp.Vector3(4, -4))

        self.assertEqual(len(fs), 9)

        # Middle fragment is 10 x 10 and covered by block
        self.check_stats(fs[4], a_eps=1000, a_mu=1000, nonlin=3000, susc=3000, cond=3000)

        for i in [0, 1, 2, 5, 8]:
            # Air
            self.check_stats(fs[i], 0, 0, 0, 0, 0)

        # 10 x 8 fragment, covered by block
        self.assertEqual(fs[3].box.low.x, -5)
        self.assertEqual(fs[3].box.low.y, -13)
        self.assertEqual(fs[3].box.high.x, 5)
        self.assertEqual(fs[3].box.high.y, -5)
        self.check_stats(fs[3], a_eps=800, a_mu=800, nonlin=2400, susc=2400, cond=2400)

        # 8 x 10 fragment, covered by block
        self.assertEqual(fs[7].box.low.x, 5)
        self.assertEqual(fs[7].box.low.y, -5)
        self.assertEqual(fs[7].box.high.x, 13)
        self.assertEqual(fs[7].box.high.y, 5)
        self.check_stats(fs[7], a_eps=800, a_mu=800, nonlin=2400, susc=2400, cond=2400)

        # 8 x 8 fragment covered by block
        self.assertEqual(fs[6].box.low.x, 5)
        self.assertEqual(fs[6].box.low.y, -13)
        self.assertEqual(fs[6].box.high.x, 13)
        self.assertEqual(fs[6].box.high.y, -5)
        self.check_stats(fs[6], a_eps=640, a_mu=640, nonlin=1920, susc=1920, cond=1920)

    def test_2d_dft_fields(self):
        # A 30 x 30 cell, with a 10 x 10 block in the middle, split into 9 10 x 10 fragments.

        # dft_fields covering 20 by 20 area in center of cell. Test with volume, and center/size
        cmpts = [mp.Ex, mp.Ey, mp.Ez]
        dft_fields_size_center = make_dft_vecs(fldc=mp.Vector3(), flds=mp.Vector3(20, 20), fld_cmp=cmpts)
        dft_fields_where = make_dft_vecs(fldw=mp.Volume(mp.Vector3(), mp.Vector3(20, 20)), fld_cmp=cmpts)

        for dft_vec in [dft_fields_size_center, dft_fields_where]:
            fs = self.get_fragment_stats(mp.Vector3(10, 10), mp.Vector3(30, 30), 2, dft_vecs=dft_vec)

            # Middle fragment is fully covered
            self.assertEqual(fs[4].num_dft_pixels, 15000)

            # 4 corners are 1/4 covered
            for i in [0, 2, 6, 8]:
                self.assertEqual(fs[i].num_dft_pixels, 3750)

            # The rest are half covered
            for i in [1, 3, 5, 7]:
                self.assertEqual(fs[i].num_dft_pixels, 7500)

    def test_3d(self):
        # A 30 x 30 x 30 cell with a 10 x 10 x 10 block placed at the center, split
        # into 27 10 x 10 x 10 fragments

        # flux covering lower-front-left fragment, near2far covering lower-middle-left,
        # force covering lower-back-left
        dft_vecs = make_dft_vecs(
            [mp.FluxRegion(mp.Vector3(-10, -10, -10), size=mp.Vector3(10, 10, 10))],
            [mp.Near2FarRegion(mp.Vector3(-10, -10, 0), size=mp.Vector3(10, 10, 10))],
            [mp.ForceRegion(mp.Vector3(-10, -10, 10), mp.X, size=mp.Vector3(10, 10, 10))]
        )
        fs = self.get_fragment_stats(mp.Vector3(10, 10, 10), mp.Vector3(30, 30, 30), 3, dft_vecs=dft_vecs)

        self.assertEqual(len(fs), 27)

        # All fragments besides the middle one have no geometry, only default_material
        for i in range(27):
            if i == 13:
                continue
            self.check_stats(fs[i], 0, 0, 0, 0, 0)

        # Middle fragments contains entire block
        idx = 13
        self.check_stats(fs[idx], a_eps=10000, a_mu=10000, nonlin=30000, susc=30000, cond=30000)

        # Check DFT regions
        for i in range(3, 27):
            self.assertEqual(fs[i].num_dft_pixels, 0)

        self.assertEqual(fs[0].num_dft_pixels, 200000)
        self.assertEqual(fs[1].num_dft_pixels, 280000)
        self.assertEqual(fs[2].num_dft_pixels, 540000)

    def test_3d_with_overlap(self):
        # A 30 x 30 x 30 cell with a 20 x 20 x 20 block placed at the center, split
        # into 27 10 x 10 x 10 fragments

        mat = mp.Medium(E_susceptibilities=[mp.DrudeSusceptibility()])
        fs = self.get_fragment_stats(mp.Vector3(20, 20, 20), mp.Vector3(30, 30, 30), 3, def_mat=mat)

        self.assertEqual(len(fs), 27)

        # Middle fragment contains entire block
        idx = 13
        self.check_stats(fs[idx], a_eps=10000, a_mu=10000, nonlin=30000, susc=30000, cond=30000)

        # Six fragments adjacent to the middle fragment faces will be half covered by the block,
        # and half covered by default_material 'mat'
        for i in [4, 10, 12, 14, 16, 22]:
            self.check_stats(fs[i], a_eps=5000, a_mu=5000, nonlin=15000, susc=20000, cond=15000)

        # The corners will be 1/8 covered by the block and 7/8 covered by default_material 'mat'
        for i in [0, 2, 6, 8, 18, 20, 24, 26]:
            self.check_stats(fs[i], a_eps=1250, a_mu=1250, nonlin=3750, susc=12500, cond=3750)

        # The rest will be 1/4 covered by the block and 3/4 covered by default_material 'mat'
        for i in [1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25]:
            self.check_stats(fs[i], a_eps=2500, a_mu=2500, nonlin=7500, susc=15000, cond=7500)

    def test_cyl(self):
        # A 30 x 30 cell, with a 10 x 10 block in the middle, split into 9 10 x 10 fragments.

        # flux covering top-left fragment, near2far covering top-middle, force covering top-right
        dft_vecs = make_dft_vecs(
            [mp.FluxRegion(mp.Vector3(-10, z=10), size=mp.Vector3(10, z=10))],
            [mp.Near2FarRegion(mp.Vector3(0, z=10), size=mp.Vector3(10, z=10))],
            [mp.ForceRegion(mp.Vector3(10, z=10), mp.X, size=mp.Vector3(10, z=10))]
        )
        fs = self.get_fragment_stats(mp.Vector3(10, 0, 10), mp.Vector3(30, 0, 30),
                                     mp.CYLINDRICAL, dft_vecs=dft_vecs)

        self.assertEqual(len(fs), 9)

        # Check fragment boxes
        self.assertEqual(fs[0].box.low.x, -15)
        self.assertEqual(fs[0].box.low.z, -15)
        self.assertEqual(fs[0].box.high.x, -5)
        self.assertEqual(fs[0].box.high.z, -5)

        self.assertEqual(fs[1].box.low.x, -15)
        self.assertEqual(fs[1].box.low.z, -5)
        self.assertEqual(fs[1].box.high.x, -5)
        self.assertEqual(fs[1].box.high.z, 5)

        self.assertEqual(fs[2].box.low.x, -15)
        self.assertEqual(fs[2].box.low.z, 5)
        self.assertEqual(fs[2].box.high.x, -5)
        self.assertEqual(fs[2].box.high.z, 15)

        # All fragments besides the middle one have no geometry, only default_material
        for i in [0, 1, 2, 3, 5, 6, 7, 8]:
            self.check_stats(fs[i], 0, 0, 0, 0, 0)

        # Middle fragment contains entire block
        idx = 4
        self.check_stats(fs[idx], a_eps=1000, a_mu=1000, nonlin=3000, susc=3000, cond=3000)

        # Check DFT regions
        for i in [0, 1, 3, 4, 6, 7]:
            self.assertEqual(fs[i].num_dft_pixels, 0)

        self.assertEqual(fs[2].num_dft_pixels, 20000)
        self.assertEqual(fs[5].num_dft_pixels, 28000)
        self.assertEqual(fs[8].num_dft_pixels, 54000)


if __name__ == '__main__':
    unittest.main()
