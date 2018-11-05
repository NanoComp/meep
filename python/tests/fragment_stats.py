from __future__ import division

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
                           def_mat=mp.air, sym=[], geom=None, pml=[]):
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

        if geom is None:
            geom = [mp.Block(size=block_size, center=box_center, material=mat)]
        sim = mp.Simulation(cell_size=cell_size, resolution=10, geometry=geom, dimensions=dims,
                            default_material=def_mat, symmetries=sym, boundary_layers=pml)

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

    def _test_1d(self, sym, pml=[]):
        # A z=30 cell, split into three fragments of size 10 each, with a block
        # covering the middle fragment.

        # flux covering first fragment, near2far covering second, force covering third
        dft_vecs = make_dft_vecs(
            [mp.FluxRegion(mp.Vector3(z=-10), size=mp.Vector3(z=10))],
            [mp.Near2FarRegion(mp.Vector3(), size=mp.Vector3(z=10))],
            [mp.ForceRegion(mp.Vector3(z=10), direction=mp.X, size=mp.Vector3(z=10))]
        )

        fs = self.get_fragment_stats(mp.Vector3(z=10), mp.Vector3(z=30), 1, dft_vecs=dft_vecs, sym=sym, pml=pml)

        self.assertEqual(len(fs), 3)

        # First and last fragments have no geometry, only default_material
        for i in [0, 2]:
            self.check_stats(fs[i], 0, 0, 0, 0, 0)

        # Second fragment contains entire block
        sym_factor = 2 if sym else 1
        self.check_stats(fs[1],
                         a_eps=100 / sym_factor,
                         a_mu=100 / sym_factor,
                         nonlin=300 / sym_factor,
                         susc=300 / sym_factor,
                         cond=300 / sym_factor)

        # Check DFT regions
        self.assertEqual(fs[0].num_dft_pixels, 8224)
        self.assertEqual(fs[1].num_dft_pixels, 11792)
        self.assertEqual(fs[2].num_dft_pixels, 21824)

        self.fs = fs

    def test_1d(self):
        self._test_1d([])

    def test_1d_with_symmetry(self):
        self._test_1d([mp.Mirror(mp.X)])

    def test_1d_with_pml(self):
        self._test_1d([], pml=[mp.PML(1)])

        for i in range(3):
            self.assertEqual(self.fs[i].num_2d_pml_pixels, 0)
            self.assertEqual(self.fs[i].num_3d_pml_pixels, 0)

        self.assertEqual(self.fs[0].num_1d_pml_pixels, 10)
        self.assertEqual(self.fs[1].num_1d_pml_pixels, 0)
        self.assertEqual(self.fs[2].num_1d_pml_pixels, 10)

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
        self.assertEqual(fs[0].num_dft_pixels, 6560)
        self.assertEqual(fs[1].num_dft_pixels, 4160)
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

        self.assertEqual(fs[0].num_dft_pixels, 4000)
        self.assertEqual(fs[1].num_dft_pixels, 80)
        self.assertEqual(fs[2].num_dft_pixels, 0)

        # Same test with volume instead of center and size
        dft_vecs = make_dft_vecs(fldw=mp.Volume(mp.Vector3(z=-10), mp.Vector3(z=10)), fld_cmp=[mp.X, mp.Y])
        fs = self.get_fragment_stats(mp.Vector3(z=10), mp.Vector3(z=30), 1, dft_vecs=dft_vecs)

        self.assertEqual(fs[0].num_dft_pixels, 4000)
        self.assertEqual(fs[1].num_dft_pixels, 80)
        self.assertEqual(fs[2].num_dft_pixels, 0)

    def _test_2d(self, sym, pml=[]):
        # A 30 x 30 cell, with a 10 x 10 block in the middle, split into 9 10 x 10 fragments.

        # flux covering top-left fragment, near2far covering top-middle, force covering top-right
        dft_vecs = make_dft_vecs(
            [mp.FluxRegion(mp.Vector3(-10, 10), size=mp.Vector3(10, 10))],
            [mp.Near2FarRegion(mp.Vector3(0, 10), size=mp.Vector3(10, 10))],
            [mp.ForceRegion(mp.Vector3(10, 10), direction=mp.X, size=mp.Vector3(10, 10))]
        )
        fs = self.get_fragment_stats(mp.Vector3(10, 10), mp.Vector3(30, 30), 2,
                                     dft_vecs=dft_vecs, sym=sym, pml=pml)

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
        sym_factor = 4 if sym else 1
        self.check_stats(fs[idx],
                         a_eps=10000 / sym_factor,
                         a_mu=10000 / sym_factor,
                         nonlin=30000 / sym_factor,
                         susc=30000 / sym_factor,
                         cond=30000 / sym_factor)

        # Check DFT regions
        for i in [0, 3, 6]:
            self.assertEqual(fs[i].num_dft_pixels, 0)

        self.assertEqual(fs[1].num_dft_pixels, 8224)
        self.assertEqual(fs[4].num_dft_pixels, 11792)
        self.assertEqual(fs[7].num_dft_pixels, 21824)

        self.assertEqual(fs[2].num_dft_pixels, 411200)
        self.assertEqual(fs[5].num_dft_pixels, 589600)
        self.assertEqual(fs[8].num_dft_pixels, 1091200)

        self.fs = fs

    def test_2d(self):
        self._test_2d([])

    def test_2d_with_symmetry(self):
        self._test_2d([mp.Mirror(mp.X), mp.Mirror(mp.Y)])

    def test_2d_with_pml_all_sides(self):
        self._test_2d([], pml=[mp.PML(1, mp.Y), mp.PML(2, mp.X, mp.Low), mp.PML(3, mp.X, mp.High)])

        # Center fragment has no PML pixels
        self.assertEqual(self.fs[4].num_1d_pml_pixels, 0)
        self.assertEqual(self.fs[4].num_2d_pml_pixels, 0)
        self.assertEqual(self.fs[4].num_3d_pml_pixels, 0)

        for i in range(len(self.fs)):
            # No regions where 3 PMLs overlap
            self.assertEqual(self.fs[i].num_3d_pml_pixels, 0)

        for i in [1, 3, 5, 7]:
            # No regions where 2 PMLs overlap
            self.assertEqual(self.fs[i].num_2d_pml_pixels, 0)

        for i in [0, 2]:
            # Lower left and top left
            self.assertEqual(self.fs[i].num_1d_pml_pixels, 2600)
            self.assertEqual(self.fs[i].num_2d_pml_pixels, 200)

        for i in [6, 8]:
            # Lower right and top right
            self.assertEqual(self.fs[i].num_1d_pml_pixels, 3400)
            self.assertEqual(self.fs[i].num_2d_pml_pixels, 300)

        for i in [3, 5]:
            # bottom center, top center
            self.assertEqual(self.fs[i].num_1d_pml_pixels, 1000)

        # Right center
        self.assertEqual(self.fs[7].num_1d_pml_pixels, 3000)
        # Left center
        self.assertEqual(self.fs[1].num_1d_pml_pixels, 2000)

    def test_2d_with_overlap(self):
        # A 30 x 30 cell, with a 20 x 20 block in the middle, split into 9 10 x 10 fragments.

        mat = mp.Medium(H_susceptibilities=[mp.DrudeSusceptibility()])
        fs = self.get_fragment_stats(mp.Vector3(20, 20), mp.Vector3(30, 30), 2, def_mat=mat)

        self.assertEqual(len(fs), 9)

        # Middle fragment contains entire block
        idx = 4
        self.check_stats(fs[idx], a_eps=10000, a_mu=10000, nonlin=30000, susc=30000, cond=30000)

        # Top-middle, bottom-middle, left-middle, and right-middle fragments are half
        # covered by the block, and half covered by default_material 'mat'.
        for i in [1, 3, 5, 7]:
            self.check_stats(fs[i], a_eps=5000, a_mu=5000, nonlin=15000, susc=20000, cond=15000)

        # The four corner fragments are quarter-filled by the block, and 3/4 filled by
        # default_material 'mat'
        for i in [0, 2, 6, 8]:
            self.check_stats(fs[i], a_eps=2500, a_mu=2500, nonlin=7500, susc=15000, cond=7500)

    def test_2d_with_partial_fragments_and_shifted_center(self):
        # A 26 x 26 cell with a 18 x 18 Block in the lower right corner

        fs = self.get_fragment_stats(mp.Vector3(18, 18), mp.Vector3(26, 26), 2, mp.Vector3(4, -4))

        self.assertEqual(len(fs), 9)

        # Middle fragment is 10 x 10 and covered by block
        self.check_stats(fs[4], a_eps=10000, a_mu=10000, nonlin=30000, susc=30000, cond=30000)

        for i in [0, 1, 2, 5, 8]:
            # Air
            self.check_stats(fs[i], 0, 0, 0, 0, 0)

        # 10 x 8 fragment, covered by block
        self.assertEqual(fs[3].box.low.x, -5)
        self.assertEqual(fs[3].box.low.y, -13)
        self.assertEqual(fs[3].box.high.x, 5)
        self.assertEqual(fs[3].box.high.y, -5)
        self.check_stats(fs[3], a_eps=8000, a_mu=8000, nonlin=24000, susc=24000, cond=24000)

        # 8 x 10 fragment, covered by block
        self.assertEqual(fs[7].box.low.x, 5)
        self.assertEqual(fs[7].box.low.y, -5)
        self.assertEqual(fs[7].box.high.x, 13)
        self.assertEqual(fs[7].box.high.y, 5)
        self.check_stats(fs[7], a_eps=8000, a_mu=8000, nonlin=24000, susc=24000, cond=24000)

        # 8 x 8 fragment covered by block
        self.assertEqual(fs[6].box.low.x, 5)
        self.assertEqual(fs[6].box.low.y, -13)
        self.assertEqual(fs[6].box.high.x, 13)
        self.assertEqual(fs[6].box.high.y, -5)
        self.check_stats(fs[6], a_eps=6400, a_mu=6400, nonlin=19200, susc=19200, cond=19200)

    def test_2d_dft_fields(self):
        # A 30 x 30 cell, with a 10 x 10 block in the middle, split into 9 10 x 10 fragments.

        # dft_fields covering 20 by 20 area in center of cell. Test with volume, and center/size
        cmpts = [mp.Ex, mp.Ey, mp.Ez]
        dft_fields_size_center = make_dft_vecs(fldc=mp.Vector3(), flds=mp.Vector3(20, 20), fld_cmp=cmpts)
        dft_fields_where = make_dft_vecs(fldw=mp.Volume(mp.Vector3(), mp.Vector3(20, 20)), fld_cmp=cmpts)

        for dft_vec in [dft_fields_size_center, dft_fields_where]:
            fs = self.get_fragment_stats(mp.Vector3(10, 10), mp.Vector3(30, 30), 2, dft_vecs=dft_vec)

            # Middle fragment is fully covered
            self.assertEqual(fs[4].num_dft_pixels, 300000)

            # 4 corners are 1/4 covered
            for i in [0, 2, 6, 8]:
                self.assertEqual(fs[i].num_dft_pixels, 75000)

            # The rest are half covered
            for i in [1, 3, 5, 7]:
                self.assertEqual(fs[i].num_dft_pixels, 150000)

    def _test_3d(self, sym, pml=[]):
        # A 30 x 30 x 30 cell with a 10 x 10 x 10 block placed at the center, split
        # into 27 10 x 10 x 10 fragments

        # flux covering lower-front-left fragment, near2far covering lower-middle-left,
        # force covering lower-back-left
        dft_vecs = make_dft_vecs(
            [mp.FluxRegion(mp.Vector3(-10, -10, -10), size=mp.Vector3(10, 10, 10))],
            [mp.Near2FarRegion(mp.Vector3(-10, -10, 0), size=mp.Vector3(10, 10, 10))],
            [mp.ForceRegion(mp.Vector3(-10, -10, 10), direction=mp.X, size=mp.Vector3(10, 10, 10))]
        )
        fs = self.get_fragment_stats(mp.Vector3(10, 10, 10), mp.Vector3(30, 30, 30), 3,
                                     dft_vecs=dft_vecs, sym=sym, pml=pml)

        self.assertEqual(len(fs), 27)

        # All fragments besides the middle one have no geometry, only default_material
        for i in range(27):
            if i == 13:
                continue
            self.check_stats(fs[i], 0, 0, 0, 0, 0)

        # Middle fragments contains entire block
        idx = 13
        sym_factor = 8 if sym else 1
        self.check_stats(fs[idx],
                         a_eps=1000000 / sym_factor,
                         a_mu=1000000 / sym_factor,
                         nonlin=3000000 / sym_factor,
                         susc=3000000 / sym_factor,
                         cond=3000000 / sym_factor)

        # Check DFT regions
        self.assertEqual(fs[0].num_dft_pixels, 20560000)
        self.assertEqual(fs[1].num_dft_pixels, 29480000)
        self.assertEqual(fs[2].num_dft_pixels, 54560000)

        self.fs = fs

    def test_3d(self):
        self._test_3d([])

    def test_3d_with_symmetry(self):
        self._test_3d([mp.Mirror(mp.X), mp.Mirror(mp.Y), mp.Mirror(mp.Z)])

    def test_3d_with_pml(self):
        self._test_3d([], pml=[mp.PML(1, mp.Y, mp.High), mp.PML(2, mp.Y, mp.Low), mp.PML(3, mp.X),
                               mp.PML(1, mp.Z, mp.High), mp.PML(2, mp.Z, mp.Low)])

        # bottom left near
        self.assertEqual(self.fs[0].num_1d_pml_pixels, 416000)
        self.assertEqual(self.fs[0].num_2d_pml_pixels, 124000)
        self.assertEqual(self.fs[0].num_3d_pml_pixels, 12000)

    def test_3d_with_overlap(self):
        # A 30 x 30 x 30 cell with a 20 x 20 x 20 block placed at the center, split
        # into 27 10 x 10 x 10 fragments

        mat = mp.Medium(E_susceptibilities=[mp.DrudeSusceptibility()])
        fs = self.get_fragment_stats(mp.Vector3(20, 20, 20), mp.Vector3(30, 30, 30), 3, def_mat=mat)

        self.assertEqual(len(fs), 27)

        # Middle fragment contains entire block
        idx = 13
        self.check_stats(fs[idx], a_eps=1000000, a_mu=1000000, nonlin=3000000, susc=3000000, cond=3000000)

        # Six fragments adjacent to the middle fragment faces will be half covered by the block,
        # and half covered by default_material 'mat'
        for i in [4, 10, 12, 14, 16, 22]:
            self.check_stats(fs[i], a_eps=500000, a_mu=500000, nonlin=1500000, susc=2000000, cond=1500000)

        # The corners will be 1/8 covered by the block and 7/8 covered by default_material 'mat'
        for i in [0, 2, 6, 8, 18, 20, 24, 26]:
            self.check_stats(fs[i], a_eps=125000, a_mu=125000, nonlin=375000, susc=1250000, cond=375000)

        # The rest will be 1/4 covered by the block and 3/4 covered by default_material 'mat'
        for i in [1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25]:
            self.check_stats(fs[i], a_eps=250000, a_mu=250000, nonlin=750000, susc=1500000, cond=750000)

    def test_cyl(self):
        # A 30 x 30 cell, with a 10 x 10 block in the middle, split into 9 10 x 10 fragments.

        # flux covering top-left fragment, near2far covering top-middle, force covering top-right
        dft_vecs = make_dft_vecs(
            [mp.FluxRegion(mp.Vector3(-10, z=10), size=mp.Vector3(10, z=10))],
            [mp.Near2FarRegion(mp.Vector3(0, z=10), size=mp.Vector3(10, z=10))],
            [mp.ForceRegion(mp.Vector3(10, z=10), direction=mp.X, size=mp.Vector3(10, z=10))]
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
        self.check_stats(fs[idx], a_eps=10000, a_mu=10000, nonlin=30000, susc=30000, cond=30000)

        # Check DFT regions
        for i in [0, 3, 6]:
            self.assertEqual(fs[i].num_dft_pixels, 0)

        self.assertEqual(fs[1].num_dft_pixels, 8224)
        self.assertEqual(fs[4].num_dft_pixels, 11792)
        self.assertEqual(fs[7].num_dft_pixels, 21824)

        self.assertEqual(fs[2].num_dft_pixels, 411200)
        self.assertEqual(fs[5].num_dft_pixels, 589600)
        self.assertEqual(fs[8].num_dft_pixels, 1091200)

    def test_no_geometry(self):
        mat = mp.Medium(
            epsilon=12,
            epsilon_offdiag=mp.Vector3(x=1),
            mu_offdiag=mp.Vector3(x=20),
            E_chi2_diag=mp.Vector3(1, 1),
            H_chi3_diag=mp.Vector3(x=1),
            E_susceptibilities=[mp.LorentzianSusceptibility(), mp.NoisyLorentzianSusceptibility()],
            H_susceptibilities=[mp.DrudeSusceptibility()],
            D_conductivity_diag=mp.Vector3(y=1),
            B_conductivity_diag=mp.Vector3(x=1, z=1)
        )
        fs = self.get_fragment_stats(mp.Vector3(), mp.Vector3(10, 10), 2, def_mat=mat, geom=[])

        self.assertEqual(len(fs), 1)
        self.check_stats(fs[0], a_eps=10000, a_mu=10000, nonlin=30000, susc=30000, cond=30000)

    def test_1d_cell_smaller_than_minimum_fragment_size(self):
        fs = self.get_fragment_stats(mp.Vector3(z=1), mp.Vector3(z=1), 1)
        self.assertEqual(len(fs), 1)
        stats = fs[0]
        self.assertEqual(stats.box.low.z, -0.5)
        self.assertEqual(stats.box.high.z, 0.5)
        self.assertEqual(stats.num_pixels_in_box, 10)

    def test_2d_cell_smaller_than_minimum_fragment_size(self):
        fs = self.get_fragment_stats(mp.Vector3(1, 1), mp.Vector3(1, 1), 2)
        self.assertEqual(len(fs), 1)
        stats = fs[0]
        self.assertEqual(stats.box.low.x, -0.5)
        self.assertEqual(stats.box.low.y, -0.5)
        self.assertEqual(stats.box.high.x, 0.5)
        self.assertEqual(stats.box.high.y, 0.5)
        self.assertEqual(stats.num_pixels_in_box, 100)

    def test_3d_cell_smaller_than_minimum_fragment_size(self):
        fs = self.get_fragment_stats(mp.Vector3(1, 1, 1), mp.Vector3(1, 1, 1), 3)
        self.assertEqual(len(fs), 1)
        stats = fs[0]
        self.assertEqual(stats.box.low.x, -0.5)
        self.assertEqual(stats.box.low.y, -0.5)
        self.assertEqual(stats.box.low.z, -0.5)
        self.assertEqual(stats.box.high.x, 0.5)
        self.assertEqual(stats.box.high.y, 0.5)
        self.assertEqual(stats.box.high.z, 0.5)
        self.assertEqual(stats.num_pixels_in_box, 1000)


class TestPMLToVolList(unittest.TestCase):

    def make_sim(self, cell, res, pml, dims):
        sim = mp.Simulation(cell_size=cell, resolution=res, boundary_layers=pml, dimensions=dims)
        sim._create_grid_volume(False)
        return sim

    def check1d(self, vol, expected_min, expected_max):
        min_vec = vol.get_min_corner()
        max_vec = vol.get_max_corner()
        min_v3 = mp.Vector3(z=min_vec.z())
        max_v3 = mp.Vector3(z=max_vec.z())
        self.assertEqual(mp.Vector3(z=expected_min), min_v3)
        self.assertEqual(mp.Vector3(z=expected_max), max_v3)

    def check2d(self, vol, expected_min, expected_max):
        min_vec = vol.get_min_corner()
        max_vec = vol.get_max_corner()
        min_v3 = mp.Vector3(min_vec.x(), min_vec.y())
        max_v3 = mp.Vector3(max_vec.x(), max_vec.y())
        self.assertEqual(expected_min, min_v3)
        self.assertEqual(expected_max, max_v3)

    def checkcyl(self, vol, expected_min, expected_max):
        min_vec = vol.get_min_corner()
        max_vec = vol.get_max_corner()
        min_v3 = mp.Vector3(min_vec.r(), 0, min_vec.z())
        max_v3 = mp.Vector3(max_vec.r(), 0, max_vec.z())
        self.assertEqual(expected_min, min_v3)
        self.assertEqual(expected_max, max_v3)

    def check3d(self, vol, expected_min, expected_max):
        min_vec = vol.get_min_corner()
        max_vec = vol.get_max_corner()
        min_v3 = mp.Vector3(min_vec.x(), min_vec.y(), min_vec.z())
        max_v3 = mp.Vector3(max_vec.x(), max_vec.y(), max_vec.z())
        self.assertEqual(expected_min, min_v3)
        self.assertEqual(expected_max, max_v3)

    def test_1d_all_sides(self):
        sim = self.make_sim(mp.Vector3(z=10), 10, [mp.PML(1)], 1)
        v1, v2, v3 = sim._pml_to_vol_list()

        self.assertFalse(v2)
        self.assertFalse(v3)
        self.assertEqual(len(v1), 2)
        self.check1d(v1[0], 4, 5)
        self.check1d(v1[1], -5, -4)

    def test_1d_high_side(self):
        sim = self.make_sim(mp.Vector3(z=10), 10, [mp.PML(1, side=mp.High)], 1)
        v1, v2, v3 = sim._pml_to_vol_list()

        self.assertFalse(v2)
        self.assertFalse(v3)
        self.assertEqual(len(v1), 1)
        self.check1d(v1[0], 4, 5)

    def test_1d_two_sides_different_thickness(self):
        sim = self.make_sim(mp.Vector3(z=10), 10, [mp.PML(1, side=mp.High), mp.PML(2, side=mp.Low)], 1)
        v1, v2, v3 = sim._pml_to_vol_list()

        self.assertFalse(v2)
        self.assertFalse(v3)
        self.assertEqual(len(v1), 2)
        self.check1d(v1[0], 4, 5)
        self.check1d(v1[1], -5, -3)

    def test_2d_all_directions_all_sides(self):
        sim = self.make_sim(mp.Vector3(10, 10), 10, [mp.PML(1)], 2)
        v1, v2, v3 = sim._pml_to_vol_list()

        self.assertFalse(v3)
        self.assertEqual(len(v1), 4)
        self.assertEqual(len(v2), 4)

        # No overlap
        self.check2d(v1[0], mp.Vector3(-4, 4), mp.Vector3(4, 5))
        self.check2d(v1[1], mp.Vector3(-4, -5), mp.Vector3(4, -4))
        self.check2d(v1[2], mp.Vector3(-5, -4), mp.Vector3(-4, 4))
        self.check2d(v1[3], mp.Vector3(4, -4), mp.Vector3(5, 4))

        # Two PMLs overlap
        self.check2d(v2[0], mp.Vector3(-5, 4), mp.Vector3(-4, 5))
        self.check2d(v2[1], mp.Vector3(4, 4), mp.Vector3(5, 5))
        self.check2d(v2[2], mp.Vector3(-5, -5), mp.Vector3(-4, -4))
        self.check2d(v2[3], mp.Vector3(4, -5), mp.Vector3(5, -4))

    def test_2d_all_sides_different_thickness_in_X(self):
        # Thickness 1 on top and bottom, 3 on right, 2 on left
        pmls = [
            mp.PML(thickness=1, direction=mp.Y),
            mp.PML(thickness=3, direction=mp.X, side=mp.High),
            mp.PML(thickness=2, direction=mp.X, side=mp.Low)
        ]
        sim = self.make_sim(mp.Vector3(10, 10), 10, pmls, 2)
        v1, v2, v3 = sim._pml_to_vol_list()

        self.assertFalse(v3)
        self.assertEqual(len(v1), 4)
        self.assertEqual(len(v2), 4)

        # No overlap
        self.check2d(v1[0], mp.Vector3(-3, 4), mp.Vector3(2, 5))
        self.check2d(v1[1], mp.Vector3(-3, -5), mp.Vector3(2, -4))
        self.check2d(v1[2], mp.Vector3(-5, -4), mp.Vector3(-3, 4))
        self.check2d(v1[3], mp.Vector3(2, -4), mp.Vector3(5, 4))

        # Two PMLs overlap
        self.check2d(v2[0], mp.Vector3(-5, 4), mp.Vector3(-3, 5))
        self.check2d(v2[1], mp.Vector3(2, 4), mp.Vector3(5, 5))
        self.check2d(v2[2], mp.Vector3(-5, -5), mp.Vector3(-3, -4))
        self.check2d(v2[3], mp.Vector3(2, -5), mp.Vector3(5, -4))

    def test_2d_three_sides_different_thickness(self):
        # Thickness 3 on top, 2 on left, 1 on right
        pmls = [
            mp.PML(3, mp.Y, mp.High),
            mp.PML(2, mp.X, mp.Low),
            mp.PML(1, mp.X, mp.High),
        ]
        sim = self.make_sim(mp.Vector3(10, 10), 10, pmls, 2)
        v1, v2, v3 = sim._pml_to_vol_list()

        self.assertFalse(v3)
        self.assertEqual(len(v1), 3)
        self.assertEqual(len(v2), 2)

        # No overlap
        self.check2d(v1[0], mp.Vector3(-3, 2), mp.Vector3(4, 5))
        self.check2d(v1[1], mp.Vector3(-5, -5), mp.Vector3(-3, 2))
        self.check2d(v1[2], mp.Vector3(4, -5), mp.Vector3(5, 2))

        # Two PMLs overlap
        self.check2d(v2[0], mp.Vector3(-5, 2), mp.Vector3(-3, 5))
        self.check2d(v2[1], mp.Vector3(4, 2), mp.Vector3(5, 5))

    def test_2d_two_sides(self):
        sim = self.make_sim(mp.Vector3(10, 10), 10, [mp.PML(1, mp.X)], 2)
        v1, v2, v3 = sim._pml_to_vol_list()

        self.assertFalse(v2)
        self.assertFalse(v3)
        self.assertEqual(len(v1), 2)
        self.check2d(v1[0], mp.Vector3(-5, -5), mp.Vector3(-4, 5))
        self.check2d(v1[1], mp.Vector3(4, -5), mp.Vector3(5, 5))

    def test_3d_all_directions_all_sides(self):
        sim = self.make_sim(mp.Vector3(10, 10, 10), 10, [mp.PML(1)], 3)
        v1, v2, v3 = sim._pml_to_vol_list()

        self.assertEqual(len(v1), 6)
        self.assertEqual(len(v2), 12)
        self.assertEqual(len(v3), 8)

        # No overlapping regions (cube faces)
        # top
        self.check3d(v1[0], mp.Vector3(-4, 4, -4), mp.Vector3(4, 5, 4))
        # bottom
        self.check3d(v1[1], mp.Vector3(-4, -5, -4), mp.Vector3(4, -4, 4))
        # left
        self.check3d(v1[2], mp.Vector3(-5, -4, -4), mp.Vector3(-4, 4, 4))
        # right
        self.check3d(v1[3], mp.Vector3(4, -4, -4), mp.Vector3(5, 4, 4))
        # near
        self.check3d(v1[4], mp.Vector3(-4, -4, -5), mp.Vector3(4, 4, -4))
        # far
        self.check3d(v1[5], mp.Vector3(-4, -4, 4), mp.Vector3(4, 4, 5))

        # Two PMLs overlap (cube edges)
        # top left
        self.check3d(v2[0], mp.Vector3(-5, 4, -4), mp.Vector3(-4, 5, 4))
        # top right
        self.check3d(v2[1], mp.Vector3(4, 4, -4), mp.Vector3(5, 5, 4))
        # top near
        self.check3d(v2[2], mp.Vector3(-4, 4, -5), mp.Vector3(4, 5, -4))
        # top far
        self.check3d(v2[3], mp.Vector3(-4, 4, 4), mp.Vector3(4, 5, 5))
        # bottom left
        self.check3d(v2[4], mp.Vector3(-5, -5, -4), mp.Vector3(-4, -4, 4))
        # bottom right
        self.check3d(v2[5], mp.Vector3(4, -5, -4), mp.Vector3(5, -4, 4))
        # bottom near
        self.check3d(v2[6], mp.Vector3(-4, -5, -5), mp.Vector3(4, -4, -4))
        # bottom far
        self.check3d(v2[7], mp.Vector3(-4, -5, 4), mp.Vector3(4, -4, 5))
        # near left
        self.check3d(v2[8], mp.Vector3(-5, -4, -5), mp.Vector3(-4, 4, -4))
        # near right
        self.check3d(v2[9], mp.Vector3(4, -4, -5), mp.Vector3(5, 4, -4))
        # far left
        self.check3d(v2[10], mp.Vector3(-5, -4, 4), mp.Vector3(-4, 4, 5))
        # far right
        self.check3d(v2[11], mp.Vector3(4, -4, 4), mp.Vector3(5, 4, 5))

        # Three PMLs overlap (cube corners)
        # top left near
        self.check3d(v3[0], mp.Vector3(-5, 4, -5), mp.Vector3(-4, 5, -4))
        # top right near
        self.check3d(v3[1], mp.Vector3(4, 4, -5), mp.Vector3(5, 5, -4))
        # top left far
        self.check3d(v3[2], mp.Vector3(-5, 4, 4), mp.Vector3(-4, 5, 5))
        # top right far
        self.check3d(v3[3], mp.Vector3(4, 4, 4), mp.Vector3(5, 5, 5))
        # bottom left near
        self.check3d(v3[4], mp.Vector3(-5, -5, -5), mp.Vector3(-4, -4, -4))
        # bottom right near
        self.check3d(v3[5], mp.Vector3(4, -5, -5), mp.Vector3(5, -4, -4))
        # bottom left far
        self.check3d(v3[6], mp.Vector3(-5, -5, 4), mp.Vector3(-4, -4, 5))
        # bottom right far
        self.check3d(v3[7], mp.Vector3(4, -5, 4), mp.Vector3(5, -4, 5))

    def test_cylindrical_all_directions_all_sides(self):
        sim = self.make_sim(mp.Vector3(10, 0, 10), 10, [mp.PML(1)], mp.CYLINDRICAL)
        v1, v2, v3 = sim._pml_to_vol_list()

        self.assertFalse(v3)
        self.assertEqual(len(v1), 4)
        self.assertEqual(len(v2), 4)

        # No overlap
        self.checkcyl(v1[0], mp.Vector3(-4, 0, 4), mp.Vector3(4, 0, 5))
        self.checkcyl(v1[1], mp.Vector3(-4, 0, -5), mp.Vector3(4, 0, -4))
        self.checkcyl(v1[2], mp.Vector3(-5, 0, -4), mp.Vector3(-4, 0, 4))
        self.checkcyl(v1[3], mp.Vector3(4, 0, -4), mp.Vector3(5, 0, 4))

        # Two PMLs overlap
        self.checkcyl(v2[0], mp.Vector3(-5, 0, 4), mp.Vector3(-4, 0, 5))
        self.checkcyl(v2[1], mp.Vector3(4, 0, 4), mp.Vector3(5, 0, 5))
        self.checkcyl(v2[2], mp.Vector3(-5, 0, -5), mp.Vector3(-4, 0, -4))
        self.checkcyl(v2[3], mp.Vector3(4, 0, -5), mp.Vector3(5, 0, -4))


if __name__ == '__main__':
    unittest.main()
