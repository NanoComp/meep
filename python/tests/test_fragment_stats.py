import unittest

import meep as mp


def make_dft_vecs(
    flx_reg=None,
    n2f_reg=None,
    frc_reg=None,
    fldc=None,
    flds=None,
    fldw=None,
    fld_cmp=None,
):
    return {
        "flux_regions": flx_reg,
        "n2f_regions": n2f_reg,
        "force_regions": frc_reg,
        "fields_center": fldc,
        "fields_size": flds,
        "fields_where": fldw,
        "fields_components": fld_cmp,
    }


def make_sim(cell, res, pml, dims, create_gv=True, k_point=False):
    sim = mp.Simulation(
        cell_size=cell,
        resolution=res,
        boundary_layers=pml,
        dimensions=dims,
        k_point=k_point,
    )
    if create_gv:
        sim._create_grid_volume(False)
    return sim


class TestFragmentStats(unittest.TestCase):
    def check_stats(self, fragment, a_eps, a_mu, nonlin, susc, cond):
        self.assertEqual(fragment.num_anisotropic_eps_pixels, a_eps)
        self.assertEqual(fragment.num_anisotropic_mu_pixels, a_mu)
        self.assertEqual(fragment.num_nonlinear_pixels, nonlin)
        self.assertEqual(fragment.num_susceptibility_pixels, susc)
        self.assertEqual(fragment.num_nonzero_conductivity_pixels, cond)

    def get_fragment_stats(
        self,
        block_size,
        cell_size,
        dims,
        box_center=mp.Vector3(),
        dft_vecs=None,
        def_mat=mp.air,
        sym=[],
        geom=None,
        pml=[],
    ):
        mat = mp.Medium(
            epsilon=12,
            epsilon_offdiag=mp.Vector3(z=1),
            mu_offdiag=mp.Vector3(x=20),
            E_chi2_diag=mp.Vector3(1, 1),
            H_chi3_diag=mp.Vector3(z=1),
            E_susceptibilities=[
                mp.LorentzianSusceptibility(),
                mp.NoisyLorentzianSusceptibility(),
            ],
            H_susceptibilities=[mp.DrudeSusceptibility()],
            D_conductivity_diag=mp.Vector3(y=1),
            B_conductivity_diag=mp.Vector3(x=1, z=1),
        )

        if geom is None:
            geom = [mp.Block(size=block_size, center=box_center, material=mat)]
        sim = mp.Simulation(
            cell_size=cell_size,
            resolution=10,
            geometry=geom,
            dimensions=dims,
            default_material=def_mat,
            symmetries=sym,
            boundary_layers=pml,
        )

        if dft_vecs:
            if dft_vecs["flux_regions"]:
                sim.add_flux(1, 0.5, 5, *dft_vecs["flux_regions"])
            if dft_vecs["n2f_regions"]:
                sim.add_near2far(1, 0.5, 7, *dft_vecs["n2f_regions"])
            if dft_vecs["force_regions"]:
                sim.add_force(1, 0.5, 9, *dft_vecs["force_regions"])
            if dft_vecs["fields_components"]:
                sim.add_dft_fields(
                    dft_vecs["fields_components"],
                    0,
                    1,
                    5,
                    where=dft_vecs["fields_where"],
                    center=dft_vecs["fields_center"],
                    size=dft_vecs["fields_size"],
                )

        gv = sim._create_grid_volume(False)
        return sim._compute_fragment_stats(gv)

    def _test_1d(self, sym, pml=[]):
        # A z=30 cell, with a size 10 block in the middle.

        # flux covering first 10 units, near2far covering second 10, and force covering third
        dft_vecs = make_dft_vecs(
            [mp.FluxRegion(mp.Vector3(z=-10), size=mp.Vector3(z=10))],
            [mp.Near2FarRegion(mp.Vector3(), size=mp.Vector3(z=10))],
            [mp.ForceRegion(mp.Vector3(z=10), direction=mp.X, size=mp.Vector3(z=10))],
        )

        fs = self.get_fragment_stats(
            mp.Vector3(z=10), mp.Vector3(z=30), 1, dft_vecs=dft_vecs, sym=sym, pml=pml
        )

        sym_factor = 2 if sym else 1
        self.check_stats(
            fs,
            a_eps=300 / sym_factor,
            a_mu=300 / sym_factor,
            nonlin=300 / sym_factor,
            susc=300 / sym_factor,
            cond=300 / sym_factor,
        )

        # Check DFT regions
        self.assertEqual(fs.num_dft_pixels, 40800)
        self.fs = fs

    def test_1d(self):
        self._test_1d([])

    def test_1d_with_symmetry(self):
        self._test_1d([mp.Mirror(mp.X)])

    def test_1d_with_pml(self):
        self._test_1d([], pml=[mp.PML(1)])
        self.assertEqual(self.fs.num_2d_pml_pixels, 0)
        self.assertEqual(self.fs.num_3d_pml_pixels, 0)
        self.assertEqual(self.fs.num_1d_pml_pixels, 20)

    def test_1d_with_overlap(self):
        # A z=30 cell, with a block covering the middle 20 units.
        mat = mp.Medium(H_susceptibilities=[mp.DrudeSusceptibility()])
        fs = self.get_fragment_stats(mp.Vector3(z=20), mp.Vector3(z=30), 1, def_mat=mat)
        self.check_stats(fs, a_eps=300, a_mu=300, nonlin=600, susc=700, cond=600)

    def test_1d_with_partial_fragment(self):
        # A cell with z=26, with a 16 unit block in the center

        # dft_flux with 2 volumes, 1 covering the first 10 units and one covering
        # half of the second 10
        dft_vecs = make_dft_vecs(
            flx_reg=[
                mp.FluxRegion(mp.Vector3(z=-9), mp.Vector3(z=8)),
                mp.FluxRegion(mp.Vector3(z=-2.5), mp.Vector3(z=5)),
            ]
        )
        fs = self.get_fragment_stats(
            mp.Vector3(z=16), mp.Vector3(z=26), 1, dft_vecs=dft_vecs
        )

        self.check_stats(fs, a_eps=260, a_mu=260, nonlin=480, susc=480, cond=480)
        # Check dft stats
        self.assertEqual(fs.num_dft_pixels, 10400)

    def test_1d_dft_fields(self):
        # A z=30 cell with a block covering the middle 10 units.

        # dft_fields covering first 10 units
        dft_vecs = make_dft_vecs(
            fldc=mp.Vector3(z=-10), flds=mp.Vector3(z=10), fld_cmp=[mp.X, mp.Y]
        )
        fs = self.get_fragment_stats(
            mp.Vector3(z=10), mp.Vector3(z=30), 1, dft_vecs=dft_vecs
        )
        self.assertEqual(fs.num_dft_pixels, 4000)

        # Same test with volume instead of center and size
        dft_vecs = make_dft_vecs(
            fldw=mp.Volume(mp.Vector3(z=-10), mp.Vector3(z=10)), fld_cmp=[mp.X, mp.Y]
        )
        fs = self.get_fragment_stats(
            mp.Vector3(z=10), mp.Vector3(z=30), 1, dft_vecs=dft_vecs
        )
        self.assertEqual(fs.num_dft_pixels, 4000)

    def _test_2d(self, sym, pml=[]):
        # A 30 x 30 cell, with a 10 x 10 block in the middle

        # flux covering top-left 10x10, near2far covering top-middle 10x10, force covering top-right
        dft_vecs = make_dft_vecs(
            [mp.FluxRegion(mp.Vector3(-10, 10), size=mp.Vector3(10, 10))],
            [mp.Near2FarRegion(mp.Vector3(0, 10), size=mp.Vector3(10, 10))],
            [
                mp.ForceRegion(
                    mp.Vector3(10, 10), direction=mp.X, size=mp.Vector3(10, 10)
                )
            ],
        )
        fs = self.get_fragment_stats(
            mp.Vector3(10, 10),
            mp.Vector3(30, 30),
            2,
            dft_vecs=dft_vecs,
            sym=sym,
            pml=pml,
        )

        # Check fragment boxes
        self.assertEqual(fs.box.low.x, -15)
        self.assertEqual(fs.box.low.y, -15)
        self.assertEqual(fs.box.high.x, 15)
        self.assertEqual(fs.box.high.y, 15)

        # Middle fragment contains entire block
        sym_factor = 4 if sym else 1
        self.check_stats(
            fs,
            a_eps=90000 / sym_factor,
            a_mu=90000 / sym_factor,
            nonlin=30000 / sym_factor,
            susc=30000 / sym_factor,
            cond=30000 / sym_factor,
        )

        # Check DFT regions
        self.assertEqual(fs.num_dft_pixels, 2040000)
        self.fs = fs

    def test_2d(self):
        self._test_2d([])

    def test_2d_with_symmetry(self):
        self._test_2d([mp.Mirror(mp.X), mp.Mirror(mp.Y)])

    def test_2d_with_pml_all_sides(self):
        self._test_2d(
            [], pml=[mp.PML(1, mp.Y), mp.PML(2, mp.X, mp.Low), mp.PML(3, mp.X, mp.High)]
        )
        self.assertEqual(self.fs.num_1d_pml_pixels, 19000)
        self.assertEqual(self.fs.num_2d_pml_pixels, 1000)
        self.assertEqual(self.fs.num_3d_pml_pixels, 0)

    def test_2d_with_absorbers(self):
        fs = self.get_fragment_stats(
            mp.Vector3(10, 10), mp.Vector3(30, 30), 2, geom=[], pml=[mp.Absorber(1)]
        )
        self.assertEqual(fs.num_1d_pml_pixels, 0)
        self.assertEqual(fs.num_2d_pml_pixels, 0)
        self.assertEqual(fs.num_3d_pml_pixels, 0)
        self.assertEqual(fs.num_nonzero_conductivity_pixels, 11600)

    def test_2d_dft_fields(self):
        # A 30 x 30 cell, with a 10 x 10 block in the middle

        # dft_fields covering 20 by 20 area in center of cell. Test with volume, and center/size
        cmpts = [mp.Ex, mp.Ey, mp.Ez]
        dft_fields_size_center = make_dft_vecs(
            fldc=mp.Vector3(), flds=mp.Vector3(20, 20), fld_cmp=cmpts
        )
        dft_fields_where = make_dft_vecs(
            fldw=mp.Volume(mp.Vector3(), mp.Vector3(20, 20)), fld_cmp=cmpts
        )

        for dft_vec in [dft_fields_size_center, dft_fields_where]:
            fs = self.get_fragment_stats(
                mp.Vector3(10, 10), mp.Vector3(30, 30), 2, dft_vecs=dft_vec
            )
            self.assertEqual(fs.num_dft_pixels, 300000 + 4 * 75000 + 4 * 150000)

    def test_2d_pml_and_absorber(self):
        blayers = [
            mp.PML(1, mp.Y, mp.High),
            mp.PML(2, mp.Y, mp.Low),
            mp.Absorber(1, mp.X, mp.High),
            mp.Absorber(3, mp.X, mp.Low),
        ]
        fs = self.get_fragment_stats(
            mp.Vector3(), mp.Vector3(30, 30), 2, pml=blayers, geom=[]
        )
        self.assertEqual(fs.num_nonzero_conductivity_pixels, 12000)
        self.assertEqual(fs.num_1d_pml_pixels, 9000)
        self.assertEqual(fs.num_2d_pml_pixels, 0)
        self.assertEqual(fs.num_3d_pml_pixels, 0)

    def _test_3d(self, sym, pml=[]):
        # A 30 x 30 x 30 cell with a 10 x 10 x 10 block placed at the center

        # flux covering lower-front-left 10x10x10, near2far covering lower-middle-left,
        # force covering lower-back-left
        dft_vecs = make_dft_vecs(
            [mp.FluxRegion(mp.Vector3(-10, -10, -10), size=mp.Vector3(10, 10, 10))],
            [mp.Near2FarRegion(mp.Vector3(-10, -10, 0), size=mp.Vector3(10, 10, 10))],
            [
                mp.ForceRegion(
                    mp.Vector3(-10, -10, 10),
                    direction=mp.X,
                    size=mp.Vector3(10, 10, 10),
                )
            ],
        )
        fs = self.get_fragment_stats(
            mp.Vector3(10, 10, 10),
            mp.Vector3(30, 30, 30),
            3,
            dft_vecs=dft_vecs,
            sym=sym,
            pml=pml,
        )

        sym_factor = 8 if sym else 1
        self.check_stats(
            fs,
            a_eps=27000000 / sym_factor,
            a_mu=27000000 / sym_factor,
            nonlin=3000000 / sym_factor,
            susc=3000000 / sym_factor,
            cond=3000000 / sym_factor,
        )

        # Check DFT regions
        self.assertEqual(fs.num_dft_pixels, 102000000)
        self.fs = fs

    def test_3d(self):
        self._test_3d([])

    def test_3d_with_symmetry(self):
        self._test_3d([mp.Mirror(mp.X), mp.Mirror(mp.Y), mp.Mirror(mp.Z)])

    def test_3d_with_pml(self):
        self._test_3d(
            [],
            pml=[
                mp.PML(1, mp.Y, mp.High),
                mp.PML(2, mp.Y, mp.Low),
                mp.PML(3, mp.X),
                mp.PML(1, mp.Z, mp.High),
                mp.PML(2, mp.Z, mp.Low),
            ],
        )

        self.assertEqual(self.fs.num_1d_pml_pixels, 8262000)
        self.assertEqual(self.fs.num_2d_pml_pixels, 1188000)
        self.assertEqual(self.fs.num_3d_pml_pixels, 54000)

    def test_3d_with_absorbers(self):
        fs = self.get_fragment_stats(
            mp.Vector3(), mp.Vector3(30, 30, 30), 3, geom=[], pml=[mp.Absorber(1)]
        )
        self.assertEqual(fs.num_1d_pml_pixels, 0)
        self.assertEqual(fs.num_2d_pml_pixels, 0)
        self.assertEqual(fs.num_3d_pml_pixels, 0)
        self.assertEqual(fs.num_nonzero_conductivity_pixels, 5048000)

    def test_cyl(self):
        # A 30 x 30 cell, with a 10 x 10 block in the middle

        # flux covering top-left fragment, near2far covering top-middle, force covering top-right
        dft_vecs = make_dft_vecs(
            [mp.FluxRegion(mp.Vector3(-10, z=10), size=mp.Vector3(10, z=10))],
            [mp.Near2FarRegion(mp.Vector3(0, z=10), size=mp.Vector3(10, z=10))],
            [
                mp.ForceRegion(
                    mp.Vector3(10, z=10), direction=mp.X, size=mp.Vector3(10, z=10)
                )
            ],
        )
        fs = self.get_fragment_stats(
            mp.Vector3(10, 0, 10),
            mp.Vector3(30, 0, 30),
            mp.CYLINDRICAL,
            dft_vecs=dft_vecs,
        )
        self.assertEqual(fs.box.low.x, -15)
        self.assertEqual(fs.box.low.z, -15)
        self.assertEqual(fs.box.high.x, 15)
        self.assertEqual(fs.box.high.z, 15)
        self.check_stats(
            fs, a_eps=90000, a_mu=90000, nonlin=30000, susc=30000, cond=30000
        )
        self.assertEqual(fs.num_dft_pixels, 2040000)

    def test_no_geometry(self):
        mat = mp.Medium(
            epsilon=12,
            epsilon_offdiag=mp.Vector3(x=1),
            mu_offdiag=mp.Vector3(x=20),
            E_chi2_diag=mp.Vector3(1, 1),
            H_chi3_diag=mp.Vector3(x=1),
            E_susceptibilities=[
                mp.LorentzianSusceptibility(),
                mp.NoisyLorentzianSusceptibility(),
            ],
            H_susceptibilities=[mp.DrudeSusceptibility()],
            D_conductivity_diag=mp.Vector3(y=1),
            B_conductivity_diag=mp.Vector3(x=1, z=1),
        )
        fs = self.get_fragment_stats(
            mp.Vector3(), mp.Vector3(10, 10), 2, def_mat=mat, geom=[]
        )
        self.check_stats(
            fs, a_eps=10000, a_mu=10000, nonlin=30000, susc=30000, cond=30000
        )

    def test_1d_cell_smaller_than_minimum_fragment_size(self):
        fs = self.get_fragment_stats(mp.Vector3(z=1), mp.Vector3(z=1), 1)
        self.assertEqual(fs.box.low.z, -0.5)
        self.assertEqual(fs.box.high.z, 0.5)
        self.assertEqual(fs.num_pixels_in_box, 10)

    def test_2d_cell_smaller_than_minimum_fragment_size(self):
        fs = self.get_fragment_stats(mp.Vector3(1, 1), mp.Vector3(1, 1), 2)
        self.assertEqual(fs.box.low.x, -0.5)
        self.assertEqual(fs.box.low.y, -0.5)
        self.assertEqual(fs.box.high.x, 0.5)
        self.assertEqual(fs.box.high.y, 0.5)
        self.assertEqual(fs.num_pixels_in_box, 100)

    def test_3d_cell_smaller_than_minimum_fragment_size(self):
        fs = self.get_fragment_stats(mp.Vector3(1, 1, 1), mp.Vector3(1, 1, 1), 3)
        self.assertEqual(fs.box.low.x, -0.5)
        self.assertEqual(fs.box.low.y, -0.5)
        self.assertEqual(fs.box.low.z, -0.5)
        self.assertEqual(fs.box.high.x, 0.5)
        self.assertEqual(fs.box.high.y, 0.5)
        self.assertEqual(fs.box.high.z, 0.5)
        self.assertEqual(fs.num_pixels_in_box, 1000)


class TestPMLToVolList(unittest.TestCase):
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
        sim = make_sim(mp.Vector3(z=10), 10, [mp.PML(1)], 1)
        v1, v2, v3 = sim._boundary_layers_to_vol_list(sim.boundary_layers)

        self.assertFalse(v2)
        self.assertFalse(v3)
        self.assertEqual(len(v1), 2)
        self.check1d(v1[0], 4, 5)
        self.check1d(v1[1], -5, -4)

    def test_1d_high_side(self):
        sim = make_sim(mp.Vector3(z=10), 10, [mp.PML(1, side=mp.High)], 1)
        v1, v2, v3 = sim._boundary_layers_to_vol_list(sim.boundary_layers)

        self.assertFalse(v2)
        self.assertFalse(v3)
        self.assertEqual(len(v1), 1)
        self.check1d(v1[0], 4, 5)

    def test_1d_two_sides_different_thickness(self):
        sim = make_sim(
            mp.Vector3(z=10), 10, [mp.PML(1, side=mp.High), mp.PML(2, side=mp.Low)], 1
        )
        v1, v2, v3 = sim._boundary_layers_to_vol_list(sim.boundary_layers)

        self.assertFalse(v2)
        self.assertFalse(v3)
        self.assertEqual(len(v1), 2)
        self.check1d(v1[0], 4, 5)
        self.check1d(v1[1], -5, -3)

    def test_2d_all_directions_all_sides(self):
        sim = make_sim(mp.Vector3(10, 10), 10, [mp.PML(1)], 2)
        v1, v2, v3 = sim._boundary_layers_to_vol_list(sim.boundary_layers)

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
            mp.PML(thickness=2, direction=mp.X, side=mp.Low),
        ]
        sim = make_sim(mp.Vector3(10, 10), 10, pmls, 2)
        v1, v2, v3 = sim._boundary_layers_to_vol_list(sim.boundary_layers)

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
        sim = make_sim(mp.Vector3(10, 10), 10, pmls, 2)
        v1, v2, v3 = sim._boundary_layers_to_vol_list(sim.boundary_layers)

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
        sim = make_sim(mp.Vector3(10, 10), 10, [mp.PML(1, mp.X)], 2)
        v1, v2, v3 = sim._boundary_layers_to_vol_list(sim.boundary_layers)

        self.assertFalse(v2)
        self.assertFalse(v3)
        self.assertEqual(len(v1), 2)
        self.check2d(v1[0], mp.Vector3(-5, -5), mp.Vector3(-4, 5))
        self.check2d(v1[1], mp.Vector3(4, -5), mp.Vector3(5, 5))

    def test_3d_all_directions_all_sides(self):
        sim = make_sim(mp.Vector3(10, 10, 10), 10, [mp.PML(1)], 3)
        v1, v2, v3 = sim._boundary_layers_to_vol_list(sim.boundary_layers)

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

    def test_3d_X_direction_only(self):
        sim = make_sim(mp.Vector3(10, 10, 10), 10, [mp.PML(1, mp.X)], 3)
        v1, v2, v3 = sim._boundary_layers_to_vol_list(sim.boundary_layers)

        self.assertEqual(len(v1), 2)
        self.assertEqual(len(v2), 0)
        self.assertEqual(len(v3), 0)

        # left
        self.check3d(v1[0], mp.Vector3(-5, -5, -5), mp.Vector3(-4, 5, 5))
        # right
        self.check3d(v1[1], mp.Vector3(4, -5, -5), mp.Vector3(5, 5, 5))

    def test_cylindrical_all_directions_all_sides(self):
        sim = make_sim(mp.Vector3(10, 0, 10), 10, [mp.PML(1)], mp.CYLINDRICAL)
        v1, v2, v3 = sim._boundary_layers_to_vol_list(sim.boundary_layers)

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


@unittest.skipIf(mp.count_processors() != 2, "MPI specific test")
class TestChunkCommunicationArea(unittest.TestCase):
    def test_2d_periodic(self):
        sim = make_sim(mp.Vector3(10, 6), 10, [mp.PML(1)], 2, k_point=mp.Vector3(1, 1))
        max_comm = sim.get_max_chunk_communication_area()
        avg_comm = sim.get_avg_chunk_communication_area()
        self.assertEqual(max_comm, 220)
        self.assertEqual(avg_comm, 220)

    def test_3d_periodic(self):
        sim = make_sim(
            mp.Vector3(10, 8, 6), 10, [mp.PML(1)], 3, k_point=mp.Vector3(1, 1, 1)
        )
        max_comm = sim.get_max_chunk_communication_area()
        avg_comm = sim.get_avg_chunk_communication_area()
        self.assertEqual(max_comm, 2360)
        self.assertEqual(avg_comm, 2360)

    def test_2d(self):
        sim = make_sim(mp.Vector3(10, 6), 10, [mp.PML(1)], 2)
        max_comm = sim.get_max_chunk_communication_area()
        avg_comm = sim.get_avg_chunk_communication_area()
        self.assertEqual(max_comm, 60)
        self.assertEqual(avg_comm, 60)

    def test_3d(self):
        sim = make_sim(mp.Vector3(10, 8, 6), 10, [mp.PML(1)], 3)
        max_comm = sim.get_max_chunk_communication_area()
        avg_comm = sim.get_avg_chunk_communication_area()
        self.assertEqual(max_comm, 480)
        self.assertEqual(avg_comm, 480)


if __name__ == "__main__":
    unittest.main()
