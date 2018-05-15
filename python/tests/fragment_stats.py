import unittest
import meep as mp


class TestFragmentStats(unittest.TestCase):

    def test_1d_stats(self):
        mat = mp.Medium(epsilon=12, epsilon_offdiag=mp.Vector3(z=1))
        geom = [mp.Block(size=mp.Vector3(z=10), material=mat)]
        sim = mp.Simulation(cell_size=mp.Vector3(z=30), resolution=10, geometry=geom, dimensions=1)
        gv = sim._create_grid_volume(False)
        fs = mp.compute_fragment_stats(geom, gv, sim.cell_size, sim.default_material,
                                       sim.subpixel_tol, sim.subpixel_maxeval, sim.ensure_periodicity)

        self.assertEqual(len(fs), 3)

        self.assertEqual(fs[0].num_anisotropic_eps_pixels, 0)
        self.assertEqual(fs[0].num_anisotropic_mu_pixels, 0)
        self.assertEqual(fs[0].num_nonlinear_pixels, 0)
        self.assertEqual(fs[0].num_susceptibility_pixels, 0)
        self.assertEqual(fs[0].num_nonzero_conductivity_pixels, 0)

        self.assertEqual(fs[1].num_anisotropic_eps_pixels, 100)
        self.assertEqual(fs[1].num_anisotropic_mu_pixels, 0)
        self.assertEqual(fs[1].num_nonlinear_pixels, 0)
        self.assertEqual(fs[1].num_susceptibility_pixels, 0)
        self.assertEqual(fs[1].num_nonzero_conductivity_pixels, 0)

        self.assertEqual(fs[2].num_anisotropic_eps_pixels, 0)
        self.assertEqual(fs[2].num_anisotropic_mu_pixels, 0)
        self.assertEqual(fs[2].num_nonlinear_pixels, 0)
        self.assertEqual(fs[2].num_susceptibility_pixels, 0)
        self.assertEqual(fs[2].num_nonzero_conductivity_pixels, 0)


if __name__ == '__main__':
    unittest.main()
