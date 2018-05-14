import unittest
import meep as mp


class TestFragmentStats(unittest.TestCase):

    def test_1d_stats(self):
        mat = mp.Medium(epsilon=12)
        geom = [mp.Block(size=mp.Vector3(z=3), material=mat)]
        sim = mp.Simulation(cell_size=mp.Vector3(z=5), resolution=10, geometry=geom, dimensions=1)
        gv = sim._create_grid_volume(False)
        fs = mp.compute_fragment_stats(geom, gv, sim.subpixel_tol, sim.subpixel_maxeval)

        self.assertEqual(len(fs), 5)
        

if __name__ == '__main__':
    unittest.main()
