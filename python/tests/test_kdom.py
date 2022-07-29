import math
import unittest

import meep as mp


class TestKdom(unittest.TestCase):
    def run_kdom(self, theta, num_band):

        resolution = 20  # pixels/um

        sx = 5
        sy = 10
        cell_size = mp.Vector3(sx, sy, 0)

        fcen = 1  # center frequency (wavelength = 1 um)
        ng = 1.5
        glass = mp.Medium(index=ng)

        # angle of incident planewave; CCW about Y axis, 0 degrees along +X axis
        theta_in = math.radians(theta)

        # k (in source medium) with correct length (plane of incidence: XY)
        k = mp.Vector3(math.cos(theta_in), math.sin(theta_in), 0).scale(fcen * ng)

        symmetries = []
        eig_parity = mp.ODD_Z
        if theta_in == 0:
            k = mp.Vector3(0, 0, 0)
            eig_parity += mp.EVEN_Y
            symmetries = [mp.Mirror(mp.Y)]

        sim = mp.Simulation(
            resolution=resolution,
            cell_size=cell_size,
            k_point=k,
            symmetries=symmetries,
            default_material=glass,
        )

        sim.init_sim()

        EigenmodeData = sim.get_eigenmode(
            fcen,
            mp.X,
            mp.Volume(center=mp.Vector3(0.3 * sx, 0, 0), size=mp.Vector3(0, sy, 0)),
            num_band,
            k,
            parity=eig_parity,
        )
        kdom = EigenmodeData.kdom

        self.assertAlmostEqual(k.y, kdom.y, places=15)

    def test_kdom(self):
        self.run_kdom(10.7, 6)
        self.run_kdom(22.9, 12)


if __name__ == "__main__":
    unittest.main()
