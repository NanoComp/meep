import unittest

import meep as mp


class TestEigfreq(unittest.TestCase):
    @unittest.skipIf(
        mp.is_single_precision(), "double-precision floating point specific test"
    )
    def test_eigfreq(self):
        w = 1.2  # width of waveguide
        r = 0.36  # radius of holes
        d = 1.4  # defect spacing (ordinary spacing = 1)
        N = 3  # number of holes on either side of defect
        sy = 6  # size of cell in y direction (perpendicular to wvg.)
        pad = 2  # padding between last hole and PML edge
        dpml = 1  # PML thickness
        sx = 2 * (pad + dpml + N) + d - 1  # size of cell in x direction

        geometry = [
            mp.Block(size=mp.Vector3(mp.inf, w, mp.inf), material=mp.Medium(epsilon=13))
        ]
        for i in range(N):
            geometry.append(mp.Cylinder(r, center=mp.Vector3(d / 2 + i)))
            geometry.append(mp.Cylinder(r, center=mp.Vector3(-(d / 2 + i))))

        fcen = 0.25
        df = 0.2
        src = [
            mp.Source(
                mp.GaussianSource(fcen, fwidth=df),
                component=mp.Hz,
                center=mp.Vector3(0),
                size=mp.Vector3(0, 0),
            )
        ]

        sim = mp.Simulation(
            cell_size=mp.Vector3(sx, sy),
            force_complex_fields=True,
            geometry=geometry,
            boundary_layers=[mp.PML(1.0)],
            sources=src,
            symmetries=[mp.Mirror(mp.X, phase=-1), mp.Mirror(mp.Y, phase=-1)],
            resolution=20,
        )
        sim.init_sim()
        eigfreq = sim.solve_eigfreq(tol=1e-6)

        self.assertAlmostEqual(eigfreq.real, 0.23445413142440263, places=5)
        self.assertAlmostEqual(eigfreq.imag, -0.0003147775697388, places=5)


if __name__ == "__main__":
    unittest.main()
