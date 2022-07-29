import os
import unittest

import h5py
from utils import ApproxComparisonTestCase

import meep as mp


class TestCavityFarfield(ApproxComparisonTestCase):

    data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "data"))

    def run_test(self, nfreqs):
        eps = 13
        w = 1.2
        r = 0.36
        d = 1.4
        N = 3
        sy = 6
        pad = 2
        dpml = 1
        sx = 2 * (pad + dpml + N) + d - 1

        cell = mp.Vector3(sx, sy, 0)

        geometry = [
            mp.Block(
                center=mp.Vector3(),
                size=mp.Vector3(mp.inf, w, mp.inf),
                material=mp.Medium(epsilon=eps),
            )
        ]

        for i in range(N):
            geometry.append(mp.Cylinder(r, center=mp.Vector3(d / 2 + i)))
            geometry.append(mp.Cylinder(r, center=mp.Vector3(d / -2 - i)))

        pml_layers = mp.PML(dpml)
        resolution = 10

        fcen = 0.25
        df = 0.2

        sources = mp.Source(
            src=mp.GaussianSource(fcen, fwidth=df), component=mp.Hz, center=mp.Vector3()
        )

        symmetries = [mp.Mirror(mp.Y, phase=-1), mp.Mirror(mp.X, phase=-1)]

        d1 = 0.2

        sim = mp.Simulation(
            cell_size=cell,
            geometry=geometry,
            sources=[sources],
            symmetries=symmetries,
            boundary_layers=[pml_layers],
            resolution=resolution,
        )

        nearfield = sim.add_near2far(
            fcen,
            0.1,
            nfreqs,
            mp.Near2FarRegion(
                mp.Vector3(0, 0.5 * w + d1), size=mp.Vector3(2 * dpml - sx)
            ),
            mp.Near2FarRegion(
                mp.Vector3(-0.5 * sx + dpml, 0.5 * w + 0.5 * d1),
                size=mp.Vector3(0, d1),
                weight=-1.0,
            ),
            mp.Near2FarRegion(
                mp.Vector3(0.5 * sx - dpml, 0.5 * w + 0.5 * d1), size=mp.Vector3(0, d1)
            ),
            decimation_factor=1,
        )
        sim.run(until=200)
        d2 = 20
        h = 4
        vol = mp.Volume(
            mp.Vector3(0, (0.5 * w) + d2 + (0.5 * h)), size=mp.Vector3(sx - 2 * dpml, h)
        )
        result = sim.get_farfields(nearfield, resolution, where=vol)
        fname = "cavity-farfield.h5" if nfreqs == 1 else "cavity-farfield-4-freqs.h5"
        ref_file = os.path.join(self.data_dir, fname)

        with h5py.File(ref_file, "r") as f:
            # Get reference data into memory
            ref_ex = mp.complexarray(f["ex.r"][()], f["ex.i"][()])
            ref_ey = mp.complexarray(f["ey.r"][()], f["ey.i"][()])
            ref_ez = mp.complexarray(f["ez.r"][()], f["ez.i"][()])
            ref_hx = mp.complexarray(f["hx.r"][()], f["hx.i"][()])
            ref_hy = mp.complexarray(f["hy.r"][()], f["hy.i"][()])
            ref_hz = mp.complexarray(f["hz.r"][()], f["hz.i"][()])

            tol = 1e-5 if mp.is_single_precision() else 1e-7
            self.assertClose(ref_ex, result["Ex"], epsilon=tol)
            self.assertClose(ref_ey, result["Ey"], epsilon=tol)
            self.assertClose(ref_ez, result["Ez"], epsilon=tol)
            self.assertClose(ref_hx, result["Hx"], epsilon=tol)
            self.assertClose(ref_hy, result["Hy"], epsilon=tol)
            self.assertClose(ref_hz, result["Hz"], epsilon=tol)

    def test_cavity_farfield(self):
        self.run_test(nfreqs=1)

    def test_cavity_farfield_four_freqs(self):
        self.run_test(nfreqs=4)


if __name__ == "__main__":
    unittest.main()
