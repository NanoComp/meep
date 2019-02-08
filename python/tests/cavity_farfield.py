import os
import unittest
import h5py
import numpy as np
import meep as mp


class TestCavityFarfield(unittest.TestCase):

    data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))

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

        geometry = [mp.Block(center=mp.Vector3(), size=mp.Vector3(mp.inf, w, mp.inf),
                            material=mp.Medium(epsilon=eps))]

        for i in range(N):
            geometry.append(mp.Cylinder(r, center=mp.Vector3(d / 2 + i)))
            geometry.append(mp.Cylinder(r, center=mp.Vector3(d / -2 - i)))

        pml_layers = mp.PML(dpml)
        resolution = 10

        fcen = 0.25
        df = 0.2

        sources = mp.Source(src=mp.GaussianSource(fcen, fwidth=df), component=mp.Hz, center=mp.Vector3())

        symmetries = [mp.Mirror(mp.Y, phase=-1), mp.Mirror(mp.X, phase=-1)]

        d1 = 0.2

        sim = mp.Simulation(cell_size=cell,
                            geometry=geometry,
                            sources=[sources],
                            symmetries=symmetries,
                            boundary_layers=[pml_layers],
                            resolution=resolution)

        nearfield = sim.add_near2far(
            fcen, 0, nfreqs,
            mp.Near2FarRegion(mp.Vector3(0, 0.5 * w + d1), size=mp.Vector3(2 * dpml - sx)),
            mp.Near2FarRegion(mp.Vector3(-0.5 * sx + dpml, 0.5 * w + 0.5 * d1), size=mp.Vector3(0, d1), weight=-1.0),
            mp.Near2FarRegion(mp.Vector3(0.5 * sx - dpml, 0.5 * w + 0.5 * d1), size=mp.Vector3(0, d1))
        )
        sim.run(until=200)
        d2 = 20
        h = 4
        vol = mp.Volume(mp.Vector3(0, (0.5 * w) + d2 + (0.5 * h)), size=mp.Vector3(sx - 2 * dpml, h))
        result = sim.get_farfields(nearfield, resolution, where=vol)
        fname = 'cavity-farfield.h5' if nfreqs == 1 else 'cavity-farfield-4-freqs.h5'
        ref_file = os.path.join(self.data_dir, fname)

        with h5py.File(ref_file, 'r') as f:
            # Get reference data into memory
            ref_ex = mp.complexarray(f['ex.r'].value, f['ex.i'].value)
            ref_ey = mp.complexarray(f['ey.r'].value, f['ey.i'].value)
            ref_ez = mp.complexarray(f['ez.r'].value, f['ez.i'].value)
            ref_hx = mp.complexarray(f['hx.r'].value, f['hx.i'].value)
            ref_hy = mp.complexarray(f['hy.r'].value, f['hy.i'].value)
            ref_hz = mp.complexarray(f['hz.r'].value, f['hz.i'].value)

            np.testing.assert_allclose(ref_ex, result['Ex'])
            np.testing.assert_allclose(ref_ey, result['Ey'])
            np.testing.assert_allclose(ref_ez, result['Ez'])
            np.testing.assert_allclose(ref_hx, result['Hx'])
            np.testing.assert_allclose(ref_hy, result['Hy'])
            np.testing.assert_allclose(ref_hz, result['Hz'])

    def test_cavity_farfield(self):
        self.run_test(nfreqs=1)

    def test_cavity_farfield_four_freqs(self):
        self.run_test(nfreqs=4)


if __name__ == '__main__':
    unittest.main()
