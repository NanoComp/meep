import math
import unittest

import numpy as np

import meep as mp


class TestGaussianBeamSource(unittest.TestCase):
    def gaussian_beam(self, rot_angle, beamfunc=mp.GaussianBeamSource):

        s = 14
        resolution = 25
        dpml = 2

        cell_size = mp.Vector3(s, s)
        boundary_layers = [mp.PML(thickness=dpml)]

        beam_x0 = mp.Vector3(0, 7.0)  # beam focus (relative to source center)
        beam_kdir = mp.Vector3(0, 1, 0).rotate(
            mp.Vector3(0, 0, 1), math.radians(rot_angle)
        )  # beam propagation direction
        beam_w0 = 0.8  # beam waist radius
        beam_E0 = mp.Vector3(0, 0, 1)
        beam_x0 = beam_x0.rotate(mp.Vector3(0, 0, 1), math.radians(rot_angle))
        fcen = 1
        src_x = 0
        src_y = -0.5 * s + dpml + 1.0
        sources = [
            beamfunc(
                src=mp.GaussianSource(fcen, fwidth=0.2 * fcen),
                center=mp.Vector3(src_x, src_y),
                size=mp.Vector3(s),
                beam_x0=beam_x0,
                beam_kdir=beam_kdir,
                beam_w0=beam_w0,
                beam_E0=beam_E0,
            )
        ]

        sim = mp.Simulation(
            resolution=resolution,
            cell_size=cell_size,
            boundary_layers=boundary_layers,
            sources=sources,
        )

        cell_dft_fields = sim.add_dft_fields(
            [mp.Ez],
            fcen,
            0,
            1,
            center=mp.Vector3(),
            size=mp.Vector3(s - 2 * dpml, s - 2 * dpml),
        )

        sim.run(until_after_sources=50)

        Ez_cell = sim.get_dft_array(cell_dft_fields, mp.Ez, 0)
        [x, y, z, w] = sim.get_array_metadata(dft_cell=cell_dft_fields)

        tol = 0.05
        idx_x = np.nonzero(
            (np.squeeze(x) > (src_x + beam_x0.x - tol))
            & (np.squeeze(x) < (src_x + beam_x0.x + tol))
        )
        idx_y = np.nonzero(
            (np.squeeze(y) > (src_y + beam_x0.y - tol))
            & (np.squeeze(y) < (src_y + beam_x0.y + tol))
        )

        Ez_beam_x0 = Ez_cell[np.squeeze(idx_x)[0], np.squeeze(idx_y)[0]]

        frac = np.abs(Ez_beam_x0) ** 2 / np.amax(np.abs(Ez_cell) ** 2)
        print(
            f"ratio of the Gaussian beam energy at the focus over the maximum beam energy for the entire cell: {frac}"
        )

        self.assertGreater(frac, 0.98)

    def test_gaussian_beam(self):
        self.gaussian_beam(-40, mp.GaussianBeam2DSource)
        self.gaussian_beam(-40, mp.GaussianBeam3DSource)


if __name__ == "__main__":
    unittest.main()
