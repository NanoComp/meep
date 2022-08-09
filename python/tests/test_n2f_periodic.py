import math
import unittest

import numpy as np
from numpy import linalg as LA

import meep as mp


class TestNear2FarPeriodicBoundaries(unittest.TestCase):
    def test_nea2far_periodic(self):
        dpml = 1.0  # PML thickness
        dsub = 3.0  # substrate thickness
        dpad = 20.0  # padding between grating and PML
        gp = 10.0  # grating period
        gh = 0.5  # grating height
        gdc = 0.5  # grating duty cycle

        sx = dpml + dsub + gh + dpad + dpml
        sy = gp

        pml_layers = [mp.PML(thickness=dpml, direction=mp.X)]

        wvl = 0.5  # center wavelength
        fcen = 1 / wvl  # center frequency
        df = 0.05 * fcen  # frequency width

        src_pt = mp.Vector3(-0.5 * sx + dpml + 0.5 * dsub)
        sources = [
            mp.Source(
                mp.GaussianSource(fcen, fwidth=df),
                component=mp.Ez,
                center=src_pt,
                size=mp.Vector3(y=sy),
            )
        ]

        glass = mp.Medium(index=1.5)
        geometry = [
            mp.Block(
                material=glass,
                size=mp.Vector3(dpml + dsub, mp.inf, mp.inf),
                center=mp.Vector3(-0.5 * sx + 0.5 * (dpml + dsub)),
            ),
            mp.Block(
                material=glass,
                size=mp.Vector3(gh, gdc * gp, mp.inf),
                center=mp.Vector3(-0.5 * sx + dpml + dsub + 0.5 * gh),
            ),
        ]

        k_point = mp.Vector3(0, 0, 0)

        symmetries = [mp.Mirror(mp.Y)]

        n2f_pt = mp.Vector3(-0.5 * sx + dpml + dsub + gh + 1.0)
        dft_pt = mp.Vector3(0.5 * sx - dpml)

        res = [20, 25, 30]
        norm = np.empty(3)

        for j in range(3):
            sim = mp.Simulation(
                resolution=res[j],
                cell_size=mp.Vector3(sx, sy),
                boundary_layers=pml_layers,
                geometry=geometry,
                k_point=k_point,
                sources=sources,
                symmetries=symmetries,
            )

            n2f_obj = sim.add_near2far(
                fcen,
                0,
                1,
                mp.Near2FarRegion(center=n2f_pt, size=mp.Vector3(y=sy)),
                nperiods=10,
            )
            dft_obj = sim.add_dft_fields(
                [mp.Ez], fcen, 0, 1, center=dft_pt, size=mp.Vector3(y=sy)
            )

            sim.run(until_after_sources=300)

            n2f_Ez = sim.get_farfields(
                n2f_obj, res[j], center=dft_pt, size=mp.Vector3(y=sy)
            )
            dft_Ez = sim.get_dft_array(dft_obj, mp.Ez, 0)

            norm[j] = LA.norm(n2f_Ez["Ez"] - dft_Ez[1:-1])
            print(f"norm:, {res[j]}, {norm[j]:.5f}")
            sim.reset_meep()

        self.assertGreater(norm[0], norm[1])
        self.assertGreater(norm[1], norm[2])


if __name__ == "__main__":
    unittest.main()
