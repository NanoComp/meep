import unittest

import meep as mp


class TestPhysical(unittest.TestCase):
    def test_physical(self):

        a = 10.0
        ymax = 3.0
        xmax = 8.0
        dx = 2.0
        w = 0.30

        cell_size = mp.Vector3(xmax, ymax)
        pml_layers = [mp.PML(ymax / 3.0)]

        sources = [
            mp.Source(
                src=mp.ContinuousSource(w),
                component=mp.Ez,
                center=mp.Vector3(-dx),
                size=mp.Vector3(),
            )
        ]

        sim = mp.Simulation(
            cell_size=cell_size,
            resolution=a,
            boundary_layers=pml_layers,
            sources=sources,
            force_complex_fields=True,
        )
        sim.init_sim()
        sim.solve_cw(tol=1e-5 if mp.is_single_precision() else 1e-6)

        p1 = mp.Vector3()
        p2 = mp.Vector3(dx)

        amp1 = sim.get_field_point(mp.Ez, p1)
        amp2 = sim.get_field_point(mp.Ez, p2)

        ratio = abs(amp1) / abs(amp2)
        ratio = ratio**2  # in 2d, decay is ~1/sqrt(r), so square to get 1/r

        fail_fmt = "Failed: amp1 = ({}, {}), amp2 = ({}, {})\nabs(amp1/amp2){} = {}, too far from 2.0"
        fail_msg = fail_fmt.format(amp1.real, amp1, amp2.real, amp2, "^2", ratio)

        self.assertTrue(ratio <= 2.12 and ratio >= 1.88, fail_msg)


if __name__ == "__main__":
    unittest.main()
