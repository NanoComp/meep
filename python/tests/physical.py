# Copyright (C) 2005-2017 Massachusetts Institute of Technology
#
#  This program is free software you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation either version 2, or (at your option)
#  any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program if not, write to the Free Software Foundation,
#  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
from __future__ import division

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

        sources = [mp.Source(src=mp.ContinuousSource(w), component=mp.Ez,
                             center=mp.Vector3(-dx), size=mp.Vector3())]

        sim = mp.Simulation(cell_size=cell_size,
                            resolution=a,
                            boundary_layers=pml_layers,
                            sources=sources,
                            force_complex_fields=True)
        sim.init_fields()
        sim.solve_cw(tol=1e-6)

        p1 = mp.Vector3()
        p2 = mp.Vector3(dx)

        amp1 = sim.get_field_point(mp.Ez, p1)
        amp2 = sim.get_field_point(mp.Ez, p2)

        ratio = abs(amp1) / abs(amp2)
        ratio = ratio ** 2  # in 2d, decay is ~1/sqrt(r), so square to get 1/r

        fail_fmt = "Failed: amp1 = ({}, {}), amp2 = ({}, {})\nabs(amp1/amp2){} = {}, too far from 2.0"
        fail_msg = fail_fmt.format(amp1.real, amp1, amp2.real, amp2, "^2", ratio)

        self.assertTrue(ratio <= 2.12 and ratio >= 1.88, fail_msg)


if __name__ == '__main__':
    unittest.main()
