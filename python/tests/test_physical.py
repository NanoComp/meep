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

import os
import sys
import unittest

mod_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(mod_dir, '..'))
import meep as mp


def one(vec):
    return 1.0


class TestPhysical(unittest.TestCase):

    a = 10.0
    ymax = 3.0
    dx = 2.0

    def radiating_base(self, sq_ratio, solve_cw=True):

        w = 0.30

        s = mp.structure(self.gv, one, mp.pml(self.ymax / 3.0))

        f = mp.fields(s)

        src = mp.continuous_src_time(w)
        f.add_point_source(mp.Ez, src, self.pnt_src_vec)

        # let the source reach steady state
        if solve_cw:
            f.solve_cw(1e-6)
        else:
            while f.time() < 400:
                f.step()

        # amp1 and amp2 are of type complex
        amp1 = f.get_field(mp.Ez, self.p1)
        amp2 = f.get_field(mp.Ez, self.p2)

        ratio = abs(amp1) / abs(amp2)
        if self.gv.dim == mp.D2:
            ratio = ratio ** 2  # in 2d, decay is ~1/sqrt(r), so square to get 1/r

        print("Ratio is {} from ({} {}) and ({} {})".format(
            ratio, amp1.real, amp1, amp2.real, amp2
        ))

        fail_fmt = "Failed: amp1 = ({}, {}), amp2 = ({}, {})\nabs(amp1/amp2){} = {}, too far from 2.0"
        fail_msg = fail_fmt.format(amp1.real, amp1, amp2.real, amp2, "^2" if sq_ratio else "", ratio)

        self.assertTrue(ratio <= 2.12 and ratio >= 1.88, fail_msg)

    def test_radiating_2d(self):

        xmax = 8.0

        # grid_volume
        self.gv = mp.voltwo(xmax, self.ymax, self.a)

        self.pnt_src_vec = mp.vec(xmax / 2 - self.dx, self.ymax / 2)
        self.p1 = mp.vec(xmax / 2 + 0 * self.dx, self.ymax / 2)
        self.p2 = mp.vec(xmax / 2 + 1 * self.dx, self.ymax / 2)

        self.radiating_base(True)

    def test_radiating_3d(self):

        xmax = 7.0

        # grid_volume
        self.gv = mp.vol3d(xmax, self.ymax, self.ymax, self.a)

        self.pnt_src_vec = mp.vec(xmax / 2.0 - self.dx, self.ymax / 2.0, self.ymax / 2.0)
        self.p1 = mp.vec(xmax / 2.0 + 0 * self.dx, self.ymax / 2.0, self.ymax / 2.0)
        self.p2 = mp.vec(xmax / 2.0 + 1 * self.dx, self.ymax / 2.0, self.ymax / 2.0)

        self.radiating_base(False)


if __name__ == '__main__':
    unittest.main()
