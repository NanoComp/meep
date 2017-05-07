# Copyright (C) 2005-2017 Massachusetts Institute of Technology
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2, or (at your option)
#  any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.


import meep as mp


class PyCallback(mp.Callback):

    def __init__(self):
        mp.Callback.__init__(self)

    def run(self, vec):
        return 1


def radiating_2d(xmax):
    a = 10.0
    ymax = 3.0

    # grid_volume
    gv = mp.voltwo(xmax, ymax, a)

    one_cb = mp.Caller()
    one_cb.setCallback(PyCallback().__disown__())
    s = mp.structure(gv, one_cb, mp.pml(ymax / 3))

    f = mp.fields(s)
    w = 0.30
    dx = 2.0
    src = mp.continuous_src_time(w)
    f.add_point_source(mp.Ez, src, mp.vec(xmax / 2 - dx, ymax / 2))

    p1 = mp.vec(xmax / 2 + 0 * dx, ymax / 2)
    p2 = mp.vec(xmax / 2 + 1 * dx, ymax / 2)

    # let the source reach steady state
    f.solve_cw(1e-6)
    # while (f.time() < 400)
    #    f.step();

    # amp1 and amp2 are of type complex<double>
    amp1 = f.get_field(mp.Ez, p1)
    amp2 = f.get_field(mp.Ez, p2)

    ratio = pow(abs(amp1) / abs(amp2), 2.0)
    print("Ratio is {} from ({} {}) and ({} {})".format(
        ratio, amp1.real, amp1, amp2.real, amp2
    ))

    fail_fmt = "Failed: amp1 = ({}, {}), amp2 = ({}, {})\nabs(amp1/amp2)^2 = {}, too far from 2.0"
    assert ratio <= 2.12 and ratio >= 1.88, fail_fmt.format(amp1.real, amp1, amp2.real, amp2, ratio)

    return 1


def attempt(name, allright):
    if allright:
        print("Passed {}".format(name))
    else:
        mp.abort("Failed {}!".format(name))


if __name__ == '__main__':
    print("Trying out some physical tests...")
    attempt("radiating source should decay spatially as 1/sqrt(r) in 2D.", radiating_2d(8.0))
