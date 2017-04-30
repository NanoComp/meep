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


def one(vec):
    return 1.0


def radiating_2d(xmax):
    a = 10.0
    ymax = 3.0

    # grid_volume
    gv = mp.voltwo(xmax, ymax, a)
    s = mp.structure(gv, one, mp.pml(ymax / 3))

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
        ratio, float(amp1), amp1, float(amp2), amp2
    ))

    fail_fmt = "Failed: amp1 = ({}, {}), amp2 = ({}, {})\nabs(amp1/amp2)^2 = {}, too far from 2.0"
    assert ratio <= 2.12 and ratio >= 1.88, fail_fmt.format(float(amp1), amp1, float(amp2), amp2, ratio)

if __name__ == '__main__':
    radiating_2d(8.0)
