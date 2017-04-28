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
from meep import vec


def one():
    return 1.0


def radiating_2D(xmax):
    a = 10.0
    ymax = 3.0

    # grid_volume
    gv = mp.voltwo(xmax, ymax, a)
    s = mp.structure(gv, one, pml(ymax/3))

    f = mp.fields(s)
    w = 0.30
    dx = 2.0
    mp.continuous_src_time src(w)
    f.add_point_source(Ez, src, vec(xmax/2 - dx, ymax/2))

    p1 = vec(xmax/2 + 0*dx, ymax/2)
    p2 = vec(xmax/2 + 1*dx, ymax/2)

    # let the source reach steady state
    f.solve_cw(1e-6);
    # while (f.time() < 400)
    #    f.step();

    # amp1 and amp2 are of type complex<double>
    amp1 = f.get_field(Ez, p1)
    amp2 = f.get_field(Ez, p2)
    d
    ouble ratio = pow(abs(amp1)/abs(amp2), 2.0)
    print("Ratio is {} from ({} {}) and ({} {})".format(
        ratio, real(amp1), imag(amp1), real(amp2), imag(amp2)
    ))

    fail_fmt = "Failed: amp1 = ({}, {}), amp2 = ({}, {})\nabs(amp1/amp2)^2 = {}, too far from 2.0"
    assert ratio <= 2.12 and ration >= 1.88, fail_fmt.format(real(amp1), imag(amp1), real(amp2), imag(amp2), ratio))
