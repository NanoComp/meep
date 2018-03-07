from __future__ import division

import os
import sys

import h5py

import meep as mp
from meep import mpb

examples_dir = os.path.realpath(os.path.dirname(__file__))
sys.path.insert(0, examples_dir)

from mpb_tri_rods import main as run_tri_rods

# run_tri_rods()

Rin = mp.Matrix(
    mp.Vector3(1, 0, 0),
    mp.Vector3(0, 1, 0),
    mp.Vector3(0, 0, 1)
)

with h5py.File('mpb_data_analysis-epsilon.h5', 'a') as f:

    if 'lattice vectors' in f.keys():
        R = f['lattice vectors'].value

        if R.shape == (3, 3):
            Rin.c1 = mp.Vector3(*R[0])
            Rin.c2 = mp.Vector3(*R[1])
            Rin.c3 = mp.Vector3(*R[2])

    if 'Bloch wavevector' in f.keys():
        kvector = f['Bloch wavevector'].value
    else:
        kvector = None

    if 'copies' in f.keys():
        copies = f['copies'].value
        Rin.c1 = Rin.c1.scale(copies[0])
        Rin.c2 = Rin.c2.scale(copies[1])
        Rin.c3 = Rin.c3.scale(copies[2])

        if kvector:
            kvector *= copies
    else:
        copies = None

    Rout = mp.Matrix(Rin.c1, Rin.c2, Rin.c3)

    if rectify:
        # Orthogonalize the output lattice vectors.  If have_ve is true, then the
        # first new lattice vector should be in the direction of the ve unit
        # vector; otherwise, the first new lattice vector is the first original
        # lattice vector. Note that we do this in such a way as to preserve the
        # volume of the unit cell, and so that our first vector (in the direction
        # of ve) smoothly interpolates between the original lattice vectors.

        if have_ve:
            ve = ve.unit()
        else:
            ve = Rout.c1.unit()

        # First, compute c1 in the direction of ve by smoothly interpolating the
        # old c1/c2/c3 (formula is slightly tricky)
        V = Rout.c1.cross(Rout.c2).dot(Rout.c3)
        Rout.c2 = Rout.c2 - Rout.c1
        Rout.c3 = Rout.c3 - Rout.c1
        Rout.c1 = ve.scale(V / Rout.c2.cross(Rout.c3).ve)

        # Now, orthogonalize c2 and c3
        Rout.c2 = Rout.c2 - ve.scale(ve.dot(Rout.c2))
        Rout.c3 = Rout.c3 - ve.scale(ve.dot(Rout.c3))
        Rout.c3 = Rout.c3 - Rout.c2.scale(Rout.c2.dot(Rout.c3) /
                                          Rout.c2.dot(Rout.c2))

        cart_map.c1 = Rout.c1.unit()
        cart_map.c2 = Rout.c2.unit()
        cart_map.c3 = Rout.c3.unit()
        cart_map = cart_map.inverse()
