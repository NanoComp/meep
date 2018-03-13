from __future__ import division

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import meep as mp
from meep import mpb

examples_dir = os.path.realpath(os.path.dirname(__file__))
sys.path.insert(0, examples_dir)

# Run mpb_tri_rods.py and get epsilon as a numpy array
from mpb_tri_rods import ms

efields = []


# Band function to get the efield at each band
def get_efields(ms, band):
    efields.append(ms.get_efield(band))


ms.run_tm(mpb.output_at_kpoint(mp.Vector3(1 / -3, 1 / 3), mpb.fix_efield_phase,
          get_efields))

md = mpb.MPBData(ms, rectify=True, resolution=32, periods=3)

for i, f in enumerate(efields):
    # Get just the z component of the efields
    # TODO: Get the dimensions instead of hard coding
    f = np.reshape(f[2::3], (32, 32))
    efields[i] = md.convert(f)

ms.run_te()

data = ms.get_epsilon()
plt.imshow(data.T, interpolation='spline36', cmap='binary')
plt.axis('off')
plt.show()

md = mpb.MPBData(ms, rectify=True, resolution=32, periods=3)
rectangular_data = md.convert(data)
plt.imshow(rectangular_data.T, interpolation='spline36', cmap='binary')
plt.axis('off')
plt.show()

for f in efields:
    plt.contour(rectangular_data.T, cmap='binary')
    plt.imshow(np.real(f).T, interpolation='spline36', cmap='RdBu', alpha=0.9)
    plt.axis('off')
    plt.show()
