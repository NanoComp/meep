from __future__ import division

import os
import sys
from meep import mpb

examples_dir = os.path.realpath(os.path.dirname(__file__))
sys.path.insert(0, examples_dir)

# Run mpb_tri_rods.py to produce the epsilon.h5 file
from mpb_tri_rods import main as run_tri_rods
run_tri_rods()

h5_fname = 'mpb_data_analysis-epsilon.h5:data'
md = mpb.MPBData(rectify=True, m=3, resolution=32, verbose=True)
md.run(h5_fname)
