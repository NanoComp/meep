import os
import sys

import matplotlib.pyplot as plt
import numpy as np

import meep as mp
from meep import mpb

examples_dir = os.path.realpath(os.path.dirname(__file__))
sys.path.insert(0, examples_dir)


def tri_rods():
    # Import the ModeSolver defined in the mpb_tri_rods.py example
    from mpb_tri_rods import ms as tr_ms

    efields = []

    # Band function to collect the efields
    def get_efields(tr_ms, band):
        efields.append(tr_ms.get_efield(band))

    tr_ms.run_tm(
        mpb.output_at_kpoint(
            mp.Vector3(1 / -3, 1 / 3), mpb.fix_efield_phase, get_efields
        )
    )

    # Create an MPBData instance to transform the efields
    md = mpb.MPBData(rectify=True, resolution=32, periods=3)

    converted = []
    for f in efields:
        # Get just the z component of the efields
        f = f[..., 0, 2]
        converted.append(md.convert(f))

    tr_ms.run_te()

    eps = tr_ms.get_epsilon()
    plt.imshow(eps.T, interpolation="spline36", cmap="binary")
    plt.axis("off")
    plt.show()

    md = mpb.MPBData(rectify=True, resolution=32, periods=3)
    rectangular_data = md.convert(eps)
    plt.imshow(rectangular_data.T, interpolation="spline36", cmap="binary")
    plt.axis("off")
    plt.show()

    for i, f in enumerate(converted):
        plt.subplot(331 + i)
        plt.contour(rectangular_data.T, cmap="binary")
        plt.imshow(np.real(f).T, interpolation="spline36", cmap="RdBu", alpha=0.9)
        plt.axis("off")

    plt.show()


def diamond():
    # Import the ModeSolver from the mpb_diamond.py example
    from mpb_diamond import ms as d_ms

    dpwr = []

    def get_dpwr(ms, band):
        dpwr.append(ms.get_dpwr(band))

    d_ms.run(mpb.output_at_kpoint(mp.Vector3(0, 0.625, 0.375), get_dpwr))

    md = mpb.MPBData(rectify=True, periods=2, resolution=32)
    converted_dpwr = [md.convert(d) for d in dpwr]

    # TODO: Plot


if __name__ == "__main__":
    tri_rods()
    diamond()
