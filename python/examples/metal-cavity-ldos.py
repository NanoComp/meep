from __future__ import division

import math
import meep as mp
import argparse

def main(args):
    resolution = 200
    sxy = 2
    dpml = 1
    sxy = sxy + 2 * dpml
    cell = mp.Vector3(sxy, sxy, 0)

    pml_layers = [mp.PML(dpml)]
    a = 1
    t = 0.1
    geometry = [mp.Block(mp.Vector3(a + 2 * t, a + 2 * t, mp.inf), material=mp.metal),
                mp.Block(mp.Vector3(a, a, mp.inf), material=mp.air)]

    w = args.w
    if w > 0:
        geometry.append(mp.Block(center=mp.Vector3(a / 2), size=mp.Vector3(2 * t, w, mp.inf),
                                 material=mp.air))

        fcen = math.sqrt(0.5) / a
        df = 0.2
        sources = [mp.Source(src=mp.GaussianSource(fcen, fwidth=df), component=mp.Ez,
                             center=mp.Vector3())]

        symmetries = [mp.Mirror(mp.Y)]

        Th = args.Th

        sim = mp.Simulation(cell_size=cell,
                            geometry=geometry,
                            boundary_layers=pml_layers,
                            sources=sources,
                            symmetries=symmetries,
                            resolution=resolution)

        h = mp.Harminv(mp.Ez, mp.Vector3(), fcen, df)
        sim.run(mp.after_sources(h), until_after_sources=Th)

        m = h.modes[0]
        f = m.freq
        Q = m.Q
        Vmode = 0.25 * a * a
        print("ldos0:, {}".format(Q / Vmode / (2 * math.pi * f * math.pi * 0.5)))

        sim.reset_meep()
        T = 2 * Q * (1 / f)
        sim.run(mp.dft_ldos(f, 0, 1), until_after_sources=T)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-Th', type=float, default=500, help='additional time after source has turned off to accumulate Harminv data')
    parser.add_argument('-w', type=float, default=0, help='width of cavity opening')
    args = parser.parse_args()
    main(args)
