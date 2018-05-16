# 1d simulation of a plane wave propagating through a Kerr medium
# and generating the third-harmonic frequency component.

from __future__ import division

import meep as mp
import argparse

def main(args):

    sz = 100          # size of cell in z direction
    fcen = 1 / 3.0    # center frequency of source
    df = fcen / 20.0  # frequency width of source
    amp = args.amp    # amplitude of source
    k = 10**args.logk # Kerr susceptibility
    dpml = 1.0        # PML thickness

    # We'll use an explicitly 1d simulation.  Setting dimensions=1 will actually
    # result in faster execution than just using two no-size dimensions.  However,
    # in this case Meep requires us to use E in the x direction (and H in y),
    # and our one no-size dimension must be z.
    dimensions = 1
    cell = mp.Vector3(0, 0, sz)
    pml_layers = mp.PML(dpml)
    resolution = 20

    # to put the same material in all space, we can just set the default material
    # and pass it to the Simulation constructor
    default_material = mp.Medium(index=1, chi3=k)

    sources = mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ex,
                        center=mp.Vector3(0, 0, -0.5*sz + dpml), amplitude=amp)

    # frequency range for flux calculation
    nfreq = 400
    fmin = fcen / 2.0
    fmax = fcen * 4

    sim = mp.Simulation(cell_size=cell,
                        geometry=[],
                        sources=[sources],
                        boundary_layers=[pml_layers],
                        default_material=default_material,
                        resolution=resolution,
                        dimensions=dimensions)

    # sim.add_flux(0.5 * (fmin + fmax), fmax - fmin, nfreq,
    #              mp.FluxRegion(mp.Vector3(0, 0, 0.5*sz - dpml - 0.5)))
    sim.add_flux(fcen, 0, 1, mp.FluxRegion(mp.Vector3(0, 0, 0.5*sz - dpml - 0.5)))
    sim.add_flux(3 * fcen, 0, 1, mp.FluxRegion(mp.Vector3(0, 0, 0.5*sz - dpml - 0.5)))

    sim.run(until_after_sources=mp.stop_when_fields_decayed(
        50, mp.Ex, mp.Vector3(0, 0, 0.5*sz - dpml - 0.5), 1e-6))

    trans1, trans3 = sim.dft_objects

    # sim.display_fluxes(trans)
    print("harmonics:, {}, {}, {}, {}".format(k, amp, mp.get_fluxes(trans1)[0], mp.get_fluxes(trans3)[0]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-amp', type=float, default=1.0, help='amplitude of source')
    parser.add_argument('-logk', type=float, default=0, help='logarithm of Kerr susceptibility')
    args = parser.parse_args()
    main(args)
