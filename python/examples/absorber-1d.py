from __future__ import division

import argparse
import math
import meep as mp


def main(args):

    resolution = 20
    cell_size = mp.Vector3(z=10)
    dimensions = 1

    # conversion factor fro eV to 1/um
    eV_um_scale = 1 / 1.23984193

    # Al, from Rakic et al., Applied Optics, vol. 32, p. 5274 (1998)
    Al_eps_inf = 1
    Al_plasma_frq = 14.98 * eV_um_scale

    Al_f0 = 0.523
    Al_frq0 = 1e-10
    Al_gam0 = 0.047 * eV_um_scale
    Al_sig0 = (Al_f0 * math.sqrt(Al_plasma_frq)) / (math.sqrt(Al_frq0))

    Al_f1 = 0.050
    Al_frq1 = 1.544 * eV_um_scale  # 803 nm
    Al_gam1 = 0.312 * eV_um_scale
    Al_sig1 = (Al_f1 * math.sqrt(Al_plasma_frq)) / (math.sqrt(Al_frq1))

    E_susceptibilities = [
        mp.DrudeSusceptibility(frequency=Al_frq0, gamma=Al_gam0, sigma=Al_sig0),
        mp.LorentzianSusceptibility(frequency=Al_frq1, gamma=Al_gam1, sigma=Al_sig1)
    ]

    Al = mp.Medium(epsilon=Al_eps_inf, E_susceptibilities=E_susceptibilities)

    pml_layers = [mp.PML(1, direction=mp.Z) if args.pml else mp.Absorber(1, direction=mp.Z)]

    sources = [mp.Source(src=mp.GaussianSource(1 / 0.803, fwidth=0.1), center=mp.Vector3(),
                         component=mp.Ex)]

    def print_stuff(sim):
        print("ex:, {}, {}".format(sim.meep_time(), sim.get_field_point(mp.Ex, mp.Vector3())))

    sim = mp.Simulation(cell_size=cell_size,
                        resolution=resolution,
                        dimensions=dimensions,
                        default_material=Al,
                        boundary_layers=pml_layers,
                        sources=sources)

    sim.run(mp.at_every(10, print_stuff),
            until_after_sources=mp.stop_when_fields_decayed(50, mp.Ex, mp.Vector3(), 1e-6))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-pml', action='store_true', default=False, help='Use PML as boundary layer')
    args = parser.parse_args()
    main(args)
