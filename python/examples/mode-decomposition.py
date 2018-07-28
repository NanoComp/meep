import meep as mp
import math
import argparse

def main(args):

    resolution = args.res

    w1 = 1            # width of waveguide 1
    w2 = 2            # width of waveguide 2
    Lw = 10           # length of waveguide 1 and 2
    Lt = args.Lt      # taper length

    Si = mp.Medium(epsilon=12.0)

    dair = 3.0
    dpml = 5.0

    sx = dpml+Lw+Lt+Lw+dpml
    sy = dpml+dair+w2+dair+dpml
    cell_size = mp.Vector3(sx,sy,0)

    prism_x = sx+1
    half_w1 = 0.5*w1
    half_w2 = 0.5*w2
    half_Lt = 0.5*Lt

    if Lt > 0:
        vertices = [mp.Vector3(-prism_x, half_w1),
                    mp.Vector3(-half_Lt, half_w1),
                    mp.Vector3(half_Lt, half_w2),
                    mp.Vector3(prism_x, half_w2),
                    mp.Vector3(prism_x, -half_w2),
                    mp.Vector3(half_Lt, -half_w2),
                    mp.Vector3(-half_Lt, -half_w1),
                    mp.Vector3(-prism_x, -half_w1)]
    else:
        vertices = [mp.Vector3(-prism_x, half_w1),
                    mp.Vector3(prism_x, half_w1),
                    mp.Vector3(prism_x, -half_w1),
                    mp.Vector3(-prism_x, -half_w1)]

    geometry = [mp.Prism(vertices, height=mp.inf, material=Si)]

    boundary_layers = [mp.PML(dpml)]

    # mode wavelength
    lcen = 6.67

    # mode frequency
    fcen = 1/lcen

    sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen, fwidth=0.2*fcen),
                                  component=mp.Ez,
                                  size=mp.Vector3(0,sy-2*dpml,0),
                                  center=mp.Vector3(-0.5*sx+dpml+0.2*Lw,0,0),
                                  eig_match_freq=True,
                                  eig_parity=mp.ODD_Z+mp.EVEN_Y)]

    symmetries=[mp.Mirror(mp.Y)]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=boundary_layers,
                        geometry=geometry,
                        sources=sources,
                        symmetries=symmetries)

    xm = -0.5*sx+dpml+0.5*Lw  # x-coordinate of monitor
    mode_monitor = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(xm,0,0), size=mp.Vector3(0,sy-2*dpml,0)))

    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(xm,0,0), 1e-9))

    coeffs, vgrp, kpoints = sim.get_eigenmode_coefficients(mode_monitor, [1], eig_parity=mp.ODD_Z+mp.EVEN_Y)

    print("mode:, {}, {:.8f}, {:.8f}".format(Lt,abs(coeffs[0,0,0])**2,abs(coeffs[0,0,1])**2))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-Lt', type=float, default=3.0, help='taper length (default: 3.0)')
    parser.add_argument('-res', type=int, default=60, help='resolution (default: 60)')
    args = parser.parse_args()
    main(args)
