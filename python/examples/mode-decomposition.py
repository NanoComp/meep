
import meep as mp
import math
import argparse
import numpy as np

def main(args):

    resolution = args.res

    w1 = 1            # width of waveguide 1
    w2 = args.w2      # width of waveguide 2
    Lw = 10           # length of waveguide 1 and 2
    Lt = args.Lt      # taper length

    Si = mp.Medium(epsilon=12.0)

    dair = 3.0
    dpml = 3.0

    sx = dpml + Lw + Lt + Lw + dpml
    sy = dpml + dair + w2 + dair + dpml
    cell_size = mp.Vector3(sx, sy, 0)

    geometry = [ mp.Block(material=Si, center=mp.Vector3(0,0,0), size=mp.Vector3(mp.inf,w1,mp.inf)),
                 mp.Block(material=Si, center=mp.Vector3(0.5*sx-0.5*(Lt+Lw+dpml),0,0), size=mp.Vector3(Lt+Lw+dpml,w2,mp.inf)) ]

    # form linear taper

    hh = w2
    ww = 2*Lt

    # taper angle (CCW, relative to +X axis)
    rot_theta = math.atan(0.5*(w2-w1)/Lt)

    pvec = mp.Vector3(-0.5*sx+dpml+Lw,0.5*w1,0)
    cvec = mp.Vector3(-0.5*sx+dpml+Lw+0.5*ww,0.5*hh+0.5*w1,0)
    rvec = cvec-pvec
    rrvec = rvec.rotate(mp.Vector3(0,0,1), rot_theta)

    geometry.append(mp.Block(material=mp.air, center=pvec+rrvec, size=mp.Vector3(ww,hh,mp.inf),
                             e1=mp.Vector3(1,0,0).rotate(mp.Vector3(0,0,1),rot_theta),
                             e2=mp.Vector3(0,1,0).rotate(mp.Vector3(0,0,1),rot_theta),
                             e3=mp.Vector3(0,0,1)))

    pvec = mp.Vector3(-0.5*sx+dpml+Lw,-0.5*w1,0)
    cvec = mp.Vector3(-0.5*sx+dpml+Lw+0.5*ww,-(0.5*hh+0.5*w1),0)
    rvec = cvec-pvec
    rrvec = rvec.rotate(mp.Vector3(0,0,1),-rot_theta)

    geometry.append(mp.Block(material=mp.air, center=pvec+rrvec, size=mp.Vector3(ww,hh,mp.inf),
                             e1=mp.Vector3(1,0,0).rotate(mp.Vector3(0,0,1),-rot_theta),
                             e2=mp.Vector3(0,1,0).rotate(mp.Vector3(0,0,1),-rot_theta),
                             e3=mp.Vector3(0,0,1)))

    boundary_layers = [ mp.PML(dpml) ]

    # mode frequency
    fcen = 0.15

    sources = [ mp.EigenModeSource(src=mp.GaussianSource(fcen, fwidth=0.5*fcen),
                                   component=mp.Ez,
                                   size=mp.Vector3(0,sy-2*dpml,0),
                                   center=mp.Vector3(-0.5*sx+dpml,0,0),
                                   eig_match_freq=True,
                                   eig_parity=mp.ODD_Z,
                                   eig_kpoint=mp.Vector3(0.4,0,0),
                                   eig_resolution=32) ]

    symmetries = [ mp.Mirror(mp.Y,+1) ]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=boundary_layers,
                        geometry=geometry,
                        sources=sources,
                        symmetries=symmetries)

    xm = -0.5*sx + dpml + 0.5*Lw  # x-coordinate of monitor
    mflux = sim.add_eigenmode(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(xm,0), size=mp.Vector3(0,sy-2*dpml)))

    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(xm,0), 1e-10))

    bands = [1]  # indices of modes for which to compute expansion coefficients
    alpha = sim.get_eigenmode_coefficients(mflux, bands)

    alpha0Plus  = alpha[2*0 + 0]  # coefficient of forward-traveling fundamental mode
    alpha0Minus = alpha[2*0 + 1]  # coefficient of backward-traveling fundamental mode

    print("refl:, {}, {:.8f}".format(Lt, abs(alpha0Minus)**2))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-Lt',  type=float, default=3.0, help='taper length (default: 3.0)')
    parser.add_argument('-w2',  type=float, default=2.0, help='width of outgoing waveguide (default: 2.0)')
    parser.add_argument('-res', type=int, default=20, help='resolution (default: 20)')
    args = parser.parse_args()
    main(args)
