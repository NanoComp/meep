# -*- coding: utf-8 -*-

import meep as mp
import argparse

def main(args):
    resolution = 30   # pixels/μm
    
    Si = mp.Medium(index=3.45)

    dpml = 1.0
    pml_layers = [mp.PML(dpml)]
    
    sx = 5
    sy = 3
    cell = mp.Vector3(sx+2*dpml,sy+2*dpml,0)

    a = 1.0     # waveguide width
    s = args.s  # waveguide separation distance

    geometry = [mp.Block(center=mp.Vector3(-0.5*(s+a)),
                         size=mp.Vector3(a,a,mp.inf),
                         material=Si),
                mp.Block(center=mp.Vector3(0.5*(s+a)),
                         size=mp.Vector3(a,a,mp.inf),
                         material=Si)]

    xodd = args.xodd
    symmetries = [mp.Mirror(mp.X, phase=-1.0 if xodd else 1.0),
                  mp.Mirror(mp.Y, phase=-1.0)]

    k_point = mp.Vector3(z=0.5)

    fcen = 0.22
    df = 0.06
    sources = [mp.Source(src=mp.GaussianSource(fcen, fwidth=df),
                         component=mp.Ey,
                         center=mp.Vector3(-0.5*(s+a)),
                         size=mp.Vector3(a,a)),
               mp.Source(src=mp.GaussianSource(fcen, fwidth=df),
                         component=mp.Ey,
                         center=mp.Vector3(0.5*(s+a)),
                         size=mp.Vector3(a,a),
                         amplitude=-1.0 if xodd else 1.0)]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell,
                        boundary_layers=pml_layers,
                        geometry=geometry,
                        symmetries=symmetries,
                        k_point=k_point,
                        sources=sources)

    h = mp.Harminv(mp.Ey, mp.Vector3(0.5*(s+a)), fcen, df)

    sim.run(mp.after_sources(h), until_after_sources=200)

    f = h.modes[0].freq
    print("freq:, {}, {}".format(s, f))

    sim.reset_meep()

    eig_sources = [mp.EigenModeSource(src=mp.GaussianSource(f, fwidth=df),
                                      size=mp.Vector3(a,a),
                                      center=mp.Vector3(-0.5*(s+a)),
                                      eig_kpoint=k_point,
                                      eig_match_freq=True,
                                      eig_parity=mp.ODD_Y),
                   mp.EigenModeSource(src=mp.GaussianSource(f, fwidth=df),
                                      size=mp.Vector3(a,a),
                                      center=mp.Vector3(0.5*(s+a)),
                                      eig_kpoint=k_point,
                                      eig_match_freq=True,
                                      eig_parity=mp.ODD_Y,
                                      amplitude=-1.0 if xodd else 1.0)]

    sim.change_sources(eig_sources)

    flux_reg = mp.FluxRegion(direction=mp.Z, center=mp.Vector3(), size=mp.Vector3(1.2*(2*a+s),1.2*a))
    wvg_flux = sim.add_flux(f, 0, 1, flux_reg)

    force_reg1 = mp.ForceRegion(mp.Vector3(0.5*s), direction=mp.X, weight=1.0, size=mp.Vector3(y=a))
    force_reg2 = mp.ForceRegion(mp.Vector3(0.5*s+a), direction=mp.X, weight=-1.0, size=mp.Vector3(y=a))
    wvg_force = sim.add_force(f, 0, 1, force_reg1, force_reg2)

    sim.run(until_after_sources=5000)

    sim.display_fluxes(wvg_flux)
    sim.display_forces(wvg_force)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-xodd', default=False, action='store_true', help='odd mirror symmetry in X direction? (default: False)')
    parser.add_argument('-s', type=float, default=1.0, help='waveguide separation distance (default: 1.0 μm)')
    args = parser.parse_args()
    main(args)
