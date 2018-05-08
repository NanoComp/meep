import meep as mp
import argparse

import sys
sys.path.insert(0, '/home/oskooi/install/meep5/python/examples')
from materials_library import *

def main(args):

    resolution = args.res

    sz = args.sz
    cell_size = mp.Vector3(0,0,sz)
    
    lambda_min = 0.4
    lambda_max = 0.8
    fmax = 1/lambda_min
    fmin = 1/lambda_max
    fcen = 0.5*(fmax+fmin)
    df = fmax-fmin
    
    dpml = 1.0
    pml_layers = [ mp.PML(dpml) ]
    
    sources = [ mp.Source(mp.GaussianSource(fcen,fwidth=df), component=mp.Ex, center=mp.Vector3(0,0,-0.5*sz+dpml)) ]

    if args.empty:
        geometry = []
    else:
        geometry = [ mp.Block(mp.Vector3(mp.inf,mp.inf,0.5*sz), center=mp.Vector3(0,0,0.25*sz), material=fused_quartz) ]

    sim = mp.Simulation(cell_size=cell_size,
                        boundary_layers=pml_layers,
                        geometry=geometry,
                        sources=sources,
                        dimensions=1,
                        k_point=mp.Vector3(0,0,0),
                        resolution=resolution)

    refl_fr = mp.FluxRegion(center=mp.Vector3(0,0,-0.25*sz))

    nfreq = 50
    refl = sim.add_flux(fcen, df, nfreq, refl_fr)

    if not args.empty:
        sim.load_minus_flux('refl-flux', refl)

    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ex, mp.Vector3(0,0,-0.5*sz+dpml), 1e-9))

    if args.empty:
        sim.save_flux('refl-flux', refl)

    sim.display_fluxes(refl)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-empty', action='store_true', default=False, help="empty? (default: False)")
    parser.add_argument('-res', type=int, default=500, help='resolution (default: 500 pixels/um)')
    parser.add_argument('-sz', type=float, default=10, help='cell size (default: 10 um)')
    args = parser.parse_args()
    main(args)
