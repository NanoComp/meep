import meep as mp
import argparse

from meep import mpb
mpb.verbosity(0)

parser = argparse.ArgumentParser()
parser.add_argument('-res', type=int,
                    default=50,
                    help='resolution (default: 50 pixels/um)')
parser.add_argument('-norm',
                    action='store_true',
                    default=False,
                    help='normalization of incident fields? (default: false)')
args = parser.parse_args()

dpml = 1.0

silicon = mp.Medium(epsilon=12)

fcen = 1/1.55
df = 0.2*fcen

sxy = 5.0

cell_size = mp.Vector3(sxy,sxy,0)

boundary_layers = [mp.PML(thickness=dpml)]

eig_parity = mp.EVEN_Y+mp.ODD_Z

symmetries = [mp.Mirror(mp.Y)]

w = 0.6 # waveguide width

geometry = [mp.Block(material=silicon,
                     center=mp.Vector3(),
                     size=mp.Vector3(mp.inf,w,mp.inf))]

sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen,fwidth=df),
                              center=mp.Vector3(-0.5*sxy+dpml,0),
                              size=mp.Vector3(0,sxy,0),
                              eig_band=1,
                              eig_parity=eig_parity,
                              eig_match_freq=True)]

sim = mp.Simulation(resolution=args.res,
                    cell_size=cell_size,
                    boundary_layers=boundary_layers,
                    sources=sources,
                    geometry=geometry,
                    symmetries=symmetries)

mode = sim.add_flux(fcen, 0, 1,
                    mp.ModeRegion(center=mp.Vector3(0.5*sxy-dpml),size=mp.Vector3(0,sxy,0)))

sim.run(until_after_sources=20)

input_data = sim.get_flux_data(mode)
input_flux = mp.get_fluxes(mode)[0]

sim.reset_meep()

sim = mp.Simulation(resolution=args.res,
                    cell_size=cell_size,
                    boundary_layers=boundary_layers,
                    sources=sources,
                    geometry=geometry,
                    symmetries=symmetries)

mode = sim.add_flux(fcen, 0, 1,
                    mp.ModeRegion(center=mp.Vector3(0.5*sxy-dpml),size=mp.Vector3(0,sxy,0)))

if args.norm:
    sim.load_minus_flux_data(mode, input_data)

sim.run(until_after_sources=20)

# mode coefficients
coeff_minus = sim.get_eigenmode_coefficients(mode,[1],eig_parity).alpha[0,0,1]

# S parameters
S11 = abs(coeff_minus)**2/input_flux

print("S11:, {}, {}".format(args.res,S11))
