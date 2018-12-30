from __future__ import division
import math
import meep as mp
import numpy as np
import argparse
import sys

##################################################
# compute modal volume from frequency-domain fields
# (if dft_cell is non-null) or time-domain fields
# (otherwise).
##################################################
def get_modal_volume(sim, box=None, dft_cell=None, nf=0):

    Exyz=[mp.Ex, mp.Ey, mp.Ez]

    if dft_cell is None:
      (Ex,Ey,Ez) = [sim.get_array(vol=box, component=c, cmplx=True) for c in Exyz]
      (X,Y,Z,W)  = sim.get_array_metadata(vol=box)
      Eps        = sim.get_array(vol=box, component=mp.Dielectric)
    else:
      (Ex,Ey,Ez) = [sim.get_dft_array(dft_cell, c, nf) for c in Exyz]
      (X,Y,Z,W)  = sim.get_dft_array_metadata(dft_cell=dft_cell)
      # slightly annoying: we need an epsilon array with empty dimensions collapsed,
      #  something not currently provided by any C++ function; for now just
      #  create it via brute-force python loop, although this could get slow.
      Eps=np.zeros(0)
      for x in X:
        for y in Y:
          for z in Z:
            Eps=np.append(Eps,sim.get_epsilon_point(mp.Vector3(x,y,z)))
      Eps=np.reshape(Eps,np.shape(W))

    EpsE2 = np.real(Eps*(np.conj(Ex)*Ex + np.conj(Ey)*Ey + np.conj(Ez)*Ez))
    num   = np.sum(W*EpsE2)
    denom = np.max(EpsE2)

    return num/denom if denom!=0.0 else 0.0

##################################################
# utility routine for color-highlighted output 
##################################################
Highlight1, Highlight2, ErrorText, PlainText='', '', '', ''
try:
    from colorama import init, Fore, Back, Style
    init()
    Highlight1 = Style.BRIGHT + Fore.WHITE + Back.CYAN
    Highlight2 = Style.BRIGHT + Fore.WHITE + Back.MAGENTA
    ErrorText  = Style.BRIGHT + Fore.RED   + Back.WHITE
    PlainText  = Style.RESET_ALL
except ImportError:
    pass

##################################################
# parse arguments
##################################################
parser = argparse.ArgumentParser()
parser.add_argument('--res',          type=float, default=20,    help='resolution')
parser.add_argument('--w',            type=float, default=0,     help='width of cavity opening')
parser.add_argument('--use_symmetry',             default=False, action='store_true')

args         = parser.parse_args()
resolution   = args.res
w            = args.w
use_symmetry = args.use_symmetry

sym_str = "True" if use_symmetry else "False"
print "use_symmetries is " + sym_str

##################################################
# define geometry
##################################################
sxy    = 2
dpml   = 1
sxy    = sxy + 2 * dpml
cell   = mp.Vector3(sxy, sxy, 0)
origin = mp.Vector3()

a = 1       # cavity width
t = 0.1     # cavity wall thickness

L_inner=a
cavity_inner_size = mp.Vector3(L_inner, L_inner, mp.inf)

L_outer=a+2*t
cavity_outer_size = mp.Vector3(L_outer, L_outer, mp.inf)

geometry = [ mp.Block( cavity_outer_size, material=mp.metal),
             mp.Block( cavity_inner_size, material=mp.air)
           ]

if w > 0:
    geometry.append(mp.Block( center=mp.Vector3(a/2),
                              size=mp.Vector3(2 * t, w, mp.inf),
                              material=mp.air
                            )
                   )

##################################################
# add sources
##################################################
fcen = math.sqrt(0.5) / a
df = 0.2
source_time_profile = mp.GaussianSource(fcen, fwidth=df)
sources = [mp.Source(src=source_time_profile, component=mp.Ez, center=origin) ]

symmetries = [mp.Mirror(mp.Y)] if use_symmetry else []

pml_layers = [mp.PML(dpml)]
sim = mp.Simulation(cell_size=cell, geometry=geometry,
                    boundary_layers=pml_layers, sources=sources,
                    symmetries=symmetries, resolution=resolution)

# mv_box is the region within which we compute modal volume
mv_box_size = mp.Vector3(L_inner, L_inner, mp.inf)
mv_box      = mp.Volume(center=origin, size=mv_box_size)
dft_cell    = sim.add_dft_fields([mp.Ex, mp.Ey, mp.Ez], fcen-df, fcen+df, 3, where=mv_box)

# Timestep until the sources are finished, pausing at fixed intervals to 
# compare the modal volume within mv_box as computed from time-domain fields 
# (a) by the built-in mode-volume computation in libmeep,
# (b) by the python routine above based on array metadata
source_end_time       = source_time_profile.swigobj.last_time_max() 
measurement_interval  = source_end_time/10.0
next_measurement_time = sim.round_time() + 2*measurement_interval # skip first interval
td_values=0
td_errors=0
RelTol=0.1 # only look for errors greater than 10%
while sim.round_time() < source_end_time:
    sim.run(until=next_measurement_time)
    next_measurement_time+=measurement_interval
    tdmv_metadata = get_modal_volume(sim, box=mv_box)   # 'time-domain model volume' computed via method (b)
    tdmv_meep     = sim.modal_volume_in_box(box=mv_box) # 'time-domain modal volume' computed via method (a)
    TestFailed = (abs(tdmv_metadata - tdmv_meep) > RelTol*abs(tdmv_meep))
    td_errors += 1 if TestFailed else 0
    sys.stdout.write(Highlight1)
    sys.stdout.write("t%i (%.3e):   mv(metadata)=%.4e   mv(meep)=%.4e   " % (td_values, sim.round_time(), tdmv_metadata, tdmv_meep) )
    sys.stdout.write(ErrorText + "FAIL" if TestFailed else "PASS")
    sys.stdout.write(PlainText + "\n")
    td_values+=1

# compute modal volume from frequency-domain fields
sys.stdout.write("\n**\n** Frequency-domain modal volumes:\n**\n\n")
for nf in range(0,3):
    freq=dft_cell.freq_min + nf*dft_cell.dfreq
    fdmv=get_modal_volume(sim, dft_cell=dft_cell, nf=nf)
    sys.stdout.write(Highlight1)
    sys.stdout.write("f%i (%.3e):   mv(metadata)=%.4e" % (nf, freq, fdmv));
    sys.stdout.write(PlainText + "\n")

if td_errors>0:
    sys.stdout.write(ErrorText + "Test FAILED with %i/%i errors in time-domain modal volumes" % (td_errors,td_values))
    sys.stdout.write(PlainText + "\n")
else:
    sys.stdout.write("Test PASSED with 0/%i errors\n" % (td_values) )

#sys.exit(td_errors)
