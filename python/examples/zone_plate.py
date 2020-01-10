import meep as mp
import numpy as np
import math
import matplotlib.pyplot as plt

resolution = 25             # pixels/μm

dpml = 1.0                  # PML thickness
dsub = 2.0                  # substrate thickness
dpad = 2.0                  # padding betweeen zone plate and PML
zh = 0.5                    # zone-plate height
zN = 25                     # number of zones (odd zones: π phase shift, even zones: none)
focal_length = 200          # focal length of zone plate

spot_length = 100           # far-field line length
ff_npts = 100               # number of far-field points to compute along spot_length

pml_layers = [mp.PML(thickness=dpml)]

wvl_cen = 0.5
frq_cen = 1/wvl_cen
dfrq = 0.2*frq_cen

## ref: eq. 7 of http://zoneplate.lbl.gov/theory
r = [math.sqrt(n*wvl_cen*(focal_length+n*wvl_cen/4)) for n in range(1,zN+1)]

sr = r[-1]+dpad+dpml
sz = dpml+dsub+zh+dpad+dpml
cell_size = mp.Vector3(sr,0,sz)

sources = [mp.Source(mp.GaussianSource(frq_cen,fwidth=dfrq,is_integrated=True),
                     component=mp.Er,
                     center=mp.Vector3(0.5*sr,0,-0.5*sz+dpml),
                     size=mp.Vector3(sr)),
           mp.Source(mp.GaussianSource(frq_cen,fwidth=dfrq,is_integrated=True),
                     component=mp.Ep,
                     center=mp.Vector3(0.5*sr,0,-0.5*sz+dpml),
                     size=mp.Vector3(sr),
                     amplitude=-1j)]

glass = mp.Medium(index=1.5)

geometry = [mp.Block(material=glass,
                     size=mp.Vector3(sr,0,dpml+dsub),
                     center=mp.Vector3(0.5*sr,0,-0.5*sz+0.5*(dpml+dsub)))]

for n in range(zN-1,-1,-1):
    geometry.append(mp.Block(material=glass if n % 2 == 0 else mp.vacuum,
                             size=mp.Vector3(r[n],0,zh),
                             center=mp.Vector3(0.5*r[n],0,-0.5*sz+dpml+dsub+0.5*zh)))

sim = mp.Simulation(cell_size=cell_size,
                    boundary_layers=pml_layers,
                    resolution=resolution,
                    sources=sources,
                    geometry=geometry,
                    dimensions=mp.CYLINDRICAL,
                    m=-1)

n2f_obj = sim.add_near2far(frq_cen, 0, 1, mp.Near2FarRegion(center=mp.Vector3(0.5*sr,0,0.5*sz-dpml),size=mp.Vector3(sr)))

sim.run(until_after_sources=100)

ff = sim.get_farfields(n2f_obj, ff_npts/spot_length, center=mp.Vector3(z=-0.5*sz+dpml+dsub+zh+focal_length),size=mp.Vector3(z=spot_length))

energy_density = np.absolute(ff['Ex'])**2+np.absolute(ff['Ey'])**2+np.absolute(ff['Ez'])**2

z = np.linspace(focal_length-0.5*spot_length,focal_length+0.5*spot_length,ff_npts)

if mp.am_master():
    plt.figure(dpi=200)
    plt.semilogy(z,energy_density,'bo-')
    plt.grid(True,axis="y",which="both",ls="-")
    plt.xlabel('z coordinate (μm)')
    plt.ylabel(r'energy density of far-field electric fields, |E|$^2$')
    plt.tight_layout()
    plt.savefig("zone_plate_farfields.png")
