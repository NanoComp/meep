---
# GDSII Import
---

This tutorial demonstrates how to set up a simulation based on importing a GDSII file. The example involves a silicon directional coupler. These devices are used as components in [photonic integrated circuits](https://en.wikipedia.org/wiki/Photonic_integrated_circuit) to split and combine an input signal. For more information on directional couplers, see Section 4.1 of [Silicon Photonics Design](https://www.amazon.com/Silicon-Photonics-Design-Devices-Systems/dp/1107085454) by Chrostowski and Hochberg.

---
## GDSII File
---

The directional coupler geometry we will investigate is described by the GDSII file [`coupler.gds`](https://github.com/stevengj/meep/blob/master/python/examples/examples/coupler.gds). A snapshot of this file viewed using the free and open-source tool [KLayout](https://www.klayout.de/) is shown below. The figure labels have been added in post processing. The design consists of two
identical waveguides which are adiabatically tapered to lie close together such that their modes couple evanescently. An input pulse from Port 1 is split in two and exits through Ports 3 and 4. The focus is to find the separation distance (d) which maximizes power in Port 4 at the telecommunication wavelength of 1.55 μm. More generally, it is possible to have two additional degrees of freedom: (1) the length of the straight waveguide sections in the two regions where the two waveguides are closest and farthest away from each other, and (2) the length of the tapered section (the taper profile is often described by a hyperbolic tangent (tanh) function).

<center>
![](../images/klayout_schematic.png)
</center>

The GDSII file is based on the [SiEPIC EBeam PDK](https://github.com/lukasc-ubc/SiEPIC_EBeam_PDK) with four major modifications:

+ the computational cell is centered at the origin of the XY plane and defined as a separate layer (0)

+ the straight waveguide sections have been made perfectly flat

+ the lower and upper branches of the coupler have been placed on separate layers (31 and 32); this facilitates varying the separation between the branches

+ the eigenmode source and flux monitors are each defined on separate layers (1-5)

As an alternative, the volume regions of the source and flux monitors could have been specified in the simulation script rather than as part of the GDSII file.

---
## Simulation Script
---

The simulation script is [`coupler.py`](https://github.com/stevengj/meep/blob/master/python/examples/coupler.py).

```python
import meep as mp
import argparse

resolution = 25 # pixels/um

gdsII_file = 'coupler.gds'
CELL_LAYER = 0
PORT1_LAYER = 1
PORT2_LAYER = 2
PORT3_LAYER = 3
PORT4_LAYER = 4
SOURCE_LAYER = 5
UPPER_BRANCH_LAYER = 31
LOWER_BRANCH_LAYER = 32
default_d = 0.3

t_oxide = 1.0
t_Si = 0.22
t_air = 0.78

eps_oxide = 2.25
oxide = mp.Medium(epsilon=eps_oxide)
eps_silicon = 12
silicon=mp.Medium(epsilon=eps_silicon)

dpml = 1

lcen = 1.55
fcen = 1/lcen
df = 0.2*fcen

# extract center and size of meep::volume
def get_center_and_size(v):
    rmin = v.get_min_corner()
    rmax = v.get_max_corner()
    v3rmin = mp.Vector3(rmin.x(),rmin.y(),rmin.z())
    v3rmax = mp.Vector3(rmax.x(),rmax.y(),rmax.z())
    if v.dim<mp.D3:
      v3rmin.z = v3rmax.z=0
    if v.dim<mp.D2:
      v3rmin.y = v3rmax.y=0
    center = 0.5*(v3rmin+v3rmax)
    size = v3rmax-v3rmin
    return (center,size)

def main(args):
    d  = args.d

    cell_thickness = dpml+t_oxide+t_Si+t_air+dpml
    cell_zmax = 0.5*cell_thickness if args.three_d else 0
    cell_zmin = -0.5*cell_thickness if args.three_d else 0
    si_zmin = 0
    si_zmax = t_Si if args.three_d else 0

    # read cell size, volumes for source region and flux monitors,
    # and coupler geometry from GDSII file
    upper_branch = mp.get_GDSII_prisms(silicon, gdsII_file, UPPER_BRANCH_LAYER, si_zmin, si_zmax)
    lower_branch = mp.get_GDSII_prisms(silicon, gdsII_file, LOWER_BRANCH_LAYER, si_zmin, si_zmax)

    (cell_center,cell_size) = get_center_and_size(mp.get_GDSII_volume(gdsII_file,CELL_LAYER, cell_zmin, cell_zmax))
    (p1_center,p1_size) = get_center_and_size(mp.get_GDSII_volume(gdsII_file,PORT1_LAYER, si_zmin, si_zmax))
    (p2_center,p2_size) = get_center_and_size(mp.get_GDSII_volume(gdsII_file,PORT2_LAYER, si_zmin, si_zmax))
    (p3_center,p3_size) = get_center_and_size(mp.get_GDSII_volume(gdsII_file,PORT3_LAYER, si_zmin, si_zmax))
    (p4_center,p4_size) = get_center_and_size(mp.get_GDSII_volume(gdsII_file,PORT4_LAYER, si_zmin, si_zmax))
    (src_center,src_size) = get_center_and_size(mp.get_GDSII_volume(gdsII_file,SOURCE_LAYER, si_zmin, si_zmax))

    # displace upper and lower branches of coupler (as well as source and flux regions)
    if d != default_d:
        delta_y = 0.5*(d-default_d)
        delta = mp.Vector3(y=delta_y)
        p1_center += delta
        p2_center -= delta
        p3_center += delta
        p4_center -= delta
        src_center += delta
        cell_size += 2*delta
        for np in range(len(lower_branch)):
            lower_branch[np].center -= delta
            for nv in range(len(lower_branch[np].vertices)):
                lower_branch[np].vertices[nv] -= delta
        for np in range(len(upper_branch)):
            upper_branch[np].center += delta
            for nv in range(len(upper_branch[np].vertices)):
                upper_branch[np].vertices[nv] += delta

    geometry = upper_branch+lower_branch

    if args.three_d:
        oxide_center = mp.Vector3(z=-0.5*t_oxide)
        oxide_size = mp.Vector3(cell_size.x,cell_size.y,t_oxide)
        oxide_layer = [mp.Block(material=oxide, center=oxide_center, size=oxide_size)]
        geometry = geometry+oxide_layer

    sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen,fwidth=df),
                                  size=src_size,
                                  center=src_center,
                                  eig_match_freq=True)]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=[mp.PML(dpml)],
                        sources=sources,
                        geometry=geometry)

    p1_region = mp.FluxRegion(center=p1_center,size=p1_size)
    flux1 = sim.add_flux(fcen,0,1,p1_region)
    p2_region = mp.FluxRegion(center=p2_center,size=p2_size)
    flux2 = sim.add_flux(fcen,0,1,p2_region)
    p3_region = mp.FluxRegion(center=p3_center,size=p3_size)
    flux3 = sim.add_flux(fcen,0,1,p3_region)
    p4_region = mp.FluxRegion(center=p4_center,size=p4_size)
    flux4 = sim.add_flux(fcen,0,1,p4_region)

    sim.run(until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,p3_center,1e-8))

    p1_flux = mp.get_fluxes(flux1)
    p2_flux = mp.get_fluxes(flux2)
    p3_flux = mp.get_fluxes(flux3)
    p4_flux = mp.get_fluxes(flux4)

    mp.master_printf("data:, {}, {}, {}, {}".format(d,-p2_flux[0]/p1_flux[0],p3_flux[0]/p1_flux[0],p4_flux[0]/p1_flux[0]))
            
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d',type=float,default=0.1,help='branch separation (default: 0.1 um)')
    parser.add_argument('--three_d',action='store_true',default=False,help='3d calculation? (default: false)')
    args = parser.parse_args()
    main(args)
```

For a given separation distance (`d`) as an input parameter, the simulation computes the fraction of the incident power exiting Ports 2-4. Each of the eight layers of the GDSII file is converted into a simulation object: the upper and lower branches are defined as a [Prism](https://meep.readthedocs.io/en/latest/Python_User_Interface.md#prism), the rectilinear source and flux monitor regions as a [Volume](https://meep.readthedocs.io/en/latest/Python_User_Interface.md#volume) and [FluxRegion](https://meep.readthedocs.io/en/latest/Python_User_Interface.md#fluxregion). The default simulation is 2d. An optional parameter (`three_d`) converts the geometry to 3d by extruding the coupler geometry in the Z direction and adding an oxide layer below. A schematic of the 3d geometry is shown below.

We compute the coupler properties for a range of separation distances, 0.02 μm to 0.30 μm in increments of 0.02 μm, via the following commands:

```
for d in `seq 0.02 0.02 0.30`; do
    mpirun -np 2 python coupler.py -d ${d} |tee -a directional_coupler.out;
done

grep data: directional_coupler.out |cut -d , -f2- > directional_coupler.dat;
```

The results are plotted below. When the two branches are sufficiently separated (`d` > 0.2 μm), practically all of the power from the source in Port 1 is transferred to Port 3. Some of the input power is lost due to scattering in the tapered sections where the translational symmetry of the incident waveguide structure is broken. For separation distances less than approximately 0.2 μm, evanescent coupling of the modes in the two branches transfers some of the input signal into Port 4. At a `d` of 0.06 μm, the power in Port 4 is maximized. For `d` less than 0.06 μm, the evanescent coupling becomes ineffective. Note that there is no power in Port 2.

<center>
![](../images/directional_coupler_flux.png)
</center>

<center>
![](../images/coupler3D.png)
</center>