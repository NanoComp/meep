import numpy as np
import gdspy
from matplotlib import pyplot as plt
import importlib
import meep as mp

# core and cladding materials
Si   = mp.Medium(index=3.4)
SiO2 = mp.Medium(index=1.4)

# layer numbers for GDS file
RING_LAYER       = 0
SOURCE0_LAYER    = 1
SOURCE1_LAYER    = 2
MONITOR_LAYER    = 3
SIMULATION_LAYER = 4

resolution = 50         # pixels/μm
dpml       = 1          # thickness of PML
zmin       = 0          # minimum z value of simulation domain (0 for 2D)
zmax       = 0          # maximum z value of simulation domain (0 for 2D)

def create_ring_gds(radius=5,waveguideWidth=0.5):
    # Reload the library each time to prevent gds library name clashes
    importlib.reload(gdspy)

    # Draw the ring
    ringCell = gdspy.Cell("ring_r{}_w{}".format(radius,waveguideWidth))

    ringCell.add(gdspy.Round((0,0),
                             inner_radius=radius-waveguideWidth/2,
                             radius=radius+waveguideWidth/2,
                             layer=RING_LAYER))

    # Draw the first source
    ringCell.add(gdspy.Rectangle((radius-waveguideWidth,0),
                                 (radius+waveguideWidth+0.1,0),
                                 SOURCE0_LAYER))

    # Draw the second source
    ringCell.add(gdspy.Rectangle((-radius+waveguideWidth,0),
                                 (-radius-waveguideWidth,0),
                                 SOURCE1_LAYER))

    # Draw the monitor location
    ringCell.add(gdspy.Rectangle((radius-waveguideWidth,0),
                                 (radius+waveguideWidth,0),
                                 MONITOR_LAYER))

    # Draw the simulation domain
    pad = 2  # padding between waveguide and edge of PML
    ringCell.add(gdspy.Rectangle((-radius-waveguideWidth/2-pad,-radius-waveguideWidth/2-pad),
                                 (radius+waveguideWidth/2+pad,radius+waveguideWidth/2+pad),
                                 SIMULATION_LAYER))

    filename = "ring_r{}_w{}.gds".format(radius,waveguideWidth)
    gdspy.write_gds(filename, unit=1.0e-6, precision=1.0e-9)

    return filename

def run_sim(filename,wavelengthCenter=1.55,bandwidth=0.05):
    # Read in the ring structure
    geometry = mp.get_GDSII_prisms(Si,filename,RING_LAYER,-100,100)

    cell = mp.GDSII_vol(filename,SIMULATION_LAYER,zmin,zmax)

    src_vol0 = mp.GDSII_vol(filename,SOURCE0_LAYER,zmin,zmax)
    src_vol1 = mp.GDSII_vol(filename,SOURCE1_LAYER,zmin,zmax)

    mon_vol = mp.GDSII_vol(filename,MONITOR_LAYER,zmin,zmax)

    # pulse center frequency
    fcen = 1/wavelengthCenter
    # pulse frequency width
    df = bandwidth*fcen

    src = [mp.Source(mp.GaussianSource(fcen, fwidth=df),
                     component=mp.Hz,
                     center=src_vol0.center,
                     size=src_vol0.size),
           mp.Source(mp.GaussianSource(fcen, fwidth=df),
                     component=mp.Hz,
                     center=src_vol1.center,
                     size=src_vol1.size,
                     amplitude=-1)]

    sim = mp.Simulation(cell_size=cell.size,
                        geometry=geometry,
                        sources=src,
                        resolution=resolution,
                        boundary_layers=[mp.PML(dpml)],
                        default_material=SiO2)

    h = mp.Harminv(mp.Hz,mon_vol.center,fcen,df)

    sim.run(mp.after_sources(h),
            until_after_sources=100)

    sim.plot2D(fields=mp.Hz)
    plt.savefig('ring_resonator_gds_Hz.png')

    freq = np.array([1/m.freq for m in h.modes])
    Q = np.array([m.Q for m in h.modes])

    sim.reset_meep()

    return freq, Q


if __name__ == '__main__':
    filename = create_ring_gds(radius=2.0,waveguideWidth=0.5)
    freq, Q = run_sim(filename)
    print("freq: {}".format(freq))
    print("Q: {}".format(Q))
    
