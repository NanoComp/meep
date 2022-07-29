import importlib

import gdspy
import numpy as np
from matplotlib import pyplot as plt

import meep as mp

# core and cladding materials
Si = mp.Medium(index=3.4)
SiO2 = mp.Medium(index=1.4)

# layer numbers for GDS file
RING_LAYER = 0
SOURCE0_LAYER = 1
SOURCE1_LAYER = 2
MONITOR_LAYER = 3
SIMULATION_LAYER = 4

resolution = 50  # pixels/μm
dpml = 1  # thickness of PML
zmin = 0  # minimum z value of simulation domain (0 for 2D)
zmax = 0  # maximum z value of simulation domain (0 for 2D)


def create_ring_gds(radius, width):
    # Reload the library each time to prevent gds library name clashes
    importlib.reload(gdspy)

    ringCell = gdspy.Cell(f"ring_resonator_r{radius}_w{width}")

    # Draw the ring
    ringCell.add(
        gdspy.Round(
            (0, 0),
            inner_radius=radius - width / 2,
            radius=radius + width / 2,
            layer=RING_LAYER,
        )
    )

    # Draw the first source
    ringCell.add(
        gdspy.Rectangle((radius - width, 0), (radius + width, 0), SOURCE0_LAYER)
    )

    # Draw the second source
    ringCell.add(
        gdspy.Rectangle((-radius - width, 0), (-radius + width, 0), SOURCE1_LAYER)
    )

    # Draw the monitor location
    ringCell.add(
        gdspy.Rectangle((radius - width / 2, 0), (radius + width / 2, 0), MONITOR_LAYER)
    )

    # Draw the simulation domain
    pad = 2  # padding between waveguide and edge of PML
    ringCell.add(
        gdspy.Rectangle(
            (-radius - width / 2 - pad, -radius - width / 2 - pad),
            (radius + width / 2 + pad, radius + width / 2 + pad),
            SIMULATION_LAYER,
        )
    )

    filename = f"ring_r{radius}_w{width}.gds"
    gdspy.write_gds(filename, unit=1.0e-6, precision=1.0e-9)

    return filename


def find_modes(filename, wvl=1.55, bw=0.05):
    # Read in the ring structure
    geometry = mp.get_GDSII_prisms(Si, filename, RING_LAYER, -100, 100)

    cell = mp.GDSII_vol(filename, SIMULATION_LAYER, zmin, zmax)

    src_vol0 = mp.GDSII_vol(filename, SOURCE0_LAYER, zmin, zmax)
    src_vol1 = mp.GDSII_vol(filename, SOURCE1_LAYER, zmin, zmax)

    mon_vol = mp.GDSII_vol(filename, MONITOR_LAYER, zmin, zmax)

    fcen = 1 / wvl
    df = bw * fcen

    src = [
        mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Hz, volume=src_vol0),
        mp.Source(
            mp.GaussianSource(fcen, fwidth=df),
            component=mp.Hz,
            volume=src_vol1,
            amplitude=-1,
        ),
    ]

    sim = mp.Simulation(
        cell_size=cell.size,
        geometry=geometry,
        sources=src,
        resolution=resolution,
        boundary_layers=[mp.PML(dpml)],
        default_material=SiO2,
    )

    h = mp.Harminv(mp.Hz, mon_vol.center, fcen, df)

    sim.run(mp.after_sources(h), until_after_sources=100)

    plt.figure()
    sim.plot2D(fields=mp.Hz, eps_parameters={"contour": True})
    plt.savefig("ring_fields.png", bbox_inches="tight", dpi=150)

    wvl = np.array([1 / m.freq for m in h.modes])
    Q = np.array([m.Q for m in h.modes])

    sim.reset_meep()

    return wvl, Q


if __name__ == "__main__":
    filename = create_ring_gds(2.0, 0.5)
    wvls, Qs = find_modes(filename, 1.55, 0.05)
    for w, Q in zip(wvls, Qs):
        print(f"mode: {w}, {Q}")
