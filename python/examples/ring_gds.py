import gdstk
from matplotlib import pyplot as plt
import meep as mp
import numpy as np

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


def get_gds_cell(fname):
    """Returns the (single) top-level cell of the GDS file `fname`."""
    return gdstk.read_gds(fname).top_level()[0]


def get_gds_prisms(material, cell, layer, datatype=0, zmin=0.0, zmax=0.0):
    """Returns a list of `mp.Prism`s, one for each polygon on (`layer`, `datatype`)."""
    prisms = []
    for poly in cell.get_polygons(layer=layer, datatype=datatype):
        vertices = [mp.Vector3(x, y, zmin) for x, y in poly.points]
        prisms.append(
            mp.Prism(
                vertices,
                height=zmax - zmin,
                axis=mp.Vector3(0, 0, 1),
                material=material,
            )
        )
    return prisms


def get_gds_vol(cell, layer, datatype=0, zmin=0.0, zmax=0.0):
    """Returns an `mp.Volume` spanning the bounding box of (`layer`, `datatype`)."""
    polygons = cell.get_polygons(layer=layer, datatype=datatype)
    xs = [x for poly in polygons for x, y in poly.points]
    ys = [y for poly in polygons for x, y in poly.points]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    center = mp.Vector3(0.5 * (xmin + xmax), 0.5 * (ymin + ymax), 0.5 * (zmin + zmax))
    size = mp.Vector3(xmax - xmin, ymax - ymin, zmax - zmin)
    dims = 2 if (zmin == 0 and zmax == 0) else 3
    return mp.Volume(center=center, size=size, dims=dims)


def create_ring_gds(radius, width):
    lib = gdstk.Library()
    ring_cell = lib.new_cell(f"ring_resonator_r{radius}_w{width}")

    # Draw the ring
    ring_cell.add(
        gdstk.ellipse(
            (0, 0),
            radius + width / 2,
            inner_radius=radius - width / 2,
            layer=RING_LAYER,
        )
    )

    # Draw the first source
    ring_cell.add(
        gdstk.rectangle((radius - width, 0), (radius + width, 0), layer=SOURCE0_LAYER)
    )

    # Draw the second source
    ring_cell.add(
        gdstk.rectangle((-radius - width, 0), (-radius + width, 0), layer=SOURCE1_LAYER)
    )

    # Draw the monitor location
    ring_cell.add(
        gdstk.rectangle(
            (radius - width / 2, 0), (radius + width / 2, 0), layer=MONITOR_LAYER
        )
    )

    # Draw the simulation domain
    pad = 2  # padding between waveguide and edge of PML
    ring_cell.add(
        gdstk.rectangle(
            (-radius - width / 2 - pad, -radius - width / 2 - pad),
            (radius + width / 2 + pad, radius + width / 2 + pad),
            layer=SIMULATION_LAYER,
        )
    )

    filename = f"ring_r{radius}_w{width}.gds"
    lib.write_gds(filename)

    return filename


def find_modes(filename, wvl=1.55, bw=0.05):
    # Read in the ring structure using gdstk
    gds_cell = get_gds_cell(filename)

    geometry = get_gds_prisms(Si, gds_cell, RING_LAYER, zmin=-100, zmax=100)

    cell = get_gds_vol(gds_cell, SIMULATION_LAYER, zmin=zmin, zmax=zmax)

    src_vol0 = get_gds_vol(gds_cell, SOURCE0_LAYER, zmin=zmin, zmax=zmax)
    src_vol1 = get_gds_vol(gds_cell, SOURCE1_LAYER, zmin=zmin, zmax=zmax)

    mon_vol = get_gds_vol(gds_cell, MONITOR_LAYER, zmin=zmin, zmax=zmax)

    fcen = 1 / wvl
    df = bw * fcen

    src = [
        mp.Source(
            mp.GaussianSource(fcen, fwidth=df),
            component=mp.Hz,
            volume=src_vol0,
        ),
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

    fig, ax = plt.subplots()
    sim.plot2D(ax=ax, fields=mp.Hz, eps_parameters={"contour": True})
    fig.savefig("ring_fields.png", bbox_inches="tight", dpi=150)

    wvl = np.array([1 / m.freq for m in h.modes])
    Q = np.array([m.Q for m in h.modes])

    sim.reset_meep()

    return wvl, Q


if __name__ == "__main__":
    filename = create_ring_gds(2.0, 0.5)
    wvls, Qs = find_modes(filename, 1.55, 0.05)
    for w, Q in zip(wvls, Qs):
        print(f"mode: {w}, {Q}")
