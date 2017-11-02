import meep as mp


# Some parameters to describe the geometry:

# dielectric constant of waveguide
eps = 13

# width of waveguide
w = 1.2

# radius of holes
r = 0.36

# defect spacing (ordinary spacing = 1)
d = 1.4

# number of holes on either side of defect
N = 3

# The cell dimensions

# size of cell in y direction (perpendicular to wvg.)
sy = 6

# padding between last hole and PML edge
pad = 2

# PML thickness
dpml = 1

# size of cell in x direction
sx = 2 * (pad + dpml + N) + d - 1

cell = mp.Vector3(sx, sy, 0)

geometry = [mp.Block(center=mp.Vector3(), size=mp.Vector3(1e20, w, 1e20),
                     material=mp.Medium(epsilon=eps))]

for i in range(N):
    geometry.append(mp.Cylinder(r, center=mp.Vector3(d / 2 + i)))
    geometry.append(mp.Cylinder(r, center=mp.Vector3(d / -2 - i)))

pml_layers = mp.PML(dpml)
resolution = 20

# pulse center frequency
fcen = 0.25

# pulse width (in frequency)
df = 0.2

# number of frequencies at which to compute flux
nfreq = 500

sources = mp.Source(src=mp.GaussianSource(fcen, fwidth=df), component=mp.Hz, center=mp.Vector3())

symmetries = [mp.Mirror(mp.Y, phase=-1), mp.Mirror(mp.X, phase=-1)]

d1 = 0.2

sim = mp.Simulation(cell_size=cell,
                    geometry=geometry,
                    sources=[sources],
                    symmetries=symmetries,
                    boundary_layers=[pml_layers],
                    resolution=resolution)

nearfield = sim.add_near2far(
    fcen, 0, 1,
    mp.Near2FarRegion(mp.Vector3(0, 0.5 * w + d1), size=mp.Vector3(2 * dpml - sx)),
    mp.Near2FarRegion(mp.Vector3(-0.5 * sx + dpml, 0.5 * w + 0.5 * d1), size=mp.Vector3(0, d1), weight=-1.0),
    mp.Near2FarRegion(mp.Vector3(0.5 * sx - dpml, 0.5 * w + 0.5 * d1), size=mp.Vector3(0, d1))
)

sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Hz, mp.Vector3(0.12, -0.37), 1e-8))

d2 = 20
h = 4

sim.output_farfields(nearfield, "spectra-{}-{}-{}".format(d1, d2, h),
                     mp.Volume(mp.Vector3(0, (0.5 * w) + d2 + (0.5 * h)), size=mp.Vector3(sx - 2 * dpml, h)),
                     resolution)
