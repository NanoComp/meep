import meep as mp

# Example file illustrating an eigenmode source, generating a waveguide mode
# (requires recent MPB version to be installed before Meep is compiled)

cell = mp.Vector3(16, 8)

# an asymmetrical dielectric waveguide:
geometry = [
    mp.Block(
        center=mp.Vector3(),
        size=mp.Vector3(mp.inf, 1, mp.inf),
        material=mp.Medium(epsilon=12),
    ),
    mp.Block(
        center=mp.Vector3(y=0.3),
        size=mp.Vector3(mp.inf, 0.1, mp.inf),
        material=mp.Medium(),
    ),
]

# create a transparent source that excites a right-going waveguide mode
sources = [
    mp.EigenModeSource(
        src=mp.ContinuousSource(0.15),
        size=mp.Vector3(y=6),
        center=mp.Vector3(x=-5),
        component=mp.Dielectric,
        eig_parity=mp.ODD_Z,
    )
]

pml_layers = [mp.PML(1.0)]

force_complex_fields = True  # so we can get time-average flux

resolution = 10

sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    boundary_layers=pml_layers,
    force_complex_fields=force_complex_fields,
    resolution=resolution,
)

sim.run(
    mp.at_beginning(mp.output_epsilon),
    mp.at_end(mp.output_png(mp.Ez, "-a yarg -A $EPS -S3 -Zc dkbluered", rm_h5=False)),
    until=200,
)

flux1 = sim.flux_in_box(
    mp.X, mp.Volume(center=mp.Vector3(-6.0), size=mp.Vector3(1.8, 6))
)
flux2 = sim.flux_in_box(
    mp.X, mp.Volume(center=mp.Vector3(6.0), size=mp.Vector3(1.8, 6))
)

# averaged over y region of width 1.8
print(f"left-going flux = {flux1 / -1.8}")

# averaged over y region of width 1.8
print(f"right-going flux = {flux2 / 1.8}")
