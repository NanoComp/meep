import meep as mp

cell_size = mp.Vector3(6, 6, 0)

geometry1 = [
    mp.Cylinder(center=mp.Vector3(), radius=1.0, material=mp.Medium(index=3.5))
]

sim1 = mp.Simulation(cell_size=cell_size, geometry=geometry1, resolution=20)

sim1.init_sim()

geometry2 = [
    mp.Cylinder(center=mp.Vector3(1, 1), radius=1.0, material=mp.Medium(index=3.5))
]

sim2 = mp.Simulation(cell_size=cell_size, geometry=geometry2, resolution=20)

sim2.init_sim()

sim1.fields.phase_in_material(sim2.structure, 10.0)

sim1.run(
    mp.at_beginning(mp.output_epsilon), mp.at_every(0.5, mp.output_epsilon), until=10
)
