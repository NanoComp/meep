import meep as mp

wavelength = 0.85
frequency = 1/wavelength

cell = mp.Vector3(16,8,0)
geometry = [mp.Block(mp.Vector3(1,8), center=mp.Vector3(), material=mp.Medium(index=3.5))]

sources = [mp.Source(mp.ContinuosSource, frequency=frequency), component=mp.Ez, center=mp.Vector3(-6,0)]

sim = mp.Simulation(cell_size=cell, boundary_layers=[mp.PML(1.0)], geometry=geometry, sources=sources, resolution=50)

sim.run(until=200)
