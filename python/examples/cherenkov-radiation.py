## moving point charge with superluminal phase velocity in dielectric media emitting Cherenkov radiation
import meep as mp

sx = 60
sy = 60
cell_size = mp.Vector3(sx, sy, 0)

dpml = 1.0
pml_layers = [mp.PML(thickness=dpml)]

v = 0.7  # velocity of point charge

symmetries = [mp.Mirror(direction=mp.Y)]

sim = mp.Simulation(
    resolution=10,
    cell_size=cell_size,
    default_material=mp.Medium(index=1.5),
    symmetries=symmetries,
    boundary_layers=pml_layers,
)


def move_source(sim):
    sim.change_sources(
        [
            mp.Source(
                mp.ContinuousSource(frequency=1e-10),
                component=mp.Ex,
                center=mp.Vector3(-0.5 * sx + dpml + v * sim.meep_time()),
            )
        ]
    )


sim.run(
    move_source,
    mp.at_every(2, mp.output_png(mp.Hz, "-vZc dkbluered -M 1")),
    until=sx / v,
)
