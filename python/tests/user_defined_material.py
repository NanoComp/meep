import os
import meep as mp


# Material function that recreates the ellipsoid-in-cylinder configuration of
# examples/cyl-ellipsoid.py
def my_material_func(p, user_data, m):
    R1X = 0.5
    R1Y = 1.0
    R2 = 3.0

    x = p.x
    y = p.y

    # test for point inside inner ellipsoid
    if (x**2 / (R1X**2) + y**2 / (R1Y**2)) < 1.0:
        nn = 1.0
    elif (x**2 / (R2**2) + y**2 / (R2**2)) < 1.0:
        nn = 3.5
    else:
        nn = 1.0

    m.epsilon_diag.x = nn**2
    m.epsilon_diag.y = nn**2
    m.epsilon_diag.z = nn**2


def main():
    eps_ref_file = "cyl-ellipsoid-eps-ref.h5"
    eps_ref_dir = os.path.abspath(os.path.realpath(os.path.join(__file__, '..', '..', 'libmeepgeom')))
    eps_ref_path = os.path.join(eps_ref_dir, eps_ref_file)

    resolution = 100
    cell = mp.Vector3(10, 10)
    symmetries = [mp.Mirror(mp.X), mp.Mirror(mp.Y)]
    boundary_layers = [mp.PML(1.0)]

    sim = mp.Simulation(cell_size=cell,
                        resolution=resolution,
                        symmetries=symmetries,
                        boundary_layers=boundary_layers,
                        material_function=my_material_func)

    sim.init_fields()

    sim.fields.output_hdf5(mp.Dielectric, sim.fields.total_volume())

if __name__ == '__main__':
    main()
