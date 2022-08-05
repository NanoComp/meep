import meep as mp

try:
    import meep.adjoint as mpa
except:
    import adjoint as mpa

import unittest

import numpy as np
from scipy.ndimage import gaussian_filter


def compute_transmittance(matgrid_symmetry=False):
    resolution = 25

    cell_size = mp.Vector3(6, 6, 0)

    boundary_layers = [mp.PML(thickness=1.0)]

    matgrid_size = mp.Vector3(2, 2, 0)
    matgrid_resolution = 2 * resolution

    Nx, Ny = int(matgrid_size.x * matgrid_resolution), int(
        matgrid_size.y * matgrid_resolution
    )

    # ensure reproducible results
    rng = np.random.RandomState(2069588)

    w = rng.rand(Nx, Ny)
    weights = w if matgrid_symmetry else 0.5 * (w + np.fliplr(w))

    matgrid = mp.MaterialGrid(
        mp.Vector3(Nx, Ny),
        mp.air,
        mp.Medium(index=3.5),
        weights=weights,
        do_averaging=False,
        grid_type="U_MEAN",
    )

    geometry = [
        mp.Block(
            center=mp.Vector3(),
            size=mp.Vector3(mp.inf, 1.0, mp.inf),
            material=mp.Medium(index=3.5),
        ),
        mp.Block(
            center=mp.Vector3(),
            size=mp.Vector3(matgrid_size.x, matgrid_size.y, 0),
            material=matgrid,
        ),
    ]

    if matgrid_symmetry:
        geometry.append(
            mp.Block(
                center=mp.Vector3(),
                size=mp.Vector3(matgrid_size.x, matgrid_size.y, 0),
                material=matgrid,
                e2=mp.Vector3(y=-1),
            )
        )

    eig_parity = mp.ODD_Y + mp.EVEN_Z

    fcen = 0.65
    df = 0.2 * fcen
    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(-2.0, 0),
            size=mp.Vector3(0, 4.0),
            eig_parity=eig_parity,
        )
    ]

    sim = mp.Simulation(
        resolution=resolution,
        cell_size=cell_size,
        boundary_layers=boundary_layers,
        sources=sources,
        geometry=geometry,
    )

    mode_mon = sim.add_flux(
        fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(2.0, 0), size=mp.Vector3(0, 4.0))
    )

    sim.run(until_after_sources=mp.stop_when_dft_decayed())

    mode_coeff = sim.get_eigenmode_coefficients(mode_mon, [1], eig_parity).alpha[
        0, :, 0
    ][0]

    tran = np.power(np.abs(mode_coeff), 2)
    print(f'tran:, {"sym" if matgrid_symmetry else "nosym"}, {tran}')

    return tran


def compute_resonant_mode_2d(res, default_mat=False):
    cell_size = mp.Vector3(1, 1, 0)

    rad = 0.301943

    fcen = 0.3
    df = 0.2 * fcen
    sources = [
        mp.Source(
            mp.GaussianSource(fcen, fwidth=df),
            component=mp.Hz,
            center=mp.Vector3(-0.1057, 0.2094, 0),
        )
    ]

    k_point = mp.Vector3(0.3892, 0.1597, 0)

    matgrid_size = mp.Vector3(1, 1, 0)
    matgrid_resolution = 1200

    # for a fixed resolution, compute the number of grid points
    # necessary which are defined on the corners of the voxels
    Nx, Ny = int(matgrid_size.x * matgrid_resolution), int(
        matgrid_size.y * matgrid_resolution
    )

    x = np.linspace(-0.5 * matgrid_size.x, 0.5 * matgrid_size.x, Nx)
    y = np.linspace(-0.5 * matgrid_size.y, 0.5 * matgrid_size.y, Ny)
    xv, yv = np.meshgrid(x, y)
    weights = np.sqrt(np.square(xv) + np.square(yv)) < rad
    filtered_weights = gaussian_filter(weights, sigma=3.0, output=np.double)

    matgrid = mp.MaterialGrid(
        mp.Vector3(Nx, Ny),
        mp.air,
        mp.Medium(index=3.5),
        weights=filtered_weights,
        do_averaging=True,
        beta=1000,
        eta=0.5,
    )

    geometry = [
        mp.Block(
            center=mp.Vector3(),
            size=mp.Vector3(matgrid_size.x, matgrid_size.y, 0),
            material=matgrid,
        )
    ]

    sim = mp.Simulation(
        resolution=res,
        cell_size=cell_size,
        default_material=matgrid if default_mat else mp.Medium(),
        geometry=[] if default_mat else geometry,
        sources=sources,
        k_point=k_point,
    )

    h = mp.Harminv(mp.Hz, mp.Vector3(0.3718, -0.2076), fcen, df)
    sim.run(mp.after_sources(h), until_after_sources=200)

    try:
        for m in h.modes:
            print(f"harminv:, {res}, {m.freq}, {m.Q}")
        freq = h.modes[0].freq
    except:
        raise RuntimeError("No resonant modes found.")

    return freq


def compute_resonant_mode_3d(use_matgrid=True):
    resolution = 25

    wvl = 1.27
    fcen = 1 / wvl
    df = 0.02 * fcen

    nSi = 3.45
    Si = mp.Medium(index=nSi)
    nSiO2 = 1.45
    SiO2 = mp.Medium(index=nSiO2)

    s = 1.0
    cell_size = mp.Vector3(s, s, s)

    rad = 0.34  # radius of sphere

    if use_matgrid:
        matgrid_resolution = 2 * resolution
        N = int(s * matgrid_resolution)
        coord = np.linspace(-0.5 * s, 0.5 * s, N)
        xv, yv, zv = np.meshgrid(coord, coord, coord)

        weights = np.sqrt(np.square(xv) + np.square(yv) + np.square(zv)) < rad
        filtered_weights = gaussian_filter(
            weights, sigma=4 / resolution, output=np.double
        )

        matgrid = mp.MaterialGrid(
            mp.Vector3(N, N, N),
            SiO2,
            Si,
            weights=filtered_weights,
            do_averaging=True,
            beta=1000,
            eta=0.5,
        )

        geometry = [mp.Block(center=mp.Vector3(), size=cell_size, material=matgrid)]
    else:
        geometry = [mp.Sphere(center=mp.Vector3(), radius=rad, material=Si)]

    sources = [
        mp.Source(
            src=mp.GaussianSource(fcen, fwidth=df),
            size=mp.Vector3(),
            center=mp.Vector3(0.13, 0.25, 0.06),
            component=mp.Ez,
        )
    ]

    k_point = mp.Vector3(0.23, -0.17, 0.35)

    sim = mp.Simulation(
        resolution=resolution,
        cell_size=cell_size,
        sources=sources,
        default_material=SiO2,
        k_point=k_point,
        geometry=geometry,
    )

    h = mp.Harminv(mp.Ez, mp.Vector3(-0.2684, 0.1185, 0.0187), fcen, df)

    sim.run(mp.after_sources(h), until_after_sources=200)

    try:
        for m in h.modes:
            print(f"harminv:, {resolution}, {m.freq}, {m.Q}")
        freq = h.modes[0].freq
    except:
        raise RuntimeError("No resonant modes found.")

    return freq


class TestMaterialGrid(unittest.TestCase):
    def test_subpixel_smoothing(self):
        # "exact" frequency computed using MaterialGrid at resolution = 300
        freq_ref = 0.29826813873225283

        res = [25, 50]
        freq_matgrid = []
        for r in res:
            freq_matgrid.append(compute_resonant_mode_2d(r))
            # verify that the frequency of the resonant mode is
            # approximately equal to the reference value
            self.assertAlmostEqual(freq_ref, freq_matgrid[-1], 2)

        # verify that the relative error is decreasing with increasing resolution
        # and is better than linear convergence because of subpixel smoothing
        self.assertLess(
            abs(freq_matgrid[1] - freq_ref) * (res[1] / res[0]),
            abs(freq_matgrid[0] - freq_ref),
        )

        freq_matgrid_default_mat = compute_resonant_mode_2d(res[0], True)
        self.assertAlmostEqual(freq_matgrid[0], freq_matgrid_default_mat)

    def test_matgrid_3d(self):
        freq_matgrid = compute_resonant_mode_3d(True)
        freq_geomobj = compute_resonant_mode_3d(False)
        self.assertAlmostEqual(freq_matgrid, freq_geomobj, places=2)

    def test_symmetry(self):
        tran_nosym = compute_transmittance(False)
        tran_sym = compute_transmittance(True)
        self.assertAlmostEqual(tran_nosym, tran_sym, places=5)


if __name__ == "__main__":
    unittest.main()
