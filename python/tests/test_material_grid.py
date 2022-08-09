'''
Checks that a material grid works with symmetries, and that the
new subpixel smoothing feature is accurate.

TODO:
    Create a 3D test that works well with the new smoothing
'''

import meep as mp
from meep import mpb
try:
    import meep.adjoint as mpa
except:
    import adjoint as mpa

import unittest

import numpy as np
import unittest

rad = 0.301943
k_point = mp.Vector3(0.3892,0.1597,0)
Si = mp.Medium(index=3.5)
cell_size = mp.Vector3(1,1,0)

def compute_transmittance(matgrid_symmetry=False):
        resolution = 25

        cell_size = mp.Vector3(6,6,0)

        boundary_layers = [mp.PML(thickness=1.0)]

        matgrid_size = mp.Vector3(2,2,0)
        matgrid_resolution = 2*resolution

        Nx, Ny = int(matgrid_size.x*matgrid_resolution), int(matgrid_size.y*matgrid_resolution)

        # ensure reproducible results
        rng = np.random.RandomState(2069588)

        w = rng.rand(Nx,Ny)
        # for subpixel smoothing, the underlying design grid must be smoothly varying
        w = mpa.tanh_projection(rng.rand(Nx,Ny),1e5,0.5)
        w = mpa.conic_filter(w,0.25,matgrid_size.x,matgrid_size.y,matgrid_resolution)
        weights = 0.5 * (w + np.fliplr(w)) if not matgrid_symmetry else w

        matgrid = mp.MaterialGrid(mp.Vector3(Nx,Ny),
                                  mp.air,
                                  mp.Medium(index=3.5),
                                  weights=weights,
                                  beta=np.inf,
                                  grid_type='U_MEAN')

        geometry = [mp.Block(center=mp.Vector3(),
                             size=mp.Vector3(mp.inf,1.0,mp.inf),
                             material=mp.Medium(index=3.5)),
                    mp.Block(center=mp.Vector3(),
                             size=mp.Vector3(matgrid_size.x,matgrid_size.y,0),
                             material=matgrid)]

        if matgrid_symmetry:
                geometry.append(mp.Block(center=mp.Vector3(),
                                         size=mp.Vector3(matgrid_size.x,matgrid_size.y,0),
                                         material=matgrid,
                                         e2=mp.Vector3(y=-1)))

        eig_parity = mp.ODD_Y + mp.EVEN_Z

        fcen = 0.65
        df = 0.2*fcen
        sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen,fwidth=df),
                                      center=mp.Vector3(-2.0,0),
                                      size=mp.Vector3(0,4.0),
                                      eig_parity=eig_parity)]

        sim = mp.Simulation(resolution=resolution,
                            cell_size=cell_size,
                            boundary_layers=boundary_layers,
                            sources=sources,
                            geometry=geometry)

        mode_mon = sim.add_flux(fcen,
                                0,
                                1,
                                mp.FluxRegion(center=mp.Vector3(2.0,0),
                                              size=mp.Vector3(0,4.0)))

        sim.run(until_after_sources=mp.stop_when_dft_decayed())

        mode_coeff = sim.get_eigenmode_coefficients(mode_mon,[1],eig_parity).alpha[0,:,0][0]

        tran = np.power(np.abs(mode_coeff),2)
        print('tran:, {}, {}'.format("sym" if matgrid_symmetry else "nosym", tran))

        return tran

def compute_resonant_mode_2d(res,radius=rad,default_mat=False,cylinder=False,cylinder_matgrid=True,do_averaging=True):
        fcen = 0.3
        df = 0.2*fcen
        sources = [mp.Source(mp.GaussianSource(fcen,fwidth=df),
                             component=mp.Hz,
                             center=mp.Vector3(-0.1057,0.2094,0))]

        matgrid_size = mp.Vector3(1,1,0)
        matgrid_resolution = 1200

        # for a fixed resolution, compute the number of grid points
        # necessary which are defined on the corners of the voxels
        Nx, Ny = int(matgrid_size.x*matgrid_resolution), int(matgrid_size.y*matgrid_resolution)

        x = np.linspace(-0.5*matgrid_size.x,0.5*matgrid_size.x,Nx)
        y = np.linspace(-0.5*matgrid_size.y,0.5*matgrid_size.y,Ny)
        xv, yv = np.meshgrid(x,y)
        weights = mpa.make_sdf(np.sqrt(np.square(xv) + np.square(yv)) < radius)
        if cylinder:
            weights = 0.0*weights+1.0          
        
        matgrid = mp.MaterialGrid(mp.Vector3(Nx,Ny),
                                  mp.air,
                                  Si,
                                  weights=weights,
                                  beta=np.inf,
                                  eta=0.5)
        
        # Use a cylinder object, not a block object
        if cylinder:
            # within the cylinder, use a (uniform) material grid
            if cylinder_matgrid:
                geometry = [mp.Cylinder(center=mp.Vector3(),
                             radius=radius,
                             material=matgrid)]
            # within the cylinder, just use a normal medium
            else:
                geometry = [mp.Cylinder(center=mp.Vector3(),
                             radius=radius,
                             material=Si)]
        # use a block object
        else:
            geometry = [mp.Block(center=mp.Vector3(),
                             size=mp.Vector3(matgrid_size.x,matgrid_size.y,0),
                             material=matgrid)]

        sim = mp.Simulation(resolution=res,
                            eps_averaging=do_averaging,
                            cell_size=cell_size,
                            default_material=matgrid if default_mat else mp.Medium(),
                            geometry=geometry if not default_mat else [],
                            sources=sources,
                            k_point=k_point)

        h = mp.Harminv(mp.Hz, mp.Vector3(0.3718,-0.2076), fcen, df)
        sim.run(mp.after_sources(h),
                until_after_sources=300)

        try:
            for m in h.modes:
                print("harminv:, {}, {}, {}".format(res,m.freq,m.Q))
            freq = h.modes[0].freq
        except:
            raise RuntimeError("No resonant modes found.")

        return freq

def compute_resonant_mode_mpb(resolution=512):
    geometry = [mp.Cylinder(rad, material=Si)]
    geometry_lattice = mp.Lattice(cell_size)

    ms = mpb.ModeSolver(num_bands=1,
                        k_points=[mp.cartesian_to_reciprocal(k_point, geometry_lattice)],
                        geometry=geometry,
                        geometry_lattice=geometry_lattice,
                        resolution=resolution)


    ms.run_te()
    return ms.freqs[0]

class TestMaterialGrid(unittest.TestCase):
    def test_subpixel_smoothing(self):
        res = 25

        def subpixel_test_matrix(radius,do_averaging):
            freq_cylinder = compute_resonant_mode_2d(res, radius, default_mat=False, cylinder=True, cylinder_matgrid=False, do_averaging=do_averaging)
            freq_matgrid = compute_resonant_mode_2d(res, radius, default_mat=False, cylinder=False, do_averaging=do_averaging)
            freq_matgrid_cylinder = compute_resonant_mode_2d(res, radius, default_mat=False, cylinder=True, cylinder_matgrid=True, do_averaging=do_averaging)
            return [freq_cylinder,freq_matgrid,freq_matgrid_cylinder]

        # when smoothing is off, all three tests should be identical to machine precision
        no_smoothing = subpixel_test_matrix(rad,False)
        self.assertEqual(no_smoothing[0], no_smoothing[1])
        self.assertEqual(no_smoothing[1], no_smoothing[2])

        # when we slightly perturb the radius, the results should be the same as before.
        no_smoothing_perturbed = subpixel_test_matrix(rad+0.01/res,False)
        self.assertEqual(no_smoothing, no_smoothing_perturbed)

        # when smoothing is on, the simple material results should be different from the matgrid results
        smoothing = subpixel_test_matrix(rad,True)
        self.assertNotEqual(smoothing[0], smoothing[1])
        self.assertAlmostEqual(smoothing[1], smoothing[2], 3)

        # when we slighty perturb the radius, the results should all be different from before.
        smoothing_perturbed = subpixel_test_matrix(rad+0.01/res,True)
        self.assertNotEqual(smoothing, smoothing_perturbed)

        # "exact" frequency computed using MPB
        freq_ref = compute_resonant_mode_mpb()
        print(freq_ref)
        res = [50, 100]
        freq_matgrid = []
        for r in res:
            freq_matgrid.append(compute_resonant_mode_2d(r,cylinder=False))
            # verify that the frequency of the resonant mode is
            # approximately equal to the reference value
            self.assertAlmostEqual(freq_ref, freq_matgrid[-1], 2)
        print("results:    ",freq_ref,freq_matgrid)

        # verify that the relative error is decreasing with increasing resolution
        # and is better than linear convergence because of subpixel smoothing
        self.assertLess(
            abs(freq_matgrid[1] - freq_ref) * (res[1] / res[0]),
            abs(freq_matgrid[0] - freq_ref),
        )

        # ensure that a material grid as a default material works
        freq_matgrid_default_mat = compute_resonant_mode_2d(res[0], rad, True)
        self.assertAlmostEqual(freq_matgrid[0], freq_matgrid_default_mat)

    def test_symmetry(self):
        tran_nosym = compute_transmittance(False)
        tran_sym = compute_transmittance(True)
        self.assertAlmostEqual(tran_nosym, tran_sym, places=5)


if __name__ == "__main__":
    unittest.main()
