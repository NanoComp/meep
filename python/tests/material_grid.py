import meep as mp
import numpy as np
from scipy.ndimage import gaussian_filter
import unittest

def compute_resonant_mode(res,default_mat=False):
        cell_size = mp.Vector3(1,1,0)

        rad = 0.301943

        fcen = 0.3
        df = 0.2*fcen
        sources = [mp.Source(mp.GaussianSource(fcen,fwidth=df),
                             component=mp.Hz,
                             center=mp.Vector3(-0.1057,0.2094,0))]

        k_point = mp.Vector3(0.3892,0.1597,0)

        design_shape = mp.Vector3(1,1,0)
        design_region_resolution = 50
        Nx = int(design_region_resolution*design_shape.x)
        Ny = int(design_region_resolution*design_shape.y)
        x = np.linspace(-0.5*design_shape.x,0.5*design_shape.x,Nx)
        y = np.linspace(-0.5*design_shape.y,0.5*design_shape.y,Ny)
        xv, yv = np.meshgrid(x,y)
        design_params = np.sqrt(np.square(xv) + np.square(yv)) < rad
        filtered_design_params = gaussian_filter(design_params,
                                                 sigma=3.0,
                                                 output=np.double)

        matgrid = mp.MaterialGrid(mp.Vector3(Nx,Ny),
                                  mp.air,
                                  mp.Medium(index=3.5),
                                  weights=filtered_design_params,
                                  do_averaging=True,
                                  beta=1000,
                                  eta=0.5)

        geometry = [mp.Block(center=mp.Vector3(),
                             size=mp.Vector3(design_shape.x,design_shape.y,0),
                             material=matgrid)]

        sim = mp.Simulation(resolution=res,
                            cell_size=cell_size,
                            default_material=matgrid if default_mat else mp.Medium(),
                            geometry=geometry if not default_mat else [],
                            sources=sources,
                            k_point=k_point)

        h = mp.Harminv(mp.Hz, mp.Vector3(0.3718,-0.2076), fcen, df)
        sim.run(mp.after_sources(h),
                until_after_sources=200)

        try:
            for m in h.modes:
                print("harminv:, {}, {}, {}".format(res,m.freq,m.Q))
            freq = h.modes[0].freq
        except:
            print("No resonant modes found.")

        sim.reset_meep()
        return freq

class TestMaterialGrid(unittest.TestCase):

    def test_material_grid(self):
        ## reference frequency computed using MaterialGrid at resolution = 300
        freq_ref = 0.306877757638932

        res = [25, 50]
        freq_matgrid = []
        for r in res:
            freq_matgrid.append(compute_resonant_mode(r))
            ## verify that the resonant mode is approximately equivalent to
            ## the reference value
            self.assertAlmostEqual(freq_ref, freq_matgrid[-1], 2)

        ## verify that the relative error is decreasing with increasing resolution
        ## and is better than linear convergence because of subpixel smoothing
        self.assertLess(abs(freq_matgrid[1]-freq_ref)*(res[1]/res[0]),
                        abs(freq_matgrid[0]-freq_ref))

        freq_matgrid_default_mat = compute_resonant_mode(res[0], True)
        self.assertAlmostEqual(freq_matgrid[0], freq_matgrid_default_mat)

if __name__ == '__main__':
    unittest.main()
