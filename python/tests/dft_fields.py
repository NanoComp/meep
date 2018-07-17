import unittest
import h5py
import numpy as np
import meep as mp


class TestDFTFields(unittest.TestCase):

    def init(self):
        resolution = 10
        n = 3.4
        w = 1.0
        r = 1.0
        pad = 4
        self.dpml = 2
        self.sxy = 2.0 * (r + w + pad + self.dpml)
        cell = mp.Vector3(self.sxy, self.sxy)
        pml_layers = [mp.PML(self.dpml)]

        geometry = [
            mp.Cylinder(r + w, material=mp.Medium(epsilon=n * n)),
            mp.Cylinder(r, material=mp.vacuum),
        ]
        self.fcen = 0.118
        self.df = 0.1

        src = mp.GaussianSource(self.fcen, fwidth=self.df)
        sources = [mp.Source(src=src, component=mp.Ez, center=mp.Vector3(r + 0.1))]

        return mp.Simulation(
            cell_size=cell,
            resolution=resolution,
            geometry=geometry,
            sources=sources,
            boundary_layers=pml_layers,
        )

    def test_get_dft_array(self):
        sim = self.init()
        sim.init_sim()
        dft_fields = sim.add_dft_fields([mp.Ez], self.fcen, self.fcen, 1)
        fr = mp.FluxRegion(mp.Vector3(), size=mp.Vector3(self.sxy, self.sxy), direction=mp.X)
        dft_flux = sim.add_flux(self.fcen, 0, 1, fr)

        # volumes with zero thickness in x and y directions to test collapsing
        # of empty dimensions in DFT array and HDF5 output routines
        thin_x_volume = mp.volume( mp.vec(0.35*self.sxy, -0.4*self.sxy),
                                   mp.vec(0.35*self.sxy, +0.4*self.sxy));
        thin_x_flux = sim.fields.add_dft_fields([mp.Ez], thin_x_volume, self.fcen, self.fcen, 1)
        thin_y_volume = mp.volume( mp.vec(-0.5*self.sxy, 0.25*self.sxy),
                                   mp.vec(+0.5*self.sxy, 0.25*self.sxy));
        thin_y_flux = sim.fields.add_dft_flux(mp.Y, thin_y_volume, self.fcen, self.fcen, 1)

        sim.run(until_after_sources=100)

        # test proper collapsing of degenerate dimensions in HDF5 files and arrays
        thin_x_array = sim.get_dft_array(thin_x_flux, mp.Ez, 0)
        thin_y_array = sim.get_dft_array(thin_y_flux, mp.Ez, 0)
        np.testing.assert_equal(thin_x_array.ndim, 1)
        np.testing.assert_equal(thin_y_array.ndim, 1)

        sim.output_dft(thin_x_flux, 'thin-x-flux')
        sim.output_dft(thin_y_flux, 'thin-y-flux')

        with h5py.File('thin-x-flux.h5', 'r') as thin_x:
            thin_x_h5 = thin_x['ez_0.r'].value + 1j * thin_x['ez_0.i'].value

        with h5py.File('thin-y-flux.h5', 'r') as thin_y:
            thin_y_h5 = thin_y['ez_0.r'].value + 1j * thin_y['ez_0.i'].value

        np.testing.assert_allclose(thin_x_array, thin_x_h5)
        np.testing.assert_allclose(thin_y_array, thin_y_h5)

        # compare array data to HDF5 file content for fields and flux
        fields_arr = sim.get_dft_array(dft_fields, mp.Ez, 0)
        flux_arr = sim.get_dft_array(dft_flux, mp.Ez, 0)

        sim.output_dft(dft_fields, 'dft-fields')
        sim.output_dft(dft_flux, 'dft-flux')

        with h5py.File('dft-fields.h5', 'r') as fields, h5py.File('dft-flux.h5', 'r') as flux:
            exp_fields = fields['ez_0.r'].value + 1j * fields['ez_0.i'].value
            exp_flux = flux['ez_0.r'].value + 1j * flux['ez_0.i'].value

        np.testing.assert_allclose(exp_fields, fields_arr)
        np.testing.assert_allclose(exp_flux, flux_arr)

if __name__ == '__main__':
    unittest.main()
