import unittest
import meep as mp


class TestDFTFields(unittest.TestCase):

    def init(self, pulse):
        resolution = 10
        n = 3.4
        w = 1.0
        r = 1.0
        pad = 4
        dpml = 2
        sxy = 2.0 * (r + w + pad + dpml)
        cell = mp.Vector3(sxy, sxy)
        symmetries = [mp.Identity(mp.Y)]
        pml_layers = [mp.PML(dpml)]

        geometry = [
            mp.Cylinder(r + w, material=mp.Medium(epsilon=n * n)),
            mp.Cylinder(r, material=mp.vacuum),
        ]
        self.fcen = 0.118
        self.df = 0.1
        self.x0 = mp.Vector3(r + 0.1)

        if pulse:
            src = mp.GaussianSource(self.fcen, df=self.df)
        else:
            src = mp.ContinuousSource(self.fcen, df=self.df)

        sources = [mp.Source(src=src, component=mp.Ez, center=self.x0)]

        return mp.Simulation(
            cell_size=cell,
            resolution=resolution,
            geometry=geometry,
            sources=sources,
            boundary_layers=pml_layers,
            symmetries=symmetries
        )

    def test_get_dft_array(self):
        sim = self.init(True)
        components = [mp.Ex, mp.Ey, mp.Ez, mp.Hx, mp.Hy, mp.Hz]
        dft_fields = sim.add_dft_fields(components, 6, self.fcen, self.fcen, 1)
        dft_flux = sim.add_dft_flux([mp.X], self.fcen, self.fcen, 1)

        sim.run(until_after_sources=100)
        sim.output_dft(dft_flux, "dft-flux")
        sim.output_dft(dft_fields, "dft-fields")

        field_array = sim.get_dft_array(dft_flux, mp.Ez, 0)

    def test_output_dft(self):
        sim = self.init(False)
        sim.solve_cw(1e-8, 10000, 10)
        # f.output_hdf5(Ez, f.v, file)
        # f.output_hdf5(Hx, f.v, file)
        # f.output_hdf5(Hy, f.v, file)


if __name__ == '__main__':
    unittest.main()
