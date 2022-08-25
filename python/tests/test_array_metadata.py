import unittest

import numpy as np
from utils import ApproxComparisonTestCase

import meep as mp


class TestArrayMetadata(ApproxComparisonTestCase):
    def test_array_metadata(self):
        resolution = 25

        n = 3.4
        w = 1
        r = 1
        pad = 4
        dpml = 2

        sxy = 2 * (r + w + pad + dpml)
        cell_size = mp.Vector3(sxy, sxy)

        nonpml_vol = mp.Volume(
            mp.Vector3(), size=mp.Vector3(sxy - 2 * dpml, sxy - 2 * dpml)
        )

        geometry = [
            mp.Cylinder(radius=r + w, material=mp.Medium(index=n)),
            mp.Cylinder(radius=r),
        ]

        fcen = 0.118
        df = 0.08

        symmetries = [mp.Mirror(mp.X, phase=-1), mp.Mirror(mp.Y, phase=+1)]

        pml_layers = [mp.PML(dpml)]

        # CW source
        src = [
            mp.Source(mp.ContinuousSource(fcen, fwidth=df), mp.Ez, mp.Vector3(r + 0.1)),
            mp.Source(
                mp.ContinuousSource(fcen, fwidth=df),
                mp.Ez,
                mp.Vector3(-(r + 0.1)),
                amplitude=-1,
            ),
        ]

        sim = mp.Simulation(
            cell_size=cell_size,
            geometry=geometry,
            sources=src,
            resolution=resolution,
            force_complex_fields=True,
            symmetries=symmetries,
            boundary_layers=pml_layers,
        )

        sim.init_sim()
        sim.solve_cw(1e-5 if mp.is_single_precision() else 1e-6, 1000, 10)

        def electric_energy(r, ez, eps):
            return np.real(eps * np.conj(ez) * ez)

        def vec_func(r):
            return r.x**2 + 2 * r.y**2

        electric_energy_total = sim.integrate_field_function(
            [mp.Ez, mp.Dielectric], electric_energy, nonpml_vol
        )
        electric_energy_max = sim.max_abs_field_function(
            [mp.Ez, mp.Dielectric], electric_energy, nonpml_vol
        )
        vec_func_total = sim.integrate_field_function([], vec_func, nonpml_vol)
        cw_modal_volume = (electric_energy_total / electric_energy_max) * vec_func_total

        sim.reset_meep()

        # pulsed source
        src = [
            mp.Source(mp.GaussianSource(fcen, fwidth=df), mp.Ez, mp.Vector3(r + 0.1)),
            mp.Source(
                mp.GaussianSource(fcen, fwidth=df),
                mp.Ez,
                mp.Vector3(-(r + 0.1)),
                amplitude=-1,
            ),
        ]

        sim = mp.Simulation(
            cell_size=cell_size,
            geometry=geometry,
            k_point=mp.Vector3(),
            sources=src,
            resolution=resolution,
            symmetries=symmetries,
            boundary_layers=pml_layers,
        )

        dft_obj = sim.add_dft_fields([mp.Ez], fcen, 0, 1, where=nonpml_vol)
        sim.run(until_after_sources=100)

        Ez = sim.get_dft_array(dft_obj, mp.Ez, 0)
        (X, Y, Z, W) = sim.get_array_metadata(dft_cell=dft_obj)
        Eps = sim.get_array(vol=nonpml_vol, component=mp.Dielectric)
        EpsE2 = np.real(Eps * np.conj(Ez) * Ez)
        xm, ym = np.meshgrid(X, Y)
        vec_func_sum = np.sum(W * (xm**2 + 2 * ym**2))
        pulse_modal_volume = np.sum(W * EpsE2) / np.max(EpsE2) * vec_func_sum

        tol = 5e-2 if mp.is_single_precision() else 1e-2
        self.assertClose(cw_modal_volume / pulse_modal_volume, 1.0, epsilon=tol)


if __name__ == "__main__":
    unittest.main()
