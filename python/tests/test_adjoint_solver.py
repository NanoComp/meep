import meep as mp

try:
    import meep.adjoint as mpa
except:
    import adjoint as mpa

import unittest
from enum import Enum
from typing import List, Union, Tuple
import numpy as np
from autograd import numpy as npa
from autograd import tensor_jacobian_product
from utils import ApproxComparisonTestCase

MonitorObject = Enum("MonitorObject", "EIGENMODE DFT LDOS")


class TestAdjointSolver(ApproxComparisonTestCase):
    @classmethod
    def setUpClass(cls):
        cls.resolution = 30  # pixels/μm

        cls.silicon = mp.Medium(epsilon=12)
        cls.sapphire = mp.Medium(
            epsilon_diag=(10.225, 10.225, 9.95),
            epsilon_offdiag=(-0.825, -0.55 * np.sqrt(3 / 2), 0.55 * np.sqrt(3 / 2)),
        )

        cls.sxy = 5.0
        cls.cell_size = mp.Vector3(cls.sxy, cls.sxy, 0)

        cls.dpml = 1.0
        cls.pml_xy = [mp.PML(thickness=cls.dpml)]
        cls.pml_x = [mp.PML(thickness=cls.dpml, direction=mp.X)]

        cls.eig_parity = mp.EVEN_Y + mp.ODD_Z
        cls.src_cmpt = mp.Ez

        cls.design_region_size = mp.Vector3(1.5, 1.5)
        cls.design_region_resolution = int(2 * cls.resolution)
        cls.Nx = int(round(cls.design_region_size.x * cls.design_region_resolution))
        cls.Ny = int(round(cls.design_region_size.y * cls.design_region_resolution))

        # ensure reproducible results
        rng = np.random.RandomState(9861548)

        # random design region
        cls.p = 0.5 * rng.rand(cls.Nx * cls.Ny)

        # random perturbation for design region
        deps = 1e-5
        cls.dp = deps * rng.rand(cls.Nx * cls.Ny)

        cls.w = 1.0
        cls.waveguide_geometry = [
            mp.Block(
                material=cls.silicon,
                center=mp.Vector3(),
                size=mp.Vector3(mp.inf, cls.w, mp.inf),
            )
        ]

        # source center frequency and bandwidth
        cls.fcen = 1 / 1.55
        cls.df = 0.05 * cls.fcen

        # monitor frequencies
        # two cases: (1) single and (2) multi frequency
        cls.mon_frqs = [
            [cls.fcen],
            [
                cls.fcen - 0.09 * cls.df,
                cls.fcen,
                cls.fcen + 0.06 * cls.df,
            ],
        ]

        cls.mode_source = [
            mp.EigenModeSource(
                src=mp.GaussianSource(cls.fcen, fwidth=cls.df),
                center=mp.Vector3(-0.5 * cls.sxy + cls.dpml, 0),
                size=mp.Vector3(0, cls.sxy - 2 * cls.dpml),
                eig_parity=cls.eig_parity,
            )
        ]

        cls.pt_source = [
            mp.Source(
                src=mp.GaussianSource(cls.fcen, fwidth=cls.df),
                center=mp.Vector3(-0.5 * cls.sxy + cls.dpml, 0),
                size=mp.Vector3(),
                component=cls.src_cmpt,
            )
        ]

        cls.line_source = [
            mp.Source(
                src=mp.GaussianSource(cls.fcen, fwidth=cls.df),
                center=mp.Vector3(-0.85, 0),
                size=mp.Vector3(0, cls.sxy - 2 * cls.dpml),
                component=cls.src_cmpt,
            )
        ]

        cls.k_point = mp.Vector3(0.23, -0.38)

        # location of DFT monitors for reflected and transmitted fields
        cls.refl_pt = mp.Vector3(-0.5 * cls.sxy + cls.dpml + 0.5, 0)
        cls.tran_pt = mp.Vector3(0.5 * cls.sxy - cls.dpml, 0)

    def adjoint_solver(
        self,
        design_params: List[float] = None,
        mon_type: MonitorObject = None,
        frequencies: List[float] = None,
        mat2: mp.Medium = None,
        need_gradient: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray]:
        matgrid = mp.MaterialGrid(
            mp.Vector3(self.Nx, self.Ny),
            mp.air,
            self.silicon if mat2 is None else mat2,
            weights=np.ones((self.Nx, self.Ny)),
        )

        matgrid_region = mpa.DesignRegion(
            matgrid,
            volume=mp.Volume(
                center=mp.Vector3(),
                size=mp.Vector3(
                    self.design_region_size.x, self.design_region_size.y, 0
                ),
            ),
        )

        matgrid_geometry = [
            mp.Block(
                center=matgrid_region.center, size=matgrid_region.size, material=matgrid
            )
        ]

        geometry = self.waveguide_geometry + matgrid_geometry

        sim = mp.Simulation(
            resolution=self.resolution,
            cell_size=self.cell_size,
            boundary_layers=self.pml_xy,
            sources=self.mode_source,
            geometry=geometry,
        )

        if not frequencies:
            frequencies = [self.fcen]

        if mon_type.name == "EIGENMODE":
            if len(frequencies) == 1:
                if mat2 is None:
                    # the incident fields of the mode source in the
                    # straight waveguide are used as normalization
                    # of the reflectance (S11) measurement.
                    ref_sim = mp.Simulation(
                        resolution=self.resolution,
                        cell_size=self.cell_size,
                        boundary_layers=self.pml_xy,
                        sources=self.mode_source,
                        geometry=self.waveguide_geometry,
                    )
                    dft_mon = ref_sim.add_mode_monitor(
                        frequencies,
                        mp.ModeRegion(
                            center=self.refl_pt,
                            size=mp.Vector3(0, self.sxy - 2 * self.dpml, 0),
                        ),
                        yee_grid=True,
                    )
                    ref_sim.run(until_after_sources=20)
                    subtracted_dft_fields = ref_sim.get_flux_data(dft_mon)
                else:
                    subtracted_dft_fields = None

                obj_list = [
                    mpa.EigenmodeCoefficient(
                        sim,
                        mp.Volume(
                            center=self.refl_pt,
                            size=mp.Vector3(0, self.sxy - 2 * self.dpml, 0),
                        ),
                        1,
                        forward=False,
                        eig_parity=self.eig_parity,
                        subtracted_dft_fields=subtracted_dft_fields,
                    ),
                    mpa.EigenmodeCoefficient(
                        sim,
                        mp.Volume(
                            center=self.tran_pt,
                            size=mp.Vector3(0, self.sxy - 2 * self.dpml, 0),
                        ),
                        2,
                        eig_parity=self.eig_parity,
                    ),
                ]

                def J(refl_mon, tran_mon):
                    return -npa.power(npa.abs(refl_mon), 2) + npa.power(
                        npa.abs(tran_mon), 2
                    )

            else:
                obj_list = [
                    mpa.EigenmodeCoefficient(
                        sim,
                        mp.Volume(
                            center=mp.Vector3(0.5 * self.sxy - self.dpml),
                            size=mp.Vector3(0, self.sxy - 2 * self.dpml, 0),
                        ),
                        1,
                        eig_parity=self.eig_parity,
                    )
                ]

                def J(tran_mon):
                    return npa.power(npa.abs(tran_mon), 2)

        elif mon_type.name == "DFT":
            obj_list = [
                mpa.FourierFields(
                    sim,
                    mp.Volume(center=mp.Vector3(1.25), size=mp.Vector3(0.25, 1, 0)),
                    self.src_cmpt,
                )
            ]

            def J(mode_mon):
                return npa.power(npa.abs(mode_mon[:, 4, 10]), 2)

        elif mon_type.name == "LDOS":
            sim.change_sources(self.line_source)
            obj_list = [mpa.LDOS(sim)]

            def J(mode_mon):
                return npa.array(mode_mon)

        opt = mpa.OptimizationProblem(
            simulation=sim,
            objective_functions=J,
            objective_arguments=obj_list,
            design_regions=[matgrid_region],
            frequencies=frequencies,
        )

        f, dJ_du = opt([design_params], need_gradient=need_gradient)

        return f, dJ_du

    def adjoint_solver_complex_fields(
        self, design_params, frequencies=None, need_gradient=True
    ) -> Tuple[np.ndarray, np.ndarray]:
        matgrid = mp.MaterialGrid(
            mp.Vector3(self.Nx, self.Ny),
            mp.air,
            self.silicon,
            weights=np.ones((self.Nx, self.Ny)),
        )

        matgrid_region = mpa.DesignRegion(
            matgrid,
            volume=mp.Volume(
                center=mp.Vector3(),
                size=mp.Vector3(
                    self.design_region_size.x, self.design_region_size.y, 0
                ),
            ),
        )

        geometry = [
            mp.Block(
                center=matgrid_region.center, size=matgrid_region.size, material=matgrid
            )
        ]

        sim = mp.Simulation(
            resolution=self.resolution,
            cell_size=self.cell_size,
            default_material=self.silicon,
            k_point=self.k_point,
            boundary_layers=self.pml_x,
            sources=self.pt_source,
            geometry=geometry,
        )

        if not frequencies:
            frequencies = [self.fcen]

        obj_list = [
            mpa.FourierFields(
                sim, mp.Volume(center=mp.Vector3(0.9), size=mp.Vector3(0.2, 0.5)), mp.Dz
            )
        ]

        def J(dft_mon):
            return npa.power(npa.abs(dft_mon[:, 3, 9]), 2)

        opt = mpa.OptimizationProblem(
            simulation=sim,
            objective_functions=J,
            objective_arguments=obj_list,
            design_regions=[matgrid_region],
            frequencies=frequencies,
        )

        f, dJ_du = opt([design_params], need_gradient=need_gradient)

        return f, dJ_du


    def test_periodic_design(self):
        """Verifies that the constaint functions are invariant when
        the design pattern is shifted along periodic directions."""
        print("*** TESTING PERIODIC DESIGN ***")

        # shifted design patterns
        idx_row_shift, idx_col_shift = int(0.19*self.Nx), int(0.28*self.Ny)
        p = self.p.reshape(self.Nx, self.Ny)
        p_row_shift = np.vstack((p[idx_row_shift:,:], p[0:idx_row_shift,:]))
        p_col_shift = np.hstack((p[:,idx_col_shift:], p[:,0:idx_col_shift]))

        eta_e  = 0.55
        eta_d = 1-eta_e
        beta, eta = 10, 0.5
        radius, c = 0.3, 400
        places = 18
        threshold_f = lambda x: mpa.tanh_projection(x, beta, eta)

        for selected_filter in (mpa.conic_filter, mpa.cylindrical_filter, mpa.gaussian_filter):
            for periodic_axes in (0, 1, (0,1)):
                filter_f = lambda x: selected_filter(x, radius,
                                                     self.design_region_size.x,
                                                     self.design_region_size.y,
                                                     self.design_region_resolution,
                                                     periodic_axes)

                constraint_solid_original = mpa.constraint_solid(self.p, c, eta_e, filter_f, threshold_f,
                                                                 self.design_region_resolution, periodic_axes)
                
                constraint_void_original = mpa.constraint_void(self.p, c, eta_d, filter_f, threshold_f,
                                                               self.design_region_resolution, periodic_axes)

                if periodic_axes in (0, (0,1)):
                    constraint_solid_row_shift = mpa.constraint_solid(p_row_shift, c, eta_e, filter_f, threshold_f,
                                                                      self.design_region_resolution, periodic_axes)
                    self.assertAlmostEqual(constraint_solid_original, constraint_solid_row_shift, places=places)

                    constraint_void_row_shift = mpa.constraint_void(p_row_shift, c, eta_d, filter_f, threshold_f,
                                                                    self.design_region_resolution, periodic_axes)
                    self.assertAlmostEqual(constraint_void_original, constraint_void_row_shift, places=places)

                if periodic_axes in (1, (0,1)):
                    constraint_solid_col_shift = mpa.constraint_solid(p_col_shift, c, eta_e, filter_f, threshold_f,
                                                                      self.design_region_resolution, periodic_axes)
                    self.assertAlmostEqual(constraint_solid_original, constraint_solid_col_shift, places=places)

                    constraint_void_col_shift = mpa.constraint_void(p_col_shift, c, eta_d, filter_f, threshold_f,
                                                                    self.design_region_resolution, periodic_axes)
                    self.assertAlmostEqual(constraint_void_original, constraint_void_col_shift, places=places)

            print(f"PASSED: filter function = {selected_filter.__name__}")


if __name__ == "__main__":
    unittest.main()