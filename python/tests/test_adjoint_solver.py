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

    def adjoint_solver_damping(
        self,
        design_params: List[float] = None,
        frequencies: List[float] = None,
        mat2: mp.Medium = None,
        need_gradient: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray]:
        matgrid = mp.MaterialGrid(
            mp.Vector3(self.Nx, self.Ny),
            mp.air,
            self.silicon if mat2 is None else mat2,
            weights=np.ones((self.Nx, self.Ny)),
            damping=np.pi * self.fcen,
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

        def J(mode_mon):
            return npa.power(npa.abs(mode_mon), 2)

        opt = mpa.OptimizationProblem(
            simulation=sim,
            objective_functions=J,
            objective_arguments=obj_list,
            design_regions=[matgrid_region],
            frequencies=frequencies,
            minimum_run_time=150,
        )

        f, dJ_du = opt([design_params], need_gradient=need_gradient)

        return f, dJ_du

    def adjoint_solver_two_objfunc(
        self,
        design_params: List[float],
        frequencies: List[float] = None,
        need_gradient: bool = True,
    ) -> Tuple[List[np.ndarray], np.ndarray]:
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
        input_flux = np.array(mp.get_fluxes(dft_mon))

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

        matgrid_geometry = [
            mp.Block(
                center=matgrid_region.center,
                size=matgrid_region.size,
                material=matgrid,
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

        obj_list = [
            mpa.EigenmodeCoefficient(
                sim,
                mp.Volume(
                    center=self.refl_pt,
                    size=mp.Vector3(0, self.sxy - 2 * self.dpml, 0),
                ),
                1,
                forward=False,
                subtracted_dft_fields=subtracted_dft_fields,
                eig_parity=self.eig_parity,
            ),
            mpa.EigenmodeCoefficient(
                sim,
                mp.Volume(
                    center=self.tran_pt,
                    size=mp.Vector3(0, self.sxy - 2 * self.dpml, 0),
                ),
                2,
                eig_parity=mp.ODD_Z,
            ),
        ]

        def J1(refl_mon, tran_mon):
            """Reflectance into first-order mode of Port 1."""
            return npa.power(npa.abs(refl_mon), 2) / input_flux

        def J2(refl_mon, tran_mon):
            """1-transmittance into second-order mode of Port 2."""
            return 1 - (npa.power(npa.abs(tran_mon), 2) / input_flux)

        opt = mpa.OptimizationProblem(
            simulation=sim,
            objective_functions=[J1, J2],
            objective_arguments=obj_list,
            design_regions=[matgrid_region],
            frequencies=frequencies,
        )

        f0, dJ_du = opt([design_params], need_gradient=need_gradient)
        f0_reflection = f0[0]
        f0_transmission = f0[1]
        f0_merged = np.concatenate((f0_reflection, f0_transmission))

        if need_gradient:
            dJ_du_reflection = dJ_du[0]
            dJ_du_transmission = dJ_du[1]
            nf = len(frequencies)
            grad = np.zeros((self.Nx * self.Ny, 2 * nf))
            if dJ_du_reflection.ndim < 2:
                dJ_du_reflection = np.expand_dims(dJ_du_reflection, axis=1)
                dJ_du_transmission = np.expand_dims(dJ_du_transmission, axis=1)
            grad[:, :nf] = dJ_du_reflection
            grad[:, nf:] = dJ_du_transmission
            return f0_merged, grad
        else:
            return f0_merged, np.array([])

    def mapping(
        self,
        x: List[float] = None,
        filter_radius: float = None,
        eta: float = None,
        beta: float = None,
    ):
        filtered_field = mpa.conic_filter(
            x,
            filter_radius,
            self.design_region_size.x,
            self.design_region_size.y,
            self.design_region_resolution,
        )

        projected_field = mpa.tanh_projection(filtered_field, beta, eta)

        return projected_field.flatten()

    def test_DFT_fields(self):
        """Verifies that the adjoint gradient for an objective function based
        on the DFT fields agrees with the finite-difference approximation."""
        print("*** TESTING DFT OBJECTIVE ***")

        for frequencies in self.mon_frqs:
            # compute objective value and its gradient for unperturbed design
            unperturbed_val, unperturbed_grad = self.adjoint_solver(
                self.p, MonitorObject.DFT, frequencies
            )

            # compute objective value for perturbed design
            perturbed_val, _ = self.adjoint_solver(
                self.p + self.dp, MonitorObject.DFT, frequencies, need_gradient=False
            )

            # compare directional derivative
            if unperturbed_grad.ndim < 2:
                unperturbed_grad = np.expand_dims(unperturbed_grad, axis=1)
            adj_dd = (self.dp[None, :] @ unperturbed_grad).flatten()
            fnd_dd = perturbed_val - unperturbed_val
            tol = 0.03 if mp.is_single_precision() else 0.002
            self.assertClose(adj_dd, fnd_dd, epsilon=tol)
            print(
                f"PASSED: frequencies={frequencies}, "
                f"adjoint gradient={adj_dd}, finite difference={fnd_dd}"
            )

    def test_eigenmode(self):
        """Verifies that the adjoint gradient for an objective function based
        on eigenmode decomposition agrees with the finite-difference
        approximation."""
        print("*** TESTING EIGENMODE OBJECTIVE ***")

        for frequencies in self.mon_frqs:
            # compute objective value and its gradient for unperturbed design
            unperturbed_val, unperturbed_grad = self.adjoint_solver(
                self.p, MonitorObject.EIGENMODE, frequencies
            )

            # compute objective for perturbed design
            perturbed_val, _ = self.adjoint_solver(
                self.p + self.dp,
                MonitorObject.EIGENMODE,
                frequencies,
                need_gradient=False,
            )

            # compare directional derivative
            if unperturbed_grad.ndim < 2:
                unperturbed_grad = np.expand_dims(unperturbed_grad, axis=1)
            adj_dd = (self.dp[None, :] @ unperturbed_grad).flatten()
            fnd_dd = perturbed_val - unperturbed_val

            # since the eigenmode source uses the mode profile
            # of the center frequency, the results for the
            # non-center frequencies of a multifrequency simulation
            # are expected to be *less* accurate than the center frequency
            if len(frequencies) == 1 and frequencies[0] == self.fcen:
                tol = 0.002 if mp.is_single_precision() else 5e-5
            else:
                tol = 0.008 if mp.is_single_precision() else 0.002

            self.assertClose(adj_dd, fnd_dd, epsilon=tol)
            print(
                f"PASSED: frequencies={frequencies}, "
                f"adjoint gradient={adj_dd}, finite difference={fnd_dd}"
            )

    def test_ldos(self):
        """Verifies that the adjoint gradient for an objective function based
        on the local density of states (LDoS) agrees with the finite-difference
        approximation."""
        print("*** TESTING LDOS OBJECTIVE ***")

        for frequencies in self.mon_frqs:
            # compute objective value and its gradient for unperturbed design
            unperturbed_val, unperturbed_grad = self.adjoint_solver(
                self.p, MonitorObject.LDOS, frequencies
            )

            # compute objective for perturbed design
            perturbed_val, _ = self.adjoint_solver(
                self.p + self.dp, MonitorObject.LDOS, frequencies, need_gradient=False
            )

            # compare directional derivative
            if unperturbed_grad.ndim < 2:
                unperturbed_grad = np.expand_dims(unperturbed_grad, axis=1)
            adj_dd = (self.dp[None, :] @ unperturbed_grad).flatten()
            fnd_dd = perturbed_val - unperturbed_val
            tol = 0.002 if mp.is_single_precision() else 0.001
            self.assertClose(adj_dd, fnd_dd, epsilon=tol)
            print(
                f"PASSED: frequencies={frequencies}, "
                f"adjoint gradient={adj_dd}, finite difference={fnd_dd}"
            )

    def test_gradient_backpropagation(self):
        """Verifies that the adjoint gradient can be back propagated through
        a differentiable mapping function applied to the design region
        and agrees with the finite-difference approximation."""
        print("*** TESTING GRADIENT BACKPROPAGATION ***")

        # filter/thresholding parameters
        filter_radius = 0.21985
        eta = 0.49093
        beta = 4.0698

        for frequencies in self.mon_frqs:
            nfrq = len(frequencies)

            mapped_p = self.mapping(self.p, filter_radius, eta, beta)

            # compute objective value and its gradient for unperturbed design
            unperturbed_val, unperturbed_grad = self.adjoint_solver(
                mapped_p, MonitorObject.EIGENMODE, frequencies
            )

            # backpropagate the gradient using vector-Jacobian product
            if nfrq > 1:
                unperturbed_grad_backprop = np.zeros(unperturbed_grad.shape)
                for i in range(nfrq):
                    unperturbed_grad_backprop[:, i] = tensor_jacobian_product(
                        self.mapping, 0
                    )(self.p, filter_radius, eta, beta, unperturbed_grad[:, i])
            else:
                unperturbed_grad_backprop = tensor_jacobian_product(self.mapping, 0)(
                    self.p, filter_radius, eta, beta, unperturbed_grad
                )

            # compute objective for perturbed design
            perturbed_val, _ = self.adjoint_solver(
                self.mapping(self.p + self.dp, filter_radius, eta, beta),
                MonitorObject.EIGENMODE,
                frequencies,
                need_gradient=False,
            )

            # compare directional derivative
            if unperturbed_grad_backprop.ndim < 2:
                unperturbed_grad_backprop = np.expand_dims(
                    unperturbed_grad_backprop, axis=1
                )
            adj_dd = (self.dp[None, :] @ unperturbed_grad_backprop).flatten()
            fnd_dd = perturbed_val - unperturbed_val

            # since the eigenmode source uses the mode profile
            # of the center frequency, the results for the
            # non-center frequencies of a multifrequency simulation
            # are expected to be less accurate than the center frequency
            if nfrq == 1 and frequencies[0] == self.fcen:
                tol = 2.1e-4 if mp.is_single_precision() else 5e-6
            else:
                tol = 0.005 if mp.is_single_precision() else 0.002

            self.assertClose(adj_dd, fnd_dd, epsilon=tol)
            print(
                f"PASSED: frequencies={frequencies}, "
                f"adjoint gradient={adj_dd}, finite difference={fnd_dd}"
            )

    def test_complex_fields(self):
        """Verifies that the adjoint gradient for an objective function based
        on the DFT fields obtained from complex time-dependent fields agrees
        with the finite-difference approximation."""
        print("*** TESTING COMPLEX FIELDS ***")

        for frequencies in self.mon_frqs:
            # compute objective value and its gradient for unperturbed design
            unperturbed_val, unperturbed_grad = self.adjoint_solver_complex_fields(
                self.p, frequencies
            )

            # compute objective value perturbed design
            perturbed_val, _ = self.adjoint_solver_complex_fields(
                self.p + self.dp, frequencies, need_gradient=False
            )

            # compare directional derivative
            if unperturbed_grad.ndim < 2:
                unperturbed_grad = np.expand_dims(unperturbed_grad, axis=1)
            adj_dd = (self.dp[None, :] @ unperturbed_grad).flatten()
            fnd_dd = perturbed_val - unperturbed_val
            tol = 0.025 if mp.is_single_precision() else 0.001
            self.assertClose(adj_dd, fnd_dd, epsilon=tol)
            print(
                f"PASSED: frequencies={frequencies}"
                f"adjoint gradient={adj_dd}, finite difference={fnd_dd}"
            )

    def test_damping(self):
        """Verifies that the adjoint gradient for a design region with a non-zero
        conductivity agrees with the finite-difference approximation."""
        print("*** TESTING CONDUCTIVITY ***")

        for frequencies in [self.mon_frqs[1]]:
            # compute objective value and its gradient for unperturbed design
            unperturbed_val, unperturbed_grad = self.adjoint_solver_damping(
                self.p, frequencies
            )

            # compute objective value perturbed design
            perturbed_val, _ = self.adjoint_solver_damping(
                self.p + self.dp, frequencies, need_gradient=False
            )

            # compare directional derivative
            if unperturbed_grad.ndim < 2:
                unperturbed_grad = np.expand_dims(unperturbed_grad, axis=1)
            adj_dd = (self.dp[None, :] @ unperturbed_grad).flatten()
            fnd_dd = perturbed_val - unperturbed_val
            tol = 0.04 if mp.is_single_precision() else 0.01
            self.assertClose(adj_dd, fnd_dd, epsilon=tol)

            print(
                f"PASSED: frequencies={frequencies}"
                f"adjoint gradient={adj_dd}, finite difference={fnd_dd}"
            )

    def test_offdiagonal(self):
        """Verifies that the adjoint gradient for a design region involving an
        anisotropic material with non-zero off-diagonal entries of the
        permittivity tensor agrees with the finite-difference approxmiation."""
        print("*** TESTING ANISOTROPIC ε ***")
        filt = lambda x: mpa.conic_filter(
            x.reshape((self.Nx, self.Ny)),
            0.25,
            self.design_region_size.x,
            self.design_region_size.y,
            self.design_region_resolution,
        ).flatten()

        for frequencies in self.mon_frqs:
            nfrq = len(frequencies)

            # compute objective value and its gradient for unperturbed design
            unperturbed_val, unperturbed_grad = self.adjoint_solver(
                filt(self.p), MonitorObject.EIGENMODE, frequencies, self.sapphire
            )

            # backpropagate the gradient using vector-Jacobian product
            if nfrq > 1:
                unperturbed_grad_backprop = np.zeros(unperturbed_grad.shape)
                for i in range(nfrq):
                    unperturbed_grad_backprop[:, i] = tensor_jacobian_product(filt, 0)(
                        self.p, unperturbed_grad[:, i]
                    )
            else:
                unperturbed_grad_backprop = tensor_jacobian_product(filt, 0)(
                    self.p, unperturbed_grad
                )

            # compute objective value perturbed design
            perturbed_val, _ = self.adjoint_solver(
                filt(self.p + self.dp),
                MonitorObject.EIGENMODE,
                frequencies,
                self.sapphire,
                need_gradient=False,
            )

            # compare directional derivative
            if unperturbed_grad_backprop.ndim < 2:
                unperturbed_grad_backprop = np.expand_dims(
                    unperturbed_grad_backprop, axis=1
                )
            adj_dd = (self.dp[None, :] @ unperturbed_grad_backprop).flatten()
            fnd_dd = perturbed_val - unperturbed_val

            # since the eigenmode source uses the mode profile
            # of the center frequency, the results for the
            # non-center frequencies of a multifrequency simulation
            # are expected to be *less* accurate than the center frequency
            if nfrq == 1 and frequencies[0] == self.fcen:
                tol = 0.04 if mp.is_single_precision() else 0.002
            else:
                tol = 0.05 if mp.is_single_precision() else 0.005

            self.assertClose(adj_dd, fnd_dd, epsilon=tol)
            print(
                f"PASSED: frequencies={frequencies}, "
                f"adjoint gradient={adj_dd}, finite difference={fnd_dd}"
            )

    def test_two_objfunc(self):
        """Verifies that the adjoint gradients from two objective functions
        each agree with the finite-difference approximation."""
        print("*** TESTING TWO OBJECTIVE FUNCTIONS ***")

        for frequencies in self.mon_frqs:
            nfrq = len(frequencies)

            # compute objective value and its gradient for unperturbed design
            unperturbed_val, unperturbed_grad = self.adjoint_solver_two_objfunc(
                self.p, frequencies
            )

            # compute objective value for perturbed design
            perturbed_val, _ = self.adjoint_solver_two_objfunc(
                self.p + self.dp,
                frequencies,
                need_gradient=False,
            )

            # since the eigenmode source uses the mode profile
            # of the center frequency, the results for the
            # non-center frequencies of a multifrequency simulation
            # are expected to be *less* accurate than the center frequency
            if nfrq == 1 and frequencies[0] == self.fcen:
                tol = 0.05 if mp.is_single_precision() else 0.0001
            else:
                tol = 0.15 if mp.is_single_precision() else 0.001

            for m in [0, 1]:
                frq_slice = slice(0, nfrq, 1) if m == 0 else slice(nfrq, 2 * nfrq, 1)
                adj_dd = (self.dp[None, :] @ unperturbed_grad[:, frq_slice]).flatten()
                fnd_dd = perturbed_val[frq_slice] - unperturbed_val[frq_slice]
                self.assertClose(adj_dd, fnd_dd, epsilon=tol)
                print(
                    f"PASSED: frequencies={frequencies}, m={m}, "
                    f"adjoint gradient={adj_dd}, finite difference={fnd_dd}"
                )

    def test_multifreq_monitor(self):
        """Verifies that the individual adjoint gradients from a multifrequency
        eigenmode-coefficient monitor are equivalent to a single-frequency
        monitor."""
        print("*** TESTING MULTIFREQUENCY MONITOR ***")

        nfrq = 5
        frqs = np.linspace(
            self.fcen - 0.2 * self.df,
            self.fcen + 0.2 * self.df,
            nfrq,
        )
        multifrq_val, multifrq_grad = self.adjoint_solver_two_objfunc(
            self.p,
            frqs,
        )

        tol = 0.005 if mp.is_single_precision() else 0.004
        for n in range(nfrq):
            frq = frqs[n]
            singlefreq_val, singlefreq_grad = self.adjoint_solver_two_objfunc(
                self.p,
                [frq],
            )
            for m in [0, 1]:
                s = n + m * nfrq
                self.assertAlmostEqual(singlefreq_val[m], multifrq_val[s], places=6)
                self.assertClose(
                    singlefreq_grad[:, m], multifrq_grad[:, s], epsilon=tol
                )
                print(f"PASSED: frequency={frq:.5f}, m={m}.")

    def test_mode_source_bandwidth(self):
        """Verifies that the accuracy of the adjoint gradient of an
        eigenmode-coefficient monitor at a single frequency is independent
        of the bandwidth of the pulsed mode source."""
        print("*** TESTING MODE SOURCE BANDWIDTH ***")

        # objective function value for unperturbed design
        unperturbed_objf, _ = self.adjoint_solver_two_objfunc(
            self.p,
            [self.fcen],
            need_gradient=False,
        )

        # objective function value for perturbed design
        perturbed_objf, _ = self.adjoint_solver_two_objfunc(
            self.p + self.dp,
            [self.fcen],
            need_gradient=False,
        )

        # finite-difference approximation
        # for the directional derivative
        fnd_dd = perturbed_objf - unperturbed_objf

        nfw = 5
        fw = np.linspace(
            0.05,
            0.25,
            nfw,
        )

        # minimum error for directional derivative of adjoint
        # gradient relative to finite-difference approximation.
        # note: this can be reduced to 0.0001 if `decay_by=1e-12`
        # is used in the `OptimizationProblem` constructor inside
        # `adjoint_solver_two_objfunc`.
        tol = 0.05

        for n in range(nfw):
            fwidth = fw[n] * self.fcen
            self.mode_source = [
                mp.EigenModeSource(
                    src=mp.GaussianSource(self.fcen, fwidth=fwidth),
                    center=mp.Vector3(-0.5 * self.sxy + self.dpml, 0),
                    size=mp.Vector3(0, self.sxy - 2 * self.dpml),
                    eig_parity=self.eig_parity,
                )
            ]

            objf, grad = self.adjoint_solver_two_objfunc(
                self.p,
                [self.fcen],
            )

            for m in [0, 1]:
                self.assertAlmostEqual(unperturbed_objf[m], objf[m], places=6)
                adj_dd = (self.dp[None, :] @ grad[:, m]).flatten()
                rel_err = abs((fnd_dd[m] - adj_dd[0]) / fnd_dd[m])
                self.assertLessEqual(rel_err, tol)
                print(f"PASSED: fwidth={fwidth:.5f}, m={m}, err={rel_err:.10f}")


if __name__ == "__main__":
    unittest.main()
