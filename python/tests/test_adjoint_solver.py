import meep as mp

try:
    import meep.adjoint as mpa
except:
    import adjoint as mpa

import unittest
from enum import Enum

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

        cls.design_region_size = mp.Vector3(1.5, 1.5)
        cls.design_region_resolution = int(2 * cls.resolution)
        cls.Nx = int(cls.design_region_size.x * cls.design_region_resolution)
        cls.Ny = int(cls.design_region_size.y * cls.design_region_resolution)

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

        cls.fcen = 1 / 1.55
        cls.df = 0.2 * cls.fcen
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
                component=mp.Ez,
            )
        ]

        cls.line_source = [
            mp.Source(
                src=mp.GaussianSource(cls.fcen, fwidth=cls.df),
                center=mp.Vector3(-0.85, 0),
                size=mp.Vector3(0, cls.sxy - 2 * cls.dpml),
                component=mp.Ez,
            )
        ]

        cls.k_point = mp.Vector3(0.23, -0.38)

    def adjoint_solver(
        self,
        design_params,
        mon_type: MonitorObject,
        frequencies=None,
        mat2=None,
        need_gradient=True,
    ):
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

        elif mon_type.name == "DFT":
            obj_list = [
                mpa.FourierFields(
                    sim,
                    mp.Volume(center=mp.Vector3(1.25), size=mp.Vector3(0.25, 1, 0)),
                    mp.Ez,
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

        if need_gradient:
            f, dJ_du = opt([design_params])
            return f, dJ_du
        else:
            f = opt([design_params], need_gradient=False)
            return f[0]

    def adjoint_solver_complex_fields(
        self, design_params, frequencies=None, need_gradient=True
    ):
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

        if need_gradient:
            f, dJ_du = opt([design_params])
            return f, dJ_du
        else:
            f = opt([design_params], need_gradient=False)
            return f[0]

    def adjoint_solver_damping(
        self, design_params, frequencies=None, mat2=None, need_gradient=True
    ):
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

        if need_gradient:
            f, dJ_du = opt([design_params])
            return f, dJ_du
        else:
            f = opt([design_params], need_gradient=False)
            return f[0]

    def mapping(self, x, filter_radius, eta, beta):
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
        print("*** TESTING DFT OBJECTIVE ***")

        # test the single frequency and multi frequency cases
        for frequencies in [[self.fcen], [1 / 1.58, self.fcen, 1 / 1.53]]:
            # compute objective value and its gradient for unperturbed design
            unperturbed_val, unperturbed_grad = self.adjoint_solver(
                self.p, MonitorObject.DFT, frequencies
            )

            # compute objective value for perturbed design
            perturbed_val = self.adjoint_solver(
                self.p + self.dp, MonitorObject.DFT, frequencies, need_gradient=False
            )

            # compare directional derivative
            if unperturbed_grad.ndim < 2:
                unperturbed_grad = np.expand_dims(unperturbed_grad, axis=1)
            adj_dd = (self.dp[None, :] @ unperturbed_grad).flatten()
            fnd_dd = perturbed_val - unperturbed_val
            print(
                f"directional derivative:, {adj_dd} (adjoint solver), {fnd_dd} (finite difference)"
            )

            tol = 0.07 if mp.is_single_precision() else 0.006
            self.assertClose(adj_dd, fnd_dd, epsilon=tol)

    def test_eigenmode(self):
        print("*** TESTING EIGENMODE OBJECTIVE ***")

        # test the single frequency and multi frequency case
        for frequencies in [[self.fcen], [1 / 1.58, self.fcen, 1 / 1.53]]:
            # compute objective value and its gradient for unperturbed design
            unperturbed_val, unperturbed_grad = self.adjoint_solver(
                self.p, MonitorObject.EIGENMODE, frequencies
            )

            # compute objective for perturbed design
            perturbed_val = self.adjoint_solver(
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
            print(
                f"directional derivative:, {adj_dd} (adjoint solver), {fnd_dd} (finite difference)"
            )

            tol = 0.04 if mp.is_single_precision() else 0.01
            self.assertClose(adj_dd, fnd_dd, epsilon=tol)

    def test_ldos(self):
        print("*** TESTING LDOS OBJECTIVE ***")

        # test the single frequency and multi frequency cases
        for frequencies in [[self.fcen], [1 / 1.58, self.fcen, 1 / 1.53]]:
            # compute objective value and its gradient for unperturbed design
            unperturbed_val, unperturbed_grad = self.adjoint_solver(
                self.p, MonitorObject.LDOS, frequencies
            )

            # compute objective for perturbed design
            perturbed_val = self.adjoint_solver(
                self.p + self.dp, MonitorObject.LDOS, frequencies, need_gradient=False
            )

            # compare directional derivative
            if unperturbed_grad.ndim < 2:
                unperturbed_grad = np.expand_dims(unperturbed_grad, axis=1)
            adj_dd = (self.dp[None, :] @ unperturbed_grad).flatten()
            fnd_dd = perturbed_val - unperturbed_val
            print(
                f"directional derivative:, {adj_dd} (adjoint solver), {fnd_dd} (finite difference)"
            )

            tol = 0.07 if mp.is_single_precision() else 0.006
            self.assertClose(adj_dd, fnd_dd, epsilon=tol)

    def test_gradient_backpropagation(self):
        print("*** TESTING GRADIENT BACKPROPAGATION ***")

        # filter/thresholding parameters
        filter_radius = 0.21985
        eta = 0.49093
        beta = 4.0698

        for frequencies in [[self.fcen], [1 / 1.58, self.fcen, 1 / 1.53]]:
            mapped_p = self.mapping(self.p, filter_radius, eta, beta)

            # compute objective value and its gradient for unperturbed design
            unperturbed_val, unperturbed_grad = self.adjoint_solver(
                mapped_p, MonitorObject.EIGENMODE, frequencies
            )

            # backpropagate the gradient using vector-Jacobian product
            if len(frequencies) > 1:
                unperturbed_grad_backprop = np.zeros(unperturbed_grad.shape)
                for i in range(len(frequencies)):
                    unperturbed_grad_backprop[:, i] = tensor_jacobian_product(
                        self.mapping, 0
                    )(self.p, filter_radius, eta, beta, unperturbed_grad[:, i])
            else:
                unperturbed_grad_backprop = tensor_jacobian_product(self.mapping, 0)(
                    self.p, filter_radius, eta, beta, unperturbed_grad
                )

            # compute objective for perturbed design
            perturbed_val = self.adjoint_solver(
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
            print(
                f"directional derivative:, {adj_dd} (adjoint solver), {fnd_dd} (finite difference)"
            )

            tol = 0.025 if mp.is_single_precision() else 0.01
            self.assertClose(adj_dd, fnd_dd, epsilon=tol)

    def test_complex_fields(self):
        print("*** TESTING COMPLEX FIELDS ***")

        for frequencies in [[self.fcen], [1 / 1.58, self.fcen, 1 / 1.53]]:
            # compute objective value and its gradient for unperturbed design
            unperturbed_val, unperturbed_grad = self.adjoint_solver_complex_fields(
                self.p, frequencies
            )

            # compute objective value perturbed design
            perturbed_val = self.adjoint_solver_complex_fields(
                self.p + self.dp, frequencies, need_gradient=False
            )

            # compare directional derivative
            if unperturbed_grad.ndim < 2:
                unperturbed_grad = np.expand_dims(unperturbed_grad, axis=1)
            adj_dd = (self.dp[None, :] @ unperturbed_grad).flatten()
            fnd_dd = perturbed_val - unperturbed_val
            print(
                f"directional derivative:, {adj_dd} (adjoint solver), {fnd_dd} (finite difference)"
            )

            tol = 0.06 if mp.is_single_precision() else 0.01
            self.assertClose(adj_dd, fnd_dd, epsilon=tol)

    def test_damping(self):
        print("*** TESTING CONDUCTIVITY ***")

        for frequencies in [[1 / 1.58, self.fcen, 1 / 1.53]]:
            # compute objective value and its gradient for unperturbed design
            unperturbed_val, unperturbed_grad = self.adjoint_solver_damping(
                self.p, frequencies
            )

            # compute objective value perturbed design
            perturbed_val = self.adjoint_solver_damping(
                self.p + self.dp, frequencies, need_gradient=False
            )

            # compare directional derivative
            if unperturbed_grad.ndim < 2:
                unperturbed_grad = np.expand_dims(unperturbed_grad, axis=1)
            adj_dd = (self.dp[None, :] @ unperturbed_grad).flatten()
            fnd_dd = perturbed_val - unperturbed_val
            print(
                f"directional derivative:, {adj_dd} (adjoint solver), {fnd_dd} (finite difference)"
            )

            tol = 0.06 if mp.is_single_precision() else 0.04
            self.assertClose(adj_dd, fnd_dd, epsilon=tol)

    def test_offdiagonal(self):
        print("*** TESTING ANISOTROPIC ε ***")
        filt = lambda x: mpa.conic_filter(
            x.reshape((self.Nx, self.Ny)),
            0.25,
            self.design_region_size.x,
            self.design_region_size.y,
            self.design_region_resolution,
        ).flatten()

        # test the single frequency and multi frequency case
        for frequencies in [[self.fcen], [1 / 1.58, self.fcen, 1 / 1.53]]:
            # compute objective value and its gradient for unperturbed design
            unperturbed_val, unperturbed_grad = self.adjoint_solver(
                filt(self.p), MonitorObject.EIGENMODE, frequencies, self.sapphire
            )

            # backpropagate the gradient using vector-Jacobian product
            if len(frequencies) > 1:
                unperturbed_grad_backprop = np.zeros(unperturbed_grad.shape)
                for i in range(len(frequencies)):
                    unperturbed_grad_backprop[:, i] = tensor_jacobian_product(filt, 0)(
                        self.p, unperturbed_grad[:, i]
                    )
            else:
                unperturbed_grad_backprop = tensor_jacobian_product(filt, 0)(
                    self.p, unperturbed_grad
                )

            # compute objective value perturbed design
            perturbed_val = self.adjoint_solver(
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
            print(
                f"directional derivative:, {adj_dd} (adjoint solver), {fnd_dd} (finite difference)"
            )

            tol = 0.1 if mp.is_single_precision() else 0.04
            self.assertClose(adj_dd, fnd_dd, epsilon=tol)


if __name__ == "__main__":
    unittest.main()
