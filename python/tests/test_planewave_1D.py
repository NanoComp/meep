# Verifies properties of the fields of a planewave using a 1D or 3D simulation.


from typing import Tuple
import unittest

import meep as mp
import numpy as np


class TestPlanewave1D(unittest.TestCase):
    def planewave_in_vacuum(
        self,
        resolution_um: float,
        polar_rad: float,
        azimuth_rad: float,
        cell_dim: int,
        yee_grid: bool,
    ) -> None:
        """
        Verifies properties of the fields of a planewave in 1D or 3D.

        Computes the DFT fields at two different positions in the 1D grid and
        determines that these fields:
            (1) have the expected relative phase, and
            (2) are the same when obtained as a single point or a slice of an
                array over the entire cell.

        Args:
            resolution_um: resolution of the grid (number of pixels per um).
            polar_rad: polar angle of the incident planewave in [0, π].
                0 is +z axis.
            azimuth_rad: azimuth angle of the incident planewave in [0, 2π].
                Rotation around the z axis. 0 is +x axis.
            cell_dim: dimension of the cell (1 or 3).
            yee_grid: whether the DFT fields are on a centered or Yee grid.
        """
        print(
            f"Testing planewaves in vacuum using {cell_dim}D simulation and "
            f"{'yee' if yee_grid else 'centered'} grid..."
        )

        pml_um = 1.0
        air_um = 10.0
        size_z_um = pml_um + air_um + pml_um
        cell_size = mp.Vector3(0, 0, size_z_um)
        pml_layers = [mp.PML(thickness=pml_um)]

        wavelength_um = 1.0
        frequency = 1 / wavelength_um
        kx = frequency * np.sin(polar_rad) * np.cos(azimuth_rad)
        ky = frequency * np.sin(polar_rad) * np.sin(azimuth_rad)
        kz = frequency * np.cos(polar_rad)

        if cell_dim == 1 and polar_rad != 0 and azimuth_rad != 0:
            raise ValueError("An oblique planewave cannot be simulated in 1D.")

        if cell_dim == 1:
            k_point = False
            dimensions = 1
        elif cell_dim == 3:
            k_point = mp.Vector3(kx, ky, 0)
            dimensions = 3
        else:
            raise ValueError("cell_dim can only be 1 or 3.")

        src_cmpt = mp.Ex
        sources = [
            mp.Source(
                src=mp.GaussianSource(frequency, fwidth=0.1 * frequency),
                component=src_cmpt,
                center=mp.Vector3(0, 0, -0.5 * air_um),
            )
        ]

        sim = mp.Simulation(
            resolution=resolution_um,
            force_complex_fields=True,
            cell_size=cell_size,
            sources=sources,
            boundary_layers=pml_layers,
            k_point=k_point,
            dimensions=dimensions,
        )

        dft_fields = sim.add_dft_fields(
            [mp.Ex],
            frequency,
            0,
            1,
            center=mp.Vector3(0, 0, 0),
            size=mp.Vector3(0, 0, size_z_um),
            yee_grid=yee_grid,
        )

        z_pos_1 = -0.234804
        dft_fields_1 = sim.add_dft_fields(
            [mp.Ex],
            frequency,
            0,
            1,
            center=mp.Vector3(0, 0, z_pos_1),
            size=mp.Vector3(0, 0, 0),
            yee_grid=yee_grid,
        )

        z_pos_2 = 2.432973
        dft_fields_2 = sim.add_dft_fields(
            [mp.Ex],
            frequency,
            0,
            1,
            center=mp.Vector3(0, 0, z_pos_2),
            size=mp.Vector3(0, 0, 0),
            yee_grid=yee_grid,
        )

        sim.run(
            until_after_sources=mp.stop_when_fields_decayed(
                25, src_cmpt, mp.Vector3(0, 0, 0.5 * air_um), 1e-6
            )
        )

        # Check that the relative phase of the fields obtained from a slice of
        # an array of the fields over the entire cell matches the analytic
        # result.
        ex_dft = sim.get_dft_array(dft_fields, mp.Ex, 0)
        z_um = np.linspace(-0.5 * size_z_um, 0.5 * size_z_um, len(ex_dft))
        z_idx_1 = np.argmin(np.abs(z_pos_1 - z_um))
        z_idx_2 = np.argmin(np.abs(z_pos_2 - z_um))
        meep_phase = ex_dft[z_idx_2] / ex_dft[z_idx_1]
        expected_phase = np.exp(1j * 2 * np.pi * kz * (z_um[z_idx_2] - z_um[z_idx_1]))
        self.assertAlmostEqual(meep_phase, expected_phase, delta=0.05)

        # Check that the relative phase of the fields obtained from a point
        # location matches the analytic result.
        ex_dft_pos_1 = sim.get_dft_array(dft_fields_1, mp.Ex, 0)
        ex_dft_pos_2 = sim.get_dft_array(dft_fields_2, mp.Ex, 0)
        meep_phase = ex_dft_pos_2 / ex_dft_pos_1
        expected_phase = np.exp(1j * 2 * np.pi * kz * (z_pos_2 - z_pos_1))
        self.assertAlmostEqual(meep_phase, expected_phase, delta=0.05)

        # Check that the fields obtained using the two approaches match.
        self.assertAlmostEqual(ex_dft[z_idx_1], ex_dft_pos_1, delta=0.08)
        self.assertAlmostEqual(ex_dft[z_idx_2], ex_dft_pos_2, delta=0.08)

        print("PASSED.")

    def test_planewave_1D(self):
        # Case 1: normal incidence with k = (0, 0, kz).
        polar_rad = 0
        azimuth_rad = 0
        self.planewave_in_vacuum(400.0, polar_rad, azimuth_rad, 1, False)
        self.planewave_in_vacuum(200.0, polar_rad, azimuth_rad, 3, False)

        self.planewave_in_vacuum(400.0, polar_rad, azimuth_rad, 1, True)
        self.planewave_in_vacuum(200.0, polar_rad, azimuth_rad, 3, True)

        # Case 2: oblique incidence with k = (kx, ky, kz).
        polar_rad = np.deg2rad(10.3)
        azimuth_rad = np.deg2rad(5.7)
        self.planewave_in_vacuum(200.0, polar_rad, azimuth_rad, 3, True)


if __name__ == "__main__":
    unittest.main()
