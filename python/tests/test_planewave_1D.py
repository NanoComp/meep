# Verifies properties of the fields of a planewave in 1D simulated in 1D or 3D.


from typing import Tuple
import unittest

import meep as mp
import numpy as np


def complex_str(complex_val: complex) -> str:
    return f"{complex_val.real:+.6f}{complex_val.imag:+.6f}*1j"


class TestPlanewave1D(unittest.TestCase):
    def planewave_in_vacuum(
        self, resolution_um: float, cell_dim: int, yee_grid: bool
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Verifies properties of the fields of a planewave in 1D or 3D.

        Computes the DFT fields for Ex at two different positions in the grid
        and determines that these fields:
            (1) have the expected relative phase, and
            (2) are the same when obtained as a single point or an array slice.

        Args:
            resolution_um: resolution of the grid in number of pixels per um.
            cell_dim: dimension of the cell: 1 or 3.
            yee_grid: whether the DFT fields obtained on a centered or Yee grid.
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

        if cell_dim == 1:
            k_point = False
            dimensions = 1
        elif cell_dim == 3:
            k_point = mp.Vector3(0, 0, 0)
            dimensions = 3
        else:
            raise ValueError("cell_dim can only be 1 or 3.")

        src_cmpt = mp.Ex
        wavelength_um = 1.0
        frequency = 1 / wavelength_um
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
                25, src_cmpt, mp.Vector3(0, 0, 0.5 * air_um), 1e-7
            )
        )

        ex_dft = sim.get_dft_array(dft_fields, mp.Ex, 0)
        z_um = np.linspace(-0.5 * size_z_um, 0.5 * size_z_um, len(ex_dft))
        z_idx_1 = np.argmin(np.abs(z_pos_1 - z_um))
        z_idx_2 = np.argmin(np.abs(z_pos_2 - z_um))
        meep_phase = ex_dft[z_idx_2] / ex_dft[z_idx_1]
        expected_phase = np.exp(1j * 2 * np.pi * (z_um[z_idx_2] - z_um[z_idx_1]))
        self.assertAlmostEqual(meep_phase, expected_phase, delta=0.05)

        ex_dft_pos_1 = sim.get_dft_array(dft_fields_1, mp.Ex, 0)
        ex_dft_pos_2 = sim.get_dft_array(dft_fields_2, mp.Ex, 0)
        meep_phase = ex_dft_pos_2 / ex_dft_pos_1
        expected_phase = np.exp(1j * 2 * np.pi * (z_pos_2 - z_pos_1))
        self.assertAlmostEqual(meep_phase, expected_phase, delta=0.05)

        self.assertAlmostEqual(ex_dft[z_idx_1], ex_dft_pos_1, delta=0.08)
        self.assertAlmostEqual(ex_dft[z_idx_2], ex_dft_pos_2, delta=0.08)

        print("PASSED.")

    def test_planewave_1D(self):
        self.planewave_in_vacuum(400.0, 1, False)
        self.planewave_in_vacuum(200.0, 3, False)

        self.planewave_in_vacuum(400.0, 1, True)
        self.planewave_in_vacuum(200.0, 3, True)


if __name__ == "__main__":
    unittest.main()
