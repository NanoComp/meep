# dispersive_eigenmode.py - Tests the meep eigenmode features (eigenmode source,
# eigenmode decomposition, and get_eigenmode) with dispersive materials.
# TODO:
#  * check materials with off diagonal components
#  * check magnetic profiles
#  * once imaginary component is supported, check that
import os
import unittest

import h5py
import numpy as np
from utils import ApproxComparisonTestCase

import meep as mp


class TestDispersiveEigenmode(ApproxComparisonTestCase):
    # ----------------------------------------- #
    # ----------- Helper Functions ------------ #
    # ----------------------------------------- #
    # Directly cals the C++ chi1 routine
    def call_chi1(self, material, frequency):

        sim = mp.Simulation(
            cell_size=mp.Vector3(1, 1, 1), default_material=material, resolution=20
        )

        sim.init_sim()
        v3 = mp.py_v3_to_vec(sim.dimensions, mp.Vector3(0, 0, 0), sim.is_cylindrical)
        chi1inv = np.zeros((3, 3), dtype=np.complex128)
        for i, com in enumerate([mp.Ex, mp.Ey, mp.Ez]):
            for k, dir in enumerate([mp.X, mp.Y, mp.Z]):
                chi1inv[i, k] = sim.structure.get_chi1inv(com, dir, v3, frequency)
        n = np.real(np.sqrt(np.linalg.inv(chi1inv.astype(np.complex128))))

        n_actual = np.real(np.sqrt(material.epsilon(frequency).astype(np.complex128)))

        tol = 1e-6 if mp.is_single_precision() else 1e-8
        self.assertClose(n, n_actual, epsilon=tol)

    @classmethod
    def setUpClass(cls):
        cls.temp_dir = mp.make_output_directory()

    @classmethod
    def tearDownClass(cls):
        mp.delete_directory(cls.temp_dir)

    def verify_output_and_slice(self, material, frequency):
        # Since the slice routines average the diagonals, we need to do that too:
        chi1 = material.epsilon(frequency).astype(np.complex128)
        chi1inv = np.linalg.inv(chi1)
        chi1inv = np.diag(chi1inv)
        N = chi1inv.size
        n = np.sqrt(N / np.sum(chi1inv))

        sim = mp.Simulation(
            cell_size=mp.Vector3(2, 2, 2),
            default_material=material,
            resolution=20,
            eps_averaging=False,
        )
        sim.use_output_directory(self.temp_dir)
        sim.init_sim()

        # Check to make sure the get_slice routine is working with frequency
        n_slice = np.sqrt(np.max(sim.get_epsilon(frequency=frequency, snap=True)))
        self.assertAlmostEqual(n, n_slice, places=4)

        # Check to make sure h5 output is working with frequency
        filename = os.path.join(
            self.temp_dir, f"{sim.get_filename_prefix()}-eps-000000.00.h5"
        )

        mp.output_epsilon(sim, frequency=frequency)
        n_h5 = 0
        mp.all_wait()
        with h5py.File(filename, "r") as f:
            n_h5 = np.sqrt(np.max(mp.complexarray(f["eps.r"][()], f["eps.i"][()])))
        self.assertAlmostEqual(n, n_h5, places=4)

    # ----------------------------------------- #
    # ----------- Test Routines --------------- #
    # ----------------------------------------- #
    def test_chi1_routine(self):
        # Checks the newly implemented get_chi1inv routines within the
        # fields and structure classes by comparing their output to the
        # python epsilon output.

        from meep.materials import Ag, Au, LiNbO3, Si

        # Check Silicon
        w0 = Si.valid_freq_range.min
        w1 = Si.valid_freq_range.max
        self.call_chi1(Si, w0)
        self.call_chi1(Si, w1)

        # Check Silver
        w0 = Ag.valid_freq_range.min
        w1 = Ag.valid_freq_range.max
        self.call_chi1(Ag, w0)
        self.call_chi1(Ag, w1)

        # Check Gold
        w0 = Au.valid_freq_range.min
        w1 = Au.valid_freq_range.max
        self.call_chi1(Au, w0)
        self.call_chi1(Au, w1)

        # Check Lithium Niobate (X,X)
        w0 = LiNbO3.valid_freq_range.min
        w1 = LiNbO3.valid_freq_range.max
        self.call_chi1(LiNbO3, w0)
        self.call_chi1(LiNbO3, w1)

        # Now let's rotate LN
        import copy

        rotLiNbO3 = copy.deepcopy(LiNbO3)
        rotLiNbO3.rotate(mp.Vector3(1, 1, 1), np.radians(34))
        self.call_chi1(rotLiNbO3, w0)
        self.call_chi1(rotLiNbO3, w1)

    def test_get_with_dispersion(self):
        # Checks the get_array_slice and output_fields method
        # with dispersive materials.

        from meep.materials import Ag, Au, LiNbO3, Si

        # Check Silicon
        w0 = Si.valid_freq_range.min
        w1 = Si.valid_freq_range.max
        self.verify_output_and_slice(Si, w0)
        self.verify_output_and_slice(Si, w1)

        # Check Silver
        w0 = Ag.valid_freq_range.min
        w1 = Ag.valid_freq_range.max
        self.verify_output_and_slice(Ag, w0)
        self.verify_output_and_slice(Ag, w1)

        # Check Gold
        w0 = Au.valid_freq_range.min
        w1 = Au.valid_freq_range.max
        self.verify_output_and_slice(Au, w0)
        self.verify_output_and_slice(Au, w1)

        # Check Lithium Niobate
        w0 = LiNbO3.valid_freq_range.min
        w1 = LiNbO3.valid_freq_range.max
        self.verify_output_and_slice(LiNbO3, w0)
        self.verify_output_and_slice(LiNbO3, w1)

        # Now let's rotate LN
        import copy

        rotLiNbO3 = copy.deepcopy(LiNbO3)
        rotLiNbO3.rotate(mp.Vector3(1, 1, 1), np.radians(34))
        self.verify_output_and_slice(rotLiNbO3, w0)
        self.verify_output_and_slice(rotLiNbO3, w1)


if __name__ == "__main__":
    unittest.main()
