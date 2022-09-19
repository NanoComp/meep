import unittest

import numpy as np

import meep as mp

# Test that is_integrated=True source correctly generates planewaves
# for sources extending into the PML, as in this tutorial:
#      https://meep.readthedocs.io/en/latest/Perfectly_Matched_Layer/#planewave-sources-extending-into-pml
# Regression test for issue #2043.


class TestIntegratedSource(unittest.TestCase):
    def test_integrated_source(self):
        sources = [
            mp.Source(
                mp.ContinuousSource(1, is_integrated=True),
                center=mp.Vector3(-2),
                size=mp.Vector3(y=6),
                component=mp.Ez,
            )
        ]
        sim = mp.Simulation(
            resolution=20,
            cell_size=(6, 6),
            boundary_layers=[mp.PML(thickness=1)],
            sources=sources,
            k_point=mp.Vector3(),
        )
        sim.run(until=30)

        # field in mid-plane should be nearly constant,
        # so compute its normalized std. dev. and check << 1
        ez = sim.get_array(mp.Ez, center=mp.Vector3(2), size=mp.Vector3(y=6))
        std = np.std(ez) / np.sqrt(np.mean(ez**2))
        print("std = ", std)
        self.assertAlmostEqual(std, 0.0, places=4 if mp.is_single_precision() else 8)


if __name__ == "__main__":
    unittest.main()
