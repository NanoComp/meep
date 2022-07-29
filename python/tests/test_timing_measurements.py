import time
import unittest

import meep as mp
from meep import timing_measurements as timing


class TimingTest(unittest.TestCase):
    def test_timing_measurements(self):
        """Tests that timing measurements have expected names and can be updated."""
        sim = mp.Simulation(
            cell_size=mp.Vector3(2, 2, 2),
            resolution=20,
        )
        time_start = time.time()
        sim.run(until=5)
        timing_measurements = timing.MeepTimingMeasurements.new_from_simulation(sim)

        # Check for expected names after updating
        self.assertSetEqual(
            set(timing_measurements.measurement_names),
            set(timing.TIMING_MEASUREMENT_IDS.keys()),
        )
        self.assertTrue(
            timing_measurements.elapsed_time > 0
            or timing_measurements.elapsed_time == -1
        )
        self.assertGreater(timing_measurements.num_time_steps, 0)
        self.assertGreaterEqual(timing_measurements.comm_efficiency, 0)
        self.assertGreaterEqual(timing_measurements.comm_efficiency_one_to_one, 0)
        self.assertGreaterEqual(timing_measurements.comm_efficiency_all_to_all, 0)


if __name__ == "__main__":
    unittest.main()
