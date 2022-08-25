import unittest

import meep.binary_partition_utils as bpu
import numpy as np
import parameterized
from meep.chunk_balancer import ChunkBalancer
from meep.timing_measurements import TIMING_MEASUREMENT_IDS, MeepTimingMeasurements

import meep as mp


class MockSimulation(mp.Simulation):
    """Class which emulates the multi-core MPI behavior while on a single core.

    This inherits all methods from mp.Simulation but overrides two behaviors:
      1. sim.time_spent_on() will provide fake timing data where all entries are
         [0.0 ... 0.0] except for time_stepping, which will have an array of
         values where the elapsed time for each process is equal to the pixel
         volume of that process's chunk.
      2. sim.structure.get_chunk_owners will return a list of process IDs which
         would be the values that you would see if running the simulation in MPI
         mode. (Otherwise in a single-core test environment, an array of zeros
         would be returned.)
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.relative_loads = None

    def _structure_get_chunk_owners(self):
        # Hacky workaround to make this test work on single-core systems
        proc_ids = [
            leaf.proc_id for leaf in bpu.enumerate_leaf_nodes(self.chunk_layout)
        ]

        return np.array(proc_ids)

    def init_sim(self):
        super().init_sim()
        setattr(self.structure, "get_chunk_owners", self._structure_get_chunk_owners)

    def time_spent_on(self, time_sink: int):
        # We're going to pretend the amount of time spent is ~volume so that
        # chunks should converge to an equal volume
        num_processes = self.chunk_layout.numchunks()
        chunk_volumes = self.structure.get_chunk_volumes()
        chunk_owners = self.structure.get_chunk_owners()

        if chunk_volumes[0].nz() > 0:  # 3D case
            volumes = [v.nx() * v.ny() * v.nz() for v in chunk_volumes]
        else:  # 2D case
            volumes = [v.nx() * v.ny() for v in chunk_volumes]

        total_volume_by_proc = np.zeros(num_processes)
        for i, v in enumerate(volumes):
            total_volume_by_proc[chunk_owners[i]] += v

        if time_sink == TIMING_MEASUREMENT_IDS["time_stepping"]:
            times = 1.0 * total_volume_by_proc
        else:
            times = 0.0 * np.ones(num_processes)

        if self.relative_loads is not None:
            times = times * self.relative_loads

        return times


# The test sim instances are lambda functions () => MockSimulation because
# the chunk balancer tests will mutate the chunk_layout attribute.

TEST_SIM_1 = lambda: MockSimulation(
    cell_size=mp.Vector3(10.0, 5.0, 0),
    resolution=20,
    chunk_layout=mp.BinaryPartition(
        data=[(mp.X, -2.0), 0, [(mp.Y, 1.5), [(mp.X, 4.0), 1, [(mp.Y, 0.5), 4, 3]], 2]]
    ),
)
TEST_SIM_2 = lambda: MockSimulation(
    cell_size=mp.Vector3(10.0, 5.0, 3.0),
    resolution=10,
    chunk_layout=mp.BinaryPartition(
        data=[(mp.X, -2.0), 0, [(mp.Y, 1.0), [(mp.X, 3.0), 1, [(mp.Y, 0.5), 4, 3]], 2]]
    ),
)
TEST_SIM_3 = lambda: MockSimulation(
    cell_size=mp.Vector3(6.0, 4.0, 0),
    resolution=10,
    chunk_layout=mp.BinaryPartition(data=[(mp.X, -2.0), 0, [(mp.X, 2.0), 1, 2]]),
)
TEST_SIM_4 = lambda: MockSimulation(
    cell_size=mp.Vector3(6.0, 4.0, 0),
    resolution=10,
    chunk_layout=mp.BinaryPartition(
        data=[(mp.X, -2.0), 0, [(mp.X, -0.5), 1, [(mp.X, 1.0), 2, [(mp.X, 2.0), 3, 4]]]]
    ),
)

TEST_SIM_DUPLICATE_PROC_ID = lambda: MockSimulation(
    cell_size=mp.Vector3(10.0, 5.0, 0),
    resolution=10,
    chunk_layout=mp.BinaryPartition(
        data=[(mp.X, -2.0), 0, [(mp.Y, 1.5), [(mp.X, 4.0), 1, [(mp.Y, 0.5), 2, 1]], 2]]
    ),
)

TEST_CHUNK_DATA_1 = {
    "chunk_layout": mp.BinaryPartition(data=[(mp.X, -2.5), 0, [(mp.X, 2.5), 1, 2]]),
    "cell_size": mp.Vector3(6.0, 4.0, 0),
    "time_stepping": [1.0, 1.0, 1.0],
    "new_chunk_layout": mp.BinaryPartition(data=[(mp.X, -2.5), 0, [(mp.X, 2.5), 1, 2]]),
}
TEST_CHUNK_DATA_2 = {
    "chunk_layout": mp.BinaryPartition(data=[(mp.X, -2.5), 0, [(mp.X, 2.5), 1, 2]]),
    "cell_size": mp.Vector3(6.0, 4.0, 0),
    "time_stepping": [3.0 - 2.5, 2.5 + 2.5, 3.0 - 2.5],
    "new_chunk_layout": mp.BinaryPartition(data=[(mp.X, -1.0), 0, [(mp.X, 1.0), 1, 2]]),
}
TEST_CHUNK_DATA_3 = {
    "chunk_layout": mp.BinaryPartition(data=[(mp.X, 2.0), 0, 1]),
    "cell_size": mp.Vector3(6.0, 4.0, 0),
    "time_stepping": [1.0, 1.0],
    "new_chunk_layout": mp.BinaryPartition(data=[(mp.X, 2.0), 0, 1]),
}
TEST_CHUNK_DATA_4 = {
    "chunk_layout": mp.BinaryPartition(data=[(mp.X, 2.0), 0, 1]),
    "cell_size": mp.Vector3(6.0, 4.0, 0),
    "time_stepping": [5.0, 1.0],
    "new_chunk_layout": mp.BinaryPartition(data=[(mp.X, 0.0), 0, 1]),
}
TEST_CHUNK_DATA_5 = {
    "chunk_layout": mp.BinaryPartition(
        data=[(mp.X, -2.0), 0, [(mp.Y, 1.5), [(mp.X, 4.0), 1, [(mp.Y, 0.5), 4, 3]], 2]]
    ),
    "cell_size": mp.Vector3(10.0, 5.0, 0),
    "time_stepping": [1.0, 1.0, 1.0, 1.0, 1.0],
    "new_chunk_layout": mp.BinaryPartition(
        data=[(mp.X, -2.0), 0, [(mp.Y, 1.5), [(mp.X, 4.0), 1, [(mp.Y, 0.5), 4, 3]], 2]]
    ),
}
TEST_CHUNK_DATA_6 = {
    "chunk_layout": mp.BinaryPartition(
        data=[(mp.X, -2.0), 0, [(mp.Y, 1.5), [(mp.X, 4.0), 1, [(mp.Y, 0.5), 4, 3]], 2]]
    ),
    "cell_size": mp.Vector3(10.0, 5.0, 0),
    "time_stepping": [1500.0, 2400.0, 700.0, 100.0, 300.0],
    "new_chunk_layout": mp.BinaryPartition(
        data=[
            (mp.X, -3.0),
            0,
            [(mp.Y, 1.25), [(mp.X, -0.3333333333333335), 1, [(mp.Y, -0.625), 4, 3]], 2],
        ]
    ),
}


class MockSimulationTest(unittest.TestCase):
    @parameterized.parameterized.expand(
        [
            (TEST_SIM_1, [6000.0, 9600.0, 2800.0, 400.0, 1200.0]),
            (TEST_SIM_2, [45000.0, 52500.0, 31500.0, 3000.0, 18000.0]),
            (TEST_SIM_3, [400.0, 1600.0, 400.0]),
            (TEST_SIM_4, [400.0, 600.0, 600.0, 400.0, 400.0]),
        ]
    )
    def test_time_spent_on(self, test_sim_constructor, expected_stepping_times):
        test_sim = test_sim_constructor()
        test_sim.init_sim()

        for time_sink in TIMING_MEASUREMENT_IDS.values():
            if time_sink == TIMING_MEASUREMENT_IDS["time_stepping"]:
                self.assertListEqual(
                    list(test_sim.time_spent_on(time_sink)), expected_stepping_times
                )
            else:
                self.assertListEqual(
                    list(test_sim.time_spent_on(time_sink)),
                    [0.0] * len(expected_stepping_times),
                )

    @parameterized.parameterized.expand(
        [
            (TEST_SIM_1, [0, 1, 4, 3, 2]),
            (TEST_SIM_2, [0, 1, 4, 3, 2]),
            (TEST_SIM_3, [0, 1, 2]),
            (TEST_SIM_4, [0, 1, 2, 3, 4]),
        ]
    )
    def test_structure_get_chunk_owners(
        self, test_sim_constructor, expected_chunk_owners
    ):
        test_sim = test_sim_constructor()
        test_sim.init_sim()

        chunk_owners = test_sim.structure.get_chunk_owners()

        self.assertListEqual(list(chunk_owners), expected_chunk_owners)


class ChunkBalancerTest(unittest.TestCase):
    @parameterized.parameterized.expand(
        [
            (TEST_SIM_1, False),
            (TEST_SIM_DUPLICATE_PROC_ID, True),
        ]
    )
    def test_validate_sim(self, test_sim_constructor, should_raise_exception):
        test_sim = test_sim_constructor()
        test_sim.init_sim()

        chunk_balancer = ChunkBalancer()

        if should_raise_exception:
            with self.assertRaises(ValueError):
                chunk_balancer._validate_sim(test_sim)
        else:
            chunk_balancer._validate_sim(test_sim)

    @parameterized.parameterized.expand(
        [
            (TEST_SIM_1,),
            (TEST_SIM_2,),
            (TEST_SIM_3,),
            (TEST_SIM_4,),
        ]
    )
    def test_chunk_layout_improvement(self, test_sim_constructor):
        """Tests that chunk_balancer improves balance after 1 iteration."""
        test_sim = test_sim_constructor()
        test_sim.init_sim()

        old_timing_measurements = MeepTimingMeasurements.new_from_simulation(
            test_sim, -1
        )

        chunk_balancer = ChunkBalancer()

        chunk_balancer.adjust_chunk_layout(test_sim, sensitivity=1.0)

        new_timing_measurements = MeepTimingMeasurements.new_from_simulation(
            test_sim, -1
        )

        old_step_times = np.array(old_timing_measurements.measurements["time_stepping"])
        new_step_times = np.array(new_timing_measurements.measurements["time_stepping"])

        old_max_time = np.max(old_step_times)
        new_max_time = np.max(new_step_times)

        old_min_time = np.min(old_step_times)
        new_min_time = np.min(new_step_times)

        self.assertLess(new_max_time, old_max_time)
        self.assertGreater(new_min_time, old_min_time)

    @parameterized.parameterized.expand(
        [
            (TEST_SIM_1,),
            (TEST_SIM_2,),
            (TEST_SIM_3,),
            (TEST_SIM_4,),
        ]
    )
    def test_chunk_layout_convergence(self, test_sim_constructor):
        """Tests that chunk_balancer converges to load balanced state."""
        test_sim = test_sim_constructor()
        test_sim.init_sim()

        chunk_balancer = ChunkBalancer()

        num_iterations = 25

        for _ in range(num_iterations):
            chunk_balancer.adjust_chunk_layout(test_sim, sensitivity=0.5)

            new_timing_measurements = MeepTimingMeasurements.new_from_simulation(
                test_sim, -1
            )

            new_step_times = np.array(
                new_timing_measurements.measurements["time_stepping"]
            )

        # Check that new stepping times have converged to close to the mean value
        tolerance = 0.05
        mean_step_time = np.mean(new_step_times)
        self.assertTrue(np.allclose(mean_step_time, new_step_times, rtol=tolerance))

    @parameterized.parameterized.expand(
        [
            (TEST_CHUNK_DATA_1,),
            (TEST_CHUNK_DATA_2,),
            (TEST_CHUNK_DATA_3,),
            (TEST_CHUNK_DATA_4,),
            (TEST_CHUNK_DATA_5,),
            (TEST_CHUNK_DATA_6,),
        ]
    )
    def test_split_pos_adjustment(self, test_chunk_data):
        chunk_layout = test_chunk_data["chunk_layout"]
        sim = MockSimulation(
            cell_size=test_chunk_data["cell_size"],
            resolution=10,
            chunk_layout=chunk_layout,
        )
        sim.init_sim()
        chunk_volumes = sim.structure.get_chunk_volumes()
        chunk_owners = sim.structure.get_chunk_owners()

        chunk_balancer = ChunkBalancer()

        measurements = {}
        num_procs = len(test_chunk_data["time_stepping"])
        for name in TIMING_MEASUREMENT_IDS.keys():
            if name == "time_stepping":
                measurements[name] = test_chunk_data["time_stepping"]
            else:
                measurements[name] = [0] * num_procs
        timing_measurements = MeepTimingMeasurements(
            measurements, -1, None, None, None, None, None
        )

        new_chunk_layout = chunk_balancer.compute_new_chunk_layout(
            timing_measurements,
            chunk_layout,
            chunk_volumes,
            chunk_owners,
            sensitivity=1.0,
        )
        expected_chunk_layout = test_chunk_data["new_chunk_layout"]

        self.assertTrue(
            bpu.partitions_are_equal(new_chunk_layout, expected_chunk_layout)
        )


if __name__ == "__main__":
    unittest.main()
