import copy
import unittest

import meep.binary_partition_utils as bpu
import numpy as np
import parameterized

import meep as mp

PARTITION_NO_DUPLICATE_PROC_ID = mp.BinaryPartition(
    data=[(mp.X, -2.0), 0, [(mp.Y, 1.5), [(mp.X, 4.0), 1, [(mp.Y, 0.5), 4, 3]], 2]]
)
PARTITION_DUPLICATE_PROC_ID = mp.BinaryPartition(
    data=[(mp.X, -2.0), 0, [(mp.Y, 1.5), [(mp.X, 4.0), 1, [(mp.Y, 0.5), 1, 3]], 0]]
)

# Mocked chunk IDs since we are in a single-core environment
CHUNK_OWNERS_NO_DUPLICATE_PROC_ID = np.array([0, 1, 4, 3, 2])
CHUNK_OWNERS_DUPLICATE_PROC_ID = np.array([0, 1, 1, 3, 0])


class BinaryPartitionUtilsTest(unittest.TestCase):
    @parameterized.parameterized.expand(
        [
            (PARTITION_NO_DUPLICATE_PROC_ID, False),
            (PARTITION_NO_DUPLICATE_PROC_ID.right, False),
            (PARTITION_NO_DUPLICATE_PROC_ID.left, True),
            (PARTITION_DUPLICATE_PROC_ID, False),
            (PARTITION_DUPLICATE_PROC_ID.right, False),
            (PARTITION_DUPLICATE_PROC_ID.left, True),
        ]
    )
    def test_is_leaf_node(self, partition, expected_leaf_status):
        self.assertEqual(bpu.is_leaf_node(partition), expected_leaf_status)

    @parameterized.parameterized.expand(
        [
            (PARTITION_NO_DUPLICATE_PROC_ID, CHUNK_OWNERS_NO_DUPLICATE_PROC_ID),
            (PARTITION_DUPLICATE_PROC_ID, CHUNK_OWNERS_DUPLICATE_PROC_ID),
        ]
    )
    def test_enumerate_leaf_nodes(self, partition, chunk_owners):
        leaf_nodes = list(bpu.enumerate_leaf_nodes(partition))
        self.assertEqual(len(leaf_nodes), partition.numchunks())
        proc_ids = [node.proc_id for node in bpu.enumerate_leaf_nodes(partition)]
        self.assertEqual(proc_ids, list(chunk_owners))  # depth first ordering

    def test_partition_has_duplicate_proc_ids(self):
        self.assertFalse(
            bpu.partition_has_duplicate_proc_ids(PARTITION_NO_DUPLICATE_PROC_ID)
        )
        self.assertTrue(
            bpu.partition_has_duplicate_proc_ids(PARTITION_DUPLICATE_PROC_ID)
        )

    @parameterized.parameterized.expand(
        [
            (
                PARTITION_NO_DUPLICATE_PROC_ID,
                [0, 1, 2, 3, 4],
            ),
            (
                PARTITION_DUPLICATE_PROC_ID,
                [0, 1, 2, 3, 4],
            ),
        ]
    )
    def test_get_total_weight(self, partition, weights_by_proc_id):
        if not bpu.partition_has_duplicate_proc_ids(partition):
            self.assertEqual(
                bpu.get_total_weight(partition, weights_by_proc_id),
                sum(weights_by_proc_id),
            )
            self.assertEqual(
                bpu.get_total_weight(partition.right.left, weights_by_proc_id),
                weights_by_proc_id[1] + weights_by_proc_id[4] + weights_by_proc_id[3],
            )
        else:
            with self.assertRaises(ValueError):
                bpu.get_total_weight(partition, weights_by_proc_id)

    @parameterized.parameterized.expand(
        [
            (
                PARTITION_NO_DUPLICATE_PROC_ID,
                [0, 1, 2, 3, 4],
            ),
            (
                PARTITION_DUPLICATE_PROC_ID,
                [0, 1, 2, 3, 4],
            ),
        ]
    )
    def test_get_left_right_total_weights(self, partition, weights_by_proc_id):
        proc_ids = [node.proc_id for node in bpu.enumerate_leaf_nodes(partition)]
        no_duplicates = len(set(proc_ids)) == len(proc_ids)
        if no_duplicates:
            self.assertEqual(
                bpu.get_left_right_total_weights(partition, weights_by_proc_id),
                (
                    bpu.get_total_weight(partition.left, weights_by_proc_id),
                    bpu.get_total_weight(partition.right, weights_by_proc_id),
                ),
            )
        else:
            with self.assertRaises(ValueError):
                bpu.get_left_right_total_weights(partition, weights_by_proc_id)

    @parameterized.parameterized.expand(
        [
            (
                PARTITION_NO_DUPLICATE_PROC_ID,
                2,
                [1500, 2400, 300, 100, 700],
            ),
            (
                PARTITION_NO_DUPLICATE_PROC_ID,
                3,
                [15000, 24000, 3000, 1000, 7000],
            ),
        ]
    )
    def test_pixel_volume(self, partition, dims, expected_pixel_volumes):
        cell_size = (
            mp.Vector3(10.0, 5.0, 1.0) if dims == 3 else mp.Vector3(10.0, 5.0, 0.0)
        )
        sim = mp.Simulation(cell_size=cell_size, resolution=10, chunk_layout=partition)
        sim.init_sim()
        chunk_volumes = sim.structure.get_chunk_volumes()
        self.assertEqual(
            [bpu.pixel_volume(vol) for vol in chunk_volumes], expected_pixel_volumes
        )

    @parameterized.parameterized.expand(
        [
            (PARTITION_NO_DUPLICATE_PROC_ID, CHUNK_OWNERS_NO_DUPLICATE_PROC_ID, 3500),
            (PARTITION_DUPLICATE_PROC_ID, CHUNK_OWNERS_DUPLICATE_PROC_ID, 3500),
        ]
    )
    def test_get_total_volume_2d(self, partition, chunk_owners, expected_total_volume):
        sim = mp.Simulation(
            cell_size=mp.Vector3(10.0, 5.0, 0.0), resolution=10, chunk_layout=partition
        )
        sim.init_sim()

        chunk_volumes = sim.structure.get_chunk_volumes()
        total_volume = sim.cell_size[0] * sim.cell_size[1] * sim.resolution**2
        self.assertEqual(
            bpu.get_total_volume(partition, chunk_volumes, chunk_owners), total_volume
        )

        if not bpu.partition_has_duplicate_proc_ids(partition):
            self.assertEqual(
                bpu.get_total_volume(partition.right, chunk_volumes, chunk_owners),
                expected_total_volume,
            )

    @parameterized.parameterized.expand(
        [
            (PARTITION_NO_DUPLICATE_PROC_ID, CHUNK_OWNERS_NO_DUPLICATE_PROC_ID, 35000),
            (PARTITION_DUPLICATE_PROC_ID, CHUNK_OWNERS_DUPLICATE_PROC_ID, 35000),
        ]
    )
    def test_get_total_volume_3d(self, partition, chunk_owners, expected_total_volume):
        sim = mp.Simulation(
            cell_size=mp.Vector3(10.0, 5.0, 1.0), resolution=10, chunk_layout=partition
        )
        sim.init_sim()

        chunk_volumes = sim.structure.get_chunk_volumes()
        total_volume = (
            sim.cell_size[0] * sim.cell_size[1] * sim.cell_size[2] * sim.resolution**3
        )
        self.assertEqual(
            bpu.get_total_volume(partition, chunk_volumes, chunk_owners), total_volume
        )

        if not bpu.partition_has_duplicate_proc_ids(partition):
            self.assertEqual(
                bpu.get_total_volume(partition.right, chunk_volumes, chunk_owners),
                expected_total_volume,
            )

    @parameterized.parameterized.expand(
        [
            (
                PARTITION_NO_DUPLICATE_PROC_ID,
                CHUNK_OWNERS_NO_DUPLICATE_PROC_ID,
                2,
            ),
            (
                PARTITION_DUPLICATE_PROC_ID,
                CHUNK_OWNERS_DUPLICATE_PROC_ID,
                2,
            ),
            (
                PARTITION_NO_DUPLICATE_PROC_ID,
                CHUNK_OWNERS_NO_DUPLICATE_PROC_ID,
                3,
            ),
            (
                PARTITION_DUPLICATE_PROC_ID,
                CHUNK_OWNERS_DUPLICATE_PROC_ID,
                3,
            ),
        ]
    )
    def test_get_left_right_total_volumes(self, partition, chunk_owners, dims):
        cell_size = (
            mp.Vector3(10.0, 5.0, 1.0) if dims == 3 else mp.Vector3(10.0, 5.0, 0.0)
        )
        sim = mp.Simulation(cell_size=cell_size, resolution=10, chunk_layout=partition)
        sim.init_sim()

        chunk_volumes = sim.structure.get_chunk_volumes()
        self.assertEqual(
            bpu.get_left_right_total_volumes(partition, chunk_volumes, chunk_owners),
            (
                bpu.get_total_volume(partition.left, chunk_volumes, chunk_owners),
                bpu.get_total_volume(partition.right, chunk_volumes, chunk_owners),
            ),
        )

    @parameterized.parameterized.expand(
        [
            (
                PARTITION_NO_DUPLICATE_PROC_ID,
                CHUNK_OWNERS_NO_DUPLICATE_PROC_ID,
                2,
            ),
            (
                PARTITION_DUPLICATE_PROC_ID,
                CHUNK_OWNERS_DUPLICATE_PROC_ID,
                2,
            ),
            (
                PARTITION_NO_DUPLICATE_PROC_ID,
                CHUNK_OWNERS_NO_DUPLICATE_PROC_ID,
                3,
            ),
            (
                PARTITION_DUPLICATE_PROC_ID,
                CHUNK_OWNERS_DUPLICATE_PROC_ID,
                3,
            ),
        ]
    )
    def test_get_grid_volumes_in_tree(self, partition, chunk_owners, dims):
        cell_size = (
            mp.Vector3(10.0, 5.0, 1.0) if dims == 3 else mp.Vector3(10.0, 5.0, 0.0)
        )
        sim = mp.Simulation(cell_size=cell_size, resolution=10, chunk_layout=partition)
        sim.init_sim()

        chunk_volumes = sim.structure.get_chunk_volumes()
        grid_volumes_in_tree = bpu.get_grid_volumes_in_tree(
            partition, chunk_volumes, chunk_owners
        )
        self.assertEqual(set(grid_volumes_in_tree), set(chunk_volumes))

        no_duplicates = len(set(chunk_owners)) == len(chunk_owners)
        grid_volumes_in_right_expected = (
            chunk_volumes[1:] if no_duplicates else chunk_volumes
        )
        grid_volumes_in_right = bpu.get_grid_volumes_in_tree(
            partition.right, chunk_volumes, chunk_owners
        )
        self.assertEqual(
            set(grid_volumes_in_right), set(grid_volumes_in_right_expected)
        )

    @parameterized.parameterized.expand(
        [
            (
                PARTITION_NO_DUPLICATE_PROC_ID,
                CHUNK_OWNERS_NO_DUPLICATE_PROC_ID,
                2,
                {0: 1500, 1: 2400, 4: 300, 3: 100, 2: 700},
            ),
            (
                PARTITION_DUPLICATE_PROC_ID,
                CHUNK_OWNERS_DUPLICATE_PROC_ID,
                2,
                {
                    0: 2200,
                    1: 2700,
                    3: 100,
                },
            ),
            (
                PARTITION_NO_DUPLICATE_PROC_ID,
                CHUNK_OWNERS_NO_DUPLICATE_PROC_ID,
                3,
                {0: 15000, 1: 24000, 4: 3000, 3: 1000, 2: 7000},
            ),
            (
                PARTITION_DUPLICATE_PROC_ID,
                CHUNK_OWNERS_DUPLICATE_PROC_ID,
                3,
                {
                    0: 22000,
                    1: 27000,
                    3: 1000,
                },
            ),
        ]
    )
    def test_get_total_volume_per_process(
        self, partition, chunk_owners, dims, expected_volumes_per_process
    ):
        cell_size = (
            mp.Vector3(10.0, 5.0, 1.0) if dims == 3 else mp.Vector3(10.0, 5.0, 0.0)
        )
        sim = mp.Simulation(cell_size=cell_size, resolution=10, chunk_layout=partition)
        sim.init_sim()
        chunk_volumes = sim.structure.get_chunk_volumes()
        volumes_per_process = bpu.get_total_volume_per_process(
            partition, chunk_volumes, chunk_owners
        )
        self.assertEqual(volumes_per_process, expected_volumes_per_process)

    @parameterized.parameterized.expand(
        [
            (
                PARTITION_NO_DUPLICATE_PROC_ID,
                CHUNK_OWNERS_NO_DUPLICATE_PROC_ID,
                2,
                (-5.0, 5.0, -2.5, 2.5, 0.0, 0.0),
                (-2.0, 5.0, -2.5, 2.5, 0.0, 0.0),
            ),
            (
                PARTITION_DUPLICATE_PROC_ID,
                CHUNK_OWNERS_DUPLICATE_PROC_ID,
                2,
                (-5.0, 5.0, -2.5, 2.5, 0.0, 0.0),
                (-5.0, 5.0, -2.5, 2.5, 0.0, 0.0),
            ),
            (
                PARTITION_NO_DUPLICATE_PROC_ID,
                CHUNK_OWNERS_NO_DUPLICATE_PROC_ID,
                3,
                (-5.0, 5.0, -2.5, 2.5, -0.5, 0.5),
                (-2.0, 5.0, -2.5, 2.5, -0.5, 0.5),
            ),
            (
                PARTITION_DUPLICATE_PROC_ID,
                CHUNK_OWNERS_DUPLICATE_PROC_ID,
                3,
                (-5.0, 5.0, -2.5, 2.5, -0.5, 0.5),
                (-5.0, 5.0, -2.5, 2.5, -0.5, 0.5),
            ),
        ]
    )
    def test_get_box_ranges(
        self,
        partition,
        chunk_owners,
        dims,
        expected_box_ranges,
        expected_right_box_ranges,
    ):
        cell_size = (
            mp.Vector3(10.0, 5.0, 1.0) if dims == 3 else mp.Vector3(10.0, 5.0, 0.0)
        )
        sim = mp.Simulation(cell_size=cell_size, resolution=10, chunk_layout=partition)
        sim.init_sim()
        chunk_volumes = sim.structure.get_chunk_volumes()
        self.assertEqual(
            bpu.get_box_ranges(partition, chunk_volumes, chunk_owners),
            expected_box_ranges,
        )
        self.assertEqual(
            bpu.get_box_ranges(partition.right, chunk_volumes, chunk_owners),
            expected_right_box_ranges,
        )

    @parameterized.parameterized.expand(
        [
            (
                copy.deepcopy(PARTITION_NO_DUPLICATE_PROC_ID),
                copy.deepcopy(PARTITION_NO_DUPLICATE_PROC_ID),
                True,
            ),
            (
                copy.deepcopy(PARTITION_DUPLICATE_PROC_ID),
                copy.deepcopy(PARTITION_DUPLICATE_PROC_ID),
                True,
            ),
            (
                copy.deepcopy(PARTITION_NO_DUPLICATE_PROC_ID),
                copy.deepcopy(PARTITION_DUPLICATE_PROC_ID),
                False,
            ),
            (
                mp.BinaryPartition(data=[(mp.X, -2.5), 0, [(mp.X, 2.5), 1, 2]]),
                mp.BinaryPartition(data=[(mp.X, -2.5), 0, [(mp.X, 2.5), 1, 2]]),
                True,
            ),
            (
                mp.BinaryPartition(data=[(mp.X, -2.5), 0, [(mp.X, 2.4), 1, 2]]),
                mp.BinaryPartition(data=[(mp.X, -2.5), 0, [(mp.X, 2.5), 1, 2]]),
                False,
            ),
        ]
    )
    def test_partitions_are_equal(self, bp1, bp2, is_equal):
        self.assertEqual(bpu.partitions_are_equal(bp1, bp2), is_equal)


if __name__ == "__main__":
    unittest.main()
