import abc
import copy
from typing import Optional, Tuple, Union

import numpy as np
from meep.timing_measurements import MeepTimingMeasurements

import meep as mp
from meep import binary_partition_utils as bpu


class AbstractChunkBalancer(abc.ABC):
    """Abstract chunk balancer for adaptive chunk layouts in Meep simulations.

    This class defines interfaces for a chunk balancer, which adjusts chunk
    layouts to optimal load balancing. It provides two main functionalities:
      1. Generating an initial chunk layout for the first iteration of an
         optimization run, using a strategy other than the default chunk
         partitioning routine in Meep
      2. Adaptively modifying chunk layouts for subsequent iterations in an
         optimization run by incorporating timing data from the previous
         iteration(s) and resizing the chunks accordingly

    Subclasses of this class can be passed as an option into optimization runs.
    """

    def make_initial_chunk_layout(
        self, sim: mp.Simulation
    ) -> Union[mp.BinaryPartition, None]:
        """Generates an initial chunk layout based on expected compute costs.

        Args:
          sim: the meep.Simulation object to generate a chunk layout for

        Returns:
          The chunk layout to be used by the simulation, or None, in which case
          Meep will use its own default logic to compute the chunk layout (this is
          the default behavior and can be overridden in child classes).
        """
        del sim
        return None

    def adjust_chunk_layout(self, sim: mp.Simulation, **kwargs) -> None:
        """Computes a new chunk layout, applies it to sim, and resets/re-inits sim.

        This method also calls self.should_adjust_chunk_layout(sim). If the current
        chunk layout is sufficiently well-balanced, no changes will be made.
        NOTE: This is a state-changing method which may reset the simulation.

        Args:
          sim: the meep.Simulation object with the chunk layout to be adjusted
          **kwargs: extra args to be passed to self.compute_new_chunk_layout

        Raises:
          ValueError: if sim.chunk_layout includes nodes with duplicate proc_ids
          ValueError: if sim.chunk_layout has proc_ids not included in
            sim.structure.get_chunk_owners() (this could occur if number of
            processes exceeds the number of physical cores)
        """
        self._validate_sim(sim)

        if self.should_adjust_chunk_layout(sim):

            old_chunk_layout = sim.chunk_layout
            chunk_volumes = sim.structure.get_chunk_volumes()
            chunk_owners = sim.structure.get_chunk_owners()

            timing_measurements = MeepTimingMeasurements.new_from_simulation(sim, -1)

            new_chunk_layout = self.compute_new_chunk_layout(
                timing_measurements,
                old_chunk_layout,
                chunk_volumes,
                chunk_owners,
                **kwargs,
            )

            sim.reset_meep()
            sim.chunk_layout = new_chunk_layout
            sim.init_sim()

    @abc.abstractmethod
    def compute_new_chunk_layout(
        self,
        timing_measurements: MeepTimingMeasurements,
        old_partition: mp.BinaryPartition,
        chunk_volumes: Tuple[mp.grid_volume],
        chunk_owners: np.ndarray,
    ) -> mp.BinaryPartition:
        """Rebalances the partition to equalize simulation time for node's children.

        Args:
          timing_measurements: elapsed time by category from previous iteration
          old_partition: the chunk_layout from the previous iteration
          chunk_volumes: obtained from sim.structure.get_chunk_volumes()
          chunk_owners: obtained from sim.structure.get_chunk_owners()

        Returns:
          The chunk layout to be used by the next iteration, with the chunks
          resized appropriately to balance load across MPI processes.
        """

    def should_adjust_chunk_layout(self, sim: mp.Simulation) -> bool:
        """Is the current layout imbalanced enough to justify rebuilding the sim?

        Args:
          sim: the meep.Simulation object with the chunk layout to be adjusted

        Returns:
          A bool specifying whether the current chunk layout is sufficiently poorly
          load-balanced to justify the upfront cost of rebuilding the Meep
          simulation object to redefine the chunks. By default, True is returned if
          this method is not overridden in a subclass.
        """
        del sim
        return True

    def _validate_sim(self, sim: mp.Simulation):
        """Ensures that chunk layout and chunk volumes are properly formatted.

        Args:
          sim: the meep.Simulation object to validate

        Raises:
          ValueError: if sim.chunk_layout includes nodes with duplicate proc_ids
          ValueError: if sim.chunk_layout has proc_ids not included in
            sim.structure.get_chunk_owners() (this could occur if number of
            processes exceeds the number of physical cores)
        """
        # Check that chunk_layout does not contain duplicate nodes for same process
        if bpu.partition_has_duplicate_proc_ids(sim.chunk_layout):
            raise ValueError("Duplicate proc_ids found in chunk_layout!")
        # Check that proc_ids in chunk_layout are assigned properly to grid owners
        chunk_owners = sim.structure.get_chunk_owners()
        proc_ids = [node.proc_id for node in bpu.enumerate_leaf_nodes(sim.chunk_layout)]
        if set(chunk_owners) != set(proc_ids):
            raise ValueError(
                f"Processes {set(proc_ids) - set(chunk_owners)} present in chunk_layout but not grid_owners!"
            )


class ChunkBalancer(AbstractChunkBalancer):
    """A chunk balancer for adaptive chunk layouts in Meep simulations.

    This class generates initial chunk layouts using Meep's built-in scheme, and
    rebalances chunks using timing data from the previous iteration in an
    optimization run to resize chunks appropriately.

    The general idea behind this chunk balancer implementation is to compute
    'load' terms defined by the ratio of (simulation time t) / (chunk volume V).
    This implicitly should account for additional background load on the various
    MPI nodes, since if one machine is more busy than another, it will have a
    higher ratio of t/V.

    The split_pos for each node is adjusted to equalize the per-core values of
    left and right simulation times, then the method is recursively called for
    each of the node's children. A sensitivity parameter determines how much the
    split_pos should change from the previous value, with 0.0 being no change
    and 1.0 being instant change to the new computed value.
    """

    def compute_new_chunk_layout(
        self,
        timing_measurements: MeepTimingMeasurements,
        partition: mp.BinaryPartition,
        chunk_volumes: Tuple[mp.grid_volume],
        chunk_owners: np.ndarray,
        new_partition: Optional[mp.BinaryPartition] = None,
        new_partition_root: Optional[mp.BinaryPartition] = None,
        xyz_bounds: Optional[Tuple[float, float, float, float, float, float]] = None,
        sensitivity: float = 0.5,
    ) -> mp.BinaryPartition:
        """Rebalances the partition to equalize simulation time for node's children.

        Args:
          timing_measurements: elapsed time by category from previous iteration
          partition: the chunk_layout from the previous iteration
          chunk_volumes: obtained from sim.structure.get_chunk_volumes()
          chunk_owners: obtained from sim.structure.get_chunk_owners()
          new_partition: reference to the new chunk_layout object
          new_partition_root: reference to the root node of the new chunk_layout
          xyz_bounds: min and max xyz values for the new sub-partition
          sensitivity: adjusts how much the new split_pos should change

        Returns:
          The chunk layout to be used by the next iteration, with the chunks
          resized appropriately to balance sim times across MPI processes.
        """

        # if we're at root, make a copy of the partition to return as new partition
        if new_partition is None:
            new_partition = copy.deepcopy(partition)
            new_partition_root = new_partition

        # if at leaf node, no more rebalancing to do, so return the partition
        if bpu.is_leaf_node(partition):
            return partition

        if xyz_bounds is None:
            xyz_bounds = bpu.get_box_ranges(partition, chunk_volumes, chunk_owners)

        # Retrieve the relevant dimensions d_min, d_max along the split axis
        # NOTE: we call the distances d_min, d_max regardless of split direction
        if partition.split_dir == mp.X:
            d_min, d_max, _, _, _, _ = xyz_bounds
        elif partition.split_dir == mp.Y:
            _, _, d_min, d_max, _, _ = xyz_bounds
        elif partition.split_dir == mp.Z:
            _, _, _, _, d_min, d_max = xyz_bounds
        else:
            raise ValueError("Split direction must be mp.X, mp.Y, or mp.Z!")

        sim_times = list(self._compute_working_times_per_process(timing_measurements))

        n_left = partition.left.numchunks()
        n_right = partition.right.numchunks()

        t_left, t_right = bpu.get_left_right_total_weights(partition, sim_times)
        v_left, v_right = bpu.get_left_right_total_volumes(
            partition, chunk_volumes, chunk_owners
        )

        v_left_new = v_left * t_right * n_left
        v_right_new = v_right * t_left * n_right
        split_frac = v_left_new / (v_left_new + v_right_new)

        old_split_pos = partition.split_pos

        new_split_pos = (d_max - d_min) * split_frac + d_min
        new_split_pos = sensitivity * new_split_pos + (1 - sensitivity) * old_split_pos

        # Adjust the split pos on new_partition. We can't modify the old partition
        # because it could affect left/right volume ratios for subsequent calls
        new_partition.split_pos = new_split_pos

        x_min, x_max, y_min, y_max, z_min, z_max = xyz_bounds
        # Update the box ranges for the new partition after the adjustment
        if partition.split_dir == mp.X:
            left_xyz_bounds = (x_min, new_split_pos, y_min, y_max, z_min, z_max)
            right_xyz_bounds = (new_split_pos, x_max, y_min, y_max, z_min, z_max)
        elif partition.split_dir == mp.Y:
            left_xyz_bounds = (x_min, x_max, y_min, new_split_pos, z_min, z_max)
            right_xyz_bounds = (x_min, x_max, new_split_pos, y_max, z_min, z_max)
        elif partition.split_dir == mp.Z:
            left_xyz_bounds = (x_min, x_max, y_min, y_max, z_min, new_split_pos)
            right_xyz_bounds = (x_min, x_max, y_min, y_max, new_split_pos, z_max)

        self.compute_new_chunk_layout(
            timing_measurements,
            partition.left,
            chunk_volumes,
            chunk_owners,
            new_partition=new_partition.left,
            new_partition_root=new_partition_root,
            xyz_bounds=left_xyz_bounds,
            sensitivity=sensitivity,
        )
        self.compute_new_chunk_layout(
            timing_measurements,
            partition.right,
            chunk_volumes,
            chunk_owners,
            new_partition=new_partition.right,
            new_partition_root=new_partition_root,
            xyz_bounds=right_xyz_bounds,
            sensitivity=sensitivity,
        )

        return new_partition

    def _compute_working_times_per_process(
        self, timing_measurements: MeepTimingMeasurements
    ) -> np.ndarray:
        """Computes the time spent by each MPI process actively working."""

        time_sinks_to_include = [
            "time_stepping",
            "boundaries_copying",
            "field_output",
            "fourier_transform",
            "mpb",
            "near_to_farfield_transform",
        ]

        num_processes = len(timing_measurements.measurements["time_stepping"])
        working_times = np.zeros(num_processes)
        for category in time_sinks_to_include:
            working_times += np.array(timing_measurements.measurements[category])

        return working_times
