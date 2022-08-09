import warnings
from typing import Dict, Generator, List, Tuple

import numpy as onp

import meep as mp


def is_leaf_node(partition: mp.BinaryPartition) -> bool:
    """Returns True if the partition has no children.

    Args:
      partition: the BinaryPartition node

    Returns:
      A boolean indicating whether partition is a leaf node.
    """
    return partition.left is None and partition.right is None


def enumerate_leaf_nodes(
    partition: mp.BinaryPartition,
) -> Generator[mp.BinaryPartition, None, None]:
    """Enumerates all leaf nodes of a partition.

    Args:
      partition: the BinaryPartition node

    Yields:
      The leaf nodes contained within the partition.
    """
    if is_leaf_node(partition):
        yield partition
    else:
        yield from enumerate_leaf_nodes(partition.left)
        yield from enumerate_leaf_nodes(partition.right)


def partition_has_duplicate_proc_ids(partition: mp.BinaryPartition) -> bool:
    """Returns True if the partition has more than one node with the same proc_id.

    Args:
      partition: the BinaryPartition node

    Returns:
      A boolean indicating if the partition contains duplicate proc_ids.
    """
    proc_ids = [node.proc_id for node in enumerate_leaf_nodes(partition)]
    return len(set(proc_ids)) != len(proc_ids)


def get_total_weight(
    partition: mp.BinaryPartition, weights_by_proc_id: List[float]
) -> float:
    """Computes the total weights contained within a BinaryPartition subtree.

    Args:
      partition: a BinaryPartition subtree to compute the associated weights for
      weights_by_proc_id: a list of weights associated with each proc_id

    Returns:
      The sum of all weights for each proc_id encountered in the subtree.

    Raises:
      ValueError: if sim.chunk_layout includes nodes with duplicate proc_ids
    """
    if partition_has_duplicate_proc_ids(partition):
        raise ValueError("Duplicate proc_ids found in chunk_layout!")
    if partition.proc_id is not None:
        return weights_by_proc_id[partition.proc_id]
    elif partition.left is not None and partition.right is not None:
        left_weight = get_total_weight(partition.left, weights_by_proc_id)
        right_weight = get_total_weight(partition.right, weights_by_proc_id)
        return left_weight + right_weight
    else:
        raise ValueError("Partition missing proc_id or left, right attributes!")


def get_left_right_total_weights(
    partition: mp.BinaryPartition, weights_by_proc_id: List[float]
) -> Tuple[float, float]:
    """Computes the total weights contained in left and right subtrees.

    Args:
      partition: a BinaryPartition tree to compute the associated weights for
      weights_by_proc_id: a list of weights associated with each proc_id

    Returns:
      The sum of weights for each proc_id encountered in the left and right
      subtrees.

    Raises:
      ValueError: if the BinaryPartition is a leaf node or improperly formatted.
    """
    if partition.left is not None and partition.right is not None:
        left_weight = get_total_weight(partition.left, weights_by_proc_id)
        right_weight = get_total_weight(partition.right, weights_by_proc_id)
        return (left_weight, right_weight)
    else:
        raise ValueError("Partition missing left, right attributes!")


def pixel_volume(grid_volume: mp.grid_volume) -> int:
    """Computes the number of pixels contained in a 2D or 3D grid_volume object.

    NOTE: This assumes that z=0 means 2D simulation and z>0 means 3D.

    Args:
      grid_volume: a meep grid_volume object

    Returns:
      The 2D or 3D number of pixels in the grid_volume.
    """
    if grid_volume.nz() > 0:  # we're in a 3D simulation
        return grid_volume.nx() * grid_volume.ny() * grid_volume.nz()
    else:  # 2D simulation
        return grid_volume.nx() * grid_volume.ny()


def get_total_volume(
    partition: mp.BinaryPartition,
    chunk_volumes: Tuple[mp.grid_volume],
    chunk_owners: onp.ndarray,
) -> float:
    """Computes the total pixel volume in a subtree from associated chunk volumes.

    NOTE: If multiple chunks are owned by the same process, this function may
    over-count the volume, since all provided grid volumes associated with a
    given process are counted.

    Args:
      partition: a BinaryPartition subtree to compute the associated volumes for
      chunk_volumes: associated grid volumes from a simulation, obtained by
        calling sim.structure.get_chunk_volumes()
      chunk_owners: list of processes associated with each grid volume, obtained
        by calling sim.structure.get_chunk_owners()

    Returns:
      The total pixel volume occupied by all chunks owned by the partition.
    """
    my_grid_volumes = get_grid_volumes_in_tree(partition, chunk_volumes, chunk_owners)
    return sum(pixel_volume(vol) for vol in my_grid_volumes)


def get_left_right_total_volumes(
    partition: mp.BinaryPartition,
    chunk_volumes: Tuple[mp.grid_volume],
    chunk_owners: onp.ndarray,
) -> Tuple[float, float]:
    """Computes the total pixel volume in left and right subtrees.

    Args:
      partition: a BinaryPartition subtree to compute the associated volumes for
      chunk_volumes: associated grid volumes from a simulation, obtained by
        calling sim.structure.get_chunk_volumes()
      chunk_owners: list of processes associated with each grid volume, obtained
        by calling sim.structure.get_chunk_owners()

    Returns:
      A tuple for left and right subtreees, where each entry is a list of the
      total pixel volume occupied by all chunks owned by each process.

    Raises:
      ValueError: if the BinaryPartition is a leaf node or improperly formatted.
    """
    if partition.left is not None and partition.right is not None:
        left_volume = get_total_volume(partition.left, chunk_volumes, chunk_owners)
        right_volume = get_total_volume(partition.right, chunk_volumes, chunk_owners)
        return (left_volume, right_volume)
    else:
        raise ValueError("Partition missing left, right attributes!")


def get_grid_volumes_in_tree(
    partition: mp.BinaryPartition,
    chunk_volumes: Tuple[mp.grid_volume],
    chunk_owners: onp.ndarray,
) -> List[mp.grid_volume]:
    """Fetches a list of grid_volumes contained in a BinaryPartition subtree.

    NOTE: If multiple chunks are owned by the same process, this function may
    over-count the volume, since all provided grid volumes associated with a
    given process are counted.

    Args:
      partition: a BinaryPartition subtree to find grid_volumes for
      chunk_volumes: associated grid volumes from a simulation, obtained by
        calling sim.structure.get_chunk_volumes()
      chunk_owners: list of processes associated with each grid volume, obtained
        by calling sim.structure.get_chunk_owners()

    Returns:
      A list of all grid_volume objects associated with any proc_id contained in
      the partition. The list is not necessarily ordered by proc_id values.
    """
    if partition_has_duplicate_proc_ids(partition):
        warnings.warn("Partition has duplicate proc_ids, overcounting possible!")

    my_proc_ids = [node.proc_id for node in enumerate_leaf_nodes(partition)]

    return [
        chunk_volumes[i]
        for (i, owner) in enumerate(chunk_owners)
        if owner in my_proc_ids
    ]


def get_total_volume_per_process(
    partition: mp.BinaryPartition,
    chunk_volumes: Tuple[mp.grid_volume],
    chunk_owners: onp.ndarray,
) -> Dict[int, float]:
    """Computes the total pixel volume per process contained in a BinaryPartition.

    Args:
      partition: a BinaryPartition subtree to compute the associated volumes for
      chunk_volumes: associated grid volumes from a simulation, obtained by
        calling sim.structure.get_chunk_volumes()
      chunk_owners: list of processes associated with each grid volume, obtained
        by calling sim.structure.get_chunk_owners()

    Returns:
      A dictionary of the total pixel volume occupied by all chunks owned by each
      process.
    """
    volumes_per_process = {}
    leaf_nodes_in_tree = enumerate_leaf_nodes(partition)
    for leaf in leaf_nodes_in_tree:
        total_volume = sum(
            pixel_volume(chunk_volumes[i])
            for i, owner in enumerate(chunk_owners)
            if owner == leaf.proc_id
        )

        volumes_per_process[leaf.proc_id] = total_volume
    return volumes_per_process


def get_box_ranges(
    partition: mp.BinaryPartition,
    chunk_volumes: Tuple[mp.grid_volume],
    chunk_owners: onp.ndarray,
) -> Tuple[float, float, float, float, float, float]:
    """Gets the max and min x, y, z dimensions spanned by a partition.

    Args:
      partition: a BinaryPartition subtree to compute the range for
      chunk_volumes: associated grid volumes from a simulation, obtained by
        calling sim.structure.get_chunk_volumes()
      chunk_owners: list of processes associated with each grid volume, obtained
        by calling sim.structure.get_chunk_owners()

    Returns:
      A 6-tuple which enumerates:
      (
        xmin: the min x value of any grid_volume associated with the partition,
        xmax: the max x value of any grid_volume associated with the partition,
        ymin: the min y value of any grid_volume associated with the partition,
        ymax: the max y value of any grid_volume associated with the partition,
        zmin: the min z value of any grid_volume associated with the partition,
        zmax: the max z value of any grid_volume associated with the partition
      )
    """
    xmins = []
    xmaxs = []
    ymins = []
    ymaxs = []
    zmins = []
    zmaxs = []
    for leaf in enumerate_leaf_nodes(partition):
        # generate a list of all volumes owned by the process
        owned_volumes = [
            chunk_volumes[i]
            for (i, owner) in enumerate(chunk_owners)
            if owner == leaf.proc_id
        ]
        # add the corners to the lists
        for vol in owned_volumes:
            xmins.append(vol.surroundings().get_min_corner().x())
            ymins.append(vol.surroundings().get_min_corner().y())
            zmins.append(vol.surroundings().get_min_corner().z())
            xmaxs.append(vol.surroundings().get_max_corner().x())
            ymaxs.append(vol.surroundings().get_max_corner().y())
            zmaxs.append(vol.surroundings().get_max_corner().z())
    return (min(xmins), max(xmaxs), min(ymins), max(ymaxs), min(zmins), max(zmaxs))


def partitions_are_equal(bp1: mp.BinaryPartition, bp2: mp.BinaryPartition) -> bool:
    """Determines if two partitions have all nodes with identical attributes.

    Args:
      bp1: a BinaryPartition object to compare equality for
      bp2: the other BinaryPartition object to compare equality for

    Returns:
      A bool if all nodes in the partitions have equal attributes
    """
    if is_leaf_node(bp1) and is_leaf_node(bp2):
        return bp1.proc_id == bp2.proc_id
    elif (not is_leaf_node(bp1)) and (not is_leaf_node(bp2)):
        return all(
            [
                bp1.split_dir == bp2.split_dir,
                bp1.split_pos == bp2.split_pos,
                partitions_are_equal(bp1.left, bp2.left),
                partitions_are_equal(bp1.right, bp2.right),
            ]
        )
    else:
        return False
