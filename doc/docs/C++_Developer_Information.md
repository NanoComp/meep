---
# C++ Developer Information
---

An overview of Meep's inner workings is provided in [Computer Physics Communications, Vol. 181, pp. 687-702, 2010](http://ab-initio.mit.edu/~oskooi/papers/Oskooi10.pdf). This page is a supplement which provides a description of the source code.

For additional details, see [Chunks and Symmetry](Chunks_and_Symmetry.md)

[TOC]

### Data Structures and Chunks

Meep's data structures are defined in `meep.hpp`. The principal data structure element is the **chunk**. A chunk is a contiguous rectangular portion of the computational grid. For example, when Meep runs on a parallel system, each process gets one or more disjoint chunks of the grid.

There are several different types of chunks:

-   `fields` and `fields_chunks`
-   `structure` and `structure_chunks`
-   `dft` and `dft_chunks`

As an example, the `fields` class encapsulates the fields over the entire grid, and one of its members is an array of `fields_chunk` variables that divides the grid. The `fields_chunk` variables store the actual field information. Every parallel process has a nearly-identical fields variable with a nearly-identical list of chunks. Chunks on one process which have been assigned to another process do not store their fields arrays; they are just placeholders.

If a given material or field is not present in a given chunk, it need not be stored. For this reason, the PML boundary regions are separated into their own chunks, even on one processor, in order that the extra data for PML need not be stored for the whole grid.

In the future, we may implement support for different chunks with different resolution, to allow nonuniform spatial resolution.

Similarly for `structure` and `structure_chunks`, except that it is only for materials parameters such as epsilon, etc. and not for the fields.

`dft_chunk` stores accumulated Fourier-transformed fields corresponding to a given chunk.

### `grid_volume` and `volume`

The `volume` class declared in `meep/vec.hpp` represents a rectilinear region, parallel to the $xyz$ axes, in "continuous space" &mdash; i.e. the corners can be at any points, not necessarily grid points. This is used, for example, whenever you want to specify the integral of some quantity (e.g., flux, energy) in a box-like region, and Meep interpolates from the grid as necessary to give an illusion of continuity.

The `grid_volume` class declared in `meep/vec.hpp` is a box of pixels. It stores the resolution, the number of pixels in each direction, the origin, etcetera. Given a `grid_volume`, there are functions to get the `volume` corresponding to the bounding box, etcetera. There is a `grid_volume` object associated with the whole computational grid, and with each chunk in the grid. There are various tricky aspects to the `grid_volume`. One is associated with the Yee grid: it has to know about different field components stored at different points. Another is associated with the fact that boundary conditions, not only the overall grid boundaries but also boundaries between chunks, are handled by an extra layer of "not-owned" pixels around the boundaries. So each chunk's `grid_volume` has "owned" grid points that the chunk is responsible for updating, and "not-owned" grid points that are updated using the boundary conditions. Due to the Yee grid which complicates everything in FDTD, unfortunately, the set of owned and not-owned coordinates is different for each field component. The `grid_volume` class keeps track of all this.

### File Organization

The core Meep C++ simulation code (all of the physics) is located in the `src/` directory, with C++
tests in the `tests/` directory.  The module `src/meepgeom.cpp` provides a C++ interface to specify Meep
geometries in terms of a list of geometric objects (spheres, cylinders, boxes) with various
material properties (via [libctl](https://libctl.readthedocs.io)'s geometry library), and is
also used by the Python interface.

The Scheme and Python interfaces are found in the `scheme/` and `python/` directories. Both interfaces use [SWIG](http://www.swig.org/) to generate wrapper code from the C++ header files,
but also have hand-written Scheme/Python code to provide a higher-level interface.  The `libpympb/`
directory contains a Python interface to MPB (which may, in the future, be moved to the MPB repository).

The following table briefly describes the purpose of some of the source files:

| Header File         | Description                                                                                                                                                                                                                                                                                                           |
|---------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| meep/vec.hpp        | Declares geometry-related classes like vec, ivec, grid_volume, volume and related utility functions.                                                                                                                                                                                                                 |
| meep/mympi.hpp      | Declares functions for initializing the meep application, cleanup, and data exchange accounting for the presence or absence of MPI. These functions present a unified interface to the rest of the application.                                                                                                       |
| meep.hpp            | All public classes likes fields, fields_chunks, structure, structure_chunks, src_time, continuous_src_time, material_function, h5_file, polarizability_identifier etc.                                                                                                                                        |
| meep_internals.hpp | Hosts declarations for classes like polarizability, polarization, src_vol, and bandsdata. Also defines macros for frequently-used loop constructs like DOCMP that are internal to Meep implementation.                                                                                                               |
| bicgstab.hpp        | Declares functions related to an implementation of an iterative solver for non-symmetric linear operators based on a generalization of the stabilized biconjugate-gradient (BiCGSTAB) algorithm proposed by van der Vorst (and described in the book "Templates for the Solution of Linear Systems" by Barrett et al. |

The following table briefly describes what is in each .cpp file:

| Source File      | Description                                                                                                |
|------------------|------------------------------------------------------------------------------------------------------------|
| polarization.cpp | Implement member functions for the polarization and polarizability classes declared in meep_internals.hpp |
| bicgstab.cpp     | Implements the solver described against bicgstab.hpp (see above)                                           |

#### Functionality Organization

| Functionality                                                                                                  | Location                                                                                                                     |
|----------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------|
| Material dispersion                                                                                            | polarization.cpp, update_from_e.cpp, and friends.                                                                          |
| Vectors, volumes etc.                                                                                          | meep/vec.hpp, vec.cpp                                                                                                        |
| Geometric objects                                                                                              | handled by [libctl](https://github.com/NanoComp/libctl) functions in libctl's geom.c, called from the Scheme front-end (not handled by Meep) |
| Fields: initialization, cleanup, chunking, stepping-plan, (dis)affiliation with sources, polarizabilities etc. | fields.cpp                                                                                                                   |
| Structure: initialization, cleanup, chunking, material parameters, boundary conditions etc.                    | structure.cpp                                                                                                                |
| MPI interface                                                                                                  | meep/mympi.hpp, mympi.cpp                                                                                                    |

### Deprecated Interfaces

Beware that some of the interfaces in the source code and in the old manual are now deprecated, as they have been superseded by newer features and may be removed at some point.

In particular, you should probably avoid:

-   The `monitor_point` class. Just declare an array to store the fields you want, get them with `fields::get_field`, and analyze them with `do_harminv`. Or, to accumulate the DFT as you run, use the `dft_chunk` class via `fields::add_dft`.
-   Slice and EPS output. This has been superseded by HDF5 output, which is much more flexible and efficient.

---
