---
# C++ Developer Information
---

Meep Design
-----------

An overview of Meep's inner workings is summarized in [Computer Physics Communications, Vol. 181, pp. 687-702, 2010](http://ab-initio.mit.edu/~oskooi/papers/Oskooi10.pdf). This page is a supplement which provides a description of the source code.

[TOC]

### Data Structures and Chunks

Meep employs several data structures declared in `meep.hpp`. The principal data structure element is the **chunk**. A chunk is a contiguous rectangular portion of the computational grid. For example, when Meep runs on a parallel system, each process gets one or more disjoint chunks of the grid. There are several different types of chunks:

-   `fields` and `fields_chunks`
-   `structure` and `structure_chunks`
-   `dft` and `dft_chunks`

As an example, the `fields` class encapsulates the simulated fields over the entire grid, and one of its members is an array of `fields_chunk` variables that the grid is divided into. The `fields_chunk` variables are where the field is actually stored. Every parallel process has a nearly-identical fields variable with a nearly-identical list of chunks. Chunks on one process which have been assigned to another process do not store their fields arrays, and are just placeholders.

If a given material or field is not present in a given chunk, it need not be stored. For this reason, the PML boundary regions are separated into their own chunks, even on one processor, in order that the extra data for PML need not be stored for the whole grid.

In the future, we may implement support for different chunks with different resolution, to allow nonuniform spatial resolution.

Similarly for `structure` and `structure_chunks`, except that it is only for materials parameters such as epsilon, etc. and not for simulated fields.

`dft_chunk` stores accumulated Fourier-transformed fields corresponding to a given chunk.

### `grid_volume` and `volume`

The `volume` class declared in `meep/vec.hpp` represents a rectangular box-like region, parallel to the $XYZ$ axes, in "continuous space" &mdash; i.e. the corners can be at any points, not necessarily grid points. This is used, for example, whenever you want to specify the integral of some quantity (e.g., flux, energy) in a box-like region, and Meep interpolates from the grid as necessary to give an illusion of continuity.

The `grid_volume` class declared in `meep/vec.hpp` is a box of pixels. It stores the resolution, the number of pixels in each direction, the origin, etcetera. Given a `grid_volume`, there are functions to get the `volume` corresponding to the bounding box, etcetera. There is a `grid_volume` object associated with the whole computational grid, and with each chunk in the grid. There are various tricky aspects to the `grid_volume`. One is associated with the Yee grid: it has to know about different field components stored at different points. Another is associated with the fact that boundary conditions, not only the overall grid boundaries but also boundaries between chunks, are handled by an extra layer of "not-owned" pixels around the boundaries. So each chunk's `grid_volume` has "owned" grid points that the chunk is responsible for updating, and "not-owned" grid points that are updated using the boundary conditions. And thanks to the Yee grid which complicates everything in FDTD, unfortunately, the set of owned and not-owned coordinates is different for each field component. The `grid_volume` class keeps track of all this.

### File Organization

The core Meep C++ simulation code (all of the physics) is located in the `src/` directory, with C++
tests in the `tests/` directory.  The `libmeepgeom/` directory provides a C++ library to specify Meep
geometries in terms of a list of geometric objects (spheres, cylinders, boxes) with various
material properties (via [libctl](https://github.com/stevengj/libctl)'s geometry library), and is
also used by the Python interface.

The Scheme and Python interfaces are found in the `scheme/` and `python/` directories, respectively;
both of them use [SWIG](http://www.swig.org/) to generate wrapper code from the C++ header files,
but also have hand-written Scheme/Python code to provide a higher-level interface.  The `libpympb/`
directory contains a Python interface to MPB, which will eventually be moved to the MPB repository.

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
| Geometric objects                                                                                              | handled by [libctl](https://github.com/stevengj/libctl) functions in libctl's geom.c, called from the Scheme front-end (not handled by Meep) |
| Fields: initialization, cleanup, chunking, stepping-plan, (dis)affiliation with sources, polarizabilities etc. | fields.cpp                                                                                                                   |
| Structure: initialization, cleanup, chunking, material parameters, boundary conditions etc.                    | structure.cpp                                                                                                                |
| MPI interface                                                                                                  | meep/mympi.hpp, mympi.cpp                                                                                                    |

### Deprecated Interfaces

Beware that some of the interfaces in the source code and in the old manual are now deprecated, as they have been superseded by newer features and may be removed at some point.

In particular, you should probably avoid:

-   The `monitor_point` class. Just declare an array to store the fields you want, get them with `fields::get_field`, and analyze them with `do_harminv`. Or, to accumulate the DFT as you run, use the `dft_chunk` class via `fields::add_dft`.
-   Slice and EPS output. This has been superseded by HDF5 output, which is much more flexible and efficient.

---
# Python Developer Information

## Overview

The `meep` Python package consists of two types of files.

1. Files created by SWIG [SWIG](http://www.swig.org/), which makes up the "low-level" interface. 
2. Files that make up the high-level Python interface.

Graphic for python source generation

## Package Organization

Description of how the `meep` package is constructed

## File Organization

| File              | Description                                                                                                                           |
|-------------------|---------------------------------------------------------------------------------------------------------------------------------------|
| geom.py           | Python geometric objects that are analogues to the geometric objects defined in `libctl/utils/geom.c`                                 |
| mpb_data.py | Definition of `MPBData`, a Python class that contains the functionality of the `mpb-data` command line program. |
| simulation.py     | Definition of the `Simulation` class and step functions associated with running a simulation.                                         |
| solver.py | Classes and functions related to the Python interface to `MPB`.
| source.py         | Definitions of Python objects relating to sources.                                                                                    |
| meep.i            | SWIG interface file for the `meep` Python module.                                                                                     |
| mpb.i             | SWIG interface file for the `meep.mpb` Python module.                                                                                 |
| numpy.i           | Typemaps for `numpy` arrays (taken from the numpy [Github repository](https://github.com/numpy/numpy/blob/master/tools/swig/numpy.i)) |
| vec.i             | SWIG interface file for `vec.hpp`. Iincluded into `meep.i`.                                                                           |
| typemap_utils.cpp | Utility functions for writing SWIG typemaps. Included into `meep.i`.                                                                  |
