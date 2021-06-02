---
# Adjoint Solver
---

Meep contains a density-based adjoint solver for efficiently computing the gradient of an arbitrary function of the mode coefficients (S-parameters) or DFT fields with respect to $\varepsilon$ on a discrete spatial grid in a subregion of the cell (i.e., a [`MaterialGrid`](../Python_User_Interface.md#materialgrid) class object) over a broad bandwidth. Regardless of the number of degrees of freedom for the grid points, just **two** separate timestepping runs are required. The first run is the "forward" calculation to compute the objective function. The second run is the "adjoint" calculation to compute the gradient of the objective function with respect to the design variables which involves a special type of source distribution and postprocessing of the results (i.e., the DFT fields within the design region). This module is implemented entirely in Python using [autograd](https://github.com/HIPS/autograd). Additionally, the adjoint solver supports [topology optimization](https://en.wikipedia.org/wiki/Topology_optimization) by providing the functionality to wrap an optimization library around the gradient computation.

There are six Jupyter notebooks that demonstrate the main features of the adjoint solver.

- [Introduction](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/01-Introduction.ipynb)

- [Waveguide Bend Optimization](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/02-Waveguide_Bend.ipynb)

- [Filtering and Thresholding Design Parameters and Broadband Objective Functions](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/03-Filtered_Waveguide_Bend.ipynb)

- [Design of a Symmetric Broadband Splitter](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/04-Splitter.ipynb)

- [Broadband Objective Function using Epigraph Formulation](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/05-Minimax.ipynb)

- [Objective Function based on Near to Far-Field Transformation](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/Near2Far-Optimization-with-Epigraph-Formulation.ipynb)

More documentation will be available soon.
