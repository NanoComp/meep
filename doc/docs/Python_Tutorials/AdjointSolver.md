---
# Adjoint Solver
---

Meep contains an adjoint solver for efficiently computing the gradient of an objective function with respect to the permittivity on a discrete spatial grid in a subregion of the cell. Regardless of the number of degrees of freedom for the grid points, just **two** separate timestepping runs are required. The first run is the "forward" calculation to compute the objective function. The second run is the "adjoint" calculation which involves a special type of source distribution and postprocessing applied to the results. This module is implemented entirely in Python and does not involve modifications to the C++ [`libmeep`](../Chunks_and_Symmetry.md) core library. At a higher level, the module implements functionality for wrapping a numerical optimizer around the gradient computation to enable automated design optimization.

Three Jupyter notebooks that explain how to use the new adjoint solver can be found here:
* [Introduction](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/01-Introduction.ipynb)
* [Waveguide Bend Optimization](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/02-Waveguide_Bend.ipynb)
* [Filtering and Thresholding Design Parameters and Broadband Objective Functions](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/03-Filtered_Waveguide_Bend.ipynb)

More documentation will come soon.
