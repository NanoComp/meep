---
# Adjoint Solver
---

Meep contains an adjoint-solver module for efficiently computing the gradient of an arbitrary function of the mode coefficients (S-parameters), DFT fields, and far fields (using the analytic near-to-far field transformation) with respect to $\varepsilon$ on a discrete spatial grid (a [`MaterialGrid`](../Python_User_Interface.md#materialgrid) class object) at multiple frequencies over a broad bandwidth. Regardless of the number of degrees of freedom for the grid points, just **two** separate timestepping runs are required. The first run is the "forward" calculation to compute the objective function and the DFT fields of the design region. The second run is the "adjoint" calculation to compute the gradient of the objective function with respect to the design variables. The adjoint run involves a special type of source distribution used to compute the DFT fields of the design region. The gradient is computed in post processing using the DFT fields from the forward and adjoint runs. This is all done automatically, and is described in:

- A. M. Hammond, A. Oskooi, M. Chen, Z. Lin, S. G. Johnson, and S. E. Ralph, “[High-performance hybrid time/frequency-domain topology optimization for large-scale photonics inverse design](https://doi.org/10.1364/OE.442074),” *Optics Express*, vol. 30, no. 3, pp. 4467–4491 (2022).

This interface to this functionality is implemented entirely in Python using [autograd](https://github.com/HIPS/autograd) and [JAX](https://github.com/google/jax). The adjoint solver supports inverse design and [topology optimization](https://en.wikipedia.org/wiki/Topology_optimization) by providing the functionality to wrap an optimization library around the gradient computation.

There are six Jupyter notebooks that demonstrate the main features of the adjoint solver.

- [Introduction](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/01-Introduction.ipynb)

- [Waveguide Bend Optimization](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/02-Waveguide_Bend.ipynb)

- [Filtering and Thresholding Design Parameters and Broadband Objective Functions](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/03-Filtered_Waveguide_Bend.ipynb)

- [Design of a Symmetric Broadband Splitter](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/04-Splitter.ipynb)

- [Broadband Objective Function using Epigraph Formulation](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/05-Minimax.ipynb)

- [Objective Function based on Near to Far-Field Transformation](https://nbviewer.jupyter.org/github/NanoComp/meep/blob/master/python/examples/adjoint_optimization/Near2Far-Optimization-with-Epigraph-Formulation.ipynb)

More documentation will be available soon.
