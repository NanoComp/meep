# Tutorials and Examples

Meep simulations are Python scripts which involve specifying the device geometry, materials, current sources, monitor fields, and everything else necessary to set up a calculation. A Python script provides the flexibility to customize the simulation for practically any application particularly those involving parameter sweeps and optimization.

Python libraries such as NumPy, SciPy, and Matplotlib can be used to augment the simulation functionality and will also be demonstrated. Much of the functionality of the low-level C++ interface has been abstracted in Python which means that you don't need to be an experienced programmer to set up simulations. Reasonable defaults are available where necessary.

Several tutorials and examples are found here. These tutorials are meant to illustrate Meep's various features in an interactive and application-oriented manner. 

## iPython/Jupyter Notebooks

Jupyter notebooks are interactive, browser based framework for displaying text, running python code, and visualizing results. There are several ways to read these notebooks online:

### Local

The recommended method to run the tutorial notebooks is by 1.) installing `meep` via `conda`, 2.) cloning this repo to your local drive, and 3.) launch a notebook server in this directory using `jupyter notebook`.

### nbviewer

`nbviewer` is a web platform that can render interactive features found inside Jupyter notebooks openly stored on web-servers. While `nbviewer` can't run python code, it can execute stored javascript code used to animate the simulations. 

### GitHub

GitHub is able to render some of the smaller notebooks as plain text. However, they are not interactive and are often too large for GitHub.

### Tutorials

Below are summaries for each tutorial, along with the features the tutorials highlight. While there is no particular order to the tutorials, they progressively incorporate more complicated features.

#### Basics 

1. __`straight-waveguide.ipynb`__ -
A simple 2D straight waveguide tutorial that explores basic meep features like `geometry`, `sources`, and `PML` layers. The tutorial also explores basic visualization and animation features.

2. __`bent-waveguide.ipynb`__-
A followup to the 2D straight waveguide tutorial by adding a bend.

3. __`bend-flux.ipynb`__-
Using the previous bent waveguide example, this tutorial calculates the loss, transmission, and reflection that the bent waveguide undergoes.

4. __`ring.ipynb`__ -
Computes the resonant mode frequencies of a 2D ring resonator using `harminv`.

5. __`visualization.ipynb`__ -
Demonstrates various visualization and animation features.

##### Ring Resonator in Cylindrical Coordinates

##### Band Diagram, Resonant Modes, and Transmission of a Waveguide Cavity

##### Material Dispersion

##### Third Harmonic Generation

##### Near to Far Field Spectra

##### Local Density of States

##### Optical Forces

##### Gyrotropic Media

##### Multilevel-Atomic Susceptibility

##### Frequency Domain Solver

##### Eigenmode Source

##### Mode Decomposition

##### GDSII Import

##### Adjoint Solver
