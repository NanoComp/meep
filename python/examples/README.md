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

* __`straight-waveguide.ipynb`__ -
A simple 2D straight waveguide tutorial that explores basic meep features like `geometry`, `sources`, and `PML` layers. The tutorial also explores basic visualization and animation features.

* __`bent-waveguide.ipynb`__-
A followup to the 2D straight waveguide tutorial by adding a bend.

* __`bend-flux.ipynb`__-
Using the previous bent waveguide example, this tutorial calculates the loss, transmission, and reflection that the bent waveguide undergoes.

* __`ring.ipynb`__ -
Computes the resonant mode frequencies of a 2D ring resonator using `harminv`.

#### Ring Resonator in Cylindrical Coordinates

#### Band Diagram, Resonant Modes, and Transmission of a Waveguide Cavity

* __`holey-wg-cavity.ipynb`__ -
Calculates the transmission and resonant modes of a waveguide photonic crystal cavity. Demonstrates the `harminv` routines and how to estimate the $Q$ of cavities.

* __`holey-wg-bands.ipynb`__ - 
Computes the band diagram of the infinite periodic waveguide by itself with no defects in the time domain. Explores the `k_point`, `run_k_point`, and periodic boundary conditions features.

#### MPB and Band diagrams

#### Material Dispersion

* __`refl-quartz.ipynb`__ -

#### Nonlinear Optics

* __`3rd-harm-1d.ipynb`__ - 
Examines 3rd harmonic generation in a $\chi^{(3)}$ material. Explores the proper way to do 1d simulations, how to include nonlinearities in materials, and compares experimental results to theory.

#### Near to Far Field Spectra

We demonstrate Meep's near-to-far field transformation feature using four different examples. There are three steps involved in this type of calculation. First, we need to define the "near" surface(s) as a set of surfaces capturing all outgoing radiation in the desired direction(s). Second, we run the simulation using a pulsed source (or alternatively, a CW source via the frequency-domain solver) to allow Meep to accumulate the Fourier transforms on the near surface(s). Third, we have Meep compute the far fields at any desired points with the option to save the far fields to an HDF5 file.

* __`antenna-radiation.ipynb`__ - 
Computes the radiation pattern of a simple point source "antenna". Explores `add_near2far`, `add_flux`, `get_fluxes`, and `get_farfield` features.

* __`metasurface_lens.ipynb`__ -

* __`binary_grating_n2f.ipynb`__ -

* __`cavity-farfield.ipynb`__ -

#### Local Density of States

#### Optical Forces

#### Gyrotropic Media

#### Multilevel-Atomic Susceptibility

#### Frequency Domain Solver

* __`solve-cw.ipynb`__ -

#### Eigenmode Source

#### Mode Decomposition

#### GDSII Import

#### Adjoint Solver
