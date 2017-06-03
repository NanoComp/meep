---
title: Perfectly matched layer
permalink: /Perfectly_matched_layer/
---

The **perfectly matched layer** (**PML**) approach to implementing absorbing boundary conditions in FDTD codes was proposed by Berenger in 1994 (see [1](http://dx.doi.org/10.1006/jcph.1994.1159)). The approach involves surrounding the computational cell with a medium that in theory absorbs without any reflection electromagnetic waves at all frequencies and angles of incidence. Berenger showed that it was sufficient to "split" Maxwell's equations into two sets of (unphysical) equations in the absorbing layers, appropriately defined. It was later shown that a similar reflectionless absorbing medium can be constructed as a lossy anisotropic dielectric and magnetic material with "matched" impedance and electrical and magnetic conductivities.

The finite-difference implementation of PML requires the conductivities to be turned on gradually over a distance of a few grid points to avoid numerical reflections from the discontinuity. It is important when using PMLs to make the computational cell sufficiently large otherwise the promixity of the PML to evanescent modes may unnecessarily extract energy from the system and thereby perturb it.

There are two useful references on PML in Meep:

-   [Notes on Perfectly Matched Layers](http://math.mit.edu/~stevenj/18.369/pml.pdf) by S. G. Johnson: a general introduction to PML concepts
-   [Notes on the UPML implementation in Meep](http://ab-initio.mit.edu/meep/pml-meep.pdf) by S. G. Johnson: a description of the precise PML formulation that is used in Meep, which is slightly different from many PML formulations described elsewhere in order to properly handle arbitrary anisotropy; a closely related formulation (with slightly different notation) is [described in our 2011 paper](http://math.mit.edu/~stevenj/papers/OskooiJo11.pdf).

PML has some limitations. See also our [recent publication](http://www.opticsinfobase.org/abstract.cfm?URI=oe-16-15-11376) on the use of adiabatic absorbers as a workaround for cases when PML fails (most notably in periodic media like photonic crystals, where the fundamental idea behind PML breaks down), another [paper on a different case where the PML idea fails](http://math.mit.edu/~stevenj/papers/LohOs09.pdf), and a [third paper](http://math.mit.edu/~stevenj/papers/OskooiJo11.pdf) on validation of PML and errors in some proposals.

[Category:Meep](Meep.md)
