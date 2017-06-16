---
# Perfectly Matched Layer
---

The **perfectly matched layer** (**PML**) approach to implementing absorbing boundary conditions in FDTD codes was proposed by [Berenger in 1994](http://dx.doi.org/10.1006/jcph.1994.1159). The approach involves surrounding the computational cell with a medium that in theory absorbs without any reflection electromagnetic waves at all frequencies and angles of incidence. Berenger showed that it was sufficient to "split" Maxwell's equations into two sets of equations in the absorbing layers, appropriately defined. These split-field equations produce wave attenuation but are unphysical. It was later shown that a similar reflectionless absorbing medium can be constructed as a lossy anisotropic dielectric and magnetic material with "matched" impedance and electrical and magnetic conductivities. This is known as the uniaxial PML (UPML).

The finite-difference implementation of PML requires the conductivities to be turned on gradually over a distance of a few grid points to avoid numerical reflections from the discontinuity. It is important when using PMLs to make the computational cell sufficiently large otherwise the promixity of the PML to evanescent modes may unnecessarily extract energy from the system and thereby perturb it.

There are two useful references on PML in Meep:

-   [Notes on Perfectly Matched Layers](http://math.mit.edu/~stevenj/18.369/pml.pdf) by S. G. Johnson: a general introduction to PML concepts
-   [Notes on the UPML implementation in Meep](http://ab-initio.mit.edu/meep/pml-meep.pdf) by S. G. Johnson: a description of the precise PML formulation that is used in Meep, which is slightly different from many PML formulations described elsewhere in order to properly handle arbitrary anisotropy; a closely related formulation with slightly different notation is described in [J. Computational Physics, vol. 230, pp. 2369-77 (2011)](http://math.mit.edu/~stevenj/papers/OskooiJo11.pdf). This paper also describes a strategy to validate PML.

PML has some limitations. See also [Optics Express, vol. 16, pp. 11376-92 (2008)](http://www.opticsinfobase.org/abstract.cfm?URI=oe-16-15-11376) on the use of adiabatic absorbers as a workaround for cases when PML fails. This occurs most notably in periodic media like photonic crystals, where the fundamental idea behind PML breaks down. [Physical Review E, vol. 79, 065601 (2011)](http://math.mit.edu/~stevenj/papers/LohOs09.pdf) provides a different case where PML fails.
