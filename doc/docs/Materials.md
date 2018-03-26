---
# Materials
---

The material structure in Maxwell's equations is determined by the relative permittivity $\varepsilon$(**r**) and permeability $\mu$(**r**).

However, $\varepsilon$ is not only a function of position. In general, it also depends on frequency (*material dispersion*) and on the electric field **E** itself (*nonlinearity*). It may also depend on the orientation of the field (*anisotropy*). Material dispersion, in turn, is generally associated with absorption loss in the material, or possibly *gain*. All of these effects can be simulated in Meep, with certain restrictions.

Similarly for the relative permeability $\mu$(**r**), for which dispersion, nonlinearity, and anisotropy are all supported as well.

In this section, we describe the form of the equations and material properties that Meep can simulate. The actual user interface where these properties are specified is described in the [User Interface](Python_User_Interface.md). 

[TOC]

Material Dispersion
-------------------

Physically, material dispersion arises because the polarization of the material does not respond instantaneously to an applied field **E**, and this is essentially the way that it is implemented in FDTD. In particular, $\mathbf{D} = \varepsilon\mathbf{E}$ is expanded to:

$$\mathbf{D} = \varepsilon_\infty \mathbf{E} + \mathbf{P}$$

where $\varepsilon_\infty$, which [must be positive](FAQ/#why-does-my-simulation-diverge-if-0), is the *instantaneous* dielectric function (the infinite-frequency response) and **P** is the remaining frequency-dependent *polarization* density in the material. **P**, in turn, has its own time-evolution equation, and the exact form of this equation determines the frequency-dependence $\varepsilon$($\omega$). **Note:** Meep's definition of $\omega$ uses a sign convention $\exp(-i\omega t)$ for the time dependence &mdash; $\varepsilon$ formulas in engineering papers that use the opposite sign convention for $\omega$ will have a sign flip in all the imaginary terms below. If you are using parameters from the literature, you should use **positive** values of $\gamma$ and $\omega$ as-is for loss; don't be confused by the difference in $\omega$ sign convention and flip the sign of the parameters. In particular, Meep supports a Lorentzian susceptibility profile which consists of a sum of harmonic resonances plus a term for the frequency-independent electric conductivity:

<center>

$$\varepsilon(\omega,\mathbf{x}) = \left( 1 + \frac{i \cdot \sigma_D(\mathbf{x})}{\omega}  \right) \left[ \varepsilon_\infty(\mathbf{x})  + \sum_n \frac{\sigma_n(\mathbf{x}) \cdot \omega_n^2 }{\omega_n^2 - \omega^2 - i\omega\gamma_n} \right] ,$$

$= \left( 1 + \frac{i \cdot \sigma_D(\mathbf{x})}{2\pi f}  \right) \left[ \varepsilon_\infty(\mathbf{x})  + \sum_n \frac{\sigma_n(\mathbf{x}) \cdot f_n^2 }{f_n^2 - f^2 - if\gamma_n/2\pi} \right] ,$

</center>

where $\sigma_D$ is the electric conductivity, $\omega_n$ and $\gamma_n$ are user-specified constants. Actually, the numbers that one specifies are $f_n = \omega_n / 2\pi$ and $\gamma_n / 2\pi$. The $\sigma_n(\mathbf{x})$ is a user-specified function of position giving the strength of the *n*-th resonance. The $\sigma$ parameters can be anisotropic (real-symmetric) tensors, while the frequency-independent term $\varepsilon_\infty$ can be an arbitrary real-symmetric tensor as well. This corresponds to evolving **P** via the equations:

$$\mathbf{P} = \sum_n \mathbf{P}_n$$

$$\frac{d^2\mathbf{P}_n}{dt^2} + \gamma_n \frac{d\mathbf{P}_n}{dt} +  \omega_n^2 \mathbf{P}_n = \sigma_n(\mathbf{x}) \omega_n^2 \mathbf{E}$$

That is, we must store and evolve a set of auxiliary fields $\mathbf{P}_n$ along with the electric field in order to keep track of the polarization **P**. Essentially any $\varepsilon(\omega)$ could be modeled by including enough of these polarization fields &mdash; Meep allows you to specify any number of these, limited only by computer memory and time which increases with the number of polarization terms you require.

Note that the conductivity $\sigma_D$ corresponds to an imaginary part of $\varepsilon$ given by $i \varepsilon_\infty \sigma_D / \omega$. This does not include the harmonic-resonance terms. When you specify frequency in Meep units, however, you are specifying *f* without the 2$\pi$, so the imaginary part of $\varepsilon$ is $i \varepsilon_\infty \sigma_D / 2\pi f$.

Meep also supports polarizations of the [Drude](https://en.wikipedia.org/wiki/Drude_model) form, typically used for metals:

$$\frac{d^2\mathbf{P}_n}{dt^2} + \gamma_n \frac{d\mathbf{P}_n}{dt} = \sigma_n(\mathbf{x}) \omega_n^2 \mathbf{E}$$

which corresponds to a term of the following form in $\varepsilon$'s $\Sigma$<sub>*n*</sub>:

$$\frac{i \sigma_n(\mathbf{x}) \cdot \omega_n^2 }{\omega (\gamma_n- i\omega)}$$,

which is equivalent to the Lorentzian model except that the $\omega_n^2$ term has been omitted from the denominator, and asymptotes to a conductivity $\sigma_n \omega_n^2 / \gamma_n$ as $\omega\to 0$. In this case, $\omega_n^2$ is just a dimensional scale factor and has no interpretation as a resonance frequency.

Numerical Stability
-------------------

If a Lorentzian resonance $\omega_n$ is specified at too high frequency relative to the time discretization $\Delta t$, the simulation becomes unstable. Essentially, the problem is that $\mathbf{P}_n$ oscillates too fast compared with the time discretization for the discretization to work properly. If this happens, there are three workarounds: increase the resolution which increases the resolution in both space and time, decrease the Courant factor which decreases $\Delta t$ compared to $\Delta x$, or use a different model function for your dielectric response.

Roughly speaking, the $\mathbf{P}_n$ equation becomes unstable for $\omega_n \Delta t / 2 > 1$. Note that, in Meep frequency units, you specify $f_n = \omega_n/2\pi$, so this quantity should be less than $1/\pi \Delta t$. Meep will check a necessary stability criterion automatically and halt with an error message if it is violated.

Finally, overlapping dispersive materials with perfectly-matched layer (PML) absorbing boundaries may produce instabilities. A workaround is to replace the PML with an absorber.

Loss and Gain
-------------

If $\gamma$ above is nonzero, then the dielectric function $\varepsilon(\omega)$ becomes *complex*, where the imaginary part is associated with absorption loss in the material if it is positive, or gain if it is negative. Alternatively, a dissipation loss or gain may be added by a positive or negative conductivity, respectively &mdash; this is often convenient if you only care about the imaginary part of $\varepsilon$ in a narrow bandwidth, and is described in detail below.

If you look at Maxwell's equations, then $d\mathbf{P}/dt$ plays exactly the same role as a current $\mathbf{J}$. Just as $\mathbf{J} \cdot \mathbf{E}$ is the rate of change of mechanical energy (the power expended by the electric field on moving the currents), therefore, the rate at which energy is lost to absorption is given by:

<center>

absorption rate $\sim \frac{d\mathbf{P}}{dt} \cdot \mathbf{E}$

</center>

Meep can keep track of this energy for the Lorentzian polarizability terms but not for the conductivity terms. For gain, this gives the amount of energy expended in amplifying the field.

Conductivity and Complex ε
--------------------------

Often, you only care about the absorption loss in a narrow bandwidth, where you just want to set the imaginary part of $\varepsilon$ (or $\mu$) to some known experimental value, in the same way that you often just care about setting a dispersionless real $\varepsilon$ that is the correct value in your bandwidth of interest.

One approach to this problem would be allowing you to specify a constant, frequency-independent, imaginary part of $\varepsilon$, but this has the disadvantage of requiring the simulation to employ complex fields which double the memory and time requirements, and also tends to be numerically unstable. Instead, the approach in Meep is for you to set the conductivity $\sigma_D$ (or $\sigma_B$ for an imaginary part of $\mu$), chosen so that $\mathrm{Im}\, \varepsilon = \varepsilon_\infty \sigma_D / \omega$ is the correct value at your frequency $\omega$ of interest. Note that, in Meep, you specify $f = \omega/2\pi$ instead of $\mu$ for the frequency, however, so you need to include the factor of 2$\pi$ when computing the corresponding imaginary part of $\varepsilon$. Conductivities can be implemented with purely real fields, so they are not nearly as expensive as implementing a frequency-independent complex $\varepsilon$ or $\mu$.

For example, suppose you want to simulate a medium with $\varepsilon = 3.4 + 0.101i$ at a frequency 0.42 (in your Meep units), and you only care about the material in a narrow bandwidth around this frequency (i.e. you don't need to simulate the full experimental frequency-dependent permittivity). Then, in Meep, you could use `meep.Medium(epsilon=3.4, D_conductivity=2*math.pi*0.42*0.101/3.4)`; i.e. $\varepsilon_\infty = \mathrm{Re}\,\varepsilon = 3.4$ and $\sigma_D = \omega \mathrm{Im} \varepsilon / \varepsilon_\infty = (2\pi 0.42) 0.101 / 3.4$.

**Note**: the "conductivity" in Meep is slightly different from the conductivity you might find in a textbook, because for computational convenience it appears as $\sigma_D \mathbf{D}$ in our Maxwell equations rather than the more-conventional $\sigma \mathbf{E}$; this just means that our definition is different from the usual electric conductivity by a factor of $\varepsilon$. Also, just as Meep uses the dimensionless relative permittivity for $\varepsilon$, it uses nondimensionalized units of 1/*a* (where *a* is your unit of distance) for the conductivities $\sigma_{D,B}$. If you have the electric conductivity $\sigma$ in SI units and want to convert to $\sigma_D$ in Meep units, you can simply use the formula: $\sigma_D = (a/c) \sigma / \varepsilon_r \varepsilon_0$ where *a* is your unit of distance in meters, *c* is the vacuum speed of light in m/s, $\varepsilon_0$ is the SI vacuum permittivity, and $\varepsilon_r$ is the real relative permittivity.

Nonlinearity
------------

In general, $\varepsilon$ can be changed anisotropically by the **E** field itself, with:

$$\Delta\varepsilon_{ij} = \sum_{k} \chi_{ijk}^{(2)} E_k + \sum_{k\ell} \chi_{ijk\ell}^{(3)} E_k E_\ell + \cdots$$

where the *ij* is the index of the change in the 3$\times$3 $\varepsilon$ tensor and the $\chi$ terms are the nonlinear susceptibilities. The $\chi^{(2)}$ sum is the [Pockels effect](https://en.wikipedia.org/wiki/Pockels_effect) and the $\chi^{(3)}$ sum is the [Kerr effect](https://en.wikipedia.org/wiki/Kerr_effect). If the above expansion is frequency-independent, then the nonlinearity is *instantaneous*; more generally, $\Delta\varepsilon$ would depend on some average of the fields at previous times.

Meep supports instantaneous, isotropic Pockels and Kerr nonlinearities, corresponding to a frequency-independent $\chi_{ijk}^{(2)} = \chi^{(2)} \cdot \delta_{ij} \delta_{jk}$ and $\chi_{ijk\ell}^{(3)} = \chi^{(3)} \cdot \delta_{ij} \delta_{k\ell}$, respectively. Thus,

$$\mathbf{D} = \left( \varepsilon_\infty(\mathbf{x}) + \chi^{(2)}(\mathbf{x})\cdot \mathrm{diag}(\mathbf{E}) + \chi^{(3)}(\mathbf{x}) \cdot |\mathbf{E}|^2 \right) \mathbf{E} + \mathbf{P}$$

Here, "diag(**E**)" indicates the diagonal 3$\times$3 matrix with the components of **E** along the diagonal.

Normally, for nonlinear systems you will want to use real fields **E**. This is the default. However, Meep uses complex fields if you have Bloch-periodic boundary conditions with a non-zero Bloch wavevector **k**, or in cylindrical coordinates with $m \neq 0$. In the C++ interface, real fields must be explicitly specified.

For complex fields in nonlinear systems, the physical interpretration of the above equations is unclear because one cannot simply obtain the physical solution by taking the real part any more. In particular, Meep simply *defines* the meaning of the nonlinearity for complex fields as follows: the real and imaginary parts of the fields do not interact nonlinearly. That is, the above equation should be taken to hold for the real and imaginary parts (of **E** and **D**) *separately*: e.g., |**E**|<sup>2</sup> is the squared magnitude of the *real* part of **E** for when computing the real part of **D**, and conversely for the imaginary part.

For a discussion of how to relate $\chi^{(3)}$ in Meep to experimental Kerr coefficients, see [Units and Nonlinearity](Units_and_Nonlinearity.md).

Magnetic Permeability μ
-----------------------

All of the above features that are supported for the electric permittivity $\varepsilon$ are also supported for the magnetic permeability $\mu$. That is, Meep supports $\mu$ with dispersion from magnetic conductivity and Lorentzian resonances, as well as magnetic $\chi^{(2)}$ and $\chi^{(3)}$ nonlinearities. The description of these is exactly the same as above, so we won't repeat it here &mdash; just take the above descriptions and replace $\varepsilon$, **E**, **D**, and $\sigma$<sub>D</sub> by $\mu$, **H**, **B**, and $\sigma$<sub>B</sub>, respectively.

Materials Library
-----------------

A materials library containing commonly used metals in optoelectronic devices is available for [Python](https://github.com/stevengj/meep/tree/master/python/examples/materials_library.py) and [Scheme](https://github.com/stevengj/meep/tree/master/scheme/examples/materials-library.scm). The data is based on results published in [A.D. Rakic et al., Applied Optics, Vol. 37, No. 22, pp. 5271-83 (1998)](https://www.osapublishing.org/ao/abstract.cfm?uri=ao-37-22-5271) [[pdf](http://faculty.kfupm.edu.sa/EE/msunaidi/EE635%20stuff/project%202/p3.pdf)]. Experimental values of the complex refractive index of 11 metals &mdash; Ag, Au, Cu, Al, Be, Cr, Ni, Pd, Pt, Ti, W &mdash; are fit to a [Drude-Lorentzian susceptibility profile](#material-dispersion) over the broadband spectrum of approximately 0.2 to 12.4 µm. Fitting parameters for the materials are defined for a unit distance of 1 µm. For simulation models which use a *different* value for the unit distance, the predefined variable `eV_um_scale` (Python) or `eV-um-scale` (Scheme) must be scaled by *multiplying* by whatever the unit distance is, in units of µm. For example, if the unit distance is 100 nm, this would require adding the line `eV_um_scale = 0.1*eV_um_scale` after the line where [`eV_um_scale` is defined](https://github.com/stevengj/meep/blob/master/python/examples/materials_library.py#L7). This change must be made directly to the materials library file.


To import the library into a Python script requires adding the following lines:

```python
import sys
sys.path.insert(0, '/path/to/materials_library.py')
from materials_library import *
```
Then, the materials can be simply used as `geometry = [ meep.Cylinder(material=Al, ... ]`.

In Scheme, the required lines are:

```scm
(include "/path/to/materials-library.scm")
```

Note: for narrowband calculations, some of the Lorentzian susceptibility terms may be unnecessary and will contribute to consuming more computational resources than are required (due to the additional storage and time stepping of the polarization fields). Computational efficiency can be improved (without significantly affecting the accuracy of the results) by removing from the material definitions those Lorentzian suspeptibility terms which are far outside the spectral region of interest.
