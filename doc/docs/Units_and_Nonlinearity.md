---
# Units and Nonlinearity
---

For linear calculations in electromagnetism, most quantities of interest are naturally expressed as dimensionless quantities, such as the ratio of the wavelength to a given lengthscale, the transmitted or reflected power as a fraction of input power, or the lifetime in units of the optical period. Matters are more complicated when one includes nonlinear effects, however, because in this case the absolute amplitude of the electric field becomes significant. This page describes how to relate Meep's units to those of experimental quantities relevant for nonlinear problems. See also [Nonlinearities](Materials.md#nonlinearity).

Kerr Nonlinearities
-------------------

The refractive index often depends on the intensity of the electric fields. There are two conventions, briefly summarized below, which are commonly used to define these third-order (or Kerr) nonlinearities. Since Meep implements only one of these two conventions, some care is required when inputting values from the literature.

The first convention defines the effective refractive index in terms of the time-averaged (and not the instantaneous) electric fields within vacuum:

$$n =n _0 + n_2 \langle \tilde{E}(t)^2 \rangle ,$$

where $n_0$ is the zeroth-order refractive index and $n_2$ is the second-order refractive index which gives the rate at which the effective refractive index varies with the field intensity. For time-harmonic fields $\tilde{E}(t)=E(\omega)e^{-i\omega t}$, this is equivalent to $n =n _0 + 2 n_2 \vert E(\omega) \vert ^2$. Based on this convention, the relationship between $n_2$ and the Kerr susceptibility $\chi^{(3)}$ is:

$$n_2 = \frac{3\chi^{(3)}}{4n_0}.$$

Meep, however, supports the second convention which defines the effective index in terms of the time-averaged electric fields **within the Kerr medium**:

$$n =n _0 + n_2 I$$

where, in Meep units, $I = 2n_0\vert E(\omega) \vert^2$. Based on this convention, the relationship between $n_2$ and $\chi^{(3)}$ is:

$$n_2 = \frac{3\chi^{(3)}}{4n_0^2}.$$

For reference, these expressions are derived in Section 4.1 of [Nonlinear Optics (3rd edition)](https://www.amazon.com/Nonlinear-Optics-Third-Robert-Boyd/dp/0123694701) by R.W. Boyd.

### Using Experimental Kerr Coefficients

The key to using correct magnitudes in nonlinear calculations in Meep is to realize that the units are still somewhat arbitrary: only the product $n_2 I$ is significant, so as long as we get this product right we are fine. Moreover, we can choose our units of distance and our units of field strength independently and arbitrarily. So, if we are given $n_2$ in some units like &#956;m<sup>2</sup>/W we are free to use &#956;m as our unit of distance and W as our unit of power.

For example, suppose now that we are given some experimental value of $n_2$ in "real" units, say $3\times10^{–8}$&#956;m<sup>2</sup>/W (silica glass). Of course, this is very small (semiconductors can be much more nonlinear), so to compensate suppose let's plan to use an unrealistic 1 MW of power in our structure (say a waveguide). To implement this in Meep, for convenience we'll use &#956;m as our unit of distance and W as our units of power. First, we'll set $\chi^{(3)}$ from $n_2$ and $n_0$ in these units: $\chi^{(3)} = 4n_0^2 (3\times 10^{-8})/3$. Then we simply monitor the power going through our waveguide (or whatever) and change the current amplitude (note that power &#8764; $J^2$) until we get $10^6$.

To monitor the power in a structure, we can use a variety of functions. See the [User Interface](Python_User_Interface.md). One can get the power flux directly through the `flux_in_box` function. Or, one can alternatively get the intensity at a single point by calling `get_field_point` and passing `Sx` etc. for the component. One thing to be cautious about is that these return the instantaneous power or intensity, and not the time average unless you use complex-valued fields which are problematic for nonlinear systems. A better approach for obtaining the time-averaged power is via the Fourier-transformed fields using the `add_flux` routine. See the [Tutorial/Third Harmonic Generation](Python_Tutorials/Third_Harmonic_Generation/) for an example.

You may have been hoping for a simple formula: set the current to *x* to get *y* power. However, this is not feasible since the amount of power or field intensity you get from a current source [depends on the source geometry, the dielectric structure, and so on](FAQ/#how-does-the-current-amplitude-relate-to-the-resulting-field-amplitude). And a formula for the units of current is not particularly useful because usually the current source in an FDTD calculation is artifically inserted to create the field, and doesn't correspond to the current source in a physical experiment.
