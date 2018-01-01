---
# Units and Nonlinearity
---

For linear calculations in electromagnetism, most quantities of interest are naturally expressed as dimensionless quantities, such as the ratio of the wavelength to a given lengthscale, the transmitted or reflected power as a fraction of input power, or the lifetime in units of the optical period. Matters are more complicated when one includes nonlinear effects, however, because in this case the absolute amplitude of the electric field becomes significant. On this page, we discuss how to relate Meep's units to those of experimental quantities relevant for nonlinear problems. See also [Nonlinearities](Materials.md#nonlinearity).

Kerr Nonlinearities
-------------------

Meep supports instantaneous Kerr nonlinearities characterized by a susceptibility $\chi^{(3)}$, corresponding to a constitutive relation (in Meep's units):

$$\mathbf{D} = \left( \varepsilon + \chi^{(3)} \cdot |\mathbf{E}|^2 \right) \mathbf{E}$$

However, the number usually reported for the strength of the Kerr nonlinearity is the AC Kerr coefficient $n_2$, defined by the effective change in refractive index $\Delta n$ for a planewave with time-average intensity $I$ travelling through a homogeneous Kerr material:

$$\Delta n = n_2 I$$

This equation itself is somewhat subtle: it is not the actual instantaneous change in refractive index at every point. Rather, it is a sort of average change in index, and in particular is the change in effective index $\beta c/\omega$ where $\beta$ is the propagation constant.

The relationship between $n_2$ and $\chi^{(3)}$, in Meep's units, is:

$$n_2 = \frac{3\chi^{(3)}}{4n_0^2}$$

where $n_0$ is the linear refractive index $n_0=\sqrt{\varepsilon}$.   (See, for example, section 4.1 of [Nonlinear Optics (3rd edition)](https://www.amazon.com/Nonlinear-Optics-Third-Robert-Boyd/dp/0123694701) by R. W. Boyd.)

**Warning:** The optics literature uses a variety of conflicting conventions for defining the dimensionful quantities $\chi^{(3)}$ and $n_2$, which may differ from the ones used here by scale factors of $n_0$ etcetera.  For example, the Boyd book mentioned above also describes a different definition of $n_2$ in which $\Delta n = n_2 |\mathbf{E}|^2$ instead of $I$, which changes the relationship between $n_2$ and $\chi^{(3)}$ by a factor of $n_0$.  If you are transcribing experimental values of $n_2$ into Meep, you may need to convert from one convention to another!

### Using Experimental Kerr Coefficients

The key to using correct magnitudes is nonlinear calculations in Meep is to realize that the units are still somewhat arbitrary: only the *product* $n_2 I$ is significant, so as long as we get this product right we are fine. Moreover, we can choose our units of distance and our units of field strength independently and arbitrarily. So, if we are given $n_2$ in some units like &#956;m<sup>2</sup>/W we are free to use &#956;m as our unit of distance and W as our units of power.

For example, suppose now that we are given some experimental value of $n_2$ in "real" units, say $3\times10^{–8}$&#956;m<sup>2</sup>/W (silica glass). Of course, this is very small (semiconductors can be much more nonlinear), so to compensate suppose let's plan to use an unrealistic 1 MW of power in our structure (say a waveguide). To implement this in Meep, for convenience we'll use &#956;m as our unit of distance and W as our units of power. First, we'll set $\chi^{(3)}$ from $n_2$ and *n* in these units: $\chi^{(3)} = 4n^2 (3\times 10^{-8})/3$ where *n* is the linear index. Then we simply monitor the power going through our waveguide (or whatever) and change the current amplitude (note that power &#8764; $J^2$) until we get $10^6$.

To monitor the power in a structure, we can use a variety of functions. See the [User Interface](Scheme_User_Interface.md). One can get the power flux directly through the `meep-fields-flux-in-box` function. Or, one can alternatively get the intensity at a single point by calling `get-field-point` and passing `Sx` etc. for the component. One thing to be cautious about is that these return the power or intensity at one instant in time, and not the time-average unless you use complex-valued fields (which are problematic for nonlinear systems).

You may have been hoping for a simple formula: set the current to *x* to get *y* power. However, this is not feasible since the amount of power or field intensity you get from a current source depends on the source geometry, the dielectric structure, and so on. And a formula for the units of current is not terribly useful because usually the current source in an FDTD calculation is artifically inserted to create the field, and doesn't correspond to the current source in a physical experiment.
