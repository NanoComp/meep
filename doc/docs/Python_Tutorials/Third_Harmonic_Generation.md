---
# Third Harmonic Generation
---

In this example, we consider wave propagation through a simple 1d nonlinear medium with a non-zero [Kerr susceptibility $\chi^{(3)}$](https://en.wikipedia.org/wiki/Kerr_effect). See also [Materials](../Materials.md#nonlinearity) and [Units and Nonlinearity](../Units_and_Nonlinearity.md). We send in a narrow-band pulse at a frequency $\omega$, and because of the nonlinearity we also get a signal at a frequency $3\omega$. See also [3rd-harm-1d.py](https://github.com/NanoComp/meep/blob/master/python/examples/3rd-harm-1d.py).

Since this is a 1d calculation, we could implement it via a 2d cell of `Vector3(S,0,0)`, specifying periodic boundary conditions in the $y$ direction. However, this is slightly inefficient since the $y$ periodic boundaries are implemented internally via extra "ghost pixels" in the $y$ direction. Instead, Meep has special support for 1d simulations in the $z$ direction. To use this, we must explicitly set `dimensions` to `1`, and in that case we can *only* use $E_x$ (and $D_x$) and $H_y$ field components. This involves no loss of generality because of the symmetry of the problem.

First, we'll load the necessary modules and define a function which will perform all the necessary computations:

```py
from typing import Tuple, List, Union
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import meep as mp


def third_harmonic_generation(logk: float, amp: float = 1.0, nfreq: int = 10,
                              flux_spectrum: bool = True) -> Union[
                                  Tuple[List[float], List[float]],
                                  Tuple[float, float]]:
```

Next, we'll define some parameters of our simulation:

```py
sz = 100  # size of cell in z direction
fcen = 1 / 3.0  # center frequency of source
df = fcen / 20.0  # frequency width of source
amp = args.amp  # amplitude of source
k = 10**logk  # Kerr susceptibility
dpml = 1.0  # PML thickness
```

Now, to define our cell, we'll do:

```py
dimensions = 1
cell = mp.Vector3(0, 0, sz)
pml_layers = [mp.PML(dpml)]
resolution = 25
```

Note that this will only put PMLs at the $\pm z$ boundaries.

In this case, we're going to fill the entire computational cell with the nonlinear medium, so we don't need to use any objects. We can just use the special `default_material` which is ordinarily vacuum:

```py
default_material = mp.Medium(index=1, chi3=k)
```

Now, our source will be a Gaussian pulse of $J_x$ just next to the $-z$ PML layer. Since this is a nonlinear calculation, we may want to play with the amplitude of the current/field, so we set the `amplitude` property explicitly to our parameter `amp`, above.

```py
sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ex,
                     center=mp.Vector3(0, 0, -0.5*sz + dpml), amplitude=amp)]
```

We'll want the frequency spectrum at the $+z$ end of the computational cell. In a linear problem, we normally look at the spectrum over the same frequency range as our source, because other frequencies are zero. In this case, however, we will look from `fcen/2` to `4*fcen`, to be sure that we can see the third-harmonic frequency.

```py
nfreq = 400
fmin = fcen / 2.0
fmax = fcen * 4

sim = mp.Simulation(cell_size=cell,
                    sources=sources,
                    boundary_layers=pml_layers,
                    default_material=default_material,
                    resolution=resolution,
                    dimensions=dimensions)

mon_pt = mp.Vector3(0, 0, 0.5 * sz - dpml - 0.5)

if flux_spectrum:
    trans = sim.add_flux(
        0.5 * (fmin + fmax), fmax - fmin, nfreq, mp.FluxRegion(mon_pt),
    )
else:
    trans1 = sim.add_flux(
        fcen, 0, 1, mp.FluxRegion(mon_pt)
    )
    trans3 = sim.add_flux(
        3 * fcen, 0, 1, mp.FluxRegion(mon_pt)
    )
```

Note that DFT decimation is off by default whenever nonlinearities are present.

Finally, we'll run the sources, plus additional time for the field to decay at the flux plane, and output the flux spectrum:

```py
sim.run(until_after_sources=mp.stop_when_fields_decayed(
        50, mp.Ex, mon_pt, 1e-6))

if flux_spectrum:
    freqs = mp.get_flux_freqs(trans)
    trans_flux = mp.get_fluxes(trans)
    return freqs, trans_flux
else:
    print(
        f"harmonics:, {k}, {amp}, {mp.get_fluxes(trans1)[0]}, "
        f"{mp.get_fluxes(trans3)[0]}"
    )
    return mp.get_fluxes(trans1)[0], mp.get_fluxes(trans3)[0]
```

In the first part of this tutorial, we plot the transmitted power spectrum for various values of $\chi^{(3)}$. In a linear calculation, we normalize the transmission against some reference spectrum, but in this case there is no obvious normalization so we will just plot the raw data for several values of `k` (i.e. of $\chi^{(3)}$):

![](../images/3rd-harm-1d-flux.png#center)

For small values of $\chi^{(3)}$, we see a peak from our source at $\omega=\frac{1}{3}$ and another peak precisely at the third-harmonic frequency $3\omega=1$. As the $\chi^{(3)}$ gets larger, frequency-mixing *within* the peaks causes them to broaden, and finally for $\chi^{(3)}=1$ we start to see a noisy, broad-spectrum transmission due to the phenomenon of **modulation instability**. Notice also that at around $10^{-13}$ the data looks weird; this is probably due to our finite simulation time, imperfect absorbing boundaries, etcetera. We haven't attempted to analyze it in detail for this case.

In the second part of the tutorial, we investigate the dependence of the power at $\omega$ and $3\omega$ as a function of $\chi^{(3)}$ and the current amplitude. We could, of course, interpolate the flux spectrum above to get the desired frequencies, but it is easier just to add two flux regions to Meep and request exactly the desired frequency components. That is, we'll add `tran1` and `tran3` monitor before `sim.run` as implemented in the function `third_harmonic_generation` shown above.

We could print these with more `display_fluxes` lines, but it is nice to print these on a single line along with $\chi^{(3)}$ and the amplitude, so that we can eventually put them all into one table in post processing. To do this, we'll use the lower-level function `get_fluxes(trans1)`, which returns a list of the flux values, and take the first element of the list since there is only one:

```py
print("harmonics:, {}, {}, {}, {}".format(k, amp, mp.get_fluxes(trans1)[0], mp.get_fluxes(trans3)[0]))
```

Notice how we separated everything with commas, and prefixed the line with `"harmonics:"` for easy grepping later.

Finally, we specify the two different parts of the computation and generate the plots based on the results.

```py
if __name__ == "__main__":
    # Part 1: plot transmitted power spectrum for several values of χ(3).
    nfreq = 400
    logk = range(-3,1)
    tflux = np.zeros((nfreq,len(logk)))
    for i, lk in enumerate(logk):
        freqs, tflux[:, i] = third_harmonic_generation(lk,nfreq=nfreq,
                                                       flux_spectrum=True)

    fig, ax = plt.subplots()
    ax.semilogy(freqs,tflux[:,0],'bo-',label='$\chi^{(3)}$=0.001')
    ax.semilogy(freqs,tflux[:,1],'ro-',label='$\chi^{(3)}$=0.01')
    ax.semilogy(freqs,tflux[:,2],'go-',label='$\chi^{(3)}$=0.1')
    ax.semilogy(freqs,tflux[:,3],'co-',label='$\chi^{(3)}$=1')
    ax.set_xlabel('frequency')
    ax.set_ylabel('transmitted power (a.u.)')
    ax.set_xlim(0.2,1.2)
    ax.set_ylim(1e-15,1e2)
    ax.legend()
    fig.savefig('transmitted_power_vs_frequency_vary_logk.png',dpi=150,
                bbox_inches='tight')

    # Part 2: plot transmittance vs. χ(3) for frequencies ω and 3ω.
    logk = np.arange(-6.0,0.2,0.2)
    first_order = np.zeros(len(logk))
    third_order = np.zeros(len(logk))
    for i, lk in enumerate(logk):
        first_order[i], third_order[i] = third_harmonic_generation(lk,
                                                                   flux_spectrum=False)

    input_flux = first_order[0]
    fig, ax = plt.subplots()
    ax.loglog(10**logk,first_order / input_flux,'ro-',label='$\omega$')
    ax.loglog(10**logk,third_order / input_flux,'bo-',label='$3\omega$')
    ax.loglog(10**logk,(10**logk)**2,'k-',label='quadratic line')
    ax.set_xlabel('$\chi^{(3)}$')
    ax.set_ylabel('transmission / incident power')
    ax.legend()
    fig.savefig('transmittance_vs_chi3.png',dpi=150,bbox_inches='tight')
```

We divide the transmitted power at all values of $\chi^{(3)}$ by the transmitted power at $\omega$ for the smallest $\chi^{(3)}$ of $10^{-6}$ which is 225.25726603587043. We then plot the fractional transmission at $\omega$ and $3\omega$ as a function of $\chi^{(3)}$ on a log-log scale:

![](../images/3rd-harm-1d-vs-chi.png#center)

As can be shown from coupled-mode theory or, equivalently, follows from [Fermi's golden rule](https://en.wikipedia.org/wiki/Fermi's_golden_rule), the third-harmonic power must go as the *square* of $\chi^{(3)}$ as long as the nonlinearity is weak (i.e. in the first Born approximation limit, where the ω source is not depleted significantly). This is precisely what we see on the above graph, where the slope of the black line indicates an exact quadratic dependence, for comparison. Once the nonlinearity becomes strong enough, however, this approximation is no longer valid and the dependence is complicated.

Finally, we note that increasing the current amplitude by a factor of $F$ or the Kerr susceptibility $\chi^{(3)}$ by a factor $F^3$ should generate the *same* third-harmonic power in the *weak* nonlinearity approximation. And indeed, we see:

```sh
harmonics:, 0.001, 1.0, 225.2091048223644, 0.021498041565849526
```

```sh
harmonics:, 1e-06, 10.0, 22525.588597389557, 0.021791784143189268
```

which have third-harmonic powers differing by about 1% (last column).
