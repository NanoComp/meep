---
# Meep Tutorial Local density of states
---

In this example, we will demonstrate the LDOS feature of Meep by investigating the Purcell enhancement phenomena in a metallic microcavity. The LDOS, in general, has many important uses for understanding classical dipole sources, but also in many physical phenomena that can be understood semiclassically in terms of dipole currents -- for example, the spontaneous emission rate of atoms (key to fluorescence and lasing phenomena) is proportional to the LDOS. The LDOS is equivalent to the power radiated by a unit dipole, $P=\frac{1}{2}\operatorname{Re}[\textbf{E}^*\cdot\textbf{J}]$, which, alternatively, is really just a measure of how much the harmonic modes of a system overlap with the source point (proportional, as well, to the radiation resistance of a dipole antenna) and is a useful quantity in electromagnetism due to the fact that the <i>same</i> current radiates a <i>different</i> amount of power depending on the surrounding geometry. Analytically, the per-polarization LDOS is exactly proportional to the power radiated by an $\ell$-oriented point-dipole current, $p(t)$, at a given position in space. For a more mathematical treatment of the theory behind the LDOS, we refer you to the relevant discussion in [chapter 4](http://arxiv.org/abs/1301.5366) of our [book](http://www.artechhouse.com/Main/Books/Advances-in-FDTD-Computational-Electrodynamics-Pho-2011.aspx), but for now we simply give the result:

$$\operatorname{LDOS}_{\ell}(\vec{x}_0,\omega)=-\frac{2}{\pi}\varepsilon(\vec{x}_0)\frac{\operatorname{Re}[\hat{E}_{\ell}(\vec{x}_0,\omega)\hat{p}(\omega)^*]}{|\hat{p}(\omega)|^2}$$

where the $|\hat{p}(\omega)|^2$ normalization is necessary for obtaining the power exerted by a unit-amplitude dipole (assuming linear materials). In FDTD, computing the LDOS therefore is straightforward: excite a point dipole source and accumulate the Fourier transforms of the field at a given point in space to obtain the entire LDOS spectrum in a single calculation. This is implemented in the `dft-ldos` feature which is the subject of this tutorial.

A lossless localized mode yields a δ-function spike in the LDOS, whereas a <i>lossy</i> (i.e. arising from either small absorption or radiation) localized mode -- a <i>resonant cavity mode</i> -- leads to a Lorentzian peak. The large enhancement in the LDOS at the resonant peak is known as a [Purcell effect](https://en.wikipedia.org/wiki/_Purcell_effect), named after Purcell's proposal for enhancing spontaneous emission of an atom in a cavity (analogous to a microwave antenna resonating in a metal box). In this case, the resonant mode's contribution to the LDOS at $\omega^{(n)}$ can be shown to be:

$$\operatorname{resonant\ LDOS} \approx \frac{2}{\pi\omega^{(n)}} \frac{Q^{(n)}}{V^{(n)}}$$

where $Q^{(n)}=\omega^{(n)}/2\gamma^{(n)}$ is the dimensionless quality factor and $V^{(n)}$ is the modal volume. This represents another way to compute the LDOS and in this tutorial we will verify this expression by comparing it to the earlier one.

Here, we consider the simple example of a two-dimensional perfect-metal $a$x$a$ cavity of finite thickness 0.1$a$, with a small notch of width $w$ on one side that allows the modes to escape. The nice thing about this example is that in the absence of the notch, the lowest-frequency TM-polarized mode is known analytically to be $E_z^{(1)}=\frac{4}{a^2}\sin(\pi x/a)\sin(\pi y/a)$, with a frequency $\omega^{(1)}=\sqrt{2}\pi c/a$ and modal volume $V^{(1)}=a^2/4$. The notch slightly perturbs this solution, but more importantly the opening allows the confined mode to radiate out into the surrounding air, yielding a finite $Q$. For $w \ll a$, this radiative escape occurs via an evanescent (sub-cutoff) mode of the channel waveguide formed by the notch, and it follows from inspection of the evanescent decay rate $\sqrt{(\pi/\omega)^2-(\omega^{(1)})^2}/c$ that the lifetime scales asymptotically as $Q^{(1)} \sim e^{\#/\omega}$ for some coefficient \#. We will validate both this prediction and the LDOS calculations above by computing the LDOS at the center of the cavity (the point of peak $|\vec{E}|$) in two ways. First, we compute the LDOS directly from the power radiated by a dipole, Fourier-transforming the result of a pulse using the `dft-ldos` command. Second, we compute the cavity mode and its lifetime $Q$ using `harminv` and then compute the LDOS by the Purcell formula shown above. The latter technique is much more efficient for high Q (small $w$), since one must run the simulation for a very long time to directly accumulate the Fourier transform of a slowly decaying mode. The two calculations, we will demonstrate, agree (to within a discretization error), verifying the LDOS analysis above, and $Q/V$ is asymptotically linear on a semilog scale versus $1/w$ as predicted.

We'll first set up the two dimensional simulation with the metal cavity and PML absorbing boundary layers

```
 (set-param! resolution 200)
 (define-param sxy 2)
 (define-param dpml 1)
 (set! sxy (+ sxy (* 2 dpml)))
 (set! geometry-lattice (make lattice (size sxy sxy no-size)))
 (set! pml-layers (list (make pml (thickness dpml))))
 (define-param a 1)
 (define-param t 0.1)
 (set! geometry (list 
     (make block (center 0 0) (size (+ a (* 2 t)) (+ a (* 2 t)) infinity) (material metal))
     (make block (center 0 0) (size a a infinity) (material air))))
```


Next we'll create a notch opening in the cavity so that the field can radiate away

```
 (define-param w 0)
 (if (> w 0)
       (set! geometry
            (append geometry
                  (list (make block (center (/ a 2) 0) (size (* 2 t) w infinity)
                        (material air))))))
```


We can now set up the TM-polarized source in the middle of the cavity (where we will also compute the LDOS; they are always co-located). We know the mode frequency of the closed cavity analytically (of course, the frequency will shift with the size of the notch which necessitates a Gaussian pulse) and also note that in Meep, frequency is specified in units of $2\pi$

```
 (define-param fcen (/ (sqrt 0.5) a))
 (define-param df 0.2)
 (set! sources (list (make source
        (src (make gaussian-src (frequency fcen) (fwidth df))) (component Ez) (center 0 0))))
```


As both the structure and sources have a mirror symmetry in the Y direction, we can exploit this to halve the size of the computational cell

```
 (set! symmetries (list (make mirror-sym (direction Y))))
```


In the first part of the calculation, we will compute the Purcell enhancement and thus need to calculate the mode frequency and its quality factor

```
 (define-param Th 500)
 (run-sources+ Th (after-sources (harminv Ez (vector3 0) fcen df)))
 (define f (harminv-freq-re (car harminv-results)))
 (define Q (harminv-Q (car harminv-results)))
 (define Vmode (* 0.25 a a))
 (print "ldos0:, " (/ Q Vmode (* 2 pi f pi 0.5)))
```


Next, we will rerun the same simulation and instead compute the LDOS using Meep's built-in `dft-ldos` feature at a single frequency corresponding to the mode frequency.

```
 (reset-meep)
 (define-param T (* 2 Q (/ f)))
 (run-sources+ T (dft-ldos f 0 1))
```


We need to run for a sufficiently long time to ensure that the Fourier transformed fields have converged and an appropriate time interval is, due to the Fourier Uncertainty Principle, just one period of the decay (which we can determine by the $Q$ we calculated previously). Of course, the smaller the notch size becomes and the higher the corresponding Q of the mode, the longer the simulation has to run which is why the former calculation is much more efficient for slowly-decaying modes.

We run several simulations spanning a number of different notch sizes and plot the result in the following figure which shows good agreement between the two methods.


![center|LDOS for a metallic cavity.](../images/Metalcavity_ldos.png)


