---
# Cylindrical Coordinates
---

Meep supports the simulation of Maxwell's equations in [cylindrical coordinates](https://en.wikipedia.org/wiki/Cylindrical_coordinate_system) for structures that have [continuous rotational symmetry around the *z* axis](../Exploiting_Symmetry.md#cylindrical-symmetry). This reduces problems in 3d to 2d, and 2d to 1d, if there is sufficient symmetry.

[TOC]

Modes of a Ring Resonator
-------------------------

In [Tutorial/Basics/Modes of a Ring Resonator](Basics.md#modes-of-a-ring-resonator), the modes of a ring resonator were computed by performing a 2d simulation. This example involves simulating the *same* structure while [exploiting](../Exploiting_Symmetry.md) the fact that the system has *continuous* rotational symmetry, by performing the simulation in cylindrical coordinates. The simulation script is in [examples/ring-cyl.ctl](https://github.com/NanoComp/meep/blob/master/scheme/examples/ring-cyl.ctl).

The parameters of the problem are defined with exactly the same values as in the 2d simulation:

```scm
(define-param n 3.4)     ; index of waveguide
(define-param w 1)       ; width of waveguide
(define-param r 1)       ; inner radius of ring
(define-param pad 4)     ; padding between waveguide and edge of PML
(define-param dpml 2)    ; thickness of PML
```

The dimensions and size of the computational cell are defined:

```scm
(define sr (+ r w pad dpml))   ; radial size (cell is from 0 to sr)
(set! dimensions CYLINDRICAL)
(set! geometry-lattice (make lattice (size sr no-size no-size)))
```

The key thing is to set the `dimensions` parameter to `CYLINDRICAL`. This means that all vectors represent ($r$,$\phi$,$z$) coordinates instead of ($x$,$y$,$z$). The computational cell in the $r$ direction is of size `sr = r + w + pad + dpml`, and runs from `0` to `sr` (by default) rather than from `-sr/2` to `sr/2` as it would for any other dimension. Note that the $z$ size is `no-size` because it is in 2d. The $\phi$ size is also `no-size`, corresponding to the continuous rotational symmetry. A finite $\phi$ size might correspond to discrete rotational symmetry, but this is not currently supported.

In particular, in systems with continuous rotational symmetry, by an analogue of Bloch's theorem, the angular dependence of the fields can always be chosen in the form $\exp(i m \phi)$ for some integer $m$. Meep uses this fact to treat the angular dependence analytically, with $m$ given by the [input variable](../Scheme_User_Interface.md#input-variables) `m` which is set to a command-line argument that is 3 by default.

```scm
(set-param! m 3)
```

This is essentially a 1d calculation where Meep must discretize the $r$ direction only. For this reason, it will be much faster than the previous 2d calculation.

The geometry is now specified by a single `block` object &mdash; remember that this is a block in cylindrical coordinates, so that it really specifies an annular ring:

```scm
(set! geometry (list
                (make block (center (+ r (/ w 2))) (size w infinity infinity)
                      (material (make dielectric (index n))))))
(set! pml-layers (list (make pml (thickness dpml))))
(set-param! resolution 10)
```

PMLs are on "all" sides. The $z$ direction has no thickness and therefore it is automatically periodic with no PML. PML is also omitted from the boundary at $r$=0 which is handled by the analytical reflection symmetry.

The remaining inputs are almost exactly the same as in the previous 2d simulation. A single Gaussian point source is added in the $z$ direction to excite $E_z$-polarized modes, with some center frequency and width:

```scm
(define-param fcen 0.15) ; pulse center frequency
(define-param df 0.1)    ; pulse width (in frequency)
(set! sources (list
               (make source
                 (src (make gaussian-src (frequency fcen) (fwidth df)))
                 (component Ez) (center (+ r 0.1) 0 0))))
              
```

Note that this isn't really a point source, however, because of the cylindrical symmetry &mdash; it is really a ring source with $\phi$ dependence $\exp(i m \phi)$. Finally, as before, the fields are timestepped until the source has turned off, plus 200 additional time units during which [Harminv](../Scheme_User_Interface.md#harminv) is used to analyze the $E_z$ field at a given point to extract the frequencies and decay rates of the modes.

```scm
(run-sources+ 200 (after-sources (harminv Ez (vector3 (+ r 0.1)) fcen df)))
```

At the very end, one period of the fields is output to create an animation. A single field output would be a 1d dataset along the $r$ direction, so to make things more interesting `to-appended` is used to append these datasets to a single HDF5 file to obtain an $r \times t$ 2d dataset. Also `in-volume` is used to specify a larger output volume than just the computational cell: in particular, the output is from `-sr` to `sr` in the $r$ direction, where the $-r$ field values are automatically inferred from the reflection symmetry.

```scm
(run-until (/ 1 fcen) 
           (in-volume (volume (center 0) (size (* 2 sr)))
                      (at-beginning output-epsilon)
                      (to-appended "ez" 
                                   (at-every (/ 1 fcen 20) output-efield-z))))
```

The simulation is ready to be run. Recall that, in the 2d calculation, three modes were obtained in this frequency range: (1) $\omega$=0.11785 with $Q$=77 and an $m$=3 field pattern, (2) $\omega$=0.14687 with $Q$=351 and an $m$=4 field pattern, and (3) $\omega$=0.17501 with $Q$=1630 and an $m$=5 field pattern. To verify the correctness of this script, the *same* modes should be obtained with some differences due to the finite resolution, except now *three* calculations are necessary, a separate one for each value of $m$. It will still be much faster than the 2d simulation because these simulations are 1d.

In particular, the three calculations are:

```sh
unix% meep m=3 ring-cyl.ctl | grep harminv
unix% meep m=4 ring-cyl.ctl | grep harminv
unix% meep m=5 ring-cyl.ctl | grep harminv
```

giving the combined output:

```
harminv0:, frequency, imag. freq., Q, |amp|, amplitude, error
harminv0:, 0.11835455441250631, -6.907792691647415e-4, 85.66741917111612, 0.02570190626349302, -0.02402703883357199-0.00912630212448642i, 5.286949731053267e-10+0.0i
harminv0:, 0.1475578747705309, -1.938438860632441e-4, 380.61008208014414, 0.19361245519715206, 0.1447225471614173+0.12861246887677943i, 5.889273063545974e-11+0.0i
harminv0:, 0.1759448592380757, -4.900590034953583e-5, 1795.1395442502285, 0.0452479314013276, -0.014395016792255884-0.042897072017212545i, 1.6343462235932872e-10+0.0i
```

This is indeed very close to the 2d simulations: the frequencies are within 1% of the previous values. The $Q$ values (lifetimes) differ by a larger amount although they are still reasonably close.

Which is more accurate, the 2d or the cylindrical simulation? This question can be answered by increasing the resolution in both cases and seeing what they converge towards. In particular, let's focus on the $m$=4 mode. In the cylindrical case, if the resolution is doubled to 20, the mode is $\omega$=0.14748 and $Q$=384. In the 2d case, if the resolution is doubled to 20 the mode is $\omega$=0.14733 and $Q$=321. It looks like the frequencies are clearly converging together and that the cylindrical simulation is more accurate (as you might expect since it describes the $\phi$ direction analytically). But the $Q$ values seem to be getting *farther* apart &mdash; what's going on?

The problem is twofold. First, there is some signal-processing error in determining $Q$ in the 2d case, as indicated by the "error" column of the `harminv` output which is only 4e-7 for the 2d simulation vs. 6e-11 for the cylindrical case. This error can be reduced by running with a narrower bandwidth source, which excites just one mode and gives a cleaner signal, or by analyzing over a longer time than 200. Doing the former, we find that the 2d value of $Q$ at a resolution of 20 should really be $Q$=343. Second, [PML](../Perfectly_Matched_Layer.md) absorbing layers are really designed to absorb planewaves incident on flat interfaces, but here we have a *cylindrical* PML layer. Because of this, there are larger numerical reflections from the PML in the cylindrical simulation, which we can rectify by pushing the PML out to a larger radius (i.e. using a larger value of `pad`) and/or increasing the PML thickness (increasing `dpml`) so that it turns on more adiabatically. In the cylindrical simulation for `resolution = 20`, if the PML thickness is increased to `dpml = 16`, the $Q$ is 343, which is in much better agreement with the 2d calculation and if the PML thickness is increased to `dpml = 32` the $Q$ is the same 343, so it seems to be converged.

This illustrates the general principle that you need to [check several parameters to ensure that results are converged in time-domain simulations](../FAQ.md#checking-convergence): the resolution, the run time, the PML thickness, etcetera.

Finally, the field images are obtained. Since one mode per `m` is being excited here anyway, according to `harminv`, there is no real need for a narrow-band source. This will be used anyway just to remind you of the general procedure, however, e.g. for the $\omega$=0.118, $m$=3 mode:

```sh
unix% meep m=3 fcen=0.118 df=0.01 ring-cyl.ctl
unix% h5topng -S 2 -Zc dkbluered -C ring-cyl-eps-001200.00.h5 ring-cyl-ez.h5
```

Note that, because of the `to-appended` command, the `ring-cyl-ez.h5` file is a 160$\times$18 dataset corresponding to an $r \times t$ slice. Repeating this for all three modes results in the images:

<center>
$E_z$ for ω=0.118 $m$=3 mode:  
![](../images/Ring-cyl-ez-0.118.png)

$E_z$ for ω=0.148 $m$=4 mode:  
![](../images/Ring-cyl-ez-0.148.png)

$E_z$ for ω=0.176 $m$=5 mode:  
![](../images/Ring-cyl-ez-0.176.png)
</center>

Because only the $\phi$=0 slice is used, the visual distinction between $m$ values is much less than with the 2d simulation. What is apparent is that, as the frequency increases, the mode becomes more localized in the waveguide and the radiating field (seen in the $r \times t$ slice as curved waves extending outward) becomes less, as expected.

Sensitivity Analysis via Perturbation Theory
--------------------------------------------

For a given mode of the ring resonator, it is often useful to know how sensitive the resonant frequency $\omega$ is to small changes in the ring radius $r$ by computing its derivative $\partial\omega/\partial r$. The gradient is also a useful quantity for shape optimization because it can be paired with fast iterative methods such as [BFGS](https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm) to find local optima. The "brute-force" approach for computing the gradient is via a finite-difference approximation requiring *two* simulations of the (1) unperturbed [$\omega(r)$] and (2) perturbed [$\omega(r+\Delta r)$] structures. Since each simulation is potentially costly, an alternative approach based on semi analytics is to use [perturbation theory](https://en.wikipedia.org/wiki/Perturbation_theory) to obtain the gradient from the fields of the unperturbed structure. This involves a single simulation and is often more computationally efficient than the brute-force approach although some care is required to set up the calculation properly.  (More generally, [adjoint methods](https://math.mit.edu/~stevenj/18.336/adjoint.pdf) can be used to compute any number of derivatives with a single additional simulation.)

[Pertubation theory for Maxwell equations involving high index-contrast dielectric interfaces](http://math.mit.edu/~stevenj/papers/KottkeFa08.pdf) is reviewed in Chapter 2 of [Photonics Crystals: Molding the Flow of Light, 2nd Edition (2008)](http://ab-initio.mit.edu/book/). The formula (equation 30 on p.19) for the frequency shift $\Delta \omega$ resulting from the displacement of a block of $\varepsilon_1$-material towards $\varepsilon_2$-material by a distance $\Delta h$ (perpendicular to the boundary) is:

<center>

$$ \Delta\omega = -\frac{\omega}{2} \frac{ \iint d^2 \vec{r} \big[ (\varepsilon_1 - \varepsilon_2) |\vec{E}_{\parallel}(\vec{r})|^2 - \big(\frac{1}{\varepsilon_1} - \frac{1}{\varepsilon_2}\big)|\varepsilon\vec{E}_{\perp}|^2\big] \Delta h}{\int d^3\vec{r} \varepsilon(\vec{r})|\vec{E}(\vec{r})|^2} + O(\Delta h^2)$$

</center>

In this expression, $\vec{E}_{\parallel}(\vec{r})$ is the component of $\vec{E}$ that is parallel to the surface, and $\varepsilon\vec{E}_{\perp}$ is the component of $\varepsilon\vec{E}$ that is perpendicular to the surface. These two components are guaranteed to be continuous across an interface between two isotropic dielectric materials. In this demonstration, $\partial\omega/\partial r$ is computed using this formula and the results are validated by comparing with the finite-difference approximation: $[\omega(r+\Delta r)-\omega(r)]/\Delta r$.

There are three parts to the calculation: (1) find the resonant frequency of the unperturbed geometry using a broadband Gaussian source, (2) find the resonant mode profile of the unperturbed geometry using a narrowband source and from these fields compute the gradient via the perturbation-theory formula, and (3) find the resonant frequency of the perturbed geometry and from this compute the gradient using the finite-difference approximation. The perturbation is applied only to the inner and outer ring radii. The ring width is constant. There are two types of modes which are computed in separate simulations using different source polarizations: parallel ($E_z$) and perpendicular ($H_z$) to the interface.

The simulation script is in [examples/perturbation-theory.ctl](https://github.com/NanoComp/meep/blob/master/scheme/examples/perturbation-theory.ctl).

```scm
(set-param! resolution 100)  ; pixels/um

(define-param perpendicular? true)
(define src-cmpt (if perpendicular? Hz Ez))
(define fcen (if perpendicular? 0.21 0.17))  ; pulse center frequency

(define-param n 3.4)  ; index of waveguide
(define-param w 1)    ; width of waveguide
(define-param r 1)    ; inner radius of ring
(define-param pad 4)  ; padding between waveguide and edge of PML
(define-param dpml 2) ; thickness of PML
(set-param! m 5)      ; angular dependence

(set! pml-layers (list (make pml (thickness dpml))))

(define sr (+ r w pad dpml)) ; radial size (cell is from 0 to sr)
(set! dimensions CYLINDRICAL)
(set! geometry-lattice (make lattice (size sr no-size no-size)))

(set! geometry (list (make block
                       (center (+ r (/ w 2)))
                       (size w infinity infinity)
                       (material (make dielectric (index n))))))

;; find resonant frequency of unperturbed geometry using broadband source

(define df (* 0.2 fcen))  ; pulse width (in frequency)

(set! sources (list
               (make source
                 (src (make gaussian-src (frequency fcen) (fwidth df)))
                 (component src-cmpt)
                 (center (+ r 0.1)))))

(run-sources+ 100 (after-sources (harminv src-cmpt (vector3 (+ r 0.1)) fcen df)))

(define frq-unperturbed (harminv-freq-re (car harminv-results)))

(reset-meep)

;; unperturbed geometry with narrowband source centered at resonant frequency

(set! pml-layers (list (make pml (thickness dpml))))

(set! geometry-lattice (make lattice (size sr no-size no-size)))

(set! geometry (list (make block
                       (center (+ r (/ w 2)))
                       (size w infinity infinity)
                       (material (make dielectric (index n))))))

(set! fcen frq-unperturbed)
(set! df (* 0.05 fcen))

(set! sources (list
               (make source
                 (src (make gaussian-src (frequency fcen) (fwidth df)))
                 (component src-cmpt)
                 (center (+ r 0.1)))))

(run-sources+ 100)

(define deps (- 1 (* n n)))
(define deps-inv (- 1 (/ (* n n))))

(define numerator-integral 0)

(if perpendicular?
    (let ((para-integral (* deps 2 pi (- (* r (sqr (magnitude (get-field-point Ep (vector3 r))))) (* (+ r w) (sqr (magnitude (get-field-point Ep (vector3 (+ r w)))))))))
          (perp-integral (* deps-inv 2 pi (- (* (+ r w) (sqr (magnitude (get-field-point Dr (vector3 (+ r w)))))) (* r (sqr (magnitude (get-field-point Dr (vector3 r)))))))))
      (set! numerator-integral (+ para-integral perp-integral)))
    (set! numerator-integral (* deps 2 pi (- (* r (sqr (magnitude (get-field-point Ez (vector3 r))))) (* (+ r w) (sqr (magnitude (get-field-point Ez (vector3 (+ r w))))))))))

(define denominator-integral (electric-energy-in-box (volume (center (* 0.5 sr)) (size sr))))
(define perturb-theory-dw-dR (/ (* -1 frq-unperturbed numerator-integral) (* 4 denominator-integral)))

(reset-meep)

;; perturbed geometry with narrowband source

(define-param dr 0.01)

(set! pml-layers (list (make pml (thickness dpml))))

(set! geometry-lattice (make lattice (size sr no-size no-size)))

(set! geometry (list (make block
                       (center (+ r dr (/ w 2)))
                       (size w infinity infinity)
                       (material (make dielectric (index n))))))

(set! sources (list
               (make source
                 (src (make gaussian-src (frequency fcen) (fwidth df)))
                 (component src-cmpt)
                 (center (+ r dr 0.1)))))

(run-sources+ 100 (after-sources (harminv src-cmpt (vector3 (+ r 0.1)) fcen df)))

(define frq-perturbed (harminv-freq-re (car harminv-results)))

(define finite-diff-dw-dR (/ (- frq-perturbed frq-unperturbed) dr))

(print "dwdR:, " perturb-theory-dw-dR " (pert. theory), " finite-diff-dw-dR " (finite diff.)\n")
```

There are three things to note. First, there is a built-in function `electric-energy-in-box` which calculates the integral of $\vec{E}\cdot\vec{D}/2 = \varepsilon|E|^2/2$. This is exactly the integral in the denominator, divided by 2. In cylindrical coordinates $(r,\phi,z)$, the integrand is [multiplied](https://en.wikipedia.org/wiki/Cylindrical_coordinate_system#Line_and_volume_elements) by the circumference $2\pi r$, or equivalently the integral is over an annular volume. Second, for the case involving the $H_z$ source, both parallel ($E_{\parallel}=E_{\phi}$) and perpendicular ($\varepsilon E_{\perp}=D_r$) fields are present which must be included in the numerator as separate terms. Field values anywhere in the grid obtained with `get-field-point` are [automatically interpolated](../Introduction.md#the-illusion-of-continuity); i.e., no additional post-processing is necessary. Third, when comparing the perturbation-theory result to the finite-difference approximation, there are *two* convergence parameters: the resolution and $\Delta r$. In principle, to demonstrate agreement with perturbation theory, the limit of the resolution should be taken to infinity and *then* the limit of $\Delta r$ to 0. In practice, this can be obtained by doubling the resolution at a given $\Delta r$ until it is sufficiently converged, then halving $\Delta r$ and repeating.

For an $E_z$ source (parallel to the interface), at a resolution of 100 the results are -0.0854469639 (perturbation theory) and -0.0852124909 (finite difference). Doubling the resolution to 200, the results are -0.0854460732 (perturbation theory) and -0.0852115350 (finite difference). Both results have converged to at least five digits. The relative error at resolution 200 is 0.3%. The mode has a $\omega$ of 0.175 and $Q$ of 1800.

For an $H_z$ source (perpendicular to the interface), at a resolution of 100 the results are -0.0805038572 (perturbation theory) and -0.0798087331 (finite difference). Doubling the resolution to 200, the results are -0.0802028346 (perturbation theory) and -0.0798088015 (finite difference). Both results have converged to at least three digits. The relative error at resolution 200 is 0.5%. The error is larger in this case due to the presence of the discontinuous fields at the dielectric interface which are handled by [subpixel smoothing](../Subpixel_Smoothing.md). The mode has a $\omega$ of 0.208 and $Q$ of 1200.

Finally, as reference, the same calculation can be set up in Cartesian coordinates as a 2d simulation. The simulation script is in [examples/perturbation-theory-2d.ctl](https://github.com/NanoComp/meep/blob/master/python/examples/perturbation-theory-2d.ctl). The results are comparable to the cylindrical coordinate case (a 1d calculation) but the 2d simulation is much slower and less accurate at the same grid resolution.