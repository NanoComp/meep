---
# Mode Decomposition
---

This tutorial demonstrates the [mode-decomposition](../Mode_Decomposition.md) feature which is used to decompose a given mode profile via the Fourier-transformed fields into a superposition of harmonic basis modes.

[TOC]

Reflectance of a Waveguide Taper
--------------------------------

This example involves computing the reflectance of the fundamental mode of a linear waveguide taper. The structure and the simulation parameters are shown in the schematic below. We will verify that computing the reflectance, the fraction of the incident power which is reflected, using two different methods produces nearly identical results: (1) mode decomposition and (2) [Poynting flux](../Introduction.md#transmittancereflectance-spectra). Also, we will demonstrate that the scaling of the reflectance with the taper length is quadratic, consistent with analytical results from [Optics Express, Vol. 16, pp. 11376-92, 2008](http://www.opticsinfobase.org/abstract.cfm?URI=oe-16-15-11376).

<center>
![](../images/waveguide-taper.png)
</center>

The structure, which can be viewed as a [two-port network](https://en.wikipedia.org/wiki/Two-port_network), consists of a single-mode waveguide of width 1 μm (`w1`) at a wavelength of 6.67 μm and coupled to a second waveguide of width 2 μm (`w2`) via a linearly-sloped taper of variable length `Lt`. The material is silicon with constant ε=12. The geometry is defined using a single [`prism`](../Scheme_User_Interface.md#prism) object with eight vertices. PML absorbing boundaries surround the entire cell. An eigenmode current source with E<sub>z</sub> polarization is used to launch the fundamental mode. The dispersion relation (or "band diagram") of the single-mode waveguide is shown in [Tutorial/Eigenmode Source](Eigenmode_Source.md). There is an eigenmode-expansion monitor placed at the midpoint of the first waveguide. This is a line monitor which extends beyond the waveguide in order to span the entire mode profile including its evanescent tails. The Fourier-transformed fields along this line monitor are used to compute the basis coefficients of the harmonic modes. These are computed separately via the eigenmode solver [MPB](https://mpb.readthedocs.io/en/latest/). This is described in [Mode Decomposition](../Mode_Decomposition.md) where it is also shown that the squared magnitude of the mode coefficient is equivalent to the power (Poynting flux) in the given eigenmode. The ratio of the complex mode coefficients can be used to compute the [S parameters](https://en.wikipedia.org/wiki/Scattering_parameters). In this example, we are computing |S<sub>11</sub>|<sup>2</sup> which is the reflectance (shown in the line prefixed by "refl:,"). Another line monitor could have been placed in the second waveguide to compute the transmittance or |S<sub>21</sub>|<sup>2</sup> into the various guided modes (since the second waveguide is multi mode). The scattered power into the radiative modes can then be computed as 1-|S<sub>11</sub>|<sup>2</sup>-|S<sub>21</sub>|<sup>2</sup>. As usual, a normalization run is required involving a straight waveguide to compute the power in the source.

The structure has mirror symmetry in the $y$ direction which can be exploited to reduce the computation size by a factor of two. This requires that we use `add-flux` rather than `add-mode-monitor` (which is not optimized for symmetry) and specify the `eig-parity` keyword argument as `ODD-Z+EVEN-Y` in the call to `get-eigenmode-coefficients`.

The simulation script is in [examples/mode-decomposition.ctl](https://github.com/NanoComp/meep/blob/master/scheme/examples/mode-decomposition.ctl).

```py
(set-param! resolution 61)  ; pixels/μm

(define-param w1 1.0)       ; width of waveguide 1
(define-param w2 2.0)       ; width of waveguide 2
(define-param Lw 10.0)      ; length of waveguides 1 and 2
(define-param Lt 8.0)       ; length of waveguide taper

(define-param dair 3.0)     ; length of air region
(define-param dpml-x 6.0)   ; length of PML in x direction
(define-param dpml-y 2.0)   ; length of PML in y direction

(define sx (+ dpml-x Lw Lt Lw dpml-x))
(define sy (+ dpml-y dair w2 dair dpml-y))

(define cell (make lattice (size sx sy no-size)))
(set! geometry-lattice cell)

(define Si (make medium (epsilon 12.0)))

(define boundary-layers (list (make pml (direction X) (thickness dpml-x))
                              (make pml (direction Y) (thickness dpml-y))))
(set! pml-layers boundary-layers)

(define-param lcen 6.67)  ; mode wavelength
(define fcen (/ lcen))    ; mode frequency

(define eig-src (list (make eigenmode-source
                        (src (make gaussian-src (frequency fcen) (fwidth (* 0.2 fcen))))
                        (center (vector3 (+ (* -0.5 sx) dpml-x (* 0.2 Lw)) 0 0))
                        (size 0 (- sy (* 2 dpml-y)) 0)
                        (eig-band 1)
                        (eig-match-freq? true)
                        (eig-parity (+ ODD-Z EVEN-Y)))))
(set! sources eig-src)

; straight waveguide
(define sw-vertices (list (vector3 (- (* -0.5 sx) 1) (* 0.5 w1) 0)
                          (vector3 (+ (* 0.5 sx) 1) (* 0.5 w1) 0)
                          (vector3 (+ (* 0.5 sx) 1) (* -0.5 w1) 0)
                          (vector3 (- (* -0.5 sx) 1) (* -0.5 w1) 0)))

(set! geometry (list (make prism
                       (vertices sw-vertices)
                       (axis 0 0 1)
                       (center auto-center)
                       (height infinity)
                       (material Si))))

(define symm (list (make mirror-sym (direction Y))))
(set! symmetries symm)

(define mon-pt (vector3 (+ (* -0.5 sx) dpml-x (* 0.7 Lw)) 0 0))
(define flux1 (add-flux fcen 0 1 (make flux-region (center mon-pt) (size 0 (- sy (* 2 dpml-y)) 0))))

(run-sources+ (stop-when-fields-decayed 50 Ez mon-pt 1e-9))

(save-flux "flux" flux1)
(define res1 (get-eigenmode-coefficients flux1 (list 1) #:eig-parity (+ ODD-Z EVEN-Y)))
(define incident-coeffs (array-ref (list-ref res1 0) 0 0 0))
(define incident-flux (list-ref (get-fluxes flux1) 0))

(reset-meep)

(set! geometry-lattice cell)

(set! pml-layers boundary-layers)

(set! sources eig-src)

; linear taper
(define tp-vertices (list (vector3 (- (* -0.5 sx) 1) (* 0.5 w1) 0)
                          (vector3 (* -0.5 Lt) (* 0.5 w1) 0)
                          (vector3 (* 0.5 Lt) (* 0.5 w2) 0)
                          (vector3 (+ (* 0.5 sx) 1) (* 0.5 w2) 0)
                          (vector3 (+ (* 0.5 sx) 1) (* -0.5 w2) 0)
                          (vector3 (* 0.5 Lt) (* -0.5 w2) 0)
                          (vector3 (* -0.5 Lt) (* -0.5 w1) 0)
                          (vector3 (- (* -0.5 sx) 1) (* -0.5 w1) 0)))

(set! geometry (list (make prism
                       (vertices tp-vertices)
                       (axis 0 0 1)
                       (center auto-center)
                       (height infinity)
                       (material Si))))

(define flux2 (add-flux fcen 0 1 (make flux-region (center mon-pt) (size 0 (- sy (* 2 dpml-y)) 0))))
(load-minus-flux "flux" flux2)

(run-sources+ (stop-when-fields-decayed 50 Ez mon-pt 1e-9))

(define res2 (get-eigenmode-coefficients flux2 (list 1) #:eig-parity (+ ODD-Z EVEN-Y)))
(define taper-coeffs (array-ref (list-ref res2 0) 0 0 1))
(define taper-flux (list-ref (get-fluxes flux2) 0))

(print "refl:, " Lt ", " (/ (sqr (magnitude taper-coeffs)) (sqr (magnitude incident-coeffs))) ", " (/ (- taper-flux) incident-flux) "\n")
```

We compute the reflectance for five different taper lengths: 1, 2, 4, 8, and 16 μm. The Bash commands to run the simulation and extract the reflectance results from the output are:

```
for m in `seq 0 4`; do
    mpirun -np 2 meep Lt=$((2**${m})) mode-decomposition.ctl |tee -a waveguide_taper.out;
done

grep refl: waveguide_taper.out |cut -d, -f2- > waveguide_taper.dat
```

The results are plotted using the Octave/Matlab script below. The plot is shown in the accompanying figure. 

```matlab
f = load('waveguide_taper.dat');
loglog(f(:,1),f(:,2),'bo-',f(:,1),f(:,3),'ro-',f(:,1),0.005./f(:,1).^2,'k-');
xlabel("taper length Lt (um)");
ylabel("reflectance");
legend("mode decomposition","Poynting flux","quadratic reference (1/Lt^2)");
axis([0.9 20 1e-6 1e-2]);
```

<center>
![](../images/refl_coeff_vs_taper_length.png)
</center>

The reflectance values computed using the two methods are nearly identical. For reference, a line with quadratic scaling is shown in black. The reflectance of the linear waveguide taper decreases quadratically with the taper length which is consistent with the analytic theory.