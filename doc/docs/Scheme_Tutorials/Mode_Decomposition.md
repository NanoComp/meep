---
# Mode Decomposition
---

This tutorial demonstrates the [mode-decomposition](../Mode_Decomposition.md) feature which is used to decompose a given mode profile via the Fourier-transformed fields into a superposition of harmonic basis modes. Examples are provided for two kinds of modes in lossless, dielectric media: (1) localized (i.e., guided) and (2) non-localized (i.e., radiative planewave).

[TOC]

Reflectance of a Waveguide Taper
--------------------------------

This example involves computing the reflectance of the fundamental mode of a linear waveguide taper. The structure and the simulation parameters are shown in the schematic below. We will verify that computing the reflectance, the fraction of the incident power which is reflected, using two different methods produces nearly identical results: (1) mode decomposition and (2) [Poynting flux](../Introduction.md#transmittancereflectance-spectra). Also, we will demonstrate that the scaling of the reflectance with the taper length is quadratic, consistent with analytical results from [Optics Express, Vol. 16, pp. 11376-92, 2008](http://www.opticsinfobase.org/abstract.cfm?URI=oe-16-15-11376).

<center>
![](../images/waveguide-taper.png)
</center>

The structure, which can be viewed as a [two-port network](https://en.wikipedia.org/wiki/Two-port_network), consists of a single-mode waveguide of width 1 μm (`w1`) at a wavelength of 6.67 μm and coupled to a second waveguide of width 2 μm (`w2`) via a linearly-sloped taper of variable length `Lt`. The material is silicon with ε=12. The taper geometry is defined using a single [`prism`](../Scheme_User_Interface.md#prism) object with eight vertices. PML absorbing boundaries surround the entire cell. An eigenmode current source with E<sub>z</sub> polarization is used to launch the fundamental mode. The dispersion relation (or "band diagram") of the single-mode waveguide is shown in [Tutorial/Eigenmode Source](Eigenmode_Source.md). There is an eigenmode-expansion monitor placed at the midpoint of the first waveguide. This is a line monitor which extends beyond the waveguide in order to span the entire mode profile including its evanescent tails. The Fourier-transformed fields along this line monitor are used to compute the basis coefficients of the harmonic modes. These are computed separately via the eigenmode solver [MPB](https://mpb.readthedocs.io/en/latest/). This is described in [Mode Decomposition](../Mode_Decomposition.md) where it is also shown that the squared magnitude of the mode coefficient is equivalent to the power (Poynting flux) in the given eigenmode. The ratio of the complex mode coefficients can be used to compute the [S parameters](https://en.wikipedia.org/wiki/Scattering_parameters). In this example, we are computing |S<sub>11</sub>|<sup>2</sup> which is the reflectance (shown in the lines prefixed by "refl:,"). Another line monitor could have been placed in the second waveguide to compute the transmittance or |S<sub>21</sub>|<sup>2</sup> into the various guided modes (since the second waveguide is multi mode). The scattered power into the radiative modes can then be computed as 1-|S<sub>11</sub>|<sup>2</sup>-|S<sub>21</sub>|<sup>2</sup>. As usual, a normalization run is required involving a straight waveguide to compute the power in the source.

The structure has mirror symmetry in the $y$ direction which can be exploited to reduce the computation size by a factor of two. This requires that we use `add-flux` rather than `add-mode-monitor` (which is not optimized for symmetry) and specify the `eig-parity` keyword argument as `ODD-Z+EVEN-Y` in the call to `get-eigenmode-coefficients`.

The simulation script is in [examples/mode-decomposition.ctl](https://github.com/NanoComp/meep/blob/master/scheme/examples/mode-decomposition.ctl).

```scm
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

(set! symmetries symm)

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

Diffraction Spectrum of a Binary Grating
----------------------------------------

The mode-decomposition feature can also be applied to planewaves in homogeneous media with scalar permittivity/permeability (i.e., no anisotropy). This will be demonstrated in this example to compute the diffraction spectrum of a binary phase [grating](https://en.wikipedia.org/wiki/Diffraction_grating). The unit cell geometry of the grating is shown in the schematic below. The grating is periodic in the $y$ direction with periodicity `gp` and has a rectangular profile of height `gh` and duty cycle `gdc`. The grating parameters are `gh`=0.5 μm, `gdc`=0.5, and `gp`=10 μm. There is a semi-infinite substrate of thickness `dsub` adjacent to the grating. The substrate and grating are glass with a refractive index of 1.5. The surrounding is air/vacuum. Perfectly matched layers (PML) of thickness `dpml` are used in the $\pm x$ boundaries.

### Transmittance Spectra for Planewave at Normal Incidence

A pulsed planewave with E<sub>z</sub> polarization spanning wavelengths of 0.4 to 0.6 μm is normally incident on the grating from the glass substrate. The eigenmode monitor is placed in the air region. We will use mode decomposition to compute the transmittance &mdash; the ratio of the power in the $+x$ direction of the diffracted mode relative to that of the incident planewave &mdash; for the first ten diffraction orders. Two simulations are required: (1) an *empty* cell of homogeneous glass to obtain the incident power of the source, and (2) the grating structure to obtain the diffraction orders. At the end of the simulation, the wavelength, angle, and transmittance for each diffraction order are computed.

The simulation script is in [examples/binary_grating.ctl](https://github.com/NanoComp/meep/blob/master/scheme/examples/binary_grating.ctl).

<center>
![](../images/grating.png)
</center>

```scm
(set-param! resolution 60)  ; pixels/μm

(define-param dpml 1.0)     ; PML thickness
(define-param dsub 3.0)     ; substrate thickness
(define-param dpad 3.0)     ; padding between grating and PML
(define-param gp 10.0)      ; grating period
(define-param gh 0.5)       ; grating height
(define-param gdc 0.5)      ; grating duty cycle

(define sx (+ dpml dsub gh dpad dpml))
(define sy gp)

(define cell (make lattice (size sx sy no-size)))
(set! geometry-lattice cell)

(define boundary-layers (list (make pml (thickness dpml) (direction X))))
(set! pml-layers boundary-layers)

(define-param wvl-min 0.4)          ; min wavelength
(define-param wvl-max 0.6)          ; max wavelength
(define fmin (/ wvl-max))           ; min frequency
(define fmax (/ wvl-min))           ; max frequency
(define fcen (* 0.5 (+ fmin fmax))) ; pulse frequency center
(define df (- fmax fmin))           ; pulse frequency width

(define pulse-src (list (make source
                          (src (make gaussian-src (frequency fcen) (fwidth df)))
                          (component Ez)
                          (center (+ (* -0.5 sx) dpml (* 0.5 dsub)) 0 0)
                          (size 0 sy 0))))

(set! sources pulse-src)

(set! k-point (vector3 0 0 0))

(define glass (make medium (index 1.5)))
(set! default-material glass)

(define symm (list (make mirror-sym (direction Y))))
(set! symmetries symm)

(define-param nfreq 21)
(define mon-pt (vector3 (- (* 0.5 sx) dpml (* 0.5 dpad)) 0 0))
(define flux-mon (add-flux fcen df nfreq (make flux-region (center mon-pt) (size 0 sy 0))))

(run-sources+ (stop-when-fields-decayed 50 Ez mon-pt 1e-9))

(define input-flux (get-fluxes flux-mon))
(define freqs (get-flux-freqs flux-mon))

(reset-meep)

(set! geometry-lattice cell)

(set! pml-layers boundary-layers)

(set! sources pulse-src)

(set! k-point (vector3 0 0 0))

(set! default-material air)

(set! geometry (list (make block
                       (material glass)
                       (size (+ dpml dsub) infinity infinity)
                       (center (+ (* -0.5 sx) (* 0.5 (+ dpml dsub))) 0 0))
                     (make block
                       (material glass)
                       (size gh (* gdc gp) infinity)
                       (center (+ (* -0.5 sx) dpml dsub (* 0.5 gh)) 0 0))))

(set! symmetries symm)

(define mode-mon (add-flux fcen df nfreq (make flux-region (center mon-pt) (size 0 sy 0))))

(run-sources+ (stop-when-fields-decayed 50 Ez mon-pt 1e-9))

(define-param nmode 10)
(define res (get-eigenmode-coefficients mode-mon (arith-sequence 1 1 nmode) #:eig-parity (+ ODD-Z EVEN-Y)))
(define coeffs (list-ref res 0))
(define kdom (list-ref res 3))

(map (lambda (nm)
       (map (lambda (nf)
              (let ((mode-wvl (/ (list-ref freqs nf)))
                    (mode-angle (rad->deg (acos (/ (vector3-x (list-ref kdom (+ (* nm nfreq) nf))) (list-ref freqs nf)))))
                    (mode-tran (/ (sqr (magnitude (array-ref coeffs nm nf 0))) (list-ref input-flux nf))))
                (if (> nm 0) (set! mode-tran (* 0.5 mode-tran)))
                (print "grating" (number->string nm) ":, " (number->string mode-wvl) ", " (number->string mode-angle) ", " (number->string mode-tran) "\n")))
              (arith-sequence 0 1 nfreq)))
     (arith-sequence 0 1 nmode))
```

Note the use of the keyword parameter argument `#eig-parity (+ ODD-Z EVEN-Y)` in the call to `get-eigenmode-coefficients`. This is important for specifying **non-degenerate** modes in MPB since the `k-point` is (0,0,0). `ODD-Z` is for modes with E<sub>z</sub> polarization. `EVEN-Y` is necessary since each diffraction order which is based on a given k<sub>x</sub> consists of *two* modes: one going in the +y direction and the other in the -y direction. `EVEN-Y` forces MPB to compute only the +k<sub>y</sub> + -k<sub>y</sub> (cosine) mode. As a result, the total transmittance must be halved in this case to obtain the transmittance for the individual +k<sub>y</sub> or -k<sub>y</sub> mode. For `ODD-Y`, MPB will compute the +k<sub>y</sub> - -k<sub>y</sub> (sine) mode but this will have zero power because the source is even. If the $y$ parity is left out, MPB will return a random superposition of the cosine and sine modes. Alternatively, in this example an input planewave with H<sub>z</sub> instead of E<sub>z</sub> polarization can be used which requires `#eig_parity (+ EVEN-Z ODD-Y)` as well as an odd mirror symmetry plane in *y*. Finally, note the use of `add-flux` instead of `add-mode-monitor` when using symmetries.

The script is run and the output is piped to a file. The results for the diffraction orders, which are displayed in the lines prefixed by `grating`, are extracted from the output and placed into a separate file using the following commands:

```
mpirun -np 2 meep binary_grating.ctl |tee diffraction_spectra.out;
grep grating diffraction_spectra.out |cut -d, -f2- > diffraction_spectra.dat;
```

The diffraction spectrum is then plotted using the following Octave/Matlab script and shown in the figure below.

```matlab
f = dlmread('diffraction_spectrum.dat');
nfreq = 21;
nmode = 10;
wvl = reshape(f(:,1),nfreq,nmode);
ang = reshape(f(:,2),nfreq,nmode);
tran = reshape(f(:,3),nfreq,nmode);

pcolor(wvl,ang,tran);
shading flat;
c = colormap("ocean");
colormap(flipud(c));
colorbar;
caxis([0 0.4]);

xlabel("wavelength (um)");
ylabel("diffraction angle (degrees)");
title("transmittance of diffraction orders");

axis([0.4 0.6 0 30]);
set(gca, 'xtick', [0.4:0.1:0.6]);
set(gca, 'ytick', [0:5:30]);
```

Each diffraction order corresponds to a single angle. In the figure below, this angle is represented by the *lower* boundary of each labeled region. For example, the m=0 order has a diffraction angle of 0° at all wavelengths. The representation of the diffraction orders as finite angular regions is an artifact of Octave/Matlab's [pcolor](https://octave.sourceforge.io/octave/function/pcolor.html) routine. Note that only the positive diffraction orders are shown as these are equivalent to the negative orders due to the symmetry of the source and the structure. The transmittance of each diffraction order should ideally be a constant for all wavelengths. The slight wavelength dependence shown in the figure is due to numerical discretization which can be mitigated by increasing the resolution.

The diffraction orders/modes are a finite set of propagating planewaves. The wavevector k<sub>x</sub> of these modes can be computed analytically: for a frequency of ω (in c=1 units), these propagating modes are the **real** solutions of sqrt(ω²n²-(k<sub>y</sub>+2πm/Λ)²) where m is the diffraction order (an integer), Λ is the periodicity of the grating, and n is the refractive index of the propagating medium. In this example, n=1, k<sub>y</sub>=0, and Λ=10 μm. Thus, at a wavelength of 0.5 μm there are a total of 20 diffraction orders of which we only computed the first 10. The wavevector k<sub>x</sub> is used to compute the angle of the diffraction order as cos<sup>-1</sup>(k<sub>x</sub>/(ωn)). Evanescent modes, those with an imaginary k<sub>x</sub>, exist for |m|>20 but these modes carry no power. Note that currently Meep does not compute the number of propagating modes for you. If the mode number passed to `get-eigenmode-coefficients` is larger than the number of propagating modes at a given frequency/wavelength, MPB's Newton solver will fail to converge and will return zero for the mode coefficient. It is therefore a good idea to know beforehand the number of propagating modes.

<center>
![](../images/grating_diffraction_spectra.png)
</center>

In the limit where the grating periodicity is much larger than the wavelength and the size of the diffracting element (i.e., more than 10 times), as it is in this example, the [diffraction efficiency](https://en.wikipedia.org/wiki/Diffraction_efficiency) can be computed analytically using scalar theory. This is described in the OpenCourseWare [Optics course](https://ocw.mit.edu/courses/mechanical-engineering/2-71-optics-spring-2009/) in the Lecture 16 (Gratings: Amplitude and Phase, Sinusoidal and Binary) [notes](https://ocw.mit.edu/courses/mechanical-engineering/2-71-optics-spring-2009/video-lectures/lecture-16-gratings-amplitude-and-phase-sinusoidal-and-binary/MIT2_71S09_lec16.pdf) and [video](https://www.youtube.com/watch?v=JmWguqCZRxk). For a review of scalar diffraction theory, see Chapter 3 ("Analysis of Two-Dimensional Signals and Systems") of [Introduction to Fourier Optics (fourth edition)](https://www.amazon.com/Introduction-Fourier-Optics-Joseph-Goodman-ebook/dp/B076TBP48F) by J.W. Goodman. From the scalar theory, the diffraction efficiency of the binary grating is 4/(mπ)<sup>2</sup> when the phase difference between the propagating distance in the glass relative to the same distance in air is π. The phase difference/contrast is (2π/λ)(n-1)s where λ is the wavelength, n is the refractive index of the grating, and s is the propagation distance in the grating (`gh` in the script). A special feature of the binary grating is that the diffraction efficiency is 0 for all *even* orders. This is verified by the diffraction spectrum shown above.

To convert the diffraction efficiency into transmittance in the *x* direction (in order to be able to compare the scalar-theory results with those from Meep), the diffraction efficiency must be multiplied by the Fresnel transmittance from air to glass and by the cosine of the diffraction angle. We compare the analytic and simulated results at a wavelength of 0.5 μm for diffraction orders 1, 3, 5, and 7. The analytic results are 0.3886, 0.0427, 0.0151, and 0.0074. The Meep results are 0.3891, 0.04287, 0.0152, and 0.0076. This corresponds to relative errors of approximately 1.3%, 0.4%, 0.8%, and 2.1% which indicates good agreement.

We can also use the complex mode coefficients to compute the phase of the diffraction orders. This can be used to generate a phase map of the binary grating as a function of its geometric parameters. Phase maps are important for the design of subwavelength phase shifters such as those used in metalenses. In this demonstration, which is adapted from the previous example, we compute the transmittance spectra and phase map of the zeroth diffraction order (at 0°) for an E<sub>z</sub>-polarized planewave pulse spanning wavelengths of 0.4 to 0.6 μm which is normally incident on a binary grating with a periodicity of 0.35 μm and height of 0.6 μm. The duty cycle of the grating is varied from 0.1 to 0.9 in separate runs.

The simulation script is in [examples/binary_grating_phasemap.ctl](https://github.com/NanoComp/meep/blob/master/python/examples/binary_grating_phasemap.ctl).

```scm
(set-param! resolution 60)  ; pixels/μm

(define-param dpml 1.0)     ; PML thickness
(define-param dsub 3.0)     ; substrate thickness
(define-param dpad 3.0)     ; padding between grating and PML
(define-param gp 0.35)      ; grating period
(define-param gh 0.6)       ; grating height
(define-param gdc 0.5)      ; grating duty cycle

(define sx (+ dpml dsub gh dpad dpml))
(define sy gp)

(define cell (make lattice (size sx sy no-size)))
(set! geometry-lattice cell)

(define boundary-layers (list (make pml (thickness dpml) (direction X))))
(set! pml-layers boundary-layers)

(define-param wvl-min 0.4)          ; min wavelength
(define-param wvl-max 0.6)          ; max wavelength
(define fmin (/ wvl-max))           ; min frequency
(define fmax (/ wvl-min))           ; max frequency
(define fcen (* 0.5 (+ fmin fmax))) ; pulse frequency center
(define df (- fmax fmin))           ; pulse frequency width
(define-param nfreq 21)             ; number of frequency bins

(define-param odd-z? true)

(define pulse-src (list (make source
                          (src (make gaussian-src (frequency fcen) (fwidth df)))
                          (component (if odd-z? Ez Hz))
                          (center (+ (* -0.5 sx) dpml (* 0.5 dsub)) 0 0)
                          (size 0 sy 0))))

(set! sources pulse-src)

(set! k-point (vector3 0 0 0))

(define glass (make medium (index 1.5)))
(set! default-material glass)

(define symm (list (make mirror-sym (direction Y) (phase (if odd-z? +1 -1)))))
(set! symmetries symm)

(define mon-pt (vector3 (- (* 0.5 sx) dpml (* 0.5 dpad)) 0 0))
(define flux-mon (add-flux fcen df nfreq (make flux-region (center mon-pt) (size 0 sy 0))))

(run-sources+ 100)

(define input-flux (get-fluxes flux-mon))
(define freqs (get-flux-freqs flux-mon))

(reset-meep)

(set! geometry-lattice cell)

(set! pml-layers boundary-layers)

(set! sources pulse-src)

(set! k-point (vector3 0 0 0))

(set! default-material air)

(set! geometry (list (make block
                       (material glass)
                       (size (+ dpml dsub) infinity infinity)
                       (center (+ (* -0.5 sx) (* 0.5 (+ dpml dsub))) 0 0))
                     (make block
                       (material glass)
                       (size gh (* gdc gp) infinity)
                       (center (+ (* -0.5 sx) dpml dsub (* 0.5 gh)) 0 0))))

(set! symmetries symm)

(define mode-mon (add-flux fcen df nfreq (make flux-region (center mon-pt) (size 0 sy 0))))

(run-sources+ 300)

(define res (get-eigenmode-coefficients mode-mon (list 1) #:eig-parity (if odd-z? (+ ODD-Z EVEN-Y) (+ EVEN-Z ODD-Y))))
(define coeffs (list-ref res 0))

(map (lambda (nf)
       (let ((mode-wvl (/ (list-ref freqs nf)))
             (mode-tran (/ (sqr (magnitude (array-ref coeffs 0 nf 0))) (list-ref input-flux nf)))
             (mode-phase (angle (array-ref coeffs 0 nf 0))))
         (print "grating" (number->string nf) ":, " (number->string mode-wvl) ", " (number->string mode-tran) ", " (number->string mode-phase) "\n")))
     (arith-sequence 0 1 nfreq))
```

The phase of the zeroth diffraction order is simply the angle of its complex mode coefficient. Note that it is generally only the relative phase (the phase difference) between different structures that is useful. The overall mode coefficient α is multiplied by a complex number given by the source amplitude, as well as an arbitrary (but deterministic) phase choice by the mode solver MPB — but as long as you keep the current source fixed as you vary the parameters of the structure, the relative phases are meaningful.

The script is run from the shell terminal using the following Bash commands:

```
for gdc in `seq 0.1 0.1 0.9`; do
    meep odd-z?=true gp=0.35 gh=0.6 gdc=${gdc} binary_grating_phasemap.ctl |tee -a phasemap.out;
done

grep grating phasemap.out |cut -d, -f2- > phasemap.dat;
```

The results are plotted using the Octave/Matlab script below. The plot is shown in the accompanying figure.

```matlab
f = dlmread('phasemap.dat');
nfreq = 21;
ngdc = 9;

gdc = ones(nfreq,1)*[0.1:0.1:0.9];

wvl = reshape(f(:,1),nfreq,ngdc);
tran = reshape(f(:,2),nfreq,ngdc);
ang = reshape(f(:,3),nfreq,ngdc);

h1 = subplot(1,2,1);
pcolor(wvl,gdc,tran);
c1 = colormap("hot");
colormap(h1,flipud(c1));
shading interp;
colorbar;
caxis([0 1.0]);

axis([0.4 0.6 0.1 0.9]);
set(gca, 'xtick', [0.4:0.1:0.6]);
set(gca, 'ytick', [0.1:0.1:0.9]);

xlabel("wavelength (um)");
ylabel("grating duty cycle");
title("transmittance");

h2 = subplot(1,2,2);
pcolor(wvl,gdc,ang);
c2 = colormap("jet");
colormap(h2, flipud(c2));
shading interp;
colorbar;

axis([0.4 0.6 0.1 0.9]);
set(gca, 'xtick', [0.4:0.1:0.6]);
set(gca, 'ytick', [0.1:0.1:0.9]);

xlabel("wavelength (um)");
ylabel("grating duty cycle");
title("phase (radians)");
```

The figure below shows the transmittance spectra (left) and phase map (right). The transmittance is nearly unity over most of the parameter space mainly because of the subwavlength dimensions of the grating. The phase variation spans the full range of -π to +π at each wavelength but varies weakly with the duty cycle due to the relatively low index of the glass grating. Higher-index materials such as [titanium dioxide](https://en.wikipedia.org/wiki/Titanium_dioxide#Thin_films) (TiO<sub>2</sub>) generally provide more control over the phase.

<center>
![](../images/grating_phasemap.png)
</center>