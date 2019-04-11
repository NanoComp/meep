---
# Near to Far Field Spectra
---

We demonstrate Meep's [near-to-far-field transformation](../Scheme_User_Interface.md#near-to-far-field-spectra) feature using two examples. There are three steps to using the near-to-far-field feature. First, we need to define the "near" surface(s) as a set of surfaces capturing *all* outgoing radiation in the desired direction(s). Second, we run the simulation using a pulsed source (or possibly, the frequency-domain solver) to allow Meep to accumulate the Fourier transforms on the near surface(s). Third, we have Meep compute the far fields at any desired points with the option to save the far fields to an HDF5 file.

[TOC]

Radiation Pattern of an Antenna
-------------------------------

In this example, we compute the [radiation pattern](https://en.wikipedia.org/wiki/Radiation_pattern) of an antenna. This involves an electric-current point dipole source as the emitter in vacuum. The source is placed at the center of a 2d square cell surrounded by PML. The near fields are obtained on a bounding box defined along the edges of the non-PML region. The far fields are computed in two ways: along the (1) sides of a square box and (2) circumference of a circle, having a length/radius many times larger than the source wavelength and lying beyond the cell. From both the near and far fields, we will also compute the total outgoing Poynting flux and demonstrate that they are equivalent. Results will be shown for three orthogonal polarizations of the input source.

The simulation geometry is shown in the following schematic.

<center>
![](../images/Near2far_simulation_geometry.png)
</center>

In the first part of the simulation, we define the cell and sources as well as the near field and flux regions. Since we are using a pulsed source (with center wavelength of 1 μm), the fields are timestepped until they have sufficiently decayed away.

The simulation script is in [examples/antenna-radiation.ctl](https://github.com/NanoComp/meep/blob/master/scheme/examples/antenna-radiation.ctl).

```scm
(set-param! resolution 50)   ; pixels/μm

(define-param sxy 4)
(define-param dpml 1)
(set! geometry-lattice (make lattice (size (+ sxy (* 2 dpml)) (+ sxy (* 2 dpml)) no-size)))

(set! pml-layers (list (make pml (thickness dpml))))

(define-param fcen 1.0)
(define-param df 0.4)
(define-param src-cmpt Ez)
(set! sources (list (make source
                      (src (make gaussian-src (frequency fcen) (fwidth df)))
                      (center 0)
                      (component src-cmpt))))

(if (= src-cmpt Ex)
    (set! symmetries (list (make mirror-sym (direction X) (phase -1))
                           (make mirror-sym (direction Y) (phase +1)))))
(if (= src-cmpt Ey)
    (set! symmetries (list (make mirror-sym (direction X) (phase +1))
                           (make mirror-sym (direction Y) (phase -1)))))
(if (= src-cmpt Ez)
    (set! symmetries (list (make mirror-sym (direction X) (phase +1))
                           (make mirror-sym (direction Y) (phase +1)))))

(define nearfield-box
  (add-near2far fcen 0 1
		(make near2far-region (center 0 (* 0.5 sxy)) (size sxy 0))
		(make near2far-region (center 0 (* -0.5 sxy)) (size sxy 0) (weight -1))
		(make near2far-region (center (* 0.5 sxy) 0) (size 0 sxy))
		(make near2far-region (center (* -0.5 sxy) 0) (size 0 sxy) (weight -1))))

(define flux-box
  (add-flux fcen 0 1
	    (make flux-region (center 0 (* 0.5 sxy)) (size sxy 0))
	    (make flux-region (center 0 (* -0.5 sxy)) (size sxy 0) (weight -1))
	    (make flux-region (center (* 0.5 sxy) 0) (size 0 sxy))
	    (make flux-region (center (* -0.5 sxy) 0) (size 0 sxy) (weight -1))))

(run-sources+ (stop-when-fields-decayed 50 src-cmpt (vector3 0 0) 1e-8))
```

After the time stepping, the flux of the near fields is computed using `get-fluxes` and displayed:

```scm
(print "near-flux:, " (list-ref (get-fluxes flux-box) 0) "\n")
```

In the first of two cases, the flux of the far fields is computed using the `flux` routine for a square box of side length 2 mm which is 2000 times larger than the source wavelength. This requires computing the outgoing flux on each of the four sides of the box separately and summing the values. The resolution of the far fields is chosen arbitrarily as 1 point/μm. This means there are 2x10<sup>6</sup> points per side length.

```scm
(define-param r (/ 1000 fcen)) ; half side length of far-field square box OR radius of far-field circle
(define-param res-ff 1)        ; resolution of far fields (points/μm)
(define far-flux (+ (list-ref (flux nearfield-box Y (volume (center 0 r 0) (size (* 2 r) 0 0)) res-ff) 0)
                    (- (list-ref (flux nearfield-box Y (volume (center 0 (- r) 0) (size (* 2 r) 0 0)) res-ff) 0))
                    (list-ref (flux nearfield-box X (volume (center r 0 0) (size 0 (* 2 r) 0)) res-ff) 0)
                    (- (list-ref (flux nearfield-box X (volume (center (- r) 0 0) (size 0 (* 2 r) 0)) res-ff) 0))))

(print "far-flux-box:, " far-flux "\n")
```

For the second of two cases, we use the `get-farfield` routine to compute the far fields by looping over a set of 100 equally-spaced points along the circumference of a circle with radius of 1 mm. The six far field components (E$_x$, E$_y$, E$_z$, H$_x$, H$_y$, H$_z$) are displayed in separate columns as complex numbers.

```scm
(define-param npts 100)           ; number of points in [0,2*pi) range of angles
(map (lambda (n)
       (let ((ff (get-farfield nearfield-box (vector3 (* r (cos (* 2 pi (/ n npts)))) (* r (sin (* 2 pi (/ n npts)))) 0))))
	 (print "farfield:, " n ", " (* 2 pi (/ n npts)))
	 (map (lambda (m)
		(print ", " (list-ref ff m)))
	      (arith-sequence 0 1 6))
	 (print "\n")))
         (arith-sequence 0 1 npts))
```

The script is run and the output piped to a file using the following shell commands. The far field data is extracted from the output and placed in a separate file.

```sh
meep src-cmpt=Ez antenna-radiation.ctl |tee source_Jz_farfields.out
grep farfield: source_Jz_farfields.out |cut -d , -f2- > source_Jz_farfields.dat
```

From the far fields at each point $\mathbf{r}$, we compute using Matlab/Octave the outgoing or radial flux: $\sqrt{P_x^2+P_y^2}$, where P$_x$ and P$_y$ are the components of the Poynting vector $\mathbf{P}(\mathbf{r})=(P_x,P_y,P_z)=\mathrm{Re}\, \mathbf{E}(\mathbf{r})^*\times\mathbf{H}(\mathbf{r})$. Note that $P_z$ is always 0 since this is a 2d simulation. The total flux is then computed and displayed:

```matlab
d = dlmread("source_Jz_farfields.dat",",");

Ex = conj(d(:,3));
Ey = conj(d(:,4));
Ez = conj(d(:,5));

Hx = d(:,6);
Hy = d(:,7);
Hz = d(:,8);

Px = real(Ey.*Hz-Ez.*Hy);
Py = real(Ez.*Hx-Ex.*Hz);
Pr = sqrt(Px.^2+Py.^2);

r = 1000;    % radius of far-field circle
disp(sprintf("far-flux-circle:, %0.6f",sum(Pr)*2*pi*r/length(Pr)));
```

By [Poynting's theorem](https://en.wikipedia.org/wiki/Poynting%27s_theorem), the total outgoing flux obtained by integrating around a *closed* surface should be the same whether it is calculated from the near or far fields (unless there are sources or absorbers in between). The flux of the near fields for the J$_z$ source is `2.456196` and that for the far fields is `2.458030` (box) and `2.457249` (circle). The ratio of near- to far-field (circle) flux is `0.999571`. Similarly, for the J$_x$ source, the values are `1.227786` (near-field), `1.227651` (far-field box), and `1.227260` (far-field circle). The ratio of near- to far-field (circle) flux is `1.000429`. The slight differences in the flux values are due to discretization effects and will decrease as the resolution is increased.

Finally, we plot the radial flux normalized by its maximum value over the entire interval to obtain a range of values between 0 and 1. These are shown below in the linearly-scaled, polar-coordinate plots. The three figures are obtained using separate runs involving a `src-cmpt` of E$_x$, E$_y$, and E$_z$. As expected, the J$_x$ and J$_y$ sources produce [dipole](https://en.wikipedia.org/wiki/Electric_dipole_moment) radiation patterns while J$_z$ has a monopole pattern.

```matlab
angles = d(:,2);
polar(angles,Pr/max(Pr),'b-');
set(gca, 'xtick', [0 0.5 1.0]);
```

<center>
![](../images/Source_radiation_pattern.png)
</center>

Far-Field Profile of a Cavity
-----------------------------

For this demonstration, we will compute the far-field spectra of a resonant cavity mode in a holey waveguide; a structure we had explored in [Tutorial/Resonant Modes and Transmission in a Waveguide Cavity](Resonant_Modes_and_Transmission_in_a_Waveguide_Cavity.md). The script is in [examples/cavity-farfield.ctl](https://github.com/NanoComp/meep/blob/master/scheme/examples/cavity-farfield.ctl). The structure is shown at the bottom of the left image below.

![center|Schematic of the computational cell for a holey waveguide with cavity showing the location of the "near" boundary surface and the far-field region.](../images/N2ff_comp_cell.png)

To set this up, we simply remove the last portion of [examples/holey-wvg-cavity.ctl](https://github.com/NanoComp/meep/blob/master/scheme/examples/holey-wvg-cavity.ctl), beginning right after the line:

```scm
(set! symmetries
      (list (make mirror-sym (direction Y) (phase -1))
            (make mirror-sym (direction X) (phase -1))))
```

and insert the following lines:

```scm
(define-param d1 0.2)
(define nearfield
         (add-near2far fcen 0 1
           (make near2far-region (center 0 (+ (* 0.5 w) d1))
                                 (size (- sx (* 2 dpml)) 0))
           (make near2far-region (center (+ (* -0.5 sx) dpml) (+ (* 0.5 w) (* 0.5 d1)))
                                 (size 0 d1)
                                 (weight -1.0))
           (make near2far-region (center (- (* 0.5 sx) dpml) (+ (* 0.5 w) (* 0.5 d1)))
                                 (size 0 d1))))
```

We are creating a "near" bounding surface, consisting of three separate regions surrounding the cavity, that captures <i>all</i> outgoing waves in the top-half of the cell. Note that the *x*-normal surface on the left has a `weight` of -1 corresponding to the direction of the *outward normal* vector relative to the *x* direction so that the far-field spectra is correctly computed from the outgoing fields, similar to the flux and force features. The parameter `d1` is the distance between the edge of the waveguide and the bounding surface, as shown in the schematic above, and we will demonstrate that changing this parameter does not change the far-field spectra which we compute at a single frequency corresponding to the cavity mode.

We then time step the fields until, at a random point, they have sufficiently decayed away as the cell is surrounded by PMLs, and output the far-field spectra over a rectangular area that lies <i>outside</i> of the cell:

```scm
(run-sources+ (stop-when-fields-decayed 50 Hz (vector3 0.12 -0.37) 1e-8))
(define-param d2 20)
(define-param h 4)
(output-farfields nearfield
 (string-append "spectra-" (number->string d1) "-" (number->string d2) "-" (number->string h))
 (volume (center 0 (+ (* 0.5 w) d2 (* 0.5 h))) (size (- sx (* 2 dpml)) h)) resolution)
```

The first item to note is that the far-field region is located <i>outside</i> of the cell, although in principle it can be located anywhere. The second is that the far-field spectra can be interpolated onto a spatial grid that has any given resolution but in this example we used the same resolution as the simulation. Note that the simulation itself used purely real fields but the output, given its analytical nature, contains complex fields. Finally, given that the far-field spectra is derived from the Fourier-transformed fields which includes an arbitrary constant factor, we should expect an overall scale and phase difference in the results obtained using the near-to-far-field feature with those from a corresponding simulation involving the full computational volume. The key point is that the results will be qualitatively but not quantitatively identical. The data will be written out to an HDF5 file having a filename prefix with the values of the three main parameters. This file will includes the far-field spectra for all six field components, including real and imaginary parts.

We run the above modified control file and in post-processing create an image of the real and imaginary parts of H$_z$ over the far-field region which is shown in insets (a) above. For comparison, we compute the steady-state fields using a larger cell that contains within it the far-field region. This involves a continuous source and complex fields. Results are shown in figure (b) above. The difference in the relative phases among any two points within each of the two field spectra is zero, which can be confirmed numerically. Also, as would be expected, it can be shown that increasing `d1` does not change the far-field spectra as long as the results are sufficiently converged. This indicates that discretization effects are irrelevant.

In general, it is tricky to interpret the overall scale and phase of the far fields, because it is related to the scaling of the Fourier transforms of the near fields. It is simplest to use the `near2far` feature in situations where the overall scaling is irrelevant, e.g. when you are computing a ratio of fields in two simulations, or a fraction of the far field in some region, etcetera.