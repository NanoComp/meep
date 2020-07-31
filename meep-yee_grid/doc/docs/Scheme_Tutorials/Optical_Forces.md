---
# Optical Forces
---

This tutorial demonstrates Meep's ability to compute classical forces via the [Maxwell stress tensor](https://en.wikipedia.org/wiki/Maxwell_stress_tensor). The geometry consists of two identical, parallel, silicon waveguides with square cross section in vacuum. A schematic of the geometry is shown below. Due to the parallel orientation of the waveguides (propagation axis along $z$ and separation in $x$), the two modes can be chosen to be either symmetric or anti-symmetric with respect to an $x$ mirror-symmetry plane between them. As the two waveguides are brought closer and closer together, their modes couple and give rise to a gradient force that is *transverse* to the waveguide axis (i.e., in the $x$ direction). This is different from [radiation pressure](https://en.wikipedia.org/wiki/Radiation_pressure) which involves momentum exchange between photons and is longitudinal in nature (i.e., along the $z$ direction). An interesting phenomena that occurs for this coupled waveguide system is that the force can be tuned to be either attractive or repulsive depending on the relative phase of the two modes. This tutorial will demonstrate this effect.

<center>
![](../images/Waveguide_forces.png)
</center>

The gradient force $F$ on each waveguide arising from the evanescent coupling of the two waveguide modes can be computed analytically:

$$F=-\frac{1}{\omega}\frac{d\omega}{ds}\Bigg\vert_\vec{k}U$$

where $\omega$ is the mode frequency of the coupled waveguide system, $s$ is the separation distance between the parallel waveguides, $k$ is the conserved wave vector, and $U$ is the total energy of the electromagnetic fields. By convention, negative and positive values correspond to attractive and repulsive forces, respectively. For more details, see [Optics Letters, Vol. 30, pp. 3042-4, 2005](https://www.osapublishing.org/ol/abstract.cfm?uri=ol-30-22-3042). This expression has been shown to be mathematically equivalent to the Maxwell stress tensor in [Optics Express, Vol. 17, pp. 18116-35, 2009](http://www.opticsinfobase.org/oe/abstract.cfm?URI=oe-17-20-18116). We will verify this result in this tutorial. In this particular example, only the fundamental mode with odd mirror-symmetry in $y$ shows the bidirectional force.

It is convenient to [normalize the force](../Scheme_User_Interface.md#force-spectra) in order to work with [dimensionless quantities](../Introduction.md#units-in-meep). Since the total transmitted power in the waveguide per unit length is $P=v_gU/L$ where $v_g$ is the group velocity, $U$ is the total energy of the electromagnetic fields (same as before), and $L$ is the waveguide length, we focus instead on the force per unit length per unit power $(F/L)(ac/P)$ where $a$ is the waveguide width/height and $c$ is the speed of light. This dimensionless quantity can be computed in a single simulation.

The gradient force $F$ can be computed using two different methods: (1) using MPB, compute the frequency $\omega$ and group velocity $v_g$ for a given mode over a range of separation distances $s$ and then use a centered [finite-difference](https://en.wikipedia.org/wiki/Finite_difference) scheme to evaluate $F$ using the formula from above, and (2) using Meep, directly compute both the gradient force $F$ and the transmitted power $P$ over the same range of separation distances $s$. This tutorial verifies that (1) and (2) produce equivalent results.

The simulation script is in [examples/parallel-wvgs-force.ctl](https://github.com/NanoComp/meep/blob/master/scheme/examples/parallel-wvgs-force.ctl).


```scm
(set-param! resolution 40)   ; pixels/um

(define Si (make medium (index 3.45)))

(define-param dpml 1.0)
(set! pml-layers (list (make pml (thickness dpml))))

(define-param sx 5)
(define-param sy 3)
(set! geometry-lattice
      (make lattice (size (+ sx (* 2 dpml)) (+ sy (* 2 dpml)) no-size)))

(define-param a 1.0)  ; waveguide width
(define-param s 1.0)  ; waveguide separation distance
(set! geometry (list
		(make block (center (* -0.5 (+ s a)) 0)
		      (size a a infinity) (material Si))
		(make block (center (* 0.5 (+ s a)) 0)
		      (size a a infinity) (material Si))))
```

Two mirror symmetries can be used to reduce the size of the computational cell by a factor of four:

```scm
(define-param xodd? true)
(set! symmetries (list
		  (make mirror-sym (direction X) (phase (if xodd? -1 +1)))
		  (make mirror-sym (direction Y) (phase -1))))
```

Next, we set the Bloch-periodic boundary condition for the mode with wavevector π/$a$:

```scm
(set! k-point (vector3 0 0 0.5))
```

Since we do not know apriori what the mode frequency is for a given waveguide separation distance, a preliminary run is required to find this out using [`Harminv`](../Scheme_User_Interface.md#harminv) and a broadband pulsed source. Since the guided mode never decays away, the runtime is chosen arbitrarily as 200 time units after the pulsed sources have turned off.

```scm
(define-param fcen 0.22)
(define-param df 0.06)
(set! sources (list
	       (make source (src (make gaussian-src (frequency fcen) (fwidth df)))
		     (component Ey) (center (* -0.5 (+ s a)) 0) (size a a 0))
	       (make source (src (make gaussian-src (frequency fcen) (fwidth df)))
		     (component Ey) (center (* 0.5 (+ s a)) 0) (size a a 0)
		     (amplitude (if xodd? -1.0 1.0)))))

(run-sources+ 200
	      (after-sources (harminv Ey (vector3 (* 0.5 (+ s a)) 0) fcen df)))

(define f (harminv-freq-re (car harminv-results)))
(print "freq:, " s ", " f "\n")
```

Once we have determined the mode frequency, we then replace the `source` with [`eigenmode-source`](../Scheme_User_Interface.md#eigenmode-source) to perform the main simulation: compute (1) the force on each waveguide due to the mode coupling and (2) the power in the mode. An `eigenmode-source` with `eig-match-freq?` set to `false` is then used to launch the guided mode using a pulse (the `frequency` but not the `fwidth` parameter of `gaussian-src` is ignored). Alternatively, a constant-amplitude point/area source can be used to launch the mode but this is less efficient as demonstrated in [Tutorial/Eigenmode Source/Index-Guided Modes in a Ridge Waveguide](Eigenmode_Source.md#index-guided-modes-in-a-ridge-waveguide).

```scm
(reset-meep)
(change-sources! (list
                  (make eigenmode-source
                    (src (make gaussian-src (frequency f) (fwidth df)))
                    (size sx sy 0)
                    (center 0 0 0)
                    (eig-band (if xodd? 2 1))
                    (eig-kpoint k-point)
                    (eig-match-freq? false)
                    (eig-parity ODD-Y))))

(define wvg-flux (add-flux f 0 1
                           (make flux-region (direction Z)
                                 (center 0 0) (size sx sy 0))))

(define wvg-force (add-force f 0 1
			     (make force-region (direction X) (weight +1)
				   (center (* 0.5 s) 0) (size 0 sy))
			     (make force-region (direction X) (weight -1)
				   (center (+ (* 0.5 s) a) 0) (size 0 sy))))

(run-sources+ 1500)

(display-fluxes wvg-flux)
(display-forces wvg-force)
```

There are four important items to note: (1) a single flux surface is used to compute the Poynting flux in $z$ which spans the entire non-PML region of the cell. This is because in the limit of small waveguide separation distance, two separate flux surfaces for each waveguide would overlap and result in overcounting. The total power through a single flux surface need, by symmetry, only be halved in order to determine the value for a single waveguide. (2) Instead of defining a closed, four-sided "box" surrounding the waveguides, the Maxwell stress tensor is computed using just two line monitors oriented in $y$ (to obtain the force in the perpendicular $x$ direction) with `weight` values of `+1`/`-1` to correctly sum the total force. The force monitors are placed in the vacuum region adjacent to the waveguide rather than on its surface so that the fields are [second-order accurate](../Subpixel_Smoothing.md). By symmetry, the force in the $y$ direction is zero and need not be computed. (3) Since the `eig-parity` parameter of `eigenmode-source` can only be specified using the $y$ and/or $z$ directions (but *not* $x$, the waveguide separation axis), the `eig_band` parameter must be set to `1`/`2` to distinguish modes with even/odd $x$-mirror symmetry.  (4) A 2d `cell-size` in $xy$ combined with a `k-point` containing a non-zero $z$ component results in a [2d simulation (which is the default)](../2d_Cell_Special_kz.md).

In this example, the fields of the guided mode never decay away to zero. [Choosing a runtime](../FAQ.md#checking-convergence) is therefore somewhat arbitrary but requires some care. A sufficiently long runtime is necessary to obtain the [steady-state response](../FAQ.md#how-do-i-compute-the-steady-state-fields). However, an excessively long runtime will lead to large values for the Fourier-transformed fields used to compute both the flux and the Maxwell stress tensor. Large floating-point numbers may contain [roundoff errors](https://en.wikipedia.org/wiki/Round-off_error) and produce inaccurate results.

The simulation is run over the range of separation distances $s$ from 0.05 to 1.00 μm in increments of 0.05 μm. The results are compared with those from MPB. This is shown in the top figure. The two methods show good agreement.

The following figure shows the $E_y$ mode profiles at a waveguide separation distance of 0.1 μm. This figure shows the source and flux monitor (red hatches), force monitors (blue lines), and PMLs (green hatches) surrounding the cell. From the force spectra shown above, at this separation distance the anti-symmetric mode is repulsive whereas the symmetric mode is attractive.

<center>
![](../images/parallel_wvgs_s0.1.png)
</center>

The MPB simulation is in [examples/parallel-wvgs-mpb.ctl](https://github.com/NanoComp/meep/blob/master/scheme/examples/parallel-wvgs-mpb.ctl). There are important differences related to the coordinate dimensions between the MPB and Meep scripts. In the MPB script, the 2d cell is defined using the $yz$ plane, the waveguide propagation axis is $x$, and the waveguide separation axis is $y$. As a consequence, the `num-bands` parameter is always `1` since the $y$ parity of the mode can be defined explicitly (i.e., `run-yodd-zodd` vs. `run-yeven-zodd`). This is different from the Meep script since Meep requires that a 2d cell be defined in the $xy$ plane. MPB has no such requirement.

```scm
(set-param! resolution 128)  ; pixels/μm

(define Si (make dielectric (index 3.45)))

(define-param syz 10)
(set! geometry-lattice (make lattice (size no-size syz syz)))

(define-param a 1.0) ; waveguide width
(define-param s 1.0) ; waveguide separation distance

(set! geometry (list
		(make block (center 0 (* -0.5 (+ s a)) 0)
		      (size infinity a a) (material Si))
		(make block (center 0 (* 0.5 (+ s a)) 0)
		      (size infinity a a) (material Si))))

(set! k-points (list (vector3 0.5 0 0)))

(set-param! num-bands 1)
(set-param! tolerance 1e-9)

(define-param yodd? true)
(if yodd? (run-yodd-zodd) (run-yeven-zodd))

(print "data:, " s ", " (list-ref freqs 0) ", "
       (list-ref (compute-group-velocity-component (vector3 1 0 0)) 0) "\n")
```

The shell script below runs the MPB simulation for each of the two symmetry configurations over a range of waveguide separate distances and pipes the results to a file. The flux and force data are extracted from the output and placed in a separate file.

```sh
#!/bin/bash

for s in `seq 0.02 0.02 1.02`; do
    mpb s=${s} yodd?=true parallel-wvgs-mpb.ctl >> parallel-wvgs-yodd.out;
    mpb s=${s} yodd?=false parallel-wvgs-mpb.ctl >> parallel-wvgs-yeven.out;
done;

grep data: parallel-wvgs-yodd.out |cut -d , -f2- > parallel-wvgs-yodd.dat;
grep data: parallel-wvgs-yeven.out |cut -d , -f2- > parallel-wvgs-yeven.dat;
```

The following Octave/Matlab script computes the gradient force based on the simulation output and plots the results in a figure.

```matlab
odd_data = dlmread("parallel-wvgs-yodd.dat",',');
force_odd = compute_force(odd_data);
s_odd = odd_data(1:end-1,1);

even_data = dlmread("parallel-wvgs-yeven.dat",',');
force_even = compute_force(even_data);
s_even = even_data(1:end-1,1);

plot(s_odd,force_odd,'b-',s_even,force_even,'r-');
xlabel("waveguide separation s/a");
ylabel("optical force (F/L)(ac/P)");
legend("antisymmetric","symmetric");

function force = compute_force(data)
  f_avg = 0.5*(data(1:end-1,2)+data(2:end,2));
  df = data(2:end,2)-data(1:end-1,2);
  vg_avg = 0.5*(data(1:end-1,3)+data(2:end,3));
  ds = data(2,1)-data(1,1);
  force =  -1./f_avg .* df/ds .* 1./vg_avg;
  return
endfunction
```