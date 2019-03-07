---
# Optical Forces
---

This tutorial demonstrates Meep's ability to compute classical forces via the [Maxwell stress tensor](https://en.wikipedia.org/wiki/Maxwell_stress_tensor) (MST). Also demonstrated is the [eigenmode source](../Scheme_User_Interface.md#eigenmode-source). The geometry consists of two identical, parallel, silicon waveguides with square cross section in vacuum. A schematic of the geometry is shown below. Due to the parallel orientation of the waveguides, the two modes can be chosen to be either symmetric or anti-symmetric with respect to a mirror-symmetry plane between them. As the two waveguides are brought closer and closer together, their modes increasingly couple and give rise to a gradient force that is *transverse* to the waveguide axis. This is different from [radiation pressure](https://en.wikipedia.org/wiki/Radiation_pressure) which involves momentum exchange between photons and is *longitudinal* in nature. An interesting phenomena that occurs for this coupled system is that the force can be tuned to be either attractive or repulsive depending on the relative phase of the modes. This tutorial will demonstrate this effect.

<center>
![](../images/Waveguide_forces.png)
</center>

The gradient force on each waveguide arising from the evanescent coupling of the two waveguide modes can be computed analytically:

$$F=-\frac{1}{\omega}\frac{d\omega}{ds}\Bigg\vert_\vec{k}U,$$

where ω is the mode frequency of the coupled-waveguide system, $s$ is the separation distance between the parallel waveguides, $k$ is the conserved wave vector and $U$ is the total energy of the electromagnetic fields. By convention, negative and positive values correspond to attractive and repulsive forces, respectively. For more details, see [Optics Letters, Vol. 30, pp. 3042-4, 2005](https://www.osapublishing.org/ol/abstract.cfm?uri=ol-30-22-3042). This expression has been shown to be mathematically equivalent to the MST in [Optics Express, Vol. 17, pp. 18116-35, 2009](http://www.opticsinfobase.org/oe/abstract.cfm?URI=oe-17-20-18116). We will verify this result in this tutorial.

It is convenient to normalize the force in order to work with dimensionless quantities. Since the total power transmitted through the waveguide is $P=v_gU/L$ where $v_g$ is the group velocity, $L$ is the waveguide length, and $U$ is defined as before, we focus instead on the force per unit length per unit power $(F/L)(ac/P)$ where $a$ is the waveguide width and $c$ is the speed of light. This dimensionless quantity enables us to compute both the flux and the force in a single simulation.

We can compute the gradient force using two different methods and verify that they are equivalent: (1) using MPB, we compute the frequency and group velocity for a given mode over a range of separation distances and then use a [finite-difference](https://en.wikipedia.org/wiki/Finite_difference) scheme to numerically evaluate the formula from above, and (2) using Meep, we directly compute both the gradient force and the power transmitted through the waveguide for the guided mode over the same range of separation distances. In this particular example, we consider just the fundamental `ODD-Y` mode which shows the bidirectional force. The range of separation distances is from 0.02 to 1.02 μm in increments of 0.02 μm.

The simulation script is in [examples/parallel-wvgs-force.ctl](https://github.com/NanoComp/meep/blob/master/scheme/examples/parallel-wvgs-force.ctl).


```scm
(set-param! resolution 30)   ; pixels/um

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

Since we do not know apriori what the mode frequency is for a given waveguide separation distance, a preliminary run is required to find this out using [`Harminv`](../Scheme_User_Interface.md#harminv) and a broadband pulsed source. Since the propagating mode never decays away, the runtime is chosen arbitrarily as 200 time units after the pulsed sources have turned off.

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

Once we have determined the mode frequency, we then replace the `source` with [`eigenmode-source`](../Scheme_User_Interface.md#eigenmodesource) to perform the main simulation: compute (1) the force on each waveguide due to the mode coupling and (2) the power in the mode. The `eigenmode-source` invokes [MPB](https://mpb.readthedocs.io) to compute the given mode of interest. The mode profile is then imported into Meep for use as the initial source amplitude. This enables a more efficient mode excitation than simply using a point or area source with constant amplitude. For more details on the eigenmode source feature, refer to Section 4.2 ("Incident Fields and Equivalent Currents") in [Chapter 4](http://arxiv.org/abs/arXiv:1301.5366) ("Electromagnetic Wave Source Conditions") of the book [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707).

```scm
(reset-meep)
(change-sources! (list
		  (make eigenmode-source
		    (src (make gaussian-src (frequency f) (fwidth df)))
		    (size a a 0)
		    (center (* -0.5 (+ s a)) 0)
		    (eig-kpoint k-point)
		    (eig-match-freq? true)
		    (eig-parity ODD-Y))
		  (make eigenmode-source
		    (src (make gaussian-src (frequency f) (fwidth df)))
		    (size a a 0)
		    (center (* 0.5 (+ s a)) 0)
		    (eig-kpoint k-point)
		    (eig-match-freq? true)
		    (eig-parity ODD-Y)
		    (amplitude (if xodd? -1.0 1.0)))))

(define wvg-flux (add-flux f 0 1
			  (make flux-region (direction Z) (center 0 0)
				(size (* 1.2 (+ (* 2 a) s)) (* 1.2 a) 0))))
(define wvg-force (add-force f 0 1
			     (make force-region (direction X) (weight +1.0)
				   (center (* 0.5 s) 0) (size 0 a))
			     (make force-region (direction X) (weight -1.0)
				   (center (+ (* 0.5 s) a) 0) (size 0 a))))

(run-sources+ 5000)

(display-fluxes wvg-flux)
(display-forces wvg-force)
```

There are two important items to note in the script: (1) We have defined a single flux surface to compute the Poynting flux in $z$ which spans an area slightly larger than both waveguides rather than two separate flux surfaces (one for each waveguide). This is because in the limit of small separation, two flux surfaces overlap whereas the total power through a single flux surface need, by symmetry, only be halved in order to determine the value for just one of the two waveguides. (2) Instead of defining a closed, four-sided "box" surrounding the waveguides for computing the MST, we chose instead to compute the MST along just two $y$-oriented lines (to obtain the force in the $x$ direction) with different `weight` values to correctly sum the total force. By symmetry, we need not consider the force in the $y$ direction. Choosing a suitable runtime requires some care. A large runtime is necessary to obtain the steady-state response but this will also lead to large values for the discrete Fourier-transformed fields used to compute both the flux and the MST. These large values may contain [roundoff errors](https://en.wikipedia.org/wiki/Round-off_error).

We run this simulation over the range of separation distances and compare the results to those obtained from MPB. This is shown in the figure above. The two methods show good agreement.

The MPB simulation is in [examples/parallel-wvgs-mpb.ctl](https://github.com/NanoComp/meep/blob/master/scheme/examples/parallel-wvgs-mpb.ctl). Note: since MPB permits symmetries only in the $y$ and $z$ directions, the coordinate axes used in the MPB script to define the waveguide geometry are different than those in the Meep script. In MPB, the propagating axis is $x$ whereas in Meep it is $z$.

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