---
# Optical Forces
---

This tutorial demonstrates Meep's ability to compute classical forces using the [Maxwell stress tensor](https://en.wikipedia.org/wiki/Maxwell_stress_tensor) (MST) as well as the eigenmode source feature which involves the mode-solver [MPB](https://mpb.readthedocs.io). See [parallel-wvgs-force.ctl](https://github.com/stevengj/meep/blob/master/scheme/examples/parallel-wvgs-force.ctl). The geometry consists of two identical silicon waveguides with square cross section in vacuum (shown in the figure below). Due to the orientation of the waveguides, the two modes can be chosen to be either symmetric or anti-symmetric with respect to a mirror-symmetry plane between them. As the two waveguides are brought closer and closer together, their modes couple more and more and give rise to an optical gradient force that is transverse to the waveguide axis. This is different from [radiation pressure](https://en.wikipedia.org/wiki/Radiation_pressure) which involves momentum exchange between photons and is longitudinal in nature. An interesting phenomena that occurs for this system is that the force can be tuned to be either attractive or repulsive depending on the relative phase of the modes. We will demonstrate this effect in this tutorial.

The optical gradient force on each waveguide arising from the evanescent coupling of the waveguide modes can also be computed analytically:

$$F=-\frac{1}{\omega}\frac{d\omega}{ds}\Bigg\vert_\vec{k}U,$$

where ω is the eigenmode frequency of the coupled-waveguide system, $s$ is the separation distance between the parallel waveguides, $k$ is the conserved wave vector and $U$ is the total energy of the electromagnetic fields. By convention, negative and positive values correspond to attractive and repulsive forces, respectively. For more details, see [Optics Letters, vol. 30, issue 22, pp. 3042-4 (2005)](http://math.mit.edu/~stevenj/papers/PovinelliLo05.pdf). This expression has been shown to be mathematically equivalent to the MST ([Optics Express, vol. 17, issue 20, pp. 18116-135](http://www.opticsinfobase.org/oe/abstract.cfm?URI=oe-17-20-18116)). We will verify this result in this tutorial.

It is convenient to normalize the force so as to eliminate the tricky units altogether. Since the total power transmitted through the waveguide is $P=v_gU/L$ where $v_g$ is the group velocity, $L$ is the waveguide length and $U$ is defined as before, we focus instead on the force per unit length and power $(F/L)(ac/P)$ where $a$ is an arbitrary unit length and $c$ is the speed of light. This dimensionless quantity enables us to compute both the flux and the force in a single simulation.

We can compute the optical gradient force in two ways: (1) using MPB, we compute the eigenfrequency and group velocity for a given mode at different separation distances and then use a finite-difference scheme to numerically evaluate the formula from above, and (2) using Meep, we compute both the MST and the power transmitted through the waveguide for the guided mode. In this particular example, we consider just the fundamental `y-odd` mode which shows the bi-directional force.

First, let's set up the computational cell. The waveguide cross section is 2d but since we are interested in a guided mode of the structure which requires specifying the axial wavevector, this will in fact be a 3d simulation:

```scm
 (set-param! resolution 30)
 (define-param nSi 3.45)
 (define Si (make medium (index nSi)))
 (define-param dpml 1.0)
 (define-param sx 5)
 (define-param sy 3)
 (set! geometry-lattice 
     (make lattice (size (+ sx (* 2 dpml)) (+ sy (* 2 dpml)) no-size)))
 (set! pml-layers (list (make pml (thickness dpml))))
 (define-param a 1.0)  ; waveguide width
 (define-param s 1.0)  ; waveguide separation distance
 (set! geometry (list 
       (make block (center (* -0.5 (+ s a)) 0)
             (size a a infinity) (material Si))
       (make block (center (*  0.5 (+ s a)) 0)
             (size a a infinity) (material Si))))
```

Two mirror symmetries can be used to reduce the size of the computational cell by a factor of four:

```scm
 (define-param xodd? true)
 (set! symmetries (list 
        (make mirror-sym (direction X) (phase (if xodd? -1 +1)))
        (make mirror-sym (direction Y) (phase -1))))
```

Next, we set the Bloch-periodic boundary condition for the mode with wavevector π/$a$:

```scm
 (define-param beta 0.5)
 (set! k-point (vector3 0 0 beta))
```

Since we do not know apriori what the eigenfrequency is at a given separation distance, we use a broadband source to excite multiple frequencies and then use `harminv` to determine the resonant frequency. The use of Bloch-periodic boundary conditions in the $Z$ direction means that the excited mode (if any) will propagate indefinitely in time which is why we stop the simulation at 200 time units after the sources have turned off. This number is somewhat arbitrary.

```scm
 (define-param fcen 0.22)
 (define-param df 0.06)
 (set! sources (list 
         (make source (src (make gaussian-src (frequency fcen) (fwidth df))) 
                (component Ey) (center (* -0.5 (+ s a)) 0) (size a a 0))
         (make source (src (make gaussian-src (frequency fcen) (fwidth df))) 
                (component Ey) (center (* 0.5 (+ s a)) 0) (size a a 0) 
                (amplitude (if xodd? -1.0 1.0)))))
 (run-sources+ 200 
     (after-sources (harminv Ey (vector3 (* 0.5 (+ s a)) 0) fcen df)))
 (define f (harminv-freq-re (car harminv-results)))
 (print "freq:, " s ", " f "\n")
```

Once we have determined the eigenmode frequency, we then excite this mode using the `eigenmode-source` feature. Also, we compute the force on each waveguide and the power in the guided mode. The `eigenmode-mode` feature invokes [MPB](https://mpb.readthedocs.io) (which needs to be pre-installed as a library) to compute the relevant mode of interest. The steady-state field profile from MPB is imported into Meep for use as the initial source amplitude. This enables an efficient excitation of the eigenmode to a much higher degree of accuracy than would otherwise be possible had we simply used a point-dipole source. For more details, refer to Section 4.2 ("Incident Fields and Equivalent Currents") in [Chapter 4](http://arxiv.org/abs/arXiv:1301.5366) ("Electromagnetic Wave Source Conditions") of the book [Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707).

```scm
 (reset-meep)
 (change-sources! (list
    (make eigenmode-source
      (src (make gaussian-src (frequency f) (fwidth df)))
      (component Ey)
      (size a a 0)
      (center (* -0.5 (+ s a)) 0)
      (eig-kpoint k-point)
      (eig-match-freq? true)
      (eig-parity ODD-Y))
    (make eigenmode-source
      (src (make gaussian-src (frequency f) (fwidth df)))
      (component Ey)
      (size a a 0)
      (center (* 0.5 (+ s a)) 0)
      (eig-kpoint k-point)
      (eig-match-freq? true)
      (eig-parity ODD-Y)
      (amplitude (if xodd? -1.0 1.0)))))
 (define wvg-pwr (add-flux f 0 1
    (make flux-region (direction Z) (center 0 0) 
      (size (* 1.2 (+ (* 2 a) s)) (* 1.2 a) 0))))
 (define wvg-force (add-force f 0 1
   (make force-region (direction X) (weight +1.0) 
         (center (* 0.5 s) 0) (size 0 a))
   (make force-region (direction X) (weight -1.0) 
         (center (+ (* 0.5 s) a) 0) (size 0 a))))
 (define-param runtime 5000)
 (run-sources+ runtime)
 (display-fluxes wvg-pwr)
 (display-forces wvg-force)
```

<b>NOTE:</b> if MPB is not pre-installed, simply remove all lines pertaining to `change-sources!`

There are two important items to note: (1) We have defined a single flux surface to compute the Poynting vector along $Z$ which spans an area slightly larger than both waveguides rather than two separate flux surfaces, one for each waveguide. This is because in the limit of small separation, two flux surfaces overlap whereas the total power through a single flux surface need, by symmetry, only be halved in order to determine the value for just one of the two waveguides. (2) Instead of defining a closed, four-sided "box" surrounding the waveguides for computing the MST, we chose instead to compute the MST along two $y$-oriented sides with different `weight` values to correctly sum the total force one on each side of the waveguide. By symmetry, we need not consider the force along the $Y$ direction. Choosing a suitable `runtime` requires some care. A large `runtime` is necessary to obtain a steady-state response but this will also lead to large values for the discrete Fourier-transformed fields used to compute both the flux and the MST. These large values might may produce floating-point errors.

We run this simulation over a range of non-zero separation distances and compare the result to that obtained from MPB. This is shown in the figure below. The two methods show good agreement. The MPB data for this plot was generated using this [control file](http://ab-initio.mit.edu/~oskooi/wiki_data/parallel-wvgs-mpb.ctl) and [shell script](http://ab-initio.mit.edu/~oskooi/wiki_data/run_wvgs_mpb.sh). The plot of the MPB data was generated using this [Jupyter notebook](http://ab-initio.mit.edu/~oskooi/wiki_data/MPB_data_plot.ipynb) ([html](http://ab-initio.mit.edu/~oskooi/wiki_data/MPB_data_plot.html)).

<center>
![](../images/Waveguide_forces.png)
</center>