---
# Meep Tutorial Ring resonator in cylindrical coordinates
---

In the [Meep tutorial](Meep_Tutorial.md), we computed the [modes of a ring resonator](Meep_Tutorial#Modes_of_a_ring_resonator.md) by performing a 2d simulation. Here, we will simulate the *same* structure while [exploiting](Exploiting_symmetry_in_Meep.md) the fact that the system has *continuous* rotational symmetry, by performing the simulation in [cylindrical coordinates](Cylindrical_coordinates_in_Meep.md). See also the `ring-cyl.ctl` example file included with Meep.

The ctl file
------------

We begin, as usual, by defining the parameters of the problem, with exactly the same values as in the 2d simulation.

```
(define-param n 3.4) ; index of waveguide
(define-param w 1) ; width of waveguide
(define-param r 1) ; inner radius of ring
(define-param pad 4) ; padding between waveguide and edge of PML
(define-param dpml 2) ; thickness of PML
```


Now, we'll define the dimensions and size of the computational cell:

```
(define sr (+ r w pad dpml)) ; radial size (cell is from 0 to sr)
(set! dimensions CYLINDRICAL)
(set! geometry-lattice (make lattice (size sr no-size no-size)))
```


The key thing here was to set the `dimensions` parameter to `CYLINDRICAL`. This means that all vectors will represent $(r,\phi,z)$ coordinates instead of $(x,y,z)$. The computational cell in the $r$ direction is of size `sr` `=` `r` `+` `w` `+` `pad` `+` `dpml`, and runs from `0` to `sr` (by default) rather than from `-sr/2` to `sr/2` as it would for any other dimension. Note that our $z$ size is `no-size` because it is two-dimensional. The $\phi$ size is also `no-size`, corresponding to the continuous rotational symmetry (a finite $\phi$ size might correspond to discrete rotational symmetry, but this is not currently supported in Meep).

In particular, in systems with continuous rotational symmetry, by an analogue of Bloch's theorem, the angular dependence of the fields can always be chosen in the form $\exp(i m \phi)$ for some integer $m$. Meep uses this fact to treat the angular dependence analytically, with $m$ given by the [input variable](Meep_Reference#Input_variables.md) `m` (which we'll set to `3`, for now).

```
(set-param! m 3)
```


Thus, we are essentially performing a 1d calculation, where Meep must discretize the $r$ direction only. For this reason, it will be much faster than the previous 2d calculation.

The geometry is now specified by a single `block` object—remember that this is a block in cylindrical coordinates, so that it really specifies an annular ring:

```
(set! geometry (list
                (make block (center (+ r (/ w 2))) (size w infinity infinity)
                      (material (make dielectric (index n))))))
(set! pml-layers (list (make pml (thickness dpml))))
(set-param! resolution 10)
```


We have added PML layers on "all" sides. Meep, however, notices that the $z$ direction has no thickness and automatically makes it periodic with no PML. Meep also omits PML from the boundary at $r=0$ (which is handled by the analytical reflection symmetry).

Now, the remaining inputs are almost exactly the same as in the previous 2d simulation. We'll add a single Gaussian point source in the $z$ direction to excite TM modes, with some center frequency and width:

```
(define-param fcen 0.15) ; pulse center frequency                            
(define-param df 0.1)  ; pulse width (in frequency) 
(set! sources (list
               (make source
                 (src (make gaussian-src (frequency fcen) (fwidth df)))
                 (component Ez) (center (+ r 0.1) 0 0))))
              
```


Note that this isn't really a "point" source, however, because of the cylindrical symmetry—it is really a "ring" source with $\phi$ dependence $\exp(i m \phi)$. Finally, as before, we run until the source has turned off, plus 200 additional time units during which we use [harminv](http://ab-initio.mit.edu/wiki/index.php/harminv) to analyze the $E_z$ field at a given point to extract the frequencies and decay rates of the modes.

```
(run-sources+ 200 (after-sources (harminv Ez (vector3 (+ r 0.1)) fcen df)))
```


At the very end, we'll also output one period of the fields to make movies, etcetera. A single field output would be a 1d dataset (along the $r$ direction), so to make things more interesting we'll use `to-appended` to append these datasets to a single HDF5 file to get an $r \times t$ 2d dataset. We'll also use `in-volume` to specify a larger output volume than just the computational cell: in particular, we'll output from `-sr` to `sr` in the $r$ direction, where Meep will automatically infer the $-r$ field values from the reflection symmetry.

```
(run-until (/ 1 fcen) 
           (in-volume (volume (center 0) (size (* 2 sr)))
                      (at-beginning output-epsilon)
                      (to-appended "ez" 
                                   (at-every (/ 1 fcen 20) output-efield-z))))
```


Results
-------

Now, we are ready to run our simulation. Recall that, in the 2d calculation, we got three modes in this frequency range: one at $\omega=0.11785$ with $Q=77$ and an $m=3$ field pattern, one at $\omega=0.14687$ with $Q=351$ and an $m=4$ field pattern, and one at $\omega=0.17501$ with $Q=1630$ and an $m=5$ field pattern. We should get the *same* modes here (with some differences due to the finite resolution), except now that we will have to run *three* calculations, a separate one for each value of $m$ (it will still be much faster than before because the simulations are 1d instead of 2d).

In particular, we'll run:

```
unix% meep m=3 ring-cyl.ctl | grep harminv
unix% meep m=4 ring-cyl.ctl | grep harminv
unix% meep m=5 ring-cyl.ctl | grep harminv
```


giving the combined output:

```
harminv0:, frequency, imag. freq., Q, |amp|, amplitude, error
harminv0:, 0.11834848194079, -6.80930025762674e-4, 86.9020879261668, 0.257477542991357, -0.234862526831655-0.105519091330034i, 2.6465657298186e-10
harminv0:, 0.147555705534808, -1.91078761299536e-4, 386.112262114517, 1.93737432741834, 1.35411847722594+1.38556213652616i, 2.73521325130449e-11
harminv0:, 0.175944214054996, -4.82976799119763e-5, 1821.45616907125, 0.45258172336278, -0.107884449861492-0.439535165601237i, 1.2656772930993e-10
```


This is indeed very close to the 2d simulations: the frequencies are within 1% of the previous values. The $Q$ values (lifetimes) differ by a larger amount (although they are still reasonably close).

Which is more accurate, the 2d or the cylindrical simulation? We can answer this question by increasing the resolutions in both cases and seeing what they converge towards. In particular, let's focus on the m=4 mode. In the cylindrical case, if I double the resolution to 20 I get $\omega=0.14748$ and $Q=383$. In the 2d case, if I double the resolution to 20 I get $\omega=0.14733$ and $Q=321$. So, it looks like the frequencies are clearly converging together and that the cylindrical simulation is more accurate (as you might expect since it describes the $\phi$ direction analytically). But the $Q$ values seem to be getting *farther* apart—what's going on?

The problem is twofold. First, there is some signal-processing error in determining $Q$ in the 2d case, as indicated by the "error" column of the `harminv` output which is only `4e-7` for the 2d simulation vs. `3e-11` for the cylindrical case. We can bring this error down by running with a narrower bandwidth source, which excites just one mode and gives a cleaner signal, or by analyzing over a longer time than 200. Doing the former, we find that the 2d value of $Q$ at a resolution of 20 should really be $Q=343$. Second, [PML](Perfectly_matched_layer.md) absorbing layers are really designed to absorb planewaves incident on flat interfaces, but here we have a *cylindrical* PML layer. Because of this, there are larger numerical reflections from the PML in the cylindrical simulation, which we can rectify by pushing the PML out to a larger radius (i.e. using a larger value of `pad`) and/or increasing the PML thickness (increasing `dpml`) so that it turns on more adiabatically. In the cylindrical simulation for `resolution=20`, if we increase to `dpml=16`, we get $Q=342$, which is in much better agreement with the 2d calculation (and if we increase to `dpml=32` we get the same $Q=342$, so it seems to be converged).

This illustrates the general principle that you need to check several parameters to ensure that results are converged in time-domain simulations: the resolution, the run time, the PML thickness, etcetera.

Finally, we can get the field images. Since we only are exciting one mode per `m` here anyway, according to `harminv`, we don't really need to use a narrow-band source. We'll do so anyway just to remind you of the general procedure, however, e.g. for the $\omega=0.118$ $m=3$ mode:

```
unix% meep m=3 fcen=0.118 df=0.01 ring-cyl.ctl
unix% h5topng -S 2 -Zc dkbluered -C ring-cyl-eps-001200.00.h5 ring-cyl-ez.h5
```


Note that, because of the `to-appended` command, the `ring-cyl-ez.h5` file is a 160×18 dataset corresponding to an $r \times t$ slice. Repeating this for all three modes results in the images:

+ $E_z$ for $\omega=0.118$ $m=3$ mode:

![center](../images/Ring-cyl-ez-0.118.png)


+ $E_z$ for $\omega=0.148$ $m=4$ mode: 

![center](../images/Ring-cyl-ez-0.148.png)


+ $E_z$ for $\omega=0.176$ $m=5$ mode:

![center](../images/Ring-cyl-ez-0.176.png)

Because we are looking only at a $\phi=0$ slice, the visual distinction between $m$ values is much less than with the 2d simulation. What is apparent is that, as the frequency increases, the mode becomes more localized in the waveguide and the radiating field (seen in the $r \times t$ slice as curved waves extending outward) becomes less, as expected.

[Category:Meep examples](Meep_examples.md)
