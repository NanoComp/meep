---
# Meep Tutorial Material dispersion
---

In this example, we will perform a simulation with a **frequency-dependent dielectric** ε(ω), corresponding to **material dispersion**. (See [Dielectric materials in Meep](Materials_in_Meep.md) for more information on how material dispersion is supported in Meep.) In particular, we will model a *uniform medium* of the dispersive material; see also the `material-dispersion.ctl` file included with Meep. From the dispersion relation $\omega(k)$, we will compute the numerical ε(ω) via the formula:

$$\varepsilon(\omega) = \left( \frac{ck}{\omega} \right) ^2$$

We will then compare this with the analytical ε(ω) that we specified.

Since this is a uniform medium, our computational cell can actually be of *zero* size (i.e. one pixel), where we will use Bloch-periodic boundary conditions to specify the wavevector *k*.

```
(set! geometry-lattice (make lattice (size no-size no-size no-size)))
(set-param! resolution 20)
```


We will then fill all space with a dispersive material:

```
(set! default-material
      (make dielectric (epsilon 2.25)
            (E-susceptibilities 
             (make lorentzian-susceptibility
               (frequency 1.1) (gamma 1e-5) (sigma 0.5))
             (make lorentzian-susceptibility
               (frequency 0.5) (gamma 0.1) (sigma 2e-5))
             )))
```


corresponding to the dielectric function:

$$\varepsilon(\omega) = \varepsilon(2\pi f) = 2.25 + \frac{1.1^2 \cdot 0.5}{1.1^2 - f^2 -if \cdot 10^{-5}/2\pi} + \frac{0.5^2 \cdot 2\cdot 10^{-5}}{0.5^2 - f^2 -if \cdot 0.1 / 2\pi}$$

The real and imaginary parts of this dielectric function ε(ω) are plotted below:


![center|Real and imaginary parts of specified analytical ε(ω).](../images/Material-dispersion-eps.png)



Here, we can see that the f=1.1 resonance causes a large change in both the real and imaginary parts of ε around that frequency. In fact, there is a range of frequencies from 1.1 to 1.2161 where ε is *negative*. In this range, no propagating modes exist—it is actually a kind of photonic band gap associated with polariton resonances in a material. (For more information on the physics of such materials, see e.g. chapter 10 of *Introduction to Solid State Physics* by C. Kittel.)

On the other hand, the f;=0.5 resonance, because the `sigma` numerator is so small, causes very little change in the real part of ε. Nevertheless, it generates a clear peak in the *imaginary* part of ε, corresponding to a resonant absorption peak.

Now, we'll set up the rest of the simulation. We'll specify a broad-band TM-polarized Gaussian source, create a list of *k* wavevectors that we want to compute $\omega(k)$ over, and compute the associated frequencies by using the `run-k-points` function:

```
(define-param fcen 1.0)
(define-param df 2.0)
(set! sources (list (make source
                      (src (make gaussian-src (frequency fcen) (fwidth df)))
                      (component Ez) (center 0 0 0))))
(define-param kmin 0.3)
(define-param kmax 2.2)
(define-param k-interp 99)
(define kpts (interpolate k-interp (list (vector3 kmin) (vector3 kmax))))
(define all-freqs (run-k-points 200 kpts)) ; a list of lists of frequencies  
```


The `run-k-points` function returns a *list of lists* of frequencies—one list of (complex) frequencies for each *k* point—which we store in the `all-freqs` variable. Finally, we want to loop over this list and print out the corresponding ε via the ratio $(ck/\omega)^2$ as described above. To do this, we will use the Scheme `map` function, which applies a given function to every element of a list (or lists), and since we have a list of lists we'll actually nest two `map` functions:

```
(map (lambda (kx fs)
       (map (lambda (f)
              (print "eps:, " (real-part f) ", " (imag-part f)
                     ", " (sqr (/ kx f)) "\n"))
            fs))
     (map vector3-x kpts) all-freqs)
```


Alternatively we could just read all of the frequencies into Matlab or a spreadsheet and compute the ratios there. After running the program with

```
unix% meep material-dispersion.ctl | tee material-dispersion.out
```


we can then `grep` for the frequencies and the computed dielectric function, and plot it. First, let's plot the dispersion relation $\omega(k)$ (for the real part of ω):


![center|Band diagram for dispersive material.](../images/Material-dispersion-bands.png)



Here, the red circles are the computed points from Meep, whereas the blue line is the analytical band diagram from the specified ε(ω). As you can see, we get *two* bands at each *k*, separated by a polaritonic gap (shaded yellow). This dispersion relation can be thought of as the interaction (anti-crossing) between the light line of the ambient ε=2.25 material (dashed black line) and the horizontal line corresponding to the phonon resonance.

Similarly, the computed and analytical real parts of the dielectric function are given by:


![center|Real part of dielectric function](../images/Material-dispersion-epsre.png)



which shows excellent agreement between the analytical (blue line) and numerical (red circles) calculations. The imaginary part, however, is more subtle:


![center|Real part of dielectric function](../images/Material-dispersion-epsim.png)



Here, the blue line is the analytical calculation from above and the red circles are the numerical value from Meep—why is the agreement so poor? There is nothing wrong with Meep, and this is *not* a numerical error; the problem is simply that we are comparing apples and oranges.

The blue line is the analytical calculation of ε(ω) for a *real* frequency ω (which corresponds to solutions with a *complex* wavevector *k*), whereas Meep is computing ε at a *complex* ω for a *real* wavevector *k*. So, the correct comparison is to plug Meep's *complex* ω into the analytical formula for ε(ω), which results in the green lines on the graph (that fall almost on top of the red circles).

Why did our comparison of the *real* part of ε look so good, then? The reason is that ε(ω) at real and complex values of ω are closely related by the analytic properties of ε. In particular, because ε is an analytic function on the real-ω axis, adding a *small* imaginary part to ω as we are doing here (the losses are small for all of the computed *k* points) does not change ε by much. The change was only significant for the imaginary ε because the imaginary ε was small to begin with. [Category:Meep examples](Meep_examples.md)
