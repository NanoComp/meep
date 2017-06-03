---
# Meep Tutorial Multilevel-atomic susceptibility
---

Meep 1.4 introduced a feature to model saturable absorption/gain via multilevel-atomic susceptibility. This is based on a generalization of the [Maxwell-Bloch equations](https://en.wikipedia.org/wiki/Maxwell-Bloch_equations) which involve the interaction of a quantized system having an arbitrary number of levels with the electromagnetic fields. Meep's implementation is similar to that described in [S.-L. Chua et al](https://www.osapublishing.org/oe/abstract.cfm?uri=oe-19-2-1539) (eqns. 1-5). We will demonstrate this feature by computing the lasing thresholds of a two-level, multimode cavity in 1d similar to the example used in [A. Cerjan et al](https://www.osapublishing.org/oe/abstract.cfm?uri=oe-20-1-474) (Fig. 2).

The cavity consists of a high-index medium with a perfect-metallic mirror on one end and an abrupt termination in air on the other. We will specify an initial population density for the ground state of the two-level system. The field within the cavity is initialized to arbitrary non-zero values and a fictitious source is used to pump the cavity at a fixed rate. The fields are time stepped until reaching steady state. Near the end of the time stepping, we output the electric field at the center of the cavity and then, in post processing, compute its Fourier transform to obtain the spectra. The simulation script is as follows.

```
(set-param! resolution 1000)
(define-param ncav 1.5)                          ; cavity refractive index
(define-param Lcav 1)                            ; cavity length
(define-param dpad 1)                            ; padding thickness
(define-param dpml 1)                            ; PML thickness
(define-param sz (+ Lcav dpad dpml))
(set! geometry-lattice (make lattice (size no-size no-size sz)))
(set! dimensions 1)
(set! pml-layers (list (make pml (thickness dpml) (side High))))
(define-param freq-21 (/ 40 (* 2 pi)))           ; emission frequency  (units of 2\pia/c)
(define-param gamma-21 (/ 4 (* 2 pi)))           ; emission linewidth  (units of 2\pia/c)
(define-param sigma-21 8e-23)                    ; dipole coupling strength
(set! sigma-21 (/ sigma-21 (sqr freq-21)))
(define-param rate-21 0.005)                     ; non-radiative rate  (units of c/a)
(define-param N0 5e23)                           ; initial population density of ground state
(define-param Rp 0)                              ; pumping rate of ground to excited state
(define two-level (make medium (index ncav)
 (E-susceptibilities (make multilevel-atom (sigma 1)
  (transitions (make transition (from-level 1) (to-level 2) (pumping-rate Rp)
                                (frequency freq-21) (gamma gamma-21) (sigma sigma-21))
               (make transition (from-level 2) (to-level 1) (transition-rate rate-21)))
  (initial-populations N0)))))
(set! geometry (list (make block (center 0 0 (+ (* -0.5 sz) (* 0.5 Lcav)))
                          (size infinity infinity Lcav) (material two-level))))
(init-fields)
(meep-fields-initialize-field fields Ex 
         (lambda (p) (if (= (vector3-z p) (+ (* -0.5 sz) (* 0.5 Lcav))) 1 0)))
(define print-field (lambda () (print "field:, " (meep-time) ", "
      (real-part (get-field-point Ex (vector3 0 0 (+ (* -0.5 sz) (* 0.5 Lcav))))) "\n")))
(define-param endt 30000)
(run-until endt (after-time (- endt 250) print-field))
```


Definition of the two-level medium involves the `multilevel-atom` sub-class of the `E-susceptibilities` material type. Each radiative and non-radiative `transition` is specified separately. The atomic resonance used to drive absorption and gain is based on a damped harmonic oscillator described in [Materials in Meep](Materials_in_Meep.md) with the same parameters. Note that the `sigma` of any given transition is multiplied by the `sigma` of its sub-class definition (1 in this example). `transition-rate` defines the rate of non-radiative decay and `pumping-rate` refers to pumping of the ground to the excited state. It is important to specify the `from-level` and `to-level` parameters correctly, otherwise the results will be undefined.

The choice of these parameters requires some care. For example, choosing a pumping rate that lies far beyond threshold will cause large inversion which is not physical and produce meaningless results. The simulation time is also important when operating near the threshold of a particular mode. The fields contain relaxation oscillations and require sufficient time to reach steady state. We also need to choose a small timestep to ensure that the data is smooth and continuous. This requires a large resolution. A large resolution is also necessary to ensure stability when the strength of the source, which depends on `sigma` and `N0`, driving the polarization is large (this also applies to a linear absorber). The spectra at a pumping rate of 0.0073 is shown below. This plot was generated using this [iPython notebook](http://ab-initio.mit.edu/~oskooi/wiki_data/fourier_transform_cavity_field.ipynb).


![center|Spectra for 1d laser cavity at pumping rate (Rp) of 0.0073](../images/Multilevel_cavity_spectra.png)



There are four modes present: two are lasing while the other two are slightly below threshold. The frequency of the passive cavity modes can be computed analytically using the equation $\omega_{cav}=(m+0.5)\pi/(n_{cav}L_{cav})$ where $n_{cav}$ and $L_{cav}$ are the cavity index and length, and $m$ is an integer. The four modes in the figure correspond to $m$=17-20 which are labelled. In the continuum limit, these modes would appear as Dirac delta functions in the spectra. The discretized model, however, produces peaks with finite width. Thus, we need to integrate a fixed number of points around each peak to smooth out the modal intensity. For this simple two-level cavity, the thresholds can be computed analytically using the steady-state ab-initio laser theory (SALT) developed by Prof. A. Douglas Stone and his group at Yale. Based on the default parameters in the script, two modes, $m$=18 and 19, should begin to lase very close to the relaxation rate. We plot the variation of the modal intensity with pumping rate.


![center|Modal intensity versus pumping rate for 1d laser cavity](../images/Multilevel_modal_intensity.png)



The two modes predicted by SALT to have the lowest thresholds are indeed the first to begin lasing. Note that the slopes of each curve for the two lasing modes are decreasing with increasing pumping rate. This gain saturation occurs because the onset of lasing from additional modes means there is less gain available to the other modes. The modal intensities reach an asymptote in the limit of large pumping rates. We can convert Meep's dimensionless parameters into real units by specifying the units of the cavity length $L_{cav}$ and then multiplying the rate terms by $L_{cav}/c$.
