;; This file realizes a 1D, one-sided Fabry-Perot laser, as described in Fig. 2 of Optics Express, Vol. 20, pp. 474-88, 2012.

;; Cavity definitions
(set-param! resolution 400)
(define-param ncav 1.5)                          ; cavity refractive index
(define-param Lcav 1)                            ; cavity length
(define-param dpad 1)                            ; padding thickness
(define-param dpml 1)                            ; PML thickness
(define-param sz (+ Lcav dpad dpml))
(set! geometry-lattice (make lattice (size no-size no-size sz)))
(set! dimensions 1)
(set! pml-layers (list (make pml (thickness dpml) (side High))))

;; For defining laser properties in MEEP, the transition rates / frequencies are specified in units of 2*pi*a/c.
;; gamma-21 in MEEP is the Full-Width Half-Max, as opposed to gamma_perp, which is the HWHM in SALT.
;; Additionally, the classical coupling element sigma = 2*theta^2*omega_a/hbar, where
;; theta is the off-diagonal dipole matrix element.

;; These different conventions can cause a bit of confusion when comparing against SALT, so here we perform
;; this transformation explicitly.

(define-param omega-a 40)                        ; omega_a in SALT
(define freq-21 (/ omega-a (* 2 pi)))            ; emission frequency  (units of 2πc/a)

(define-param gamma-perp 4)                      ; HWHM in angular frequency, SALT
(define gamma-21 (/ (* 2 gamma-perp) (* 2 pi)))  ; FWHM emission linewidth in sec^-1  (units of 2πc/a)
					         ; Note that 2*pi*gamma-21 = 2*gamma_perp in SALT.

(define-param theta 1)                           ; theta, the off-diagonal dipole matrix element, in SALT
(define sigma-21 (* 2 theta theta omega-a))      ; dipole coupling strength (hbar = 1)

;; The gain medium in MEEP is allowed to have an arbitrary number of levels, and is not
;; restricted to a two-level gain medium, as it simulates the populations of every individual
;; atomic energy level.

;; If you are using a 2 level gain model, you can compare against
;; results which only simulate the atomic inversion, using the definitions
;; gamma_parallel = pumping-rate + rate-21
;; D_0 = (pumping-rate - rate-21)/(pumping-rate + rate-21) * N0

;; In fact, even if you arn't using a 2 level gain model, you can compare against an effective
;; two level model using the formula provided in Cerjan et al., Opt. Express 20, 474 (2012).

;; Here, D_0 as written above is not yet in "SALT" units. To make this conversion,
;; D_0 (SALT) = theta^2/(hbar*gamma_perp) * D_0 (as written above)

;; Finally, note the lack of 4*pi in the above conversion that is written in many published SALT papers.
;; This 4*pi comes from using Gaussian units, in which the displacement field, D = E + 4*pi*P, whereas
;; in SI units, D = eps0*E + P, which is what MEEP uses.

;; Gain medium pump and decay rates are specified in units of c/a.

(define-param rate-21 0.005)                     ; non-radiative rate  (units of c/a)
(define-param N0 37)                             ; initial population density of ground state
(define-param Rp 0.0051)                         ; pumping rate of ground to excited state
;; so for example, these parameters have D_0 (SALT) = 0.0693.

;; Make the actual medium in MEEP:
(define two-level (make medium (index ncav)
			(E-susceptibilities (make multilevel-atom (sigma-diag 1 0 0)
			  (transitions (make transition (from-level 1) (to-level 2) (pumping-rate Rp)
					     (frequency freq-21) (gamma gamma-21) (sigma sigma-21))
				       (make transition (from-level 2) (to-level 1) (transition-rate rate-21)))
			  (initial-populations N0)))))

;; Specify the cavity geometry:
(set! geometry (list (make block (center 0 0 (+ (* -0.5 sz) (* 0.5 Lcav)))
			   (size infinity infinity Lcav) (material two-level))))

(init-fields)
(meep-fields-initialize-field fields Ex
			      (lambda (p) (if (= (vector3-z p) (+ (* -0.5 sz) (* 0.5 Lcav))) 1 0)))

;; Specify the end time:
(define-param endt 7000)
;; Note that the total number of time steps run is endt*resolution*2. This is the origin of the extra
;; factor of 2 in the definition of dt in fieldfft_meep.m.

(define print-field (lambda () (print "field:, " (meep-time) ", "
	       (real-part (get-field-point Ex (vector3 0 0 (+ (* -0.5 sz) Lcav (* 0.5 dpad))))) "\n")))

(run-until endt (after-time (- endt 250) print-field))
