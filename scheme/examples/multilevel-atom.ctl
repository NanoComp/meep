;; Cavity definitions
(set-param! resolution 500)
(define-param ncav 1.5)                          ; cavity refractive index
(define-param Lcav 1)                            ; cavity length
(define-param dpad 0.1)                            ; padding thickness
(define-param dpml 0.1)                            ; PML thickness
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
(define freq-21 (/ omega-a (* 2 pi)))            ; emission frequency  (units of 2\pia/c)

(define-param gamma-perp 4)                      ; HWHM in angular frequency, SALT
(define gamma-21 (/ (* 2 gamma-perp) (* 2 pi)))  ; FWHM emission linewidth in sec^-1 (units of 2\pia/c)
					         ; Note that 2*pi*gamma-21 = 2*gamma_perp in SALT.

(define-param theta 1)                           ; theta, the off-diagonal dipole matrix element, in SALT
(define sigma-21 (* 2 theta theta omega-a))      ; dipole coupling strength (hbar = 1)

;; The gain medium in MEEP is allowed to have an arbitrary number of levels, and is not
;; restricted to a two-level gain medium, as it simulates the populations of every individual
;; atomic energy level. To compare against results which only simulate the atomic inversion, use
;; a two-level gain medium, and note that:
;; gamma_parallel = pumping-rate + rate-21
;; D_0 = (pumping-rate - rate-21)/(pumping-rate + rate-21) * N0

;; Note that D_0 here is not yet in "SALT" units. To make this conversion,
;; D_0 (SALT) = 4*pi*theta^2/(hbar*gamma_perp) * D_0 (as written above)

;; Gain medium pump and decay rates are specified in units of c/a.

(define-param rate-21 0.005)                         ; non-radiative rate  (units of c/a)
(define-param N0 10)                               ; initial population density of ground state
(define-param Rp 1.00514)                          ; pumping rate of ground to excited state
;; so for example, these parameters have D_0 (SALT) = 0.43.

;; Make the actual medium in MEEP:
(define two-level (make medium (index ncav)
			(E-susceptibilities (make multilevel-atom (sigma 1)
			  (transitions (make transition (from-level 1) (to-level 2) (pumping-rate Rp)
					     (frequency freq-21) (gamma gamma-21) (sigma sigma-21))
				       (make transition (from-level 2) (to-level 1) (transition-rate rate-21)))
			  (initial-populations N0)))))

;; Specify the cavity geometry:
(set! geometry (list (make block (center 0 0 (+ (* -0.5 sz) (* 0.5 Lcav)))
			   (size infinity infinity Lcav) (material two-level))))

;; Initialize the fields, has to be non-zero, doesn't really matter what.
(init-fields)
(meep-fields-initialize-field fields Ex
			      (lambda (p) (if (= (vector3-z p) (+ (* -0.5 sz) (* 0.5 Lcav))) 1 0)))
;			      (lambda (p) (if (<= (vector3-z p) (+ (* -0.5 sz) Lcav)) 1 0)))

;; Specify the end time:
(define-param endt 5000)

;; Run it:
(define print-field (lambda () (print "field:, " (meep-time) ", "
	       (real-part (get-field-point Ex (vector3 0 0 (+ (* -0.5 sz) (* 0.5 Lcav))))) "\n")))
;	       (real-part (get-field-point Ex (vector3 0 0 (+ (* -0.5 sz) Lcav (* 0.5 dpad))))) "\n")))
(run-until endt (after-time (- endt 250) print-field))
