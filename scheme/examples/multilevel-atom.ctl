(set-param! resolution 1000)
(define-param ncav 1.5)                          ; cavity refractive index
(define-param Lcav 1)                            ; cavity length
(define-param dpad 1)                            ; padding thickness
(define-param dpml 1)                            ; PML thickness
(define-param sz (+ Lcav dpad dpml))
(set! geometry-lattice (make lattice (size no-size no-size sz)))
(set! dimensions 1)
(set! pml-layers (list (make pml (thickness dpml) (side High))))

(define-param freq-21 (/ 40 (* 2 pi)))           ; emission frequency  (units of 2\pia/c)
(define-param gamma-21 (/ 4 (* 2 pi)))           ; emission linewidth  (units of 2\pia/c)
(define-param sigma-21 8e-23)                    ; dipole coupling strength
(set! sigma-21 (/ sigma-21 (sqr freq-21)))
(define-param rate-21 0.005)                     ; non-radiative rate  (units of c/a)
(define-param N0 5e23)                           ; initial population density of ground state
(define-param Rp 0)                              ; pumping rate of ground to excited state

(define two-level (make medium (index ncav)
	(E-susceptibilities (make multilevel-atom (sigma 1)
	  (transitions (make transition (from-level 1) (to-level 2) (pumping-rate Rp)
				        (frequency freq-21) (gamma gamma-21) (sigma sigma-21))
		       (make transition (from-level 2) (to-level 1) (transition-rate rate-21)))
	  (initial-populations N0)))))

(set! geometry (list (make block (center 0 0 (+ (* -0.5 sz) (* 0.5 Lcav)))
			   (size infinity infinity Lcav) (material two-level))))

(init-fields)
(meep-fields-initialize-field fields Ex
      (lambda (p) (if (= (vector3-z p) (+ (* -0.5 sz) (* 0.5 Lcav))) 1 0)))
(define print-field (lambda () (print "field:, " (meep-time) ", "
      (real-part (get-field-point Ex (vector3 0 0 (+ (* -0.5 sz) (* 0.5 Lcav))))) "\n")))
(define-param endt 30000)
(run-until endt (after-time (- endt 250) print-field))
