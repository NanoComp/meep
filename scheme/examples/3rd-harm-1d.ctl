; 1d simulation of a plane wave propagating through a Kerr medium
; and generating the third-harmonic frequency component.

(define-param sz 100) ; size of cell in z direction
(define-param fcen (/ 1 3)) ; center frequency of source
(define-param df (/ fcen 20)) ; frequency width of source
(define-param amp 1.0) ; amplitude of source
(define-param k 1e-2) ; Kerr susceptibility

(define-param dpml 1.0) ; PML layer thickness

; We'll use an explicitly 1d simulation.  Setting dimensions=1 will actually
; result in faster execution than just using two no-size dimensions.  However,
; in this case Meep requires us to use E in the x direction (and H in y),
; and our one no-size dimension must be z.
(set-param! dimensions 1)
(set! geometry-lattice (make lattice (size no-size no-size sz)))

; to put the same material in all space, we can just set the default material
(set! default-material (make dielectric (index 1) (chi3 k)))

(set! pml-layers (list (make pml (thickness dpml))))
(set-param! resolution 20)

(set! sources (list
               (make source
                 (src (make gaussian-src (frequency fcen) (fwidth df)))
                 (component Ex)
                 (center 0 0 (+ (* -0.5 sz) dpml))
		 (amplitude amp))))

; frequency range for flux calculation
(define-param nfreq 400)
(define-param fmin (/ fcen 2))
(define-param fmax (* fcen 4))

(define trans ; transmitted flux
  (add-flux (* 0.5 (+ fmin fmax)) (- fmax fmin) nfreq
	    (make flux-region (center 0 0 (- (* 0.5 sz) dpml 0.5)))))

; also compute a ''single'' flux point at fcen and 3*fcen
(define trans1 (add-flux fcen 0 1 (make flux-region
				    (center 0 0 (- (* 0.5 sz) dpml 0.5)))))
(define trans3 (add-flux (* 3 fcen) 0 1 (make flux-region
					(center 0 0 (- (* 0.5 sz) dpml 0.5)))))

(run-sources+
 (stop-when-fields-decayed 50 Ex
			   (vector3 0 0 (- (* 0.5 sz) dpml 0.5))
			   1e-6))

(display-fluxes trans)

(print "harmonics:, " k ", " amp ", "
       (first (get-fluxes trans1)) ", " (first (get-fluxes trans3)) "\n")
