; Meep Tutorial: Hz-polarized transmission and reflection through a cavity
; formed by a periodic sequence of holes in a dielectric waveguide,
; with a defect formed by a larger spacing between one pair of holes.

; This structure is based on one analyzed in:
;    S. Fan, J. N. Winn, A. Devenyi, J. C. Chen, R. D. Meade, and
;    J. D. Joannopoulos, "Guided and defect modes in periodic dielectric
;    waveguides," J. Opt. Soc. Am. B, 12 (7), 1267-1272 (1995).

; Some parameters to describe the geometry:
(define-param eps 13) ; dielectric constant of waveguide
(define-param w 1.2) ; width of waveguide
(define-param r 0.36) ; radius of holes
(define-param d 1.4) ; defect spacing (ordinary spacing = 1)
(define-param N 3) ; number of holes on either side of defect

; The cell dimensions
(define-param sy 6) ; size of cell in y direction (perpendicular to wvg.)
(define-param pad 2) ; padding between last hole and PML edge
(define-param dpml 1) ; PML thickness

(define sx (+ (* 2 (+ pad dpml N)) d -1)) ; size of cell in x direction
(set! geometry-lattice (make lattice (size sx sy no-size)))

(set! geometry
      (append ; combine lists of objects:
       (list (make block (center 0 0) (size infinity w infinity)
		   (material (make dielectric (epsilon eps)))))
       (geometric-object-duplicates (vector3 1 0) 0 (- N 1)
	(make cylinder (center (/ d 2) 0) (radius r) (height infinity)
	      (material air)))
       (geometric-object-duplicates (vector3 -1 0) 0 (- N 1)
	(make cylinder (center (/ d -2) 0) (radius r) (height infinity)
	      (material air)))))

(set! pml-layers (list (make pml (thickness dpml))))
(set-param! resolution 20)

(define-param fcen 0.25) ; pulse center frequency
(define-param df 0.2)  ; pulse width (in frequency)

(define-param nfreq 500) ; number of frequencies at which to compute flux

; false = transmission spectrum, true = resonant modes:
(define-param compute-mode? false)

(if compute-mode?
    (begin
      (set! sources (list
		     (make source
		       (src (make gaussian-src (frequency fcen) (fwidth df)))
		       (component Hz) (center 0 0))))

      (set! symmetries
	    (list (make mirror-sym (direction Y) (phase -1))
		  (make mirror-sym (direction X) (phase -1))))

      (run-sources+ 400
		    (at-beginning output-epsilon)
		    (after-sources (harminv Hz (vector3 0) fcen df)))
      (run-until (/ 1 fcen) (at-every (/ 1 fcen 20) output-hfield-z))
      )
    (begin
      (set! sources (list
		     (make source
		       (src (make gaussian-src (frequency fcen) (fwidth df)))
		       (component Ey)
		       (center (+ dpml (* -0.5 sx)) 0)
		       (size 0 w))))

      (set! symmetries (list (make mirror-sym (direction Y) (phase -1))))

      (let ((trans ; transmitted flux
	      (add-flux fcen df nfreq
		        (make flux-region
		          (center (- (* 0.5 sx) dpml 0.5) 0) (size 0 (* w 2))))))

          (run-sources+ (stop-when-fields-decayed
		        50 Ey
		        (vector3 (- (* 0.5 sx) dpml 0.5) 0)
		        1e-3)
		        (at-beginning output-epsilon)
		        (during-sources
                        (in-volume (volume (center 0 0) (size sx 0))
                                   (to-appended "hz-slice" (at-every 0.4 output-hfield-z)))))

          (display-fluxes trans) ; print out the flux spectrum
      )))
