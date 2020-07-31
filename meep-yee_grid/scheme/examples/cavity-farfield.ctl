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

(set! sources (list
	       (make source
		 (src (make gaussian-src (frequency fcen) (fwidth df)))
		 (component Hz) (center 0 0))))

(set! symmetries
      (list (make mirror-sym (direction Y) (phase -1))
	    (make mirror-sym (direction X) (phase -1))))

(define-param d1 0.2)
(define nearfield
  (add-near2far fcen 0 1
    (make near2far-region (center 0 (+ (* 0.5 w) d1))
	  (size (- sx (* 2 dpml)) 0))
    (make near2far-region (center (+ (* -0.5 sx) dpml) (+ (* 0.5 w) (* 0.5 d1)))
	  (size 0 d1)
	  (weight -1.0))
    (make near2far-region (center (- (* 0.5 sx) dpml) (+ (* 0.5 w) (* 0.5 d1)))
	  (size 0 d1))))

(run-sources+ (stop-when-fields-decayed 50 Hz (vector3 0.12 -0.37) 1e-8))
(define-param d2 20)
(define-param h 4)
(output-farfields nearfield
  (string-append "spectra-" (number->string d1) "-" (number->string d2) "-" (number->string h))
  (volume (center 0 (+ (* 0.5 w) d2 (* 0.5 h))) (size (- sx (* 2 dpml)) h)) resolution)
