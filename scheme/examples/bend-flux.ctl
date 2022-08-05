; From the Meep tutorial: transmission around a 90-degree waveguide
; bend in 2d.

(define-param sx 16) ; size of cell in X direction
(define-param sy 32) ; size of cell in Y direction
(set! geometry-lattice (make lattice (size sx sy no-size)))

(define-param pad 4) ; padding distance between waveguide and cell edge
(define-param w 1) ; width of waveguide

(define wvg-ycen (* -0.5 (- sy w (* 2 pad)))) ; y center of horiz. wvg
(define wvg-xcen (* 0.5 (- sx w (* 2 pad)))) ; x center of vert. wvg

(define-param no-bend? false) ; if true, have straight waveguide, not bend

(set! geometry
      (if no-bend?
	  (list
	   (make block
	     (center 0 wvg-ycen)
	     (size infinity w infinity)
	     (material (make dielectric (epsilon 12)))))
	  (list
	   (make block
	     (center (* -0.5 pad) wvg-ycen)
	     (size (- sx pad) w infinity)
	     (material (make dielectric (epsilon 12))))
	   (make block
	     (center wvg-xcen (* 0.5 pad))
	     (size w (- sy pad) infinity)
	     (material (make dielectric (epsilon 12)))))))

(define-param fcen 0.15) ; pulse center frequency
(define-param df 0.1)  ; pulse width (in frequency)
(set! sources (list
	       (make source
		 (src (make gaussian-src (frequency fcen) (fwidth df)))
		 (component Ez)
		 (center (+ 1 (* -0.5 sx)) wvg-ycen)
		 (size 0 w))))

(set! pml-layers (list (make pml (thickness 1.0))))
(set-param! resolution 10)

(define-param nfreq 100) ; number of frequencies at which to compute flux
(define trans ; transmitted flux
      (add-flux fcen df nfreq
		(if no-bend?
		    (make flux-region
		      (center (- (/ sx 2) 1.5) wvg-ycen) (size 0 (* w 2)))
		    (make flux-region
		      (center wvg-xcen (- (/ sy 2) 1.5)) (size (* w 2) 0)))))
(define refl ; reflected flux
      (add-flux fcen df nfreq
		(make flux-region
		  (center (+ (* -0.5 sx) 1.5) wvg-ycen) (size 0 (* w 2)))))

; for normal run, load negated fields to subtract incident from refl. fields
(if (not no-bend?) (load-minus-flux "refl-flux" refl))

(run-sources+
 (stop-when-fields-decayed 50 Ez
			   (if no-bend?
			       (vector3 (- (/ sx 2) 1.5) wvg-ycen)
			       (vector3 wvg-xcen (- (/ sy 2) 1.5)))
			   1e-3)
 (at-beginning output-epsilon))

; for normalization run, save flux fields for refl. plane
(if no-bend? (save-flux "refl-flux" refl))

(display-fluxes trans refl)
