(set-param! resolution 40)   ; pixels/um

(define Si (make medium (index 3.45)))

(define-param dpml 1.0)
(set! pml-layers (list (make pml (thickness dpml))))

(define-param sx 5)
(define-param sy 3)
(set! geometry-lattice
      (make lattice (size (+ sx (* 2 dpml)) (+ sy (* 2 dpml)) no-size)))

(define-param a 1.0)  ; waveguide width/height
(define-param s 1.0)  ; waveguide separation distance
(set! geometry (list
		(make block (center (* -0.5 (+ s a)) 0)
		      (size a a infinity) (material Si))
		(make block (center (* 0.5 (+ s a)) 0)
		      (size a a infinity) (material Si))))

(define-param xodd? true)
(set! symmetries (list
		  (make mirror-sym (direction X) (phase (if xodd? -1 +1)))
		  (make mirror-sym (direction Y) (phase -1))))

(set! k-point (vector3 0 0 0.5))

(define-param fcen 0.22)
(define-param df 0.06)
(set! sources (list
	       (make source (src (make gaussian-src (frequency fcen) (fwidth df)))
		     (component Ey) (center (* -0.5 (+ s a)) 0) (size a a 0))
	       (make source (src (make gaussian-src (frequency fcen) (fwidth df)))
		     (component Ey) (center (* 0.5 (+ s a)) 0) (size a a 0)
		     (amplitude (if xodd? -1.0 1.0)))))

(run-sources+ 200
	      (after-sources (harminv Ey (vector3 (* 0.5 (+ s a)) 0) fcen df)))

(define f (harminv-freq-re (car harminv-results)))
(print "freq:, " s ", " f "\n")

(reset-meep)
(change-sources! (list
                  (make eigenmode-source
                    (src (make gaussian-src (frequency f) (fwidth df)))
                    (size sx sy 0)
                    (center 0 0 0)
                    (eig-band (if xodd? 2 1))
                    (eig-kpoint k-point)
                    (eig-match-freq? false)
                    (eig-parity ODD-Y))))

(define wvg-flux (add-flux f 0 1
                           (make flux-region (direction Z)
                                 (center 0 0) (size sx sy 0))))

(define wvg-force (add-force f 0 1
			     (make force-region (direction X) (weight +1)
				   (center (* 0.5 s) 0) (size 0 sy))
			     (make force-region (direction X) (weight -1)
				   (center (+ (* 0.5 s) a) 0) (size 0 sy))))

(run-sources+ 1500)

(display-fluxes wvg-flux)
(display-forces wvg-force)
