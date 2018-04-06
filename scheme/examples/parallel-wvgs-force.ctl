(set-param! resolution 30)
(define-param nSi 3.45)
(define Si (make medium (index nSi)))
(define-param dpml 1.0)
(define-param sx 5)
(define-param sy 3)
(set! geometry-lattice
      (make lattice (size (+ sx (* 2 dpml)) (+ sy (* 2 dpml)) no-size)))
(set! pml-layers (list (make pml (thickness dpml))))
(define-param a 1.0) ; waveguide width
(define-param s 1.0)  ; waveguide separation distance
(set! geometry (list
		(make block (center (* -0.5 (+ s a)) 0)
		      (size a a infinity) (material Si))
		(make block (center (*  0.5 (+ s a)) 0)
		      (size a a infinity) (material Si))))

(define-param xodd? true)
(set! symmetries (list
		  (make mirror-sym (direction X) (phase (if xodd? -1 +1)))
		  (make mirror-sym (direction Y) (phase -1))))

(define-param beta 0.5)
(set! k-point (vector3 0 0 beta))

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
		    (component Ey)
		    (size a a 0)
		    (center (* -0.5 (+ s a)) 0)
		    (eig-kpoint k-point)
		    (eig-match-freq? true)
		    (eig-parity ODD-Y))
		  (make eigenmode-source
		    (src (make gaussian-src (frequency f) (fwidth df)))
		    (component Ey)
		    (size a a 0)
		    (center (* 0.5 (+ s a)) 0)
		    (eig-kpoint k-point)
		    (eig-match-freq? true)
		    (eig-parity ODD-Y)
		    (amplitude (if xodd? -1.0 1.0)))))
(define wvg-pwr (add-flux f 0 1
			  (make flux-region (direction Z) (center 0 0)
				(size (* 1.2 (+ (* 2 a) s)) (* 1.2 a) 0))))
(define wvg-force (add-force f 0 1
			     (make force-region (direction X) (weight +1.0)
				   (center (* 0.5 s) 0) (size 0 a))
			     (make force-region (direction X) (weight -1.0)
				   (center (+ (* 0.5 s) a) 0) (size 0 a))))
(define-param runtime 5000)
(run-sources+ runtime)
(display-fluxes wvg-pwr)
(display-forces wvg-force)
