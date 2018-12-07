(set-param! resolution 100)

(set! geometry-lattice (make lattice (size 10 10 no-size)))

(set! geometry (list (make cylinder (center 0 0 0) (radius 3) (height infinity)
			   (material (make medium (index 3.5))))
		     (make ellipsoid (center 0 0 0) (size 1 2 infinity)
			   (material air))))

(set! pml-layers (list (make pml (thickness 1))))

;(define-param src-cmpt Ez) ; S-polarization: Ez / P-polarization: Hz
(define-param src-cmpt Hz)
(set! sources (list (make source (src (make gaussian-src (frequency 1) (fwidth 0.1)))
			  (center 0 0 0) (component src-cmpt))))

(if (= src-cmpt Ez)
    (set! symmetries (list (make mirror-sym (direction X))
			   (make mirror-sym (direction Y)))))

(if (= src-cmpt Hz)
    (set! symmetries (list (make mirror-sym (direction X) (phase -1))
			   (make mirror-sym (direction Y) (phase -1)))))

;(define print-stuff (lambda () (print "field:, " (get-field-point src-cmpt (vector3 4.13 3.75 0)) "\n")))
(define print-stuff (lambda () (print "t, Ez: " (meep-round-time) " " (get-field-point src-cmpt (vector3 4.13 3.75 0)) "\n")))

(run-until 23 (at-beginning output-epsilon)
	      (at-every 0.25 print-stuff)
	      (at-end print-stuff) 
              (at-end output-efield-z))

(print "stopped at meep time = " (meep-round-time) )
