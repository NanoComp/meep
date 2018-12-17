;; moving point charge with superluminal phase velocity in dielectric media emitting Cherenkov radiation

(set-param! resolution 10)

(define-param sx 60)
(define-param sy 60)
(set! geometry-lattice (make lattice (size sx sy no-size)))

(define-param dpml 1.0)
(set! pml-layers (list (make pml (thickness dpml))))

(set! default-material (make dielectric (index 1.5)))

(define-param v 0.7) ; velocity of point charge

(set! symmetries (list (make mirror-sym (direction Y))))

(run-until (/ sx v)
	   (lambda ()
	     (change-sources! (list (make source
				      (src (make continuous-src (frequency 1e-10)))
				      (component Ex)
				      (center (+ (* -0.5 sx) dpml (* v (meep-time))))))))
	              (at-every 2 (output-png Hz "-vZc dkbluered -M 1")))
