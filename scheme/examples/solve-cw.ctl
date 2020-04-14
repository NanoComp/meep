; Calculating 2d ring-resonator modes using frequency-domain solver, from the Meep tutorial.

(define-param n 3.4)
(define-param w 1)
(define-param r 1)
(define-param pad 4)
(define-param dpml 2)

(define sxy (* 2 (+ r w pad dpml)))
(set! geometry-lattice (make lattice (size sxy sxy no-size)))

(set! geometry (list
		(make cylinder (center 0 0) (height infinity)
		      (radius (+ r w)) (material (make dielectric (index n))))
		(make cylinder (center 0 0) (height infinity)
		      (radius r) (material air))))

(set! pml-layers (list (make pml (thickness dpml))))
(set-param! resolution 10)

(define-param fcen 0.118)
(set! sources (list
               (make source
                 (src (make continuous-src (frequency fcen)))
                 (component Ez)
                 (center (+ r 0.1) 0))
               (make source
                 (src (make continuous-src (frequency fcen)))
                 (component Ez)
                 (center (- (+ r 0.1)) 0)
                 (amplitude -1))))

(set! symmetries (list (make mirror-sym (direction X) (phase -1))
                       (make mirror-sym (direction Y) (phase +1))))

(set! force-complex-fields? true)

(define-param solve-cw-tol 1e-8)
(define-param solve-cw-maxiters 10000)
(define-param solve-cw-L 10)

(define (ez-real r ez) (real-part ez))

(init-fields)
(meep-fields-solve-cw fields solve-cw-tol solve-cw-maxiters solve-cw-L)
(in-volume (volume (center 0 0) (size (- sxy (* 2 dpml)) (- sxy (* 2 dpml))))
	   (output-epsilon)
	   (output-real-field-function "ez-real" (list Ez) ez-real))

(exit)
