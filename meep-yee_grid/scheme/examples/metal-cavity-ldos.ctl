(set-param! resolution 200)
(define-param sxy 2)
(define-param dpml 1)
(set! sxy (+ sxy (* 2 dpml)))
(set! geometry-lattice (make lattice (size sxy sxy no-size)))
(set! pml-layers (list (make pml (thickness dpml))))
(define-param a 1)
(define-param t 0.1)
(set! geometry (list
      (make block (center 0 0) (size (+ a (* 2 t)) (+ a (* 2 t)) infinity) (material metal))
      (make block (center 0 0) (size a a infinity) (material air))))

(define-param w 0)
(if (> w 0)
    (set! geometry
	  (append geometry
		  (list (make block (center (/ a 2) 0) (size (* 2 t) w infinity)
                        (material air))))))

(define-param fcen (/ (sqrt 0.5) a))
(define-param df 0.2)
(set! sources (list (make source
       (src (make gaussian-src (frequency fcen) (fwidth df))) (component Ez) (center 0 0))))


(set! symmetries (list (make mirror-sym (direction Y))))

(define-param Th 500)
(run-sources+ Th (after-sources (harminv Ez (vector3 0) fcen df)))
(define f (harminv-freq-re (car harminv-results)))
(define Q (harminv-Q (car harminv-results)))
(define Vmode (* 0.25 a a))
(print "ldos0:, " (/ Q Vmode (* 2 pi f pi 0.5)))

(reset-meep)
(define-param T (* 2 Q (/ f)))
(run-sources+ T (dft-ldos f 0 1))
