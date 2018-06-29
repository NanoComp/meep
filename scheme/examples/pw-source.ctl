; This example creates an approximate Ez-polarized planewave in vacuum
; propagating at a 45-degree angle, by using a couple of current sources
; with amplitude exp(ikx) corresponding to the desired planewave.

(define-param s 11) ; the size of the computational cell, not including PML
(define-param dpml 1) ; thickness of PML layers

(define sxy (+ s (* 2 dpml))) ; cell size, including PML
(set! geometry-lattice (make lattice (size sxy sxy no-size)))

(set! pml-layers (list (make pml (thickness dpml))))
(set-param! resolution 10)

; pw-amp is a function that returns the amplitude exp(ik(x+x0)) at a
; given point x.  (We need the x0 because current amplitude functions
; in Meep are defined relative to the center of the current source,
; whereas we want a fixed origin.)  Actually, it is a function of k
; and x0 that returns a function of x ...
(define (pw-amp k x0) (lambda (x)
  (exp (* 0+1i (vector3-dot k (vector3+ x x0))))))

(define-param fcen 0.8) ; pulse center frequency
(define-param df 0.02) ; turn-on bandwidth
(define-param kdir (vector3 1 1)) ; direction of k (length is irrelevant)
(define-param n 1) ; refractive index of material containing the source
(define k (vector3-scale (* 2 pi fcen n)
                         (unit-vector3 kdir))) ; k with correct length

(set! sources
      (list

       ; left
       (make source
	 (src (make continuous-src (frequency fcen) (fwidth df)))
	 (component Ez) (center (* -0.5 s) 0) (size 0 s)
	 (amp-func (pw-amp k (vector3 (* -0.5 s) 0))))

       ; bottom
       (make source
	 (src (make continuous-src (frequency fcen) (fwidth df)))
	 (component Ez) (center 0 (* -0.5 s)) (size s 0)
	 (amp-func (pw-amp k (vector3 0 (* -0.5 s)))))

       ))

(define-param T 400) ; run time
(run-until T (at-end output-efield-z))
