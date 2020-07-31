;; linear-chirped pulse planewave with higher frequencies at the front (down-chirp)

(set-param! resolution 40)

(define-param dpml 2)
(set! pml-layers (list (make pml (thickness dpml) (direction X))))

(define-param sx 40)
(define-param sy 6)

(set! geometry-lattice (make lattice (size (+ sx (* 2 dpml)) sy no-size)))

(define-param v0 1.0) ; pulse center frequency
(define-param a 0.2)  ; Gaussian envelope half-width
(define-param b -0.5) ; linear chirp rate (positive: up-chirp, negative: down-chirp)
(define-param t0 15)  ; peak time

(define chirp (lambda (t) (* (exp (* 0+1i 2 pi v0 (- t t0))) (exp (+ (* (- a) (sqr (- t t0))) (* 0+1i b (sqr (- t t0))))))))

(set! sources (list (make source
                      (src (make custom-src (src-func chirp)))
                      (center (* -0.5 sx) 0 0)
                      (size 0 sy 0)
                      (component Ez))))

(set! k-point (vector3 0 0 0))

(run-until (+ t0 50) (in-volume (volume (center 0 0 0) (size sx sy 0)) (at-every 2.7 output-efield-z)))
