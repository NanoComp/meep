(set-param! resolution 400) ; pixels/um

(define-param sz 10)
(set! geometry-lattice (make lattice (size no-size no-size sz)))
(set! dimensions 1)

(define lambda-min 0.4)
(define lambda-max 0.8)
(define fmax (/ lambda-min))
(define fmin (/ lambda-max))
(define fcen (* 0.5 (+ fmax fmin)))
(define df (- fmax fmin))

(define dpml 1.0)
(set! pml-layers (list (make pml (thickness dpml))))

(set! k-point (vector3 0 0 0))

(set! sources (list (make source (src (make gaussian-src (frequency fcen) (fwidth df))) (component Ex) (center 0 0 (+ (* -0.5 sz) dpml)))))

(define-param empty? true)

(if (not empty?)
    (set! geometry (list (make block (size infinity infinity (* 0.5 sz)) (center 0 0 (* 0.25 sz)) (material fused-quartz)))))

(define nfreq 50)
(define refl (add-flux fcen df nfreq (make flux-region (center 0 0 (* -0.25 sz)))))

(if (not empty?) (load-minus-flux "refl-flux" refl))

(run-sources+ (stop-when-fields-decayed 50 Ex (vector3 0 0 (+ (* -0.5 sz) dpml)) 1e-9))

(if empty? (save-flux "refl-flux" refl))

(display-fluxes refl)
