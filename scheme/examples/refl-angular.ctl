(set-param! resolution 200) ; pixels/um

(define-param sz 10)
(define-param dpml 1)
(set! sz (+ sz (* 2 dpml)))
(set! pml-layers (list (make pml (thickness dpml))))

(set! geometry-lattice (make lattice (size no-size no-size sz)))

(define lambda-cen 0.6)
(define fcen (/ lambda-cen))
(define df (* 0.2 fcen))

; rotation angle of source: CCW relative to y axis
(define-param theta 0)
(define theta-r (deg->rad theta))

(set! dimensions (if (= theta-r 0) 1 3))

; plane of incidence is xz
(set! k-point (vector3* fcen (vector3 (sin theta-r) 0 (cos theta-r))))

(set! sources (list (make source (src (make gaussian-src (frequency fcen) (fwidth df)))
			         (component Ex) (center 0 0 (+ (* -0.5 sz) dpml)))))

(define-param empty? true)

(if (not empty?)
    (set! geometry (list (make block (size infinity infinity (* 0.5 sz)) (center 0 0 (* 0.25 sz)) (material (make medium (index 3.5)))))))

(define refl (add-flux fcen 0 1 (make flux-region (center 0 0 (* -0.25 sz)))))

(if (not empty?) (load-minus-flux "refl-flux" refl))

(run-sources+ (stop-when-fields-decayed 50 Ex (vector3 0 0 (+ (* -0.5 sz) dpml)) 1e-9))

(if empty? (save-flux "refl-flux" refl))

(display-fluxes refl)
