(set-param! resolution 100)

(define dpml 1.0)
(define sx 10)
(set! sx (+ sx (* 2 dpml)))

(define cell (make lattice (size sx no-size no-size)))
(set! geometry-lattice cell)

(define boundary-layers (list (make pml (thickness dpml))))
(set! pml-layers boundary-layers)

(define fcen 1.0)

;; rotation angle of source: CCW around Y axis, 0 degrees along +X axis
(define-param theta 19.4)
(set! theta (deg->rad theta))

;; plane of incidence is XZ
(define k (vector3* fcen (vector3 0 0 (sin theta))))
(set! k-point k)

(set-param! kz-2d "real/imag")

(define pw-src (list (make source
                       (src (make gaussian-src (frequency fcen) (fwidth (* 0.2 fcen))))
                       (component Ey)
                       (center (+ (* -0.5 sx) dpml)))))
(set! sources pw-src)

(define refl (add-flux fcen 0 1 (make flux-region (center (* -0.25 sx)))))

(run-sources+ (stop-when-fields-decayed 50 Ey (vector3 (+ (* -0.5 sx) dpml)) 1e-9))

(define input-flux (list-ref (get-fluxes refl) 0))

(save-flux "refl-flux" refl)

(reset-meep)

(set! geometry-lattice cell)
(set! pml-layers boundary-layers)
(set! k-point k)
(set! sources pw-src)

;; add a block with n=3.5 for the air-dielectric interface
(set! geometry (list (make block
                       (center (* 0.25 sx))
                       (size (* 0.5 sx) infinity infinity)
                       (material (make medium (index 3.5))))))

(define refl (add-flux fcen 0 1 (make flux-region (center (* -0.25 sx)))))
(load-minus-flux "refl-flux" refl)

(run-sources+ (stop-when-fields-decayed 50 Ey (vector3 (+ (* -0.5 sx) dpml)) 1e-9))

(define refl-flux (list-ref (get-fluxes refl) 0))

(print "refl:, " (- (/ refl-flux input-flux)) "\n")
